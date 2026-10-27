# ==============================================================================
# 08_water_pans_analysis.R  -  Water-pan catchments and per-pan diagnostics
# ==============================================================================
# PURPOSE
#   Snap each stream-fed pan to the stream network, delineate its upslope catchment, compute seasonal catchment
#   metrics, accessibility sheds (3 / 5 / 10 / 15 km), the ranch's water gap, refugia value, and the
#   status-specific diagnostics of DESCRIPTION.md 6.6:
#     - existing_degraded pans: four failure-mode hypotheses, ranked by EO evidence,
#     - the existing_functional_under_pressure pan: the demand it absorbs,
#     - the existing_desilting pan: the PRE-INTERVENTION baseline for a 6 / 12 / 24 month re-run.
#
# INPUTS
#   data/water_pans/water_pans.csv (+ image-read centres and stand-off uncertainty from amaka_field_points.csv)
#   outputs/rasters/hydrology/_buffered/* (03), outputs/rasters/{rpi_*,pml_*}.tif, outputs/rasters/s2/*.tif,
#   data/gee_exports/{worldcover,chirps,smap,modis_lst}*.tif, outputs/tables/cells_*.csv (03-06),
#   data/boundary/context/Tsavo_Settled_Areas.shp, Tsavo_{East,West}_NP.shp
#
# OUTPUTS
#   outputs/tables/pan_catchment_metrics.csv / .xlsx   pan x season (static metrics repeated on each row)
#   outputs/tables/pan_catchment_static.csv            one row per pan
#   outputs/tables/pan_diagnostics.csv                 per-pan diagnosis (failure modes / pressure / baseline)
#   outputs/tables/pan_hypothesis_scores.csv           the four hypothesis scores per stream-fed pan + the functional control
#   outputs/tables/pan_sheds.csv, water_gap_summary.csv
#   outputs/tables/desilting_baseline_PRE_INTERVENTION_<date>.csv   (or a clear NO_DESILTING_COORDINATES marker)
#   outputs/rasters/pan_catchments.tif (one band per stream-fed pan) + pan_catchments.shp (all pans, with the method),
#       water_gap.tif (primary radius), water_gap_all_radii.tif, dist_to_functional_pan_m.tif, dist_to_any_pan_m.tif
#
# RUNTIME NOTES
#   * About 20 s. Pour points start from the image-read centre where there is one (Cfg.Snap.Img.m = 60 m snap), else
#     from the recorded coordinates with the stand-off uncertainty as snap distance (brief: 100 m; known_issues F11).
#     Snapping uses WhiteboxTools' Jenson method onto the calibrated stream network (03).
#   * Only STREAM-FED pans (Cfg.Pan.StreamFed.Types: earth dams, natural pans) are snapped and given a catchment.
#     The brief snaps every pan, but a rock-catchment reservoir or borehole is not filled by a channel: the HQ reservoir's
#     centre sits ON a Strahler-3 channel draining 13,000+ ha, which is the river's basin, not the structure's supply.
#     Those pans get a context disc (Cfg.Pan.Context.m) for site context and their catchment hypotheses are not scored
#     (known_issues F12).
#   * Each pan's catchment is delineated on its own, so a catchment includes any pan upslope of it.
#   * If water_pans.csv has no rows, a NO_PANS_YET marker is written and the script stops cleanly.
#   * Pans outside the Amaka AOI (in the neighbouring Kambanga polygon) ARE analysed and flagged inside_aoi = FALSE.
#     The water gap and everything "inside Amaka" use the AOI only.
#   * The hypothesis scores are SCREENING evidence from satellite data, not diagnoses. Each score is the normal CDF of
#     (catchment mean - ranch mean) / SD of ranch means over windows of the same size as the catchment (so a 560 ha
#     catchment is compared with 560 ha windows, not with single 21 ha cells); 0.5 = typical of the ranch. The functional
#     dam is scored by the same tests as a CONTROL: a hypothesis that the working dam also "fails" does not discriminate.
#     Wall integrity, spillway, leakage and sediment depth are invisible to EO.
#   * Refugia value uses distance to the ranch boundary and to settled areas; there is no local road layer, so the roads
#     term is omitted (roads_available = FALSE).
# ==============================================================================

if (!exists("Setup.Loaded")) source("R/00_setup.R")
whitebox::wbt_verbose(FALSE)

Run.Pan.Analysis <- function() {

  # ---- PHASE 1: Pans and inputs ---------------------------------------------------------
  log_phase("PHASE 1: pans and inputs")
  Pans.All <- read_pans()
  if (nrow(Pans.All) == 0) {
    writeLines("NO_PANS_YET: data/water_pans/water_pans.csv has no rows. Add the existing pans and re-run.", file.path(Dir.Out.Tables, "NO_PANS_YET.txt"))
    message("  water_pans.csv is empty: wrote NO_PANS_YET and stopped cleanly."); return(invisible(NULL))
  }
  Pans.All$stream_fed <- Pans.All$type %in% Cfg.Pan.StreamFed.Types
  nrm <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
  wb  <- function(name) nrm(file.path(Dir.Out.Hydro.Buf, name))
  AOI.Sf <- read_aoi(); AOI.U <- sf::st_union(sf::st_geometry(AOI.Sf)); Aoi.ha <- as.numeric(sf::st_area(AOI.U)) / 1e4
  Aoi.Vec <- terra::vect(AOI.Sf)
  Dem <- terra::rast(file.path(Dir.DEM, "amaka_dem_utm37s.tif"))
  rd  <- function(name) { r <- terra::rast(wb(name)); if (is.na(terra::crs(r, describe = TRUE)$code)) terra::crs(r) <- paste0("EPSG:", Cfg.CRS.Metric); r }
  R.Slope <- rd("slope_deg.tif"); R.Twi <- rd("twi.tif"); R.Tpi <- rd("tpi_500m.tif"); R.Acc <- rd("flow_accum_d8.tif"); R.Ord <- rd("strahler.tif")
  Cell.ha <- prod(terra::res(R.Acc)) / 1e4
  Located <- Pans.All[Pans.All$located, ]
  Pans.Sf <- sf::st_transform(sf::st_as_sf(Located, coords = c("start_lon", "start_lat"), crs = Cfg.CRS.Geo, remove = FALSE), Cfg.CRS.Metric)
  message(sprintf("  %d pans in the file, %d with coordinates (%s)", nrow(Pans.All), nrow(Located),
                  paste(sprintf("%s: %s%s", Located$pan_id, sub("existing_", "", Located$status), ifelse(Located$stream_fed, "", " [not stream-fed]")), collapse = "; ")))
  message("  without coordinates: ", ifelse(any(!Pans.All$located), paste(Pans.All$pan_id[!Pans.All$located], collapse = ", "), "none"))
  one_geom <- function(g) sf::st_sfc(sf::st_geometry(g)[[1]], crs = Cfg.CRS.Metric)      # first geometry as a length-1 sfc


  # ---- PHASE 2: Snap, delineate, polygonise ---------------------------------------------------
  log_phase("PHASE 2: catchments")
  Dir.Pan <- file.path(Dir.Out.Hydro.Buf, "pans"); dir.create(Dir.Pan, showWarnings = FALSE)
  Catch <- list(); Snap <- list(); Ws.R <- list()
  for (i in seq_len(nrow(Pans.Sf))) {
    id <- Pans.Sf$pan_id[i]
    acc.start <- terra::extract(R.Acc, terra::vect(Pans.Sf[i, ]))[[2]]       # cells draining to the START cell, before any snapping
    if (!Pans.Sf$stream_fed[i]) {                                             # rock catchment, borehole, ...: context disc only
      disc <- sf::st_buffer(sf::st_geometry(Pans.Sf[i, ]), Cfg.Pan.Context.m)
      Catch[[id]] <- list(poly = disc, area_ha = as.numeric(sf::st_area(disc)) / 1e4, touches_edge = FALSE,
                          method = sprintf("context disc %g m (not stream-fed)", Cfg.Pan.Context.m))
      Snap[[id]]  <- list(pt = Pans.Sf[i, ], moved_m = 0, snap_used_m = 0, acc = acc.start, acc.start = acc.start)
      message(sprintf("  %s: type '%s' is not stream-fed -> no snapping; %g m context disc. The start cell drains %s ha (%s cells)%s",
                      id, Pans.Sf$type[i], Cfg.Pan.Context.m, format(round(acc.start * Cell.ha), big.mark = ","), format(acc.start, big.mark = ","),
                      ifelse(acc.start >= Cfg.Stream.Threshold, ": it sits ON a mapped channel", "")))
      next
    }
    f.pt <- file.path(Dir.Pan, paste0(id, "_start.shp")); f.sn <- file.path(Dir.Pan, paste0(id, "_snapped.shp")); f.ws <- file.path(Dir.Pan, paste0(id, "_catchment.tif"))
    write_shp(Pans.Sf[i, "pan_id"], f.pt)
    snap.m <- Pans.Sf$snap_m[i]
    repeat {                                               # widen the snap distance (up to 300 m) until the point is on the stream network
      whitebox::wbt_jenson_snap_pour_points(pour_pts = nrm(f.pt), streams = wb("streams.tif"), output = nrm(f.sn), snap_dist = snap.m)
      sp <- sf::st_read(f.sn, quiet = TRUE); sf::st_crs(sp) <- Cfg.CRS.Metric
      acc <- if (nrow(sp)) terra::extract(R.Acc, terra::vect(sp))[[2]] else NA_real_
      if ((is.finite(acc) && acc >= Cfg.Stream.Threshold) || snap.m >= 300) break
      snap.m <- snap.m + 60
    }
    if (!nrow(sp)) stop(id, ": the pour point could not be snapped within 300 m of the stream network", call. = FALSE)
    moved <- as.numeric(sf::st_distance(sp, Pans.Sf[i, ]))
    whitebox::wbt_watershed(d8_pntr = wb("d8_pointer.tif"), pour_pts = nrm(f.sn), output = nrm(f.ws))
    ws <- terra::rast(f.ws); if (is.na(terra::crs(ws, describe = TRUE)$code)) terra::crs(ws) <- paste0("EPSG:", Cfg.CRS.Metric)
    ws <- terra::ifel(ws > 0, 1, NA)
    pg <- sf::st_union(sf::st_as_sf(terra::as.polygons(ws, dissolve = TRUE)))
    ring <- terra::rast(ws); terra::values(ring) <- 0; ring[1, ] <- 1; ring[nrow(ring), ] <- 1; ring[, 1] <- 1; ring[, ncol(ring)] <- 1
    Catch[[id]] <- list(poly = one_geom(sf::st_set_crs(pg, Cfg.CRS.Metric)), method = "D8 watershed from the snapped pour point",
                        area_ha = sum(!is.na(terra::values(ws))) * prod(terra::res(ws)) / 1e4,
                        touches_edge = terra::global(!is.na(ws) & ring == 1, "max", na.rm = TRUE)[[1]] == 1)
    Snap[[id]]  <- list(pt = sp, moved_m = moved, snap_used_m = snap.m, acc = acc, acc.start = acc.start)
    Ws.R[[id]]  <- ws
    message(sprintf("  %s: snapped %.0f m (snap distance %g m), %s upslope cells at the pour point, catchment %.0f ha%s",
                    id, moved, snap.m, format(acc, big.mark = ","), Catch[[id]]$area_ha, ifelse(Catch[[id]]$touches_edge, " [TOUCHES THE DEM EDGE: truncated]", "")))
  }
  Catch.Sf <- sf::st_sf(pan_id = names(Catch), method = vapply(Catch, function(x) x$method, ""), geometry = do.call(c, lapply(Catch, function(x) x$poly)))
  write_shp(Catch.Sf, file.path(Dir.Out.Rasters, "pan_catchments.shp"))
  if (length(Ws.R)) {
    Pc <- do.call(c, unname(Ws.R)); names(Pc) <- names(Ws.R)
    terra::writeRaster(Pc, file.path(Dir.Out.Rasters, "pan_catchments.tif"), overwrite = TRUE, datatype = "INT1U")
  }


  # ---- PHASE 3: Rasters used for the catchment metrics -------------------------------------------
  log_phase("PHASE 3: metric layers")
  rr <- function(...) terra::rast(file.path(Dir.Out.Rasters, ...))
  Wc <- terra::rast(gee_file("worldcover_2021_envelope"))
  Pml <- lapply(setNames(Seasons.All, Seasons.All), function(s) setNames(lapply(c("gpp", "et", "esfrac", "ecfrac", "wue"), function(v) rr(sprintf("pml_%s_%s.tif", v, s))), c("gpp", "et", "esfrac", "ecfrac", "wue")))
  clim_mean <- function(prefix, s) terra::app(gee_stack(prefix, s), "mean", na.rm = TRUE)
  Clim <- lapply(setNames(Seasons.All, Seasons.All), function(s) list(rain = clim_mean("chirps_precip_sum", s), sm_root = clim_mean("smap_rootzone_mean", s),
                                                                       sm_surf = clim_mean("smap_surface_mean", s), lst = clim_mean("modis_lst_day_mean", s)))
  S2.Bsi.Dry <- rr("s2", "s2_bsi_DRY_mean.tif"); S2.Bsi.Lr <- rr("s2", "s2_bsi_LR_mean.tif"); S2.Ndvi.Dry <- rr("s2", "s2_ndvi_DRY_mean.tif")
  Site.Rad.m <- 150          # slope and TPI "at the site" = mean within this radius of the pour point (one 30 m DSM cell is too noisy)
  lab.sterile <- paste0("shrub-dominated (ratio > ", Cfg.Sterile.Ratio, ")")

  one <- function(r, poly) zonal_stat(r, sf::st_sf(geometry = poly))[[1]]
  site_mean <- function(r, pt) terra::extract(r, terra::vect(sf::st_buffer(sf::st_geometry(pt), Site.Rad.m)), fun = "mean", na.rm = TRUE)[[2]]
  catch_metrics <- function(id) {
    poly <- Catch[[id]]$poly; pt <- Snap[[id]]$pt
    fr <- exactextractr::exact_extract(Wc, sf::st_transform(sf::st_sf(geometry = poly), terra::crs(Wc)), fun = "frac", progress = FALSE)
    g <- function(k) if (paste0("frac_", k) %in% names(fr)) fr[[paste0("frac_", k)]][1] else 0
    trees <- g(10); shrub <- g(20)
    st <- tibble::tibble(
      pan_id = id, catchment_method = Catch[[id]]$method, catch_ha = Catch[[id]]$area_ha, catch_touches_dem_edge = Catch[[id]]$touches_edge,
      snap_moved_m = Snap[[id]]$moved_m, snap_distance_used_m = Snap[[id]]$snap_used_m, pour_flow_acc_cells = Snap[[id]]$acc,
      start_cell_upslope_ha = Snap[[id]]$acc.start * Cell.ha,
      pour_strahler = terra::extract(R.Ord, terra::vect(pt))[[2]], pour_elev_m = terra::extract(Dem, terra::vect(pt))[[2]],
      site_slope_deg = site_mean(R.Slope, pt), tpi_site_m = site_mean(R.Tpi, pt),
      slope_mean_deg = one(R.Slope, poly), twi_mean = one(R.Twi, poly), tpi_mean_m = one(R.Tpi, poly),
      rpi_now = one(rr("rpi_now.tif"), poly), rpi_sen_slope = one(rr("rpi_sen_slope.tif"), poly), rpi_degrading_share = one(rr("rpi_degrading.tif"), poly),
      wc_trees = trees, wc_shrub = shrub, wc_grass = g(30), wc_bare = g(60),
      shrub_tree_ratio = ifelse(trees >= 0.01, shrub / trees, NA_real_),
      shrub_tree_class = dplyr::case_when(trees < 0.01 ~ "no trees (<1%)", shrub / trees > Cfg.Sterile.Ratio ~ lab.sterile, shrub / trees >= 1 ~ "mixed (1-5)", TRUE ~ "tree-dominated (<1)"),
      bsi_DRY = one(S2.Bsi.Dry, poly), bsi_LR = one(S2.Bsi.Lr, poly), ndvi_DRY = one(S2.Ndvi.Dry, poly))
    se <- dplyr::bind_rows(lapply(Seasons.All, function(s) tibble::tibble(
      pan_id = id, season = s, gpp_gC_m2_d = one(Pml[[s]]$gpp, poly), et_mm_d = one(Pml[[s]]$et, poly), es_frac = one(Pml[[s]]$esfrac, poly),
      ec_frac = one(Pml[[s]]$ecfrac, poly), wue_gC_per_mm = one(Pml[[s]]$wue, poly), soil_moisture_rootzone = one(Clim[[s]]$sm_root, poly),
      soil_moisture_surface = one(Clim[[s]]$sm_surf, poly), rainfall_mm = one(Clim[[s]]$rain, poly), lst_day_C = one(Clim[[s]]$lst, poly))))
    list(static = st, seasonal = se)
  }
  Res <- lapply(names(Catch), catch_metrics)
  Static <- dplyr::bind_rows(lapply(Res, `[[`, "static")); Seasonal <- dplyr::bind_rows(lapply(Res, `[[`, "seasonal"))
  Static <- dplyr::left_join(Static, Pans.All[, c("pan_id", "stream_fed", "start_source")], by = "pan_id")
  Annual.Rain <- Seasonal$rainfall_mm[Seasonal$season == "ANNUAL"]; names(Annual.Rain) <- Seasonal$pan_id[Seasonal$season == "ANNUAL"]
  Static$annual_rainfall_mm <- unname(Annual.Rain[Static$pan_id])
  Static$runoff_yield_m3_per_yr <- ifelse(Static$stream_fed, Static$catch_ha * 1e4 * (Static$annual_rainfall_mm / 1000) * Cfg.Runoff.Coef, NA_real_)   # runoff coefficient ASSUMED


  # ---- PHASE 4: Distances, sheds, water gap --------------------------------------------------------------
  log_phase("PHASE 4: distances, accessibility sheds and the water gap")
  Settled <- sf::st_union(sf::st_geometry(sf::st_read(file.path(Dir.Boundary.Ctx, "Tsavo_Settled_Areas.shp"), quiet = TRUE)))
  Np <- sf::st_union(c(sf::st_geometry(sf::st_read(file.path(Dir.Boundary.Ctx, "Tsavo_East_NP.shp"), quiet = TRUE)), sf::st_geometry(sf::st_read(file.path(Dir.Boundary.Ctx, "Tsavo_West_NP.shp"), quiet = TRUE))))
  Pt <- Pans.Sf; inside <- lengths(sf::st_intersects(Pt, AOI.U)) > 0
  d.bnd <- as.numeric(sf::st_distance(Pt, sf::st_boundary(AOI.U)))
  Dist <- tibble::tibble(pan_id = Pt$pan_id, inside_aoi = inside, dist_to_boundary_m = ifelse(inside, d.bnd, -d.bnd),
                         dist_settled_m = as.numeric(sf::st_distance(Pt, Settled)), dist_tsavo_np_m = as.numeric(sf::st_distance(Pt, Np)))
  dm <- sf::st_distance(Pt); diag(dm) <- Inf
  Dist$nearest_pan_m <- if (nrow(Pt) > 1) apply(dm, 1, min) else NA_real_
  Dist$nearest_pan_id <- if (nrow(Pt) > 1) Pt$pan_id[apply(dm, 1, which.min)] else NA_character_
  Static <- dplyr::left_join(Static, Dist, by = "pan_id")

  # refugia value: mean of the saturating distance terms that exist (no road layer locally)
  sat <- function(d, ref.km) pmin(pmax(d / (ref.km * 1000), 0), 1)
  Static$refugia_boundary <- sat(Static$dist_to_boundary_m, Cfg.Refugia.Ref.km[["boundary"]])
  Static$refugia_settled  <- sat(Static$dist_settled_m, Cfg.Refugia.Ref.km[["settled"]])
  Static$refugia_value <- rowMeans(Static[, c("refugia_boundary", "refugia_settled")]); Static$roads_available <- FALSE

  Radii <- Cfg.Shed.Radii.km
  Func <- Pt[Pt$functional, ]; Any <- Pt
  shed <- function(pts, r.km) if (nrow(pts)) sf::st_union(sf::st_buffer(sf::st_geometry(pts), r.km * 1000)) else sf::st_sfc(sf::st_polygon(), crs = Cfg.CRS.Metric)
  ha <- function(g) if (length(g) == 0 || all(sf::st_is_empty(g))) 0 else as.numeric(sum(sf::st_area(g))) / 1e4
  Sheds <- dplyr::bind_rows(lapply(seq_len(nrow(Pt)), function(i) dplyr::bind_rows(lapply(Radii, function(r) {
    b <- sf::st_buffer(sf::st_geometry(Pt[i, ]), r * 1000)
    tibble::tibble(pan_id = Pt$pan_id[i], radius_km = r, shed_ha = ha(b), shed_in_aoi_ha = ha(suppressWarnings(sf::st_intersection(b, AOI.U))),
                   share_of_aoi_pct = 100 * ha(suppressWarnings(sf::st_intersection(b, AOI.U))) / Aoi.ha)
  }))))
  Gap <- dplyr::bind_rows(lapply(Radii, function(r) {
    wf <- shed(Func, r); wa <- shed(Any, r)
    in.f <- if (nrow(Func)) ha(suppressWarnings(sf::st_intersection(wf, AOI.U))) else 0; in.a <- ha(suppressWarnings(sf::st_intersection(wa, AOI.U)))
    tibble::tibble(radius_km = r, aoi_ha = Aoi.ha, watered_by_functional_ha = in.f, gap_functional_ha = Aoi.ha - in.f, gap_functional_pct = 100 * (Aoi.ha - in.f) / Aoi.ha,
                   watered_by_any_pan_ha = in.a, gap_any_pan_ha = Aoi.ha - in.a, gap_any_pan_pct = 100 * (Aoi.ha - in.a) / Aoi.ha)
  }))
  readr::write_csv(Sheds, file.path(Dir.Out.Tables, "pan_sheds.csv")); readr::write_csv(Gap, file.path(Dir.Out.Tables, "water_gap_summary.csv"))
  print(as.data.frame(Gap |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)

  # distance rasters on the DEM grid. water_gap = 1 where an AOI cell is farther than the radius from EVERY functional pan
  d.func <- if (nrow(Func)) terra::distance(Dem, terra::vect(Func)) else Dem * 0 + 1e6
  d.any  <- terra::distance(Dem, terra::vect(Any))
  names(d.func) <- "dist_to_functional_pan_m"; names(d.any) <- "dist_to_any_pan_m"
  terra::writeRaster(d.func, file.path(Dir.Out.Rasters, "dist_to_functional_pan_m.tif"), overwrite = TRUE)
  terra::writeRaster(d.any,  file.path(Dir.Out.Rasters, "dist_to_any_pan_m.tif"), overwrite = TRUE)
  Gap.R <- terra::rast(lapply(Radii, function(r) terra::mask(terra::crop((d.func > r * 1000) * 1, Aoi.Vec), Aoi.Vec)))
  names(Gap.R) <- sprintf("gap_%gkm", Radii)
  terra::writeRaster(Gap.R, file.path(Dir.Out.Rasters, "water_gap_all_radii.tif"), overwrite = TRUE, datatype = "INT1U")
  Gap.Main <- Gap.R[[sprintf("gap_%gkm", Cfg.WaterGap.Radius.km)]]; names(Gap.Main) <- "water_gap"
  terra::writeRaster(Gap.Main, file.path(Dir.Out.Rasters, "water_gap.tif"), overwrite = TRUE, datatype = "INT1U")


  # ---- PHASE 5: Status-specific diagnostics ------------------------------------------------------------------------
  log_phase("PHASE 5: diagnostics")
  Grid <- get_analysis_grid(); Cells <- Grid$cells; Grid <- Grid$grid
  Ref <- dplyr::left_join(Cells, readr::read_csv(file.path(Dir.Out.Tables, "cells_terrain.csv"), show_col_types = FALSE), by = "cell_id") |>
    dplyr::left_join(readr::read_csv(file.path(Dir.Out.Tables, "cells_rpi.csv"), show_col_types = FALSE), by = "cell_id") |>
    dplyr::left_join(readr::read_csv(file.path(Dir.Out.Tables, "cells_pml.csv"), show_col_types = FALSE), by = "cell_id") |>
    dplyr::left_join(readr::read_csv(file.path(Dir.Out.Tables, "cells_landcover.csv"), show_col_types = FALSE), by = "cell_id") |> dplyr::filter(in_aoi)
  Cell.Grid.ha <- stats::median(Cells$area_ha)
  win <- function(catch.ha) max(3, 2 * floor(sqrt(catch.ha / Cell.Grid.ha) / 2) + 1)          # odd window (in grid cells) with the catchment's area

  # SD of window means of one metric over the ranch: how much a catchment-sized patch varies, which is the right yardstick
  Sd.Cache <- list()
  ref_sd <- function(col, w) {
    key <- paste(col, w)
    if (is.null(Sd.Cache[[key]])) {
      r <- terra::rast(Grid); v <- rep(NA_real_, terra::ncell(r)); v[Ref$cell_id] <- Ref[[col]]; terra::values(r) <- v
      n <- terra::focal(terra::ifel(is.na(r), 0, 1), w = w, fun = "sum")                       # AOI cells inside each window
      m <- terra::ifel(n >= ceiling(0.5 * w * w), terra::focal(r, w = w, fun = "mean", na.rm = TRUE), NA)
      Sd.Cache[[key]] <<- stats::sd(terra::values(m)[Ref$cell_id], na.rm = TRUE)
    }
    Sd.Cache[[key]]
  }
  ev <- function(x, col, w, direction = 1) {
    z <- (x - mean(Ref[[col]], na.rm = TRUE)) / ref_sd(col, w)
    ifelse(is.finite(z), stats::pnorm(direction * z), 0.5)
  }
  Tpi.Aoi <- terra::mask(R.Tpi, Aoi.Vec)                       # 30 m TPI inside the AOI: the yardstick for the TPI at a pan site
  Tpi.M <- terra::global(Tpi.Aoi, "mean", na.rm = TRUE)[[1]]; Tpi.S <- terra::global(Tpi.Aoi, "sd", na.rm = TRUE)[[1]]

  SD <- dplyr::left_join(Static, Pans.All[, c("pan_id", "name", "status")], by = "pan_id")
  Es.Dry <- Seasonal$es_frac[Seasonal$season == "DRY"]; names(Es.Dry) <- Seasonal$pan_id[Seasonal$season == "DRY"]
  SD$esfrac_DRY <- unname(Es.Dry[SD$pan_id])
  Ref.Func <- SD[SD$status == "existing_functional_under_pressure" & SD$stream_fed, ]
  ref.ha <- if (nrow(Ref.Func)) min(Ref.Func$catch_ha) else NA_real_

  Hyp <- dplyr::bind_rows(lapply(which(SD$stream_fed), function(i) {
    d <- SD[i, ]; w <- win(d$catch_ha)
    h1 <- if (is.finite(ref.ha)) max(0, min(1, 1 - d$catch_ha / ref.ha)) else NA_real_
    c2 <- c(bare = ev(d$bsi_DRY, "bsi_DRY", w), Es_ET = ev(d$esfrac_DRY, "esfrac_DRY", w), RPI_decline = ev(d$rpi_sen_slope, "rpi_sen_slope", w, -1))
    c3 <- c(low_RPI = ev(d$rpi_now, "rpi_now", w, -1), low_TWI = ev(d$twi_mean, "twi_mean", w, -1), steep = ev(d$slope_mean_deg, "slope_mean", w))
    c4 <- c(site_slope_off = 1 - trapezoid(d$site_slope_deg, Cfg.Terrain.Slope.Deg), TPI_ridge = stats::pnorm((d$tpi_site_m - Tpi.M) / Tpi.S))
    h2 <- mean(c2); h3 <- mean(c3); h4 <- mean(c4)
    comp <- function(x) paste0(" [parts: ", paste(sprintf("%s %.2f", names(x), x), collapse = ", "), "]")
    tibble::tibble(pan_id = d$pan_id, hypothesis = c("H1 catchment too small to fill", "H2 siltation from upslope overgrazing", "H3 flashy runoff from degraded vegetation", "H4 terrain geometry wrong for retention"),
                   score = c(h1, h2, h3, h4),
                   evidence = c(sprintf("catchment %.0f ha%s", d$catch_ha, ifelse(is.finite(ref.ha), sprintf(" = %.0f%% of the functional dam's %.0f ha; estimated yield %s m3/yr at runoff coefficient %.2f", 100 * d$catch_ha / ref.ha, ref.ha, format(round(d$runoff_yield_m3_per_yr), big.mark = ","), Cfg.Runoff.Coef), " (no functional reference)")),
                                sprintf("catchment dry-season BSI %.3f (ranch %.3f), Es/ET dry %.2f (ranch %.2f), RPI slope %+.4f/yr (ranch %+.4f)%s", d$bsi_DRY, mean(Ref$bsi_DRY, na.rm = TRUE), d$esfrac_DRY, mean(Ref$esfrac_DRY, na.rm = TRUE), d$rpi_sen_slope, mean(Ref$rpi_sen_slope, na.rm = TRUE), comp(c2)),
                                sprintf("catchment RPI %.3f (ranch %.3f), TWI %.2f (ranch %.2f), mean slope %.2f deg (ranch %.2f)%s", d$rpi_now, mean(Ref$rpi_now, na.rm = TRUE), d$twi_mean, mean(Ref$twi_mean, na.rm = TRUE), d$slope_mean_deg, mean(Ref$slope_mean, na.rm = TRUE), comp(c3)),
                                sprintf("slope within %g m of the pour point %.2f deg (viable %.1f-%.1f), TPI there %+.2f m (ridge > 0, depression < 0)%s", Site.Rad.m, d$site_slope_deg, Cfg.Terrain.Slope.Deg[2], Cfg.Terrain.Slope.Deg[3], d$tpi_site_m, comp(c4))))
  }))
  Ctl <- Hyp[Hyp$pan_id %in% Ref.Func$pan_id, ] |> dplyr::group_by(hypothesis) |> dplyr::summarise(control_score = mean(score), .groups = "drop")
  Hyp <- dplyr::left_join(Hyp, Ctl, by = "hypothesis"); Hyp$excess_over_functional <- Hyp$score - Hyp$control_score

  Diag <- dplyr::bind_rows(lapply(seq_len(nrow(Pans.All)), function(i) {
    p <- Pans.All[i, ]
    if (!p$located) return(tibble::tibble(pan_id = p$pan_id, name = p$name, status = p$status, inside_aoi = NA, diagnosis_type = "no coordinates", summary = "No coordinates: no catchment, shed or baseline can be computed. Supply coordinates (GPS or from imagery)."))
    s <- SD[SD$pan_id == p$pan_id, ]; dist <- Dist[Dist$pan_id == p$pan_id, ]
    if (p$status == "existing_degraded" && !p$stream_fed) {
      on.ch <- s$start_cell_upslope_ha >= Cfg.Stream.Threshold * Cell.ha
      return(tibble::tibble(pan_id = p$pan_id, name = p$name, status = p$status, inside_aoi = dist$inside_aoi, diagnosis_type = "not stream-fed: catchment hypotheses not scored",
                            summary = paste0(sprintf("'%s' (type '%s') is not a stream-fed dam, so the four catchment hypotheses are not scored. ", p$name, p$type),
                                             if (on.ch) sprintf("Its position sits on a mapped channel (Strahler %s, at least %s ha upslope): either the structure lies on the river or the 30 m surface model cannot separate it from the channel. ", s$pour_strahler, format(round(s$start_cell_upslope_ha), big.mark = ",")) else "",
                                             sprintf("What EO cannot see (wall, spillway, leakage, sediment depth, the outcrop that feeds it) needs a field check; %g m context metrics are in pan_catchment_metrics.csv.", Cfg.Pan.Context.m))))
    }
    if (p$status == "existing_degraded") {
      h <- Hyp[Hyp$pan_id == p$pan_id, ] |> dplyr::arrange(dplyr::desc(score))
      ctl <- function(k) if (is.finite(h$control_score[k])) sprintf("; functional dam %.2f on the same test", h$control_score[k]) else "; no functional reference"
      return(tibble::tibble(pan_id = p$pan_id, name = p$name, status = p$status, inside_aoi = dist$inside_aoi, diagnosis_type = "failure-mode hypotheses (ranked)",
                            rank1 = h$hypothesis[1], rank1_score = h$score[1], rank1_functional_control = h$control_score[1],
                            rank2 = h$hypothesis[2], rank2_score = h$score[2], rank2_functional_control = h$control_score[2],
                            rank3 = h$hypothesis[3], rank3_score = h$score[3], rank3_functional_control = h$control_score[3],
                            rank4 = h$hypothesis[4], rank4_score = h$score[4], rank4_functional_control = h$control_score[4],
                            summary = paste0("Most supported: ", h$hypothesis[1], " (", round(h$score[1], 2), ctl(1), "): ", h$evidence[1], ". Next: ", h$hypothesis[2], " (", round(h$score[2], 2), ctl(2), ").")))
    }
    if (p$status == "existing_functional_under_pressure") {
      others <- Pt[Pt$pan_id != p$pan_id, ]
      dd <- if (nrow(others)) as.numeric(sf::st_distance(Pt[Pt$pan_id == p$pan_id, ], others)) else numeric()
      n.f <- vapply(c(5, 10, 15), function(r) sum(dd[others$functional] <= r * 1000), 0)
      n.d <- vapply(c(5, 10, 15), function(r) sum(dd[!others$functional] <= r * 1000), 0)
      sh <- Sheds[Sheds$pan_id == p$pan_id, ]
      return(tibble::tibble(pan_id = p$pan_id, name = p$name, status = p$status, inside_aoi = dist$inside_aoi, diagnosis_type = "pressure (demand absorbed)",
                            other_functional_within_5km = n.f[1], other_functional_within_10km = n.f[2], other_functional_within_15km = n.f[3],
                            degraded_or_desilting_within_5km = n.d[1], degraded_or_desilting_within_10km = n.d[2], degraded_or_desilting_within_15km = n.d[3],
                            aoi_share_within_10km_pct = sh$share_of_aoi_pct[sh$radius_km == 10], aoi_share_within_15km_pct = sh$share_of_aoi_pct[sh$radius_km == 15],
                            summary = sprintf("%d other functional pans within 5 km, %d within 10 km, %d within 15 km; %d non-functional (degraded or desilting) pans with coordinates within 15 km. It reaches %.0f%% of the AOI at 10 km and %.0f%% at 15 km%s.",
                                              n.f[1], n.f[2], n.f[3], n.d[3], sh$share_of_aoi_pct[sh$radius_km == 10], sh$share_of_aoi_pct[sh$radius_km == 15],
                                              ifelse(n.f[3] == 0, ", so it is the only reliable dry-season water for that area", ""))))
    }
    tibble::tibble(pan_id = p$pan_id, name = p$name, status = p$status, inside_aoi = dist$inside_aoi, diagnosis_type = "baseline", summary = "Baseline table written (see desilting_baseline_PRE_INTERVENTION_*.csv).")
  }))
  readr::write_csv(Diag, file.path(Dir.Out.Tables, "pan_diagnostics.csv"))
  readr::write_csv(Hyp, file.path(Dir.Out.Tables, "pan_hypothesis_scores.csv"))
  print(as.data.frame(Hyp |> dplyr::mutate(score = round(score, 2), control_score = round(control_score, 2), excess_over_functional = round(excess_over_functional, 2)) |> dplyr::select(-evidence)), row.names = FALSE)
  cat("\n"); for (k in seq_len(nrow(Diag))) cat(sprintf("%s [%s]: %s\n", Diag$pan_id[k], Diag$diagnosis_type[k], Diag$summary[k]))


  # ---- PHASE 6: Tables and the desilting baseline -------------------------------------------------------------------------------
  log_phase("PHASE 6: outputs and the desilting baseline")
  Long <- dplyr::left_join(Seasonal, Static, by = "pan_id") |> dplyr::left_join(Pans.All[, c("pan_id", "name", "status", "type")], by = "pan_id") |>
    dplyr::relocate(pan_id, name, status, type, season)
  readr::write_csv(Static, file.path(Dir.Out.Tables, "pan_catchment_static.csv"))
  readr::write_csv(Long, file.path(Dir.Out.Tables, "pan_catchment_metrics.csv"))
  write_xlsx_min(Long, file.path(Dir.Out.Tables, "pan_catchment_metrics.xlsx"), sheet = "pan_catchment_metrics")

  test.pan <- Sys.getenv("AMAKA_TEST_DESILTING_PAN", "")                      # code-path test only: writes to tempdir, never to outputs/
  desilt.id <- if (nzchar(test.pan)) test.pan else Pans.All$pan_id[Pans.All$status == "existing_desilting"]
  mark <- file.path(Dir.Out.Tables, "NO_DESILTING_COORDINATES.txt")
  if (!length(desilt.id)) {
    message("  no existing_desilting pan in the file: no baseline to write")
  } else if (!all(desilt.id %in% Static$pan_id)) {
    writeLines(c("NO_DESILTING_COORDINATES: the desilting dam in water_pans.csv has no coordinates, so its pre-intervention baseline cannot be computed.",
                 "Add its longitude / latitude (GPS, or read from imagery) and re-run R/08_water_pans_analysis.R; the baseline will then be written to",
                 "outputs/tables/desilting_baseline_PRE_INTERVENTION_<date>.csv, to be re-run 6, 12 and 24 months after the works.",
                 sprintf("Pans affected: %s", paste(setdiff(desilt.id, Static$pan_id), collapse = ", "))), mark)
    message("  desilting dam has NO COORDINATES: wrote ", basename(mark), " instead of a baseline")
  } else {
    if (file.exists(mark)) file.remove(mark)
    base <- Long[Long$pan_id %in% desilt.id, ]; base$baseline_date <- format(Sys.Date()); base$note <- "PRE-INTERVENTION baseline: re-run 6, 12 and 24 months after the works"
    out <- if (nzchar(test.pan)) file.path(tempdir(), sprintf("TEST_desilting_baseline_%s.csv", format(Sys.Date()))) else file.path(Dir.Out.Tables, sprintf("desilting_baseline_PRE_INTERVENTION_%s.csv", format(Sys.Date())))
    readr::write_csv(base, out); message("  wrote ", out)
  }

  message(sprintf("  %d catchments, %d pan-season rows, %d shed rows", nrow(Static), nrow(Long), nrow(Sheds)))
  invisible(list(static = Static, seasonal = Seasonal, diag = Diag, hyp = Hyp, gap = Gap))
}

Pan.Result <- Run.Pan.Analysis()
log_phase("DONE")
