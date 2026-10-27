# ==============================================================================
# 11_pan_siting.R  -  Multi-criteria pan siting
# ==============================================================================
# PURPOSE
#   Suitability surface from Cfg.Suit.Weights (DESCRIPTION.md 6.9), then the top Cfg.Siting.N candidate sites by non-maximum
#   suppression at >= Cfg.Siting.MinSpace.m spacing, each with its catchment, yield, water-gap gain, comfort gain, woody cost,
#   distances and a one-line rationale. Existing non-functional pans get the same gap-gain numbers, so "build new" can be weighed
#   against "restore what is there".
#
# INPUTS
#   recommendation_zones_bitmask.tif (10), dist_to_functional_pan_m.tif (08), TWI / slope / flow accumulation / D8 pointer / streams (03),
#   WorldCover and Hansen (gee/08), CHIRPS (gee/03), Settled Areas, the Amaka AOI, water_pans.csv
#
# OUTPUTS
#   outputs/rasters/pan_suitability.tif (the weighted sum as defined in the brief), pan_suitability_smoothed.tif (eligible pixels, smoothed)
#   outputs/tables/candidate_pan_sites.csv, rehabilitation_gap_gain.csv, siting_sensitivity.csv
#   outputs/figures/pan_siting_map.png
#   (optional) candidate_new rows appended to data/water_pans/water_pans.csv: ONLY when Cfg.Siting.Append is TRUE, i.e. after Isaac confirms
#
# RUNTIME NOTES
#   * About a minute (the weight sweep re-picks sites 16 times). Everything is on the 30 m DEM grid. The weights are arguments of
#     Run.Pan.Siting() (defaults Cfg.Suit.Weights = the brief's v1 values).
#   * Terms (all on [0, 1]):  water_gap = mean over Cfg.Siting.Gap.Radii.km of "farther than r from every FUNCTIONAL pan" (the brief's binary
#     is_in_water_gap at four radii averaged, so 12 km from the only working pan scores 0.75, not 1; give one radius for the literal binary);
#     terrain = geometric mean of slope (5 x 5 smoothed) and flow-accumulation trapezoids (Cfg.Terrain.*) and scaled TWI, so a site must
#     sit on a drainage line of dam-sized catchment; zone_d_dry = share of the 5 km neighbourhood (inside Amaka) that is large-herbivore
#     comfort habitat (cover and TWI part of Zone D) - read literally, "inside Zone D-dry" would contradict the water-gap term, because D-dry
#     is by definition near EXISTING water; refugia = mean of distance-to-boundary (3 km) and distance-to-settled (15 km) saturations;
#     cover = shrub + tree fraction (already on [0, 1], not stretched); zone_e = 1 inside Zone E; woody_risk = mature-tree share within
#     Cfg.Woody.Shed.km, scaled to the ranch; settled = 1 - distance-to-settled / 15 km.
#   * Hard rules on top of the weighted sum: a candidate must be inside Amaka and at least Cfg.Siting.Edge.m inside its boundary (the polygon is
#     coarse), outside Zone E (the brief: Zone E is an EXCLUSION criterion) and have terrain feasibility >= Cfg.Siting.Min.Terrain. Candidates are
#     picked on the suitability smoothed over Cfg.Siting.Smooth.m. These rules leave only about 1-2% of the ranch buildable (drainage lines), which
#     is why fewer than N sites fit at the full spacing: the script then repeats at half the spacing as a tagged alternative.
#   * stability_pct = share of 16 runs with each weight halved or raised by half in which a site within half the spacing is also picked.
#     It is high mostly BECAUSE so little ground is buildable (the weights choose among the same few channel reaches), so it says which sites
#     survive a change of weights, not that they are right.
#   * Yield = catchment area x CHIRPS annual rainfall x Cfg.Runoff.Coef (an ASSUMPTION, no gauge data). Whether a pan holds water depends on
#     soils and wall design, which EO cannot see: every site needs a field visit. The two failed pans in the west are in the same water gap, so
#     rehabilitation_gap_gain.csv should be read next to the candidates.
# ==============================================================================

if (!exists("Setup.Loaded")) source("R/00_setup.R")
whitebox::wbt_verbose(FALSE)

Run.Pan.Siting <- function(weights = Cfg.Suit.Weights, n.sites = Cfg.Siting.N, min.space.m = Cfg.Siting.MinSpace.m, tag = "") {
  # `tag` is appended to the table and figure names of an alternative run (e.g. relaxed spacing); the rasters, the rehabilitation table and
  # water_pans.csv are only written by the untagged run.

  # ---- PHASE 1: Layers on the 30 m grid ----------------------------------------------------------------
  log_phase("PHASE 1: layers on the 30 m grid")
  nrm <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
  Dem <- terra::rast(file.path(Dir.DEM, "amaka_dem_utm37s.tif"))
  Aoi.Sf <- read_aoi(); Aoi.U <- sf::st_union(sf::st_geometry(Aoi.Sf)); Aoi.Mask <- terra::rasterize(terra::vect(Aoi.Sf), Dem, field = 1)
  Px.ha <- prod(terra::res(Dem)) / 1e4
  rr <- function(...) terra::rast(file.path(Dir.Out.Rasters, ...))
  Hb <- function(n) terra::rast(file.path(Dir.Out.Hydro.Buf, n))
  avg_to_dem <- function(r) terra::resample(r, Dem, method = "average")
  tf <- function(x) terra::ifel(is.na(x), 0, x)
  sat <- function(d, ref.m) terra::clamp(d / ref.m, 0, 1)
  wmean_focal <- function(x, w) {                                                  # weighted focal mean that ignores NA cells
    num <- terra::focal(terra::ifel(is.na(x), 0, x), w = w, fun = "sum", fillvalue = 0)
    den <- terra::focal(terra::ifel(is.na(x), 0, 1), w = w, fun = "sum", fillvalue = 0)
    terra::ifel(den > 0, num / den, NA)
  }
  nb_mean <- function(x, radius.m, fact = 5) {                                     # mean of x within radius.m, computed on a coarser grid
    xa <- terra::aggregate(x, fact = fact, fun = "mean", na.rm = TRUE)
    terra::resample(wmean_focal(xa, terra::focalMat(xa, radius.m, type = "circle")), x, method = "bilinear")
  }

  Slope <- Hb("slope_deg.tif"); Twi <- Hb("twi.tif"); Acc <- Hb("flow_accum_d8.tif")
  Zb <- rr("recommendation_zones_bitmask.tif"); bit <- function(k) (Zb %/% k) %% 2           # A=1 B=2 C=4 D-wet=8 D-dry=16 E=32
  Zone.E <- bit(32); Zone.Dwet <- bit(8)
  D.Func <- rr("dist_to_functional_pan_m.tif")
  Wc <- terra::rast(gee_file("worldcover_2021_envelope")); Woody <- avg_to_dem(terra::ifel(Wc == 10 | Wc == 20, 1, 0))
  Hansen <- terra::rast(gee_file("hansen_gfc_v1_11_envelope")); Mature <- avg_to_dem(terra::ifel(Hansen$treecover2000 >= Cfg.Tree.Mature.Pct, 1, 0))
  Settled.Sf <- sf::st_read(file.path(Dir.Boundary.Ctx, "Tsavo_Settled_Areas.shp"), quiet = TRUE)
  Settled.Near <- Settled.Sf[lengths(sf::st_intersects(Settled.Sf, sf::st_buffer(Aoi.U, 40000))) > 0, ]
  Dem.150 <- terra::aggregate(Dem, 5)
  Dist.Set <- if (nrow(Settled.Near)) terra::resample(terra::distance(Dem.150, terra::vect(Settled.Near)), Dem, method = "bilinear") else Dem * 0 + 1e5
  Dist.Bnd <- terra::distance(terra::ifel(is.na(Aoi.Mask), 1, NA))                              # inside the AOI: distance to the nearest pixel outside it
  message(sprintf("  AOI %s ha; nearest settled area %.1f km from the AOI; functional pans: %d",
                  format(round(sum(!is.na(terra::values(Aoi.Mask))) * Px.ha), big.mark = ","), as.numeric(sf::st_distance(Aoi.U, sf::st_union(sf::st_geometry(Settled.Sf)))) / 1000,
                  sum(read_pans()$functional)))


  # ---- PHASE 2: Terms and the suitability surface ------------------------------------------------------------
  log_phase("PHASE 2: suitability")
  Rad <- Cfg.Siting.Gap.Radii.km
  s.slope <- trapezoid(terra::focal(Slope, w = 5, fun = "mean", na.rm = TRUE), Cfg.Terrain.Slope.Deg)
  s.acc   <- trapezoid(Acc, Cfg.Terrain.Acc.Cells)
  s.twi   <- scale_ref(Twi, Aoi.Mask, min.range = 0.5, label = "TWI")
  Tm <- list(water_gap  = Reduce(`+`, lapply(Rad, function(r) terra::ifel(D.Func > r * 1000, 1, 0))) / length(Rad),
             terrain    = (s.slope * s.acc * s.twi)^(1 / 3),
             zone_d_dry = nb_mean(Zone.Dwet, Cfg.Zone.D.Dist.km * 1000),
             refugia    = (sat(Dist.Bnd, Cfg.Refugia.Ref.km[["boundary"]] * 1000) + sat(Dist.Set, Cfg.Refugia.Ref.km[["settled"]] * 1000)) / 2,
             cover      = Woody,
             zone_e     = Zone.E,
             woody_risk = scale_ref(nb_mean(Mature, Cfg.Woody.Shed.km * 1000), Aoi.Mask, label = "mature-tree share in the shed"),
             settled    = 1 - sat(Dist.Set, Cfg.Refugia.Ref.km[["settled"]] * 1000))
  for (k in names(Tm)) names(Tm[[k]]) <- k
  stopifnot(setequal(names(weights), names(Tm)))
  suit <- function(w) terra::mask(Reduce(`+`, lapply(names(w), function(k) w[[k]] * Tm[[k]])), Aoi.Mask)
  Elig <- tf(Tm$terrain >= Cfg.Siting.Min.Terrain) * tf(Zone.E == 0) * tf(Aoi.Mask) * tf(Dist.Bnd >= Cfg.Siting.Edge.m)   # buildable terrain, outside Zone E, inside Amaka and off its edge
  Smooth.W <- terra::focalMat(Dem, Cfg.Siting.Smooth.m, type = "circle")
  smooth_suit <- function(S) { Se <- terra::mask(S, terra::ifel(Elig == 1, 1, NA)); terra::mask(wmean_focal(Se, Smooth.W), Se) }
  pick_sites <- function(Ss, n, spacing) {                                                     # non-maximum suppression on the smoothed surface
    v <- terra::values(Ss)[, 1]; cell <- which(is.finite(v)); xy <- terra::xyFromCell(Ss, cell)
    df <- data.frame(cell = cell, x = xy[, 1], y = xy[, 2], score = v[cell]); out <- NULL
    while (nrow(df) && (is.null(out) || nrow(out) < n)) {
      i <- which.max(df$score); out <- rbind(out, df[i, ]); df <- df[sqrt((df$x - df$x[i])^2 + (df$y - df$y[i])^2) >= spacing, ]
    }
    out
  }
  S0 <- suit(weights); names(S0) <- "suitability"; S0s <- smooth_suit(S0); names(S0s) <- "suitability_smoothed"
  if (!nzchar(tag)) {
    terra::writeRaster(S0, file.path(Dir.Out.Rasters, "pan_suitability.tif"), overwrite = TRUE)
    terra::writeRaster(S0s, file.path(Dir.Out.Rasters, "pan_suitability_smoothed.tif"), overwrite = TRUE)
  }
  a <- terra::values(S0)[, 1]; a <- a[is.finite(a)]
  Elig.ha <- sum(terra::values(Elig)[, 1] == 1, na.rm = TRUE) * Px.ha
  message(sprintf("  suitability over the AOI: %.2f to %.2f (mean %.2f); %s ha eligible (%.1f%% of the AOI): terrain >= %.1f, outside Zone E, >= %g m inside the boundary",
                  min(a), max(a), mean(a), format(round(Elig.ha), big.mark = ","), 100 * Elig.ha / (length(a) * Px.ha), Cfg.Siting.Min.Terrain, Cfg.Siting.Edge.m))
  Sites <- pick_sites(S0s, n.sites, min.space.m)
  message(sprintf("  %d of %d candidate sites found at >= %.1f km spacing", nrow(Sites), n.sites, min.space.m / 1000))
  if (!nrow(Sites)) stop("No eligible pixel: relax Cfg.Siting.Min.Terrain or check the zone and terrain layers.", call. = FALSE)


  # ---- PHASE 3: Sensitivity to the weights ------------------------------------------------------------------------
  log_phase("PHASE 3: weight sensitivity (16 runs)")
  Runs <- list()
  for (k in names(weights)) for (f in c(0.5, 1.5)) {
    w <- weights; w[[k]] <- w[[k]] * f
    Runs[[sprintf("%s x%.1f", k, f)]] <- pick_sites(smooth_suit(suit(w)), n.sites, min.space.m)
  }
  stab <- vapply(seq_len(nrow(Sites)), function(i) 100 * mean(vapply(Runs, function(d) any(sqrt((d$x - Sites$x[i])^2 + (d$y - Sites$y[i])^2) <= min.space.m / 2), TRUE)), 0)
  Sens <- dplyr::bind_rows(lapply(names(Runs), function(n) data.frame(run = n, rank = seq_len(nrow(Runs[[n]])), x_utm = Runs[[n]]$x, y_utm = Runs[[n]]$y, smoothed_score = Runs[[n]]$score)))
  readr::write_csv(Sens, file.path(Dir.Out.Tables, sprintf("siting_sensitivity%s.csv", tag)))
  message("  stability of the baseline sites (% of 16 runs that pick a site within ", min.space.m / 2000, " km): ", paste(sprintf("C%02d %.0f%%", seq_along(stab), stab), collapse = ", "))


  # ---- PHASE 4: Attributes of each candidate ------------------------------------------------------------------------
  log_phase("PHASE 4: candidate attributes")
  Dir.Pan <- file.path(Dir.Out.Hydro.Buf, "pans"); dir.create(Dir.Pan, showWarnings = FALSE)
  Aoi.Px <- which(!is.na(terra::values(Aoi.Mask)[, 1])); Aoi.Xy <- terra::xyFromCell(Aoi.Mask, Aoi.Px)
  Aoi.Df <- terra::values(D.Func)[Aoi.Px, 1]; Aoi.Dw <- terra::values(Zone.Dwet)[Aoi.Px, 1] == 1
  Rain.Ann <- terra::app(gee_stack("chirps_precip_sum", "ANNUAL"), "mean", na.rm = TRUE)
  H.Layers <- c(Hansen$treecover2000 >= Cfg.Tree.Mature.Pct, Hansen$loss, Hansen$loss * Hansen$treecover2000 / 100); names(H.Layers) <- c("mature_tree", "loss", "loss_canopy")
  Pans <- sf::st_transform(sf::st_as_sf(dplyr::filter(read_pans(), located), coords = c("start_lon", "start_lat"), crs = Cfg.CRS.Geo, remove = FALSE), Cfg.CRS.Metric)

  ext1 <- function(r, x, y) { v <- terra::extract(r, cbind(x, y)); v[[ncol(v)]] }                  # value of a single-layer raster at one point
  gap_ha <- function(x, y, r.km) sum(sqrt((Aoi.Xy[, 1] - x)^2 + (Aoi.Xy[, 2] - y)^2) <= r.km * 1000 & Aoi.Df > r.km * 1000) * Px.ha        # AOI inside r of the site and outside every functional shed
  comfort <- function(x, y) { near <- sqrt((Aoi.Xy[, 1] - x)^2 + (Aoi.Xy[, 2] - y)^2) <= Cfg.Zone.D.Dist.km * 1000
                              c(dwet_ha = sum(near & Aoi.Dw) * Px.ha, new_ddry_ha = sum(near & Aoi.Dw & Aoi.Df > Cfg.Zone.D.Dist.km * 1000) * Px.ha) }
  woody <- function(x, y) { shed <- sf::st_transform(sf::st_sf(geometry = sf::st_buffer(sf::st_sfc(sf::st_point(c(x, y)), crs = Cfg.CRS.Metric), Cfg.Woody.Shed.km * 1000)), terra::crs(Hansen))
                            wh <- exactextractr::exact_extract(H.Layers, shed, fun = "sum", progress = FALSE) * prod(terra::res(Hansen)) / 1e4
                            c(mature_ha = wh[[1]], loss_ha = wh[[2]]) }
  catchment <- function(i, x, y) {
    f.pt <- file.path(Dir.Pan, sprintf("cand%02d_start.shp", i)); f.sn <- file.path(Dir.Pan, sprintf("cand%02d_snapped.shp", i)); f.ws <- file.path(Dir.Pan, sprintf("cand%02d_catchment.tif", i))
    write_shp(sf::st_sf(id = i, geometry = sf::st_sfc(sf::st_point(c(x, y)), crs = Cfg.CRS.Metric)), f.pt)
    whitebox::wbt_jenson_snap_pour_points(pour_pts = nrm(f.pt), streams = nrm(file.path(Dir.Out.Hydro.Buf, "streams.tif")), output = nrm(f.sn), snap_dist = 90)
    whitebox::wbt_watershed(d8_pntr = nrm(file.path(Dir.Out.Hydro.Buf, "d8_pointer.tif")), pour_pts = nrm(f.sn), output = nrm(f.ws))
    ws <- terra::rast(f.ws); if (is.na(terra::crs(ws, describe = TRUE)$code)) terra::crs(ws) <- paste0("EPSG:", Cfg.CRS.Metric)
    ws <- terra::ifel(ws > 0, 1, NA); pg <- sf::st_union(sf::st_as_sf(terra::as.polygons(ws, dissolve = TRUE))); sf::st_crs(pg) <- Cfg.CRS.Metric
    ha <- sum(!is.na(terra::values(ws))) * Px.ha; rain <- zonal_stat(Rain.Ann, sf::st_sf(geometry = pg))[[1]]
    up <- Pans$pan_id[lengths(sf::st_intersects(Pans, pg)) > 0]                    # existing pans inside this catchment = upslope of the site
    list(vals = c(catch_ha = ha, rain_mm = rain, yield_m3 = ha * 1e4 * rain / 1000 * Cfg.Runoff.Coef), upstream = paste(up, collapse = ";"))
  }

  Pt.Ll <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(Sites, coords = c("x", "y"), crs = Cfg.CRS.Metric), Cfg.CRS.Geo))
  Tm.At <- terra::extract(c(Tm$water_gap, Tm$terrain, Tm$zone_d_dry, Tm$refugia, Tm$cover, Tm$zone_e, Tm$woody_risk, Tm$settled), cbind(Sites$x, Sites$y)); Tm.At <- Tm.At[, setdiff(names(Tm.At), "ID"), drop = FALSE]
  Cand <- dplyr::bind_rows(lapply(seq_len(nrow(Sites)), function(i) {
    x <- Sites$x[i]; y <- Sites$y[i]; p <- sf::st_sfc(sf::st_point(c(x, y)), crs = Cfg.CRS.Metric)
    cat.i <- catchment(i, x, y); cm <- cat.i$vals; cf <- comfort(x, y); wd <- woody(x, y); np <- as.numeric(sf::st_distance(p, Pans)); j <- which.min(np)
    g5 <- gap_ha(x, y, 5); g10 <- gap_ha(x, y, 10)
    d.set <- as.numeric(sf::st_distance(p, sf::st_union(sf::st_geometry(Settled.Sf)))); d.bnd <- as.numeric(sf::st_distance(p, sf::st_boundary(Aoi.U)))
    tibble::tibble(rank = i, site_id = sprintf("AMK_C%02d", i), longitude = Pt.Ll[i, 1], latitude = Pt.Ll[i, 2], x_utm37s = x, y_utm37s = y,
                   suitability = ext1(S0, x, y), suitability_smoothed = Sites$score[i], stability_pct = stab[i],
                   terrain_feasibility = Tm.At$terrain[i], term_water_gap = Tm.At$water_gap[i], term_zone_d_dry = Tm.At$zone_d_dry[i], term_refugia = Tm.At$refugia[i], term_cover = Tm.At$cover[i],
                   term_woody_risk = Tm.At$woody_risk[i], term_settled = Tm.At$settled[i], in_zone_E = Tm.At$zone_e[i] == 1,
                   eligible_area_ha = Elig.ha, spacing_km = min.space.m / 1000,
                   catchment_ha = cm[["catch_ha"]], annual_rainfall_mm = cm[["rain_mm"]], annual_yield_m3 = cm[["yield_m3"]], runoff_coefficient_assumed = Cfg.Runoff.Coef,
                   existing_pans_upslope = cat.i$upstream,
                   newly_covered_gap_ha_5km = g5, newly_covered_gap_ha_10km = g10, comfort_habitat_ha_5km = cf[["dwet_ha"]], new_dry_comfort_ha_5km = cf[["new_ddry_ha"]],
                   mature_tree_ha_in_5km_shed = wd[["mature_ha"]], hansen_loss_ha_in_5km_shed = wd[["loss_ha"]], expected_loss_ha_per_decade = wd[["loss_ha"]] / 23 * 10,
                   dist_settled_km = d.set / 1000, dist_boundary_km = d.bnd / 1000,
                   nearest_pan_id = Pans$pan_id[j], nearest_pan_status = Pans$status[j], nearest_pan_km = np[j] / 1000)
  }))
  Cand$rationale <- sprintf("Brings water to %s ha of the 5 km gap (%s ha at 10 km) and %s ha of large-herbivore comfort habitat newly within 5 km of water; catchment %s ha (about %s m3/yr at an assumed runoff coefficient of %.2f); %.1f km from settled areas, %.1f km inside the boundary; %s ha of mature trees in its 5 km shed (%.1f ha lost since 2001); outside Zone E; %s%s",
                            format(round(Cand$newly_covered_gap_ha_5km), big.mark = ","), format(round(Cand$newly_covered_gap_ha_10km), big.mark = ","), format(round(Cand$new_dry_comfort_ha_5km), big.mark = ","),
                            format(round(Cand$catchment_ha), big.mark = ","), format(round(Cand$annual_yield_m3, -3), big.mark = ","), Cfg.Runoff.Coef, Cand$dist_settled_km, Cand$dist_boundary_km,
                            format(round(Cand$mature_tree_ha_in_5km_shed), big.mark = ","), Cand$hansen_loss_ha_in_5km_shed,
                            sprintf("picked in %.0f%% of 16 weight variations (only %s ha is buildable, so the weights mostly choose among the same few channels)", Cand$stability_pct, format(round(Cand$eligible_area_ha), big.mark = ",")),
                            paste0(ifelse(Cand$nearest_pan_km < 1.5, sprintf("; only %.1f km from existing pan %s (%s): compare with restoring it", Cand$nearest_pan_km, Cand$nearest_pan_id, sub("existing_", "", Cand$nearest_pan_status)), ""),
                                   ifelse(nzchar(Cand$existing_pans_upslope), sprintf("; its catchment contains existing pan(s) %s, so water kept there is water this site does not get", Cand$existing_pans_upslope), "")))
  readr::write_csv(Cand, file.path(Dir.Out.Tables, sprintf("candidate_pan_sites%s.csv", tag)))
  print(as.data.frame(Cand |> dplyr::select(rank, site_id, longitude, latitude, suitability, stability_pct, terrain_feasibility, catchment_ha, annual_yield_m3, newly_covered_gap_ha_5km, newly_covered_gap_ha_10km,
                                            new_dry_comfort_ha_5km, mature_tree_ha_in_5km_shed, dist_settled_km, nearest_pan_id, nearest_pan_km) |> dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4)))), row.names = FALSE)

  # the same gap-gain numbers for the existing pans that do not work (rehabilitation option)
  Rehab <- Pans[!Pans$functional, ]
  Rehab.Tab <- dplyr::bind_rows(lapply(seq_len(nrow(Rehab)), function(i) {
    xy <- sf::st_coordinates(Rehab[i, ]); cf <- comfort(xy[1], xy[2]); wd <- woody(xy[1], xy[2])
    tibble::tibble(pan_id = Rehab$pan_id[i], name = Rehab$name[i], status = Rehab$status[i], inside_aoi = lengths(sf::st_intersects(Rehab[i, ], Aoi.U)) > 0,
                   gap_ha_covered_5km = gap_ha(xy[1], xy[2], 5), gap_ha_covered_10km = gap_ha(xy[1], xy[2], 10), comfort_habitat_ha_5km = cf[["dwet_ha"]], new_dry_comfort_ha_5km = cf[["new_ddry_ha"]],
                   mature_tree_ha_in_5km_shed = wd[["mature_ha"]], hansen_loss_ha_in_5km_shed = wd[["loss_ha"]])
  }))
  if (!nzchar(tag)) {
    readr::write_csv(Rehab.Tab, file.path(Dir.Out.Tables, "rehabilitation_gap_gain.csv"))
    print(as.data.frame(Rehab.Tab |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)
  }


  # ---- PHASE 5: Map ----------------------------------------------------------------------------------------------------------
  log_phase("PHASE 5: siting map")
  Bb <- sf::st_bbox(c(sf::st_buffer(Aoi.U, 2500), sf::st_buffer(sf::st_geometry(Pans), 1500)))                  # the AOI and every pan with a margin
  Hs <- terra::crop(Hb("hillshade.tif"), terra::ext(Bb[["xmin"]], Bb[["xmax"]], Bb[["ymin"]], Bb[["ymax"]])); names(Hs) <- "shade"
  Ez <- terra::ifel(Zone.E == 1, 1, NA); names(Ez) <- "zone_E"
  Cand.Sf <- sf::st_as_sf(Cand, coords = c("x_utm37s", "y_utm37s"), crs = Cfg.CRS.Metric, remove = FALSE); Cand.Sf$label <- as.character(Cand$rank)
  Pans.Plot <- Pans; Pans.Plot$state <- ifelse(Pans.Plot$functional, "existing, functional", "existing, failed or degraded")
  Map.Siting <-
    tmap::tm_shape(Hs) + tmap::tm_raster(col.scale = tmap::tm_scale_continuous(values = "brewer.greys"), col.legend = tmap::tm_legend_hide()) +
    tmap::tm_shape(terra::mask(S0s, terra::ifel(Elig == 1, 1, NA))) + tmap::tm_raster(col.scale = tmap::tm_scale_continuous(values = "scico.vik"), col.legend = tmap::tm_legend(title = "Suitability\n(eligible pixels)")) +
    tmap::tm_shape(Ez) + tmap::tm_raster(col.scale = tmap::tm_scale_categorical(values = "#6a3d9a", labels = "Zone E (excluded)"), col.legend = tmap::tm_legend(title = ""), col_alpha = 0.3) +
    tmap::tm_shape(sf::st_geometry(Aoi.Sf)) + tmap::tm_borders(col = "black", lwd = 2.5) +
    tmap::tm_shape(Pans.Plot) + tmap::tm_symbols(fill = "state", fill.scale = tmap::tm_scale_categorical(values = c("existing, functional" = "#1f78b4", "existing, failed or degraded" = "#e31a1c")),
                                                 fill.legend = tmap::tm_legend(title = "Existing pans"), shape = 21, size = 0.9, col = "white") +
    tmap::tm_shape(Cand.Sf) + tmap::tm_symbols(fill = "#ffd92f", col = "black", shape = 23, size = 1.2) + tmap::tm_text("label", size = 0.9, fontface = "bold") +
    tmap::tm_title(sprintf("Candidate pan sites (yellow diamonds, ranked): %d at >= %.1f km spacing, weights v1\nSuitability shown only where a pan is buildable (%s ha: drainage lines of dam-sized catchment, outside Zone E)",
                           nrow(Cand), min.space.m / 1000, format(round(Elig.ha), big.mark = ",")), size = 0.7)
  tmap::tmap_save(Map.Siting, file.path(Dir.Out.Figures, sprintf("pan_siting_map%s.png", tag)), width = 8, height = 6.4, dpi = 200)
  message("  wrote ", sprintf("pan_siting_map%s.png", tag))


  # ---- PHASE 6: Optional append to water_pans.csv -----------------------------------------------------------------------------
  if (nzchar(tag)) {
    message("  alternative run '", tag, "': water_pans.csv not touched")
  } else if (isTRUE(Cfg.Siting.Append)) {
    Wp <- readr::read_csv(Path.Pans, show_col_types = FALSE, col_types = readr::cols(.default = readr::col_character()))
    Wp <- Wp[Wp$status != "candidate_new" | is.na(Wp$status), ]
    New <- tibble::tibble(pan_id = Cand$site_id, name = sprintf("Candidate site %d", Cand$rank), longitude = as.character(round(Cand$longitude, 6)), latitude = as.character(round(Cand$latitude, 6)),
                          type = "earth_dam", status = "candidate_new", season_functional = NA_character_, notes = sprintf("R/11 %s: rank %d, suitability %.2f, stability %.0f%%", format(Sys.Date()), Cand$rank, Cand$suitability, Cand$stability_pct))
    readr::write_csv(dplyr::bind_rows(Wp, New[, intersect(names(New), names(Wp))]), Path.Pans, na = "")
    message("  appended ", nrow(New), " candidate_new rows to water_pans.csv (Cfg.Siting.Append = TRUE)")
  } else message("  water_pans.csv NOT touched: set Cfg.Siting.Append <- TRUE only after Isaac confirms the candidates")

  invisible(list(candidates = Cand, sensitivity = Sens, rehabilitation = Rehab.Tab))
}

Siting.Result <- Run.Pan.Siting()
if (nrow(Siting.Result$candidates) < Cfg.Siting.N) {          # the ranch is only ~6.5 x 7 km: fewer than N sites fit at the full spacing (DESCRIPTION.md 10: "adjust if ranch geometry is tight")
  message("\nOnly ", nrow(Siting.Result$candidates), " of ", Cfg.Siting.N, " sites fit at ", Cfg.Siting.MinSpace.m / 1000, " km: repeating at half the spacing as a clearly tagged alternative")
  Siting.Alt <- Run.Pan.Siting(min.space.m = Cfg.Siting.MinSpace.m / 2, tag = sprintf("_spacing_%gm", Cfg.Siting.MinSpace.m / 2))
}
log_phase("DONE")
