# 10_integrated_stress.R  -  Composite stress index and recommendation zones A-E
# PURPOSE
#   Min-max scaled stress index for the ANNUAL, DRY and LR seasons (weights from Cfg.Stress.Weights) and the recommendation
#   zones of DESCRIPTION.md 6.8: A regeneration feasibility, B catchment interventions, C pan protection, D-wet / D-dry large-
#   herbivore comfort, E waterless retention buffer. Every pan also gets the woody-sacrifice counter-metric.
#
# INPUTS
#   RPI (04), PML seasonal (05), Sentinel-2 BSI (06), TWI / streams / flow accumulation (03), pan catchments, distances and
#   water gap (08), Hansen tree cover and loss (gee/08), WorldCover (gee/08)
#
# OUTPUTS
#   outputs/rasters/stress_<ANNUAL|DRY|LR>_v1_030_020_020_015_015.tif (+ plain stress_annual.tif, stress_DRY.tif, stress_LR.tif)
#   outputs/rasters/stress_amplification_DRY_minus_LR.tif
#   outputs/rasters/recommendation_zones.tif (categorical, priority order B > A > C > E > D-dry > D-wet) and
#       recommendation_zones_bitmask.tif (A=1 B=2 C=4 D-wet=8 D-dry=16 E=32, zones overlap)
#   outputs/tables/stress_summary.csv, stress_terms.csv, stress_scaling.csv, cells_stress.csv, recommendation_zone_summary.csv,
#       pan_stress_zoneC.csv (stress, Zone C test and woody-sacrifice counter-metric per pan)
#   outputs/figures/10_stress_check.png, 10_zones_check.png
#
# RUNTIME NOTES
#   * About 30 s. Everything is computed on the 30 m DEM grid. Coarse terms (RPI 460 m, PML 500 m) are REPLICATED onto it with
#     nearest-neighbour (no interpolation), fine terms (TWI 30 m, Sentinel-2 BSI 10 m) keep their detail, so within a coarse cell
#     the 30 m variation comes from TWI and BSI only. Cell tables (cells_stress.csv) are area-weighted means over the RPI-grid cells.
#   * Each term is min-max scaled between the 2nd and 98th percentile of the AOI (Cfg.Stress.Clip) and clamped, so the index
#     is RELATIVE to Amaka: it ranks places within the ranch and says nothing about absolute stress. Terms and their limits are in
#     stress_scaling.csv. A flat layer is not stretched (minimum range guard).
#   * bare_fraction = Sentinel-2 BSI (season-specific): WorldCover has no bare class on Amaka (known_issues D5). Es_frac is
#     season-specific too. Es/ET is highest in the RAINY seasons (wet soil evaporates), so the DRY Es term reads as the bare-soil
#     share of the little water that moves, not as dryness (known_issues). is_degrading_RPI is 1 in about 1% of cells.
#   * "High" / "low" cut-offs are quantiles of the AOI (Cfg.Relative.Q). Zone A = top-quartile TWI and top-quartile dry-season BSI and
#     not degrading. Zone B = Strahler order 1-2 reaches with top-quartile MEAN UPSLOPE stress (D8 mass flux of the DRY stress,
#     divided by flow accumulation) within Cfg.Zone.B.Pan.km of a pan or inside a pan catchment. Zone C = pan catchment whose mean
#     DRY stress is in the top quartile of catchment-sized windows of the ranch AND whose RPI Sen slope is negative.
#     Zone D-wet = shrub + tree cover > Cfg.Zone.Cover.Min and TWI above the ranch median; D-dry = D-wet within Cfg.Zone.D.Dist.km
#     of a FUNCTIONAL pan. Zone E = farther than Cfg.Zone.E.Dist.km from every functional pan, with that woody cover and RPI_now at or
#     above the Cfg.Zone.E.Rpi.Q quantile.
#   * The index is structural only (DESCRIPTION.md 9.9: non-structural pathways are >2x stronger in McCleery 2025) and daytime-biased
#     (9.8). Zone C never reads as "exclusion is good": the trade-off note travels with every row (McCleery 2025; Jones 2025).

if (!exists("Setup.Loaded")) source("R/00_setup.R")
whitebox::wbt_verbose(FALSE)

# PHASE 1: Layers on the 30 m grid
log_phase("PHASE 1: layers on the 30 m grid")
nrm <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
Dem <- terra::rast(file.path(Dir.DEM, "amaka_dem_utm37s.tif"))
Aoi.Sf <- read_aoi(); Aoi.Vec <- terra::vect(Aoi.Sf); Aoi.Mask <- terra::rasterize(Aoi.Vec, Dem, field = 1)
Px.ha <- prod(terra::res(Dem)) / 1e4
rep_to_dem <- function(r) terra::project(r, Dem, method = "near")                # coarse layers: replicated, not interpolated
avg_to_dem <- function(r) terra::resample(r, Dem, method = "average")           # fine layers: mean of the fine pixels
rr <- function(...) terra::rast(file.path(Dir.Out.Rasters, ...))
Hb <- function(n) terra::rast(file.path(Dir.Out.Hydro.Buf, n))
q_aoi <- function(x, p) global_q(terra::mask(x, Aoi.Mask), p)
tf <- function(x) terra::ifel(is.na(x), 0, x)                                    # logical / NA raster -> 0 / 1

Twi <- Hb("twi.tif"); Streams <- Hb("streams.tif"); Order <- Hb("strahler.tif"); Acc <- Hb("flow_accum_d8.tif")
Rpi.Now <- rep_to_dem(rr("rpi_now.tif")); Rpi.Sen <- rep_to_dem(rr("rpi_sen_slope.tif")); Rpi.Deg <- tf(rep_to_dem(rr("rpi_degrading.tif")))
Es  <- setNames(lapply(c("ANNUAL", "DRY", "LR"), function(s) rep_to_dem(rr(sprintf("pml_esfrac_%s.tif", s)))), c("ANNUAL", "DRY", "LR"))
Bsi <- setNames(lapply(c("ANNUAL", "DRY", "LR"), function(s) avg_to_dem(rr("s2", sprintf("s2_bsi_%s_mean.tif", s)))), c("ANNUAL", "DRY", "LR"))
Wc <- terra::rast(gee_file("worldcover_2021_envelope"))
Woody <- avg_to_dem(terra::ifel(Wc == 10 | Wc == 20, 1, 0)); names(Woody) <- "woody_cover"      # shrub + tree fraction
Hansen <- terra::rast(gee_file("hansen_gfc_v1_11_envelope"))
D.Func <- rr("dist_to_functional_pan_m.tif"); D.Any <- rr("dist_to_any_pan_m.tif")
In.Catch <- terra::app(rr("pan_catchments.tif"), "max", na.rm = TRUE)                              # 1 inside any stream-fed pan catchment
message(sprintf("  grid %d x %d cells at %g m; AOI = %s ha (%s pixels)", terra::nrow(Dem), terra::ncol(Dem), terra::res(Dem)[1],
                format(round(sum(!is.na(terra::values(Aoi.Mask))) * Px.ha), big.mark = ","), format(sum(!is.na(terra::values(Aoi.Mask))), big.mark = ",")))


# PHASE 2: Scaled terms and the three stress surfaces 
log_phase("PHASE 2: stress index")
Scale.Log <- tibble::tibble()
scale_term <- function(x, label, min.range = 0.01) {                             # scale_ref() against the AOI, logging the limits used
  lim <- scale_limits(x, Aoi.Mask)
  Scale.Log <<- dplyr::bind_rows(Scale.Log, tibble::tibble(term = label, lower = lim[1], upper = lim[2], span = diff(lim)))
  scale_ref(x, Aoi.Mask, min.range = min.range, label = label)
}
Wt <- Cfg.Stress.Weights; stopifnot(abs(sum(Wt) - 1) < 1e-9, identical(names(Wt), c("rpi", "degrading", "es_dry", "bare", "twi")))
Ver <- paste0("v1_", paste(sprintf("%03d", round(100 * Wt)), collapse = "_"))
T.Rpi <- 1 - scale_term(Rpi.Now, "RPI_now (1 - scaled)"); T.Twi <- 1 - scale_term(Twi, "TWI (1 - scaled)", min.range = 0.5)

Stress <- list(); Terms <- list()
for (s in c("ANNUAL", "DRY", "LR")) {
  t.es <- scale_term(Es[[s]], paste0("Es_frac_", s)); t.bare <- scale_term(Bsi[[s]], paste0("BSI_", s))
  Terms[[s]] <- c(Wt[["rpi"]] * T.Rpi, Wt[["degrading"]] * Rpi.Deg, Wt[["es_dry"]] * t.es, Wt[["bare"]] * t.bare, Wt[["twi"]] * T.Twi)
  names(Terms[[s]]) <- c("rpi_deficit", "degrading", "es_frac", "bare_bsi", "twi_deficit")
  Stress[[s]] <- sum(Terms[[s]]); names(Stress[[s]]) <- paste0("stress_", s)
  nm <- c(ANNUAL = "annual", DRY = "DRY", LR = "LR")[[s]]
  terra::writeRaster(Stress[[s]], file.path(Dir.Out.Rasters, sprintf("stress_%s_%s.tif", nm, Ver)), overwrite = TRUE)
  terra::writeRaster(Stress[[s]], file.path(Dir.Out.Rasters, sprintf("stress_%s.tif", nm)), overwrite = TRUE)
}
Amp <- Stress$DRY - Stress$LR; names(Amp) <- "stress_amplification"
terra::writeRaster(Amp, file.path(Dir.Out.Rasters, "stress_amplification_DRY_minus_LR.tif"), overwrite = TRUE)
Scale.Log <- dplyr::distinct(Scale.Log); readr::write_csv(Scale.Log, file.path(Dir.Out.Tables, "stress_scaling.csv"))

av <- function(r) { v <- terra::values(terra::mask(r, Aoi.Mask))[, 1]; v[is.finite(v)] }
Sum.Tab <- dplyr::bind_rows(lapply(names(Stress), function(s) { v <- av(Stress[[s]])
  tibble::tibble(season = s, weights_version = Ver, aoi_mean = mean(v), aoi_sd = sd(v), p10 = quantile(v, .1), p50 = median(v), p90 = quantile(v, .9), min = min(v), max = max(v),
                 cor_with_LR = cor(v, av(Stress$LR))) }))
Term.Tab <- dplyr::bind_rows(lapply(names(Terms), function(s) { m <- terra::mask(Terms[[s]], Aoi.Mask)
  tibble::tibble(season = s, term = names(m), weighted_mean = unname(terra::global(m, "mean", na.rm = TRUE)[[1]]), weighted_sd = unname(terra::global(m, "sd", na.rm = TRUE)[[1]])) })) |>
  dplyr::group_by(season) |> dplyr::mutate(share_of_mean_pct = 100 * weighted_mean / sum(weighted_mean), share_of_sd_pct = 100 * weighted_sd / sum(weighted_sd)) |> dplyr::ungroup()
readr::write_csv(Sum.Tab, file.path(Dir.Out.Tables, "stress_summary.csv")); readr::write_csv(Term.Tab, file.path(Dir.Out.Tables, "stress_terms.csv"))
print(as.data.frame(Sum.Tab |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
print(as.data.frame(Term.Tab |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
message(sprintf("  DRY vs LR stress, correlation over the AOI: %.2f; mean amplification (DRY - LR) %+.3f", cor(av(Stress$DRY), av(Stress$LR)), mean(av(Amp))))


# PHASE 3: Cell table (RPI grid) 
log_phase("PHASE 3: cell table")
Grid <- get_analysis_grid(); Cells <- Grid$cells; Cell.Poly <- cell_polygons(Grid$grid)
Cs <- zonal_stat(c(Stress$ANNUAL, Stress$DRY, Stress$LR, Amp), Cell.Poly, "mean")
Cells.Stress <- dplyr::bind_cols(tibble::tibble(cell_id = Cell.Poly$cell_id), Cs)
readr::write_csv(Cells.Stress, file.path(Dir.Out.Tables, "cells_stress.csv"))
message(sprintf("  %d cells; AOI cells' DRY stress %.3f-%.3f", nrow(Cells.Stress), min(Cells.Stress$stress_DRY[Cells$in_aoi], na.rm = TRUE), max(Cells.Stress$stress_DRY[Cells$in_aoi], na.rm = TRUE)))


# PHASE 4: Recommendation zones 
log_phase("PHASE 4: recommendation zones")
Thr <- c(twi_hi = q_aoi(Twi, Cfg.Relative.Q), twi_med = q_aoi(Twi, 0.5), bsi_dry_hi = q_aoi(Bsi$DRY, Cfg.Relative.Q), rpi_lo = q_aoi(Rpi.Now, Cfg.Zone.E.Rpi.Q))
message("  relative thresholds (AOI): ", paste(sprintf("%s = %.3f", names(Thr), Thr), collapse = ", "))
Zone <- list()
Zone$A    <- tf(Twi >= Thr[["twi_hi"]]) * tf(Bsi$DRY >= Thr[["bsi_dry_hi"]]) * (1 - Rpi.Deg)
Zone$Dwet <- tf(Woody > Cfg.Zone.Cover.Min) * tf(Twi >= Thr[["twi_med"]])
Zone$Ddry <- Zone$Dwet * tf(D.Func <= Cfg.Zone.D.Dist.km * 1000)
Zone$E    <- tf(D.Func > Cfg.Zone.E.Dist.km * 1000) * tf(Woody > Cfg.Zone.Cover.Min) * tf(Rpi.Now >= Thr[["rpi_lo"]])

# Zone B: mean UPSLOPE stress of every stream cell = D8 mass flux of the DRY stress / flow accumulation (QC: loading = 1 reproduces accumulation exactly)
f.load <- file.path(Dir.Out.Hydro.Buf, "tmp_stress_loading.tif"); f.one <- file.path(Dir.Out.Hydro.Buf, "tmp_ones.tif")
f.zero <- file.path(Dir.Out.Hydro.Buf, "tmp_zeros.tif"); f.flux <- file.path(Dir.Out.Hydro.Buf, "tmp_stress_flux.tif")
terra::writeRaster(terra::ifel(is.na(Stress$DRY), global_q(Stress$DRY, 0.5), Stress$DRY), f.load, overwrite = TRUE)
terra::writeRaster(Dem * 0 + 1, f.one, overwrite = TRUE); terra::writeRaster(Dem * 0, f.zero, overwrite = TRUE)
whitebox::wbt_d8_mass_flux(dem = nrm(file.path(Dir.Out.Hydro.Buf, "dem_breached.tif")), loading = nrm(f.load), efficiency = nrm(f.one), absorption = nrm(f.zero), output = nrm(f.flux))
Up.Stress <- terra::rast(f.flux) / Acc; names(Up.Stress) <- "upslope_stress_DRY"
terra::writeRaster(Up.Stress, file.path(Dir.Out.Rasters, "upslope_stress_DRY.tif"), overwrite = TRUE)
unlink(c(f.load, f.one, f.zero, f.flux))
Low.Order <- tf(Streams > 0) * tf(Order >= min(Cfg.Zone.B.Orders)) * tf(Order <= max(Cfg.Zone.B.Orders)) * tf(Aoi.Mask)
Near.Pan  <- tf(D.Any <= Cfg.Zone.B.Pan.km * 1000) + tf(In.Catch) > 0
Thr.B <- global_q(terra::mask(Up.Stress, terra::ifel(Low.Order == 1, 1, NA)), Cfg.Relative.Q)
Zone$B <- Low.Order * tf(Up.Stress >= Thr.B) * Near.Pan
message(sprintf("  Zone B: %d stream cells of order %s in the AOI, %d with top-quartile upslope stress (>= %.3f) near a pan",
                sum(terra::values(Low.Order)[, 1] == 1, na.rm = TRUE), paste(range(Cfg.Zone.B.Orders), collapse = "-"), sum(terra::values(Zone$B)[, 1] == 1, na.rm = TRUE), Thr.B))

# Zone C: stress and RPI trend per pan catchment (stream-fed pans) judged against catchment-sized windows of the ranch
Pc <- sf::st_read(file.path(Dir.Out.Rasters, "pan_catchments.shp"), quiet = TRUE)
Pan.Static <- readr::read_csv(file.path(Dir.Out.Tables, "pan_catchment_static.csv"), show_col_types = FALSE)
Pans.Sf <- sf::st_transform(sf::st_as_sf(dplyr::filter(read_pans(), located), coords = c("start_lon", "start_lat"), crs = Cfg.CRS.Geo, remove = FALSE), Cfg.CRS.Metric)
Cell.Grid.ha <- stats::median(Cells$area_ha)
H.Layers <- c(Hansen$treecover2000 >= Cfg.Tree.Mature.Pct, Hansen$loss, Hansen$loss * Hansen$treecover2000 / 100)
names(H.Layers) <- c("mature_tree", "loss", "loss_canopy")                       # 30 m pixel counts inside a shed x pixel area = ha
window_means <- function(v, w) {                                                  # focal means of an AOI-cell metric over w x w cell windows
  r <- terra::rast(Grid$grid); full <- rep(NA_real_, terra::ncell(r)); ids <- Cells$cell_id[Cells$in_aoi]; full[ids] <- v[match(ids, Cells.Stress$cell_id)]
  terra::values(r) <- full
  n <- terra::focal(terra::ifel(is.na(r), 0, 1), w = w, fun = "sum")
  m <- terra::ifel(n >= ceiling(0.5 * w * w), terra::focal(r, w = w, fun = "mean", na.rm = TRUE), NA)
  x <- terra::values(m)[ids, 1]; x[is.finite(x)]
}
Pan.Tab <- dplyr::bind_rows(lapply(seq_len(nrow(Pans.Sf)), function(i) {
  id <- Pans.Sf$pan_id[i]; st <- Pan.Static[Pan.Static$pan_id == id, ]; poly <- Pc[Pc$pan_id == id, ]
  stream.fed <- isTRUE(st$stream_fed)
  m <- sapply(c("ANNUAL", "DRY", "LR"), function(s) exactextractr::exact_extract(Stress[[s]], sf::st_transform(poly, terra::crs(Dem)), fun = "mean", progress = FALSE))
  w <- max(3, 2 * floor(sqrt(st$catch_ha / Cell.Grid.ha) / 2) + 1)
  pctl <- if (stream.fed) 100 * mean(window_means(Cells.Stress$stress_DRY, w) <= m[["DRY"]]) else NA_real_
  shed <- sf::st_buffer(sf::st_geometry(Pans.Sf[i, ]), Cfg.Woody.Shed.km * 1000)
  wh <- exactextractr::exact_extract(H.Layers, sf::st_transform(sf::st_sf(geometry = shed), terra::crs(Hansen)), fun = "sum", progress = FALSE)
  px30 <- prod(terra::res(Hansen)) / 1e4
  tibble::tibble(pan_id = id, name = Pans.Sf$name[i], status = Pans.Sf$status[i], stream_fed = stream.fed, area_basis = poly$method,
                 stress_ANNUAL = m[["ANNUAL"]], stress_DRY = m[["DRY"]], stress_LR = m[["LR"]], stress_amplification = m[["DRY"]] - m[["LR"]],
                 dry_stress_percentile_vs_catchment_sized_windows = pctl, rpi_sen_slope = st$rpi_sen_slope,
                 zone_C = stream.fed && is.finite(pctl) && pctl >= 100 * Cfg.Relative.Q && st$rpi_sen_slope < 0,
                 zone_C_reason = dplyr::case_when(!stream.fed ~ "not stream-fed: no catchment", pctl < 100 * Cfg.Relative.Q & st$rpi_sen_slope >= 0 ~ "stress not in the top quartile and RPI not declining",
                                                  pctl < 100 * Cfg.Relative.Q ~ "stress not in the top quartile", st$rpi_sen_slope >= 0 ~ "RPI not declining", TRUE ~ "both criteria met"),
                 woody_shed_km = Cfg.Woody.Shed.km, mature_tree_ha_in_shed = wh[[1]] * px30, hansen_loss_ha_in_shed_2001_2023 = wh[[2]] * px30,
                 hansen_loss_canopy_ha_in_shed = wh[[3]] * px30, expected_loss_ha_per_decade = wh[[2]] * px30 / 23 * 10,
                 tradeoff_note = "Excluding large herbivores can release cryptic herbivores (small mammals) that suppress tree recruitment (McCleery 2025); recruitment is bottlenecked by cryptic herbivores, grass competition and abiotic factors whatever the large-herbivore pressure (Jones 2025). Exclusion is not unambiguously positive.")
}))
readr::write_csv(Pan.Tab, file.path(Dir.Out.Tables, "pan_stress_zoneC.csv"))
print(as.data.frame(Pan.Tab |> dplyr::select(pan_id, stream_fed, stress_DRY, stress_LR, stress_amplification, dry_stress_percentile_vs_catchment_sized_windows, rpi_sen_slope, zone_C, zone_C_reason,
                                             mature_tree_ha_in_shed, hansen_loss_ha_in_shed_2001_2023, expected_loss_ha_per_decade) |> dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
C.Ids <- Pan.Tab$pan_id[Pan.Tab$zone_C]
Zone$C <- if (length(C.Ids)) tf(terra::app(rr("pan_catchments.tif")[[C.Ids]], "max", na.rm = TRUE)) else Dem * 0
message(sprintf("  Zone C pan catchments: %s", ifelse(length(C.Ids), paste(C.Ids, collapse = ", "), "none meet both criteria")))

Zone <- lapply(Zone, function(z) terra::mask(z, Aoi.Mask))
Bits <- c(A = 1, B = 2, C = 4, Dwet = 8, Ddry = 16, E = 32)
Bitmask <- sum(Reduce(c, lapply(names(Bits), function(k) Zone[[k]] * Bits[[k]]))); names(Bitmask) <- "zone_bitmask"
terra::writeRaster(Bitmask, file.path(Dir.Out.Rasters, "recommendation_zones_bitmask.tif"), overwrite = TRUE, datatype = "INT1U")
Codes <- c(A = 1, B = 2, C = 3, Dwet = 4, Ddry = 5, E = 6); Labels <- c("A regeneration feasibility", "B catchment interventions", "C pan protection", "D-wet comfort", "D-dry comfort", "E waterless buffer")
Cat <- Aoi.Mask * 0
for (k in c("Dwet", "Ddry", "E", "C", "A", "B")) Cat <- terra::ifel(Zone[[k]] == 1, Codes[[k]], Cat)            # lowest priority first, so higher priorities overwrite
Cat <- terra::ifel(Cat == 0, NA, Cat); names(Cat) <- "zone"
levels(Cat) <- data.frame(id = unname(Codes), zone = Labels)
terra::writeRaster(Cat, file.path(Dir.Out.Rasters, "recommendation_zones.tif"), overwrite = TRUE, datatype = "INT1U")

Aoi.ha <- sum(!is.na(terra::values(Aoi.Mask))) * Px.ha
Zone.Tab <- dplyr::bind_rows(lapply(names(Codes), function(k) tibble::tibble(zone = Labels[[unname(Codes[k])]], overlapping_ha = sum(terra::values(Zone[[k]])[, 1] == 1, na.rm = TRUE) * Px.ha,
                                                                              exclusive_ha = sum(terra::values(Cat)[, 1] == Codes[[k]], na.rm = TRUE) * Px.ha))) |>
  dplyr::mutate(overlapping_pct_of_aoi = 100 * overlapping_ha / Aoi.ha, exclusive_pct_of_aoi = 100 * exclusive_ha / Aoi.ha)
readr::write_csv(Zone.Tab, file.path(Dir.Out.Tables, "recommendation_zone_summary.csv"))
print(as.data.frame(Zone.Tab |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)


# PHASE 5: Check figures 
log_phase("PHASE 5: check figures")
Hs <- Hb("hillshade.tif"); names(Hs) <- "shade"
Ps <- sf::st_geometry(Pans.Sf)
stress_map <- function(r, title, lim = c(0, 1), pal = "scico.lajolla") {
  tmap::tm_shape(Hs) + tmap::tm_raster(col.scale = tmap::tm_scale_continuous(values = "brewer.greys"), col.legend = tmap::tm_legend_hide()) +
    tmap::tm_shape(terra::mask(r, Aoi.Mask)) + tmap::tm_raster(col.scale = tmap::tm_scale_continuous(values = pal, limits = lim), col.legend = tmap::tm_legend(title = "")) +
    tmap::tm_shape(sf::st_geometry(Aoi.Sf)) + tmap::tm_borders(col = "black", lwd = 2) + tmap::tm_shape(Pans.Sf) + tmap::tm_symbols(fill = "cyan", col = "black", size = 0.5) +
    tmap::tm_title(title, size = 0.8)
}
Rng <- range(c(av(Stress$DRY), av(Stress$LR)))
Fig <- tmap::tmap_arrange(stress_map(Stress$LR, "Stress, long rains", Rng), stress_map(Stress$DRY, "Stress, dry season", Rng),
                          stress_map(Amp, "Dry minus long rains", c(-max(abs(av(Amp))), max(abs(av(Amp)))), "scico.vik"), ncol = 3)
tmap::tmap_save(Fig, file.path(Dir.Out.Figures, "10_stress_check.png"), width = 15, height = 6, dpi = 150)
Zcol <- c("#e6ab02", "#1f78b4", "#e31a1c", "#b2df8a", "#33a02c", "#6a3d9a")
Map.Z <- tmap::tm_shape(Hs) + tmap::tm_raster(col.scale = tmap::tm_scale_continuous(values = "brewer.greys"), col.legend = tmap::tm_legend_hide()) +
  tmap::tm_shape(Cat) + tmap::tm_raster(col.scale = tmap::tm_scale_categorical(values = Zcol, levels = Labels), col.legend = tmap::tm_legend(title = "Zone"), col_alpha = 0.85) +
  tmap::tm_shape(sf::st_geometry(Aoi.Sf)) + tmap::tm_borders(col = "black", lwd = 2) + tmap::tm_shape(Pans.Sf) + tmap::tm_symbols(fill = "cyan", col = "black", size = 0.6) +
  tmap::tm_text("pan_id", size = 0.6, ymod = -0.9) + tmap::tm_title("Recommendation zones (overlaps resolved B > A > C > E > D-dry > D-wet)", size = 0.8)
tmap::tmap_save(Map.Z, file.path(Dir.Out.Figures, "10_zones_check.png"), width = 8, height = 8, dpi = 170)
message("  wrote 10_stress_check.png and 10_zones_check.png")

log_phase("DONE")
