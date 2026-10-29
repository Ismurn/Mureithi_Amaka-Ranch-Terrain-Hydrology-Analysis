# 12_maps_figures.R  -  Figure factory
# PURPOSE
#   Report-ready maps and charts with stable file names (embedded by analysis/Amaka_Analysis.Rmd): terrain, RPI, dry-season versus
#   wet-season panels, riparian strips, per-pan catchment maps, the water gap, hypothesis scores, stress and recommendation zones.
#   tmap for maps, ggplot2 for charts; hillshade underlay on every map (DESCRIPTION.md 6.10).
#
# INPUTS
#   everything under outputs/ (rasters, tables, shapefiles) written by R/01-R/11
#
# OUTPUTS (outputs/figures/)
#   fig_terrain.png, fig_rpi_maps.png, fig_rpi_series.png, fig_pml_dry_vs_wet.png, fig_rainfall_series.png, fig_bsi_ndvi_seasons.png,
#   fig_landcover_sources.png, fig_riparian_map.png, fig_riparian_contrast.png, fig_riparian_gradient.png, fig_water_gap.png,
#   fig_pan_<pan_id>.png (one per pan with coordinates), fig_hypotheses.png, fig_stress.png, fig_zones.png, fig_zone_areas.png
#   (the check figures 01_*, 02_*, 03_*, 09_* and pan_siting_map.png are written by their own scripts)
#
# RUNTIME NOTES
#   * About a minute. Each figure is built inside build(), so one failing figure is reported at the end instead of stopping the rest.
#   * Palettes: scico roma (RPI), davos (ET), lajolla (stress), batlow (TWI), vik (diverging, suitability) - DESCRIPTION.md 6.10.
#   * Maps are cropped to the AOI plus a margin; coarse products (RPI 463 m, PML 500 m) are drawn as the blocks they are.
#     Stress and zone maps are masked to the AOI because they are scaled against the ranch and mean nothing outside it.
#   * Dry and wet panels of a variable share one colour scale, so the seasons can be compared by eye.

if (!exists("Setup.Loaded")) source("R/00_setup.R")

# PHASE 1: Common objects 
log_phase("PHASE 1: common objects")
Aoi.Sf <- read_aoi(); Aoi.U <- sf::st_union(sf::st_geometry(Aoi.Sf)); Aoi.Vec <- terra::vect(Aoi.Sf)
Dem <- terra::rast(file.path(Dir.DEM, "amaka_dem_utm37s.tif")); Aoi.Mask <- terra::rasterize(Aoi.Vec, Dem, field = 1)
bb_of <- function(g, margin) sf::st_bbox(sf::st_buffer(sf::st_union(sf::st_geometry(g)), margin))
bb_ext <- function(bb) terra::ext(bb[["xmin"]], bb[["xmax"]], bb[["ymin"]], bb[["ymax"]])
crop_bb <- function(r, bb) terra::crop(r, bb_ext(bb))
rep_to_dem <- function(r) terra::project(r, Dem, method = "near")                  # coarse layers: replicated onto the 30 m grid, not interpolated
rr  <- function(...) terra::rast(file.path(Dir.Out.Rasters, ...))
Hb  <- function(n) terra::rast(file.path(Dir.Out.Hydro.Buf, n))
tab <- function(f) readr::read_csv(file.path(Dir.Out.Tables, f), show_col_types = FALSE)
q_aoi <- function(r, p) global_q(terra::mask(terra::crop(r, terra::ext(Aoi.Vec)), Aoi.Vec), p)       # quantile inside the AOI, on any grid

Hs.All <- Hb("hillshade.tif"); names(Hs.All) <- "shade"
Pans <- sf::st_transform(sf::st_as_sf(dplyr::filter(read_pans(), located), coords = c("start_lon", "start_lat"), crs = Cfg.CRS.Geo, remove = FALSE), Cfg.CRS.Metric)
Pans$state <- ifelse(Pans$functional, "existing, functional", ifelse(Pans$status == "existing_desilting", "existing, desilting", "existing, degraded"))
Pans$label <- sub("AMK_", "", Pans$pan_id)
Pan.Cols <- c("existing, functional" = "#1f78b4", "existing, degraded" = "#e31a1c", "existing, desilting" = "#ff7f00")
Aoi.Bb <- bb_of(Aoi.Sf, 1500); Wide.Bb <- bb_of(c(sf::st_geometry(Aoi.Sf), sf::st_geometry(Pans)), 1500)    # AOI + margin; AOI + all pans + margin
Streams <- sf::st_read(file.path(Dir.Out.Hydro.Buf, "streams_strahler_utm.shp"), quiet = TRUE)

L.Shade <- function(bb = Aoi.Bb) tmap::tm_shape(crop_bb(Hs.All, bb)) + tmap::tm_raster(col.scale = tmap::tm_scale_continuous(values = "brewer.greys"), col.legend = tmap::tm_legend_hide())
L.Aoi   <- function() tmap::tm_shape(sf::st_geometry(Aoi.Sf)) + tmap::tm_borders(col = "black", lwd = 2)
L.Pans  <- function(size = 0.7) tmap::tm_shape(Pans) + tmap::tm_symbols(fill = "state", fill.scale = tmap::tm_scale_categorical(values = Pan.Cols), fill.legend = tmap::tm_legend_hide(), shape = 21, size = size, col = "white") +
                                    tmap::tm_text("label", size = 0.6, ymod = -0.9, col = "black")
map_r <- function(r, title, pal, limits = NULL, midpoint = NA, legend = "", bb = Aoi.Bb, mask.aoi = FALSE, pans = TRUE) {
  r <- crop_bb(r, bb); if (mask.aoi) r <- terra::mask(r, Aoi.Vec)
  if (!is.null(limits)) r <- terra::clamp(r, limits[1], limits[2])
  m <- L.Shade(bb) + tmap::tm_shape(r) + tmap::tm_raster(col.scale = tmap::tm_scale_continuous(values = pal, limits = limits, midpoint = midpoint), col.legend = tmap::tm_legend(title = legend)) + L.Aoi()
  if (pans) m <- m + L.Pans()
  m + tmap::tm_title(title, size = 0.8)
}
save_tm <- function(maps, name, ncol, w, h) tmap::tmap_save(tmap::tmap_arrange(maps, ncol = ncol), file.path(Dir.Out.Figures, name), width = w, height = h, dpi = 150)
save_gg <- function(p, name, w, h) ggplot2::ggsave(file.path(Dir.Out.Figures, name), p, width = w, height = h, dpi = 150)
Fig.Log <- tibble::tibble(figure = character(), status = character())
build <- function(name, expr) {
  res <- tryCatch({ force(expr); "ok" }, error = function(e) paste("FAILED:", conditionMessage(e)))
  Fig.Log <<- dplyr::bind_rows(Fig.Log, tibble::tibble(figure = name, status = res)); message(sprintf("  %-26s %s", name, res))
}
theme_amaka <- ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), plot.title.position = "plot", legend.position = "bottom")
wrap <- function(x, width = 95) stringr::str_wrap(x, width)                         # subtitles are long: wrap them so they are not clipped


# PHASE 2: Terrain
log_phase("PHASE 2: terrain, RPI and seasonal panels")
build("fig_terrain", {
  Str <- suppressWarnings(sf::st_crop(Streams, Aoi.Bb)); Str$order <- factor(Str$order)
  m.str <- L.Shade() + tmap::tm_shape(Str) + tmap::tm_lines(col = "order", col.scale = tmap::tm_scale_categorical(values = c("#9ecae1", "#4292c6", "#2171b5", "#08306b")[seq_along(levels(Str$order))]),
                                                              col.legend = tmap::tm_legend(title = "Strahler order"), lwd = 1.6) + L.Aoi() + L.Pans() + tmap::tm_title(sprintf("Streams (%s-cell threshold)", format(Cfg.Stream.Threshold, big.mark = ",")), size = 0.8)
  Twi <- Hb("twi.tif"); Slp <- Hb("slope_deg.tif")
  save_tm(list(m.str,
               map_r(Twi, "Topographic wetness index (TWI)", "scico.batlow", c(q_aoi(Twi, .02), q_aoi(Twi, .98)), legend = "TWI"),
               map_r(Slp, "Slope (degrees)", "scico.davos", c(0, q_aoi(Slp, .98)), legend = "degrees"),
               map_r(Dem, "Elevation (m, GLO-30)", "scico.imola", c(q_aoi(Dem, .01), q_aoi(Dem, .99)), legend = "m")), "fig_terrain.png", 2, 9, 9)
})

#  RPI 
build("fig_rpi_maps", {
  Rpi <- rep_to_dem(rr("rpi_now.tif")); Sen <- rep_to_dem(rr("rpi_sen_slope.tif"))
  lim.s <- max(abs(c(q_aoi(Sen, .02), q_aoi(Sen, .98)))) * c(-1, 1)
  save_tm(list(map_r(Rpi, "RPI now (mean of 2020-2022)", "scico.roma", c(q_aoi(Rpi, .02), q_aoi(Rpi, .98)), legend = "RPI"),
               map_r(Sen, "RPI trend 2000-2022 (Sen slope per year)", "scico.vik", lim.s, midpoint = 0, legend = "per year")), "fig_rpi_maps.png", 2, 10, 5)
})
build("fig_rpi_series", {
  z <- dplyr::filter(tab("rpi_zonal_annual.csv"), variable == "rpi"); tr <- dplyr::filter(tab("rpi_trend_summary.csv"), variable == "rpi")
  b0 <- stats::median(z$mean - tr$sen_slope_per_year * z$year)
  p <- ggplot2::ggplot(z, ggplot2::aes(year, mean)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = p10, ymax = p90), fill = "grey85") + ggplot2::geom_line(colour = "grey30") +
    ggplot2::geom_point(colour = "#08306b") + ggplot2::geom_abline(intercept = b0, slope = tr$sen_slope_per_year, colour = "firebrick", linetype = "dashed") +
    ggplot2::labs(x = NULL, y = "Rangeland Productivity Index", title = "RPI over Amaka, 2000-2022",
                  subtitle = wrap(sprintf("AOI mean (line) and 10th-90th percentile of cells (band); Sen slope %+.4f per year (dashed), Mann-Kendall p = %.2f: no significant trend", tr$sen_slope_per_year, tr$mk_p), 100)) + theme_amaka
  save_gg(p, "fig_rpi_series.png", 8, 4.2)
})

#PML dry versus wet
build("fig_pml_dry_vs_wet", {
  vars <- list(et = c("ET (mm per day)", "scico.davos"), gpp = c("GPP (gC per m2 per day)", "scico.bamako"), esfrac = c("Soil evaporation share of ET (Es/ET)", "scico.lajolla"))
  maps <- list()
  for (v in names(vars)) {
    rs <- list(DRY = rr(sprintf("pml_%s_DRY.tif", v)), LR = rr(sprintf("pml_%s_LR.tif", v)))
    lim <- range(unlist(lapply(rs, function(r) terra::values(crop_bb(r, Aoi.Bb)))), na.rm = TRUE)
    for (s in c("DRY", "LR")) maps[[paste(v, s)]] <- map_r(rs[[s]], sprintf("%s: %s", vars[[v]][1], c(DRY = "dry season", LR = "long rains")[[s]]), vars[[v]][2], lim)
  }
  save_tm(maps, "fig_pml_dry_vs_wet.png", 2, 9, 12)
})
build("fig_rainfall_series", {
  d <- tab("climate_context_series.csv") |> dplyr::filter(dataset == "CHIRPS", variable == "rainfall_mm", season %in% c("LR", "SR", "DRY")) |>
    dplyr::mutate(season = factor(season, c("LR", "SR", "DRY"), c("Long rains (Mar-May)", "Short rains (Oct-Dec)", "Dry (Jan-Feb + Jun-Sep)")),
                  year_type = dplyr::case_when(drought ~ "drought year (SPI <= -1)", wet ~ "wet year (SPI >= 1)", TRUE ~ "typical"))
  p <- ggplot2::ggplot(d, ggplot2::aes(year, value, fill = year_type)) + ggplot2::geom_col() + ggplot2::facet_wrap(~ season, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("drought year (SPI <= -1)" = "#d95f02", "typical" = "grey70", "wet year (SPI >= 1)" = "#1f78b4"), name = NULL) +
    ggplot2::labs(x = NULL, y = "Rainfall (mm, CHIRPS, Amaka AOI)", title = "Seasonal rainfall, 2000-2024") + theme_amaka
  save_gg(p, "fig_rainfall_series.png", 8, 7)
})

# Sentinel-2 seasons and land cover 
log_phase("PHASE 3: vegetation, riparian, pans")
build("fig_bsi_ndvi_seasons", {
  s2 <- function(v, s) terra::rast(file.path(Dir.Out.Rasters, "s2", sprintf("s2_%s_%s_mean.tif", v, s)))
  maps <- list()
  for (v in c("bsi", "ndvi")) {
    rs <- list(DRY = s2(v, "DRY"), LR = s2(v, "LR")); lim <- c(min(sapply(rs, q_aoi, p = .02)), max(sapply(rs, q_aoi, p = .98)))
    for (s in c("DRY", "LR")) maps[[paste(v, s)]] <- map_r(rs[[s]], sprintf("%s, %s (Sentinel-2, 10 m)", c(bsi = "Bare soil index", ndvi = "NDVI")[[v]], c(DRY = "dry season", LR = "long rains")[[s]]),
                                                           c(bsi = "scico.lajolla", ndvi = "scico.bamako")[[v]], lim)
  }
  save_tm(maps, "fig_bsi_ndvi_seasons.png", 2, 9, 9)
})
build("fig_landcover_sources", {
  d <- tab("landcover_fractions.csv") |> tidyr::pivot_longer(c(trees, shrubland, grassland, bare), names_to = "class", values_to = "fraction") |>
    dplyr::mutate(class = factor(class, c("trees", "shrubland", "grassland", "bare")), source = factor(source, rev(unique(source))))
  p <- ggplot2::ggplot(d, ggplot2::aes(fraction, source, fill = class)) + ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = c(trees = "#1b7837", shrubland = "#bf812d", grassland = "#a6d96a", bare = "#d9d9d9"), name = NULL) +
    ggplot2::scale_x_continuous(labels = scales::percent) + ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, 34)) +
    ggplot2::labs(x = "Share of area", y = NULL, title = "Land cover of Amaka by source", subtitle = wrap("About 94-97% shrubland whichever source is used; no bare class in the 10 m WorldCover", 100)) + theme_amaka
  save_gg(p, "fig_landcover_sources.png", 8, 3.6)
})

# Riparian
build("fig_riparian_map", {
  Rip <- sf::st_read(file.path(Dir.Out.Rasters, "riparian_order2plus.shp"), quiet = TRUE); Str <- suppressWarnings(sf::st_crop(Streams[Streams$order >= 2, ], Aoi.Bb))
  m <- L.Shade() + tmap::tm_shape(Rip) + tmap::tm_polygons(fill = "#4292c6", fill_alpha = 0.45, col = NA) + tmap::tm_shape(Str) + tmap::tm_lines(col = "#08306b", lwd = 1.4) + L.Aoi() + L.Pans() +
    tmap::tm_title(sprintf("Riparian strips: Strahler order >= 2, %g m each side (%.1f%% of the AOI)", Cfg.Riparian.m, tab("riparian_sets.csv")$pct_of_aoi[tab("riparian_sets.csv")$set == "order2plus"]), size = 0.8)
  tmap::tmap_save(m, file.path(Dir.Out.Figures, "fig_riparian_map.png"), width = 6, height = 6, dpi = 150)
})
build("fig_riparian_contrast", {
  nice <- c(ndvi_dry = "NDVI, dry season", ndvi_lr = "NDVI, long rains", bsi_dry = "Bare soil index, dry", ndwi_dry = "NDWI, dry", ndwi_lr = "NDWI, long rains", wc_trees = "Tree fraction", wc_shrub = "Shrub fraction",
            wc_grass = "Grass fraction", twi = "TWI", treecover2000 = "Hansen tree cover 2000", tree_loss = "Hansen tree loss", rpi_now = "RPI now", rpi_sen_slope = "RPI trend", esfrac_DRY = "Es/ET, dry",
            esfrac_LR = "Es/ET, long rains", gpp_DRY = "GPP, dry", gpp_LR = "GPP, long rains", wue_ANNUAL = "WUE, annual")
  d <- tab("riparian_condition.csv") |> dplyr::filter(set == "order2plus") |> dplyr::mutate(metric = factor(nice[metric], rev(nice)), resolves = ifelse(resolves_strip, "10-30 m layer: resolves the strip", "coarse layer: strip is mostly its surroundings"))
  p <- ggplot2::ggplot(d, ggplot2::aes(std_diff_vs_interior, metric, fill = resolves)) + ggplot2::geom_col() + ggplot2::geom_vline(xintercept = 0) +
    ggplot2::scale_fill_manual(values = c("10-30 m layer: resolves the strip" = "#2171b5", "coarse layer: strip is mostly its surroundings" = "grey70"), name = NULL) +
    ggplot2::labs(x = "Riparian minus interior, in ranch standard deviations", y = NULL, title = "Riparian strips against the rest of the ranch",
                  subtitle = "Effect size, not a significance test: neighbouring pixels are not independent") + theme_amaka
  save_gg(p, "fig_riparian_contrast.png", 8, 5.5)
})
build("fig_riparian_gradient", {
  d <- tab("riparian_gradient.csv") |> dplyr::mutate(dist_class = factor(dist_class, unique(dist_class))) |>
    tidyr::pivot_longer(c(ndvi_dry, bsi_dry, wc_shrub, wc_grass), names_to = "metric", values_to = "value") |>
    dplyr::mutate(metric = factor(metric, c("ndvi_dry", "bsi_dry", "wc_shrub", "wc_grass"), c("NDVI, dry season", "Bare soil index, dry", "Shrub fraction", "Grass fraction")))
  p <- ggplot2::ggplot(d, ggplot2::aes(dist_class, value, group = 1)) + ggplot2::geom_line(colour = "grey40") + ggplot2::geom_point(colour = "#08306b", size = 2) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 2) + ggplot2::labs(x = "Distance from the nearest order >= 2 stream", y = NULL, title = "Condition with distance from the streams") + theme_amaka +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  save_gg(p, "fig_riparian_gradient.png", 8, 5)
})

# Water gap and per-pan maps
build("fig_water_gap", {
  G <- rr("water_gap_all_radii.tif"); gap <- tab("water_gap_summary.csv"); Func <- Pans[Pans$functional, ]
  mk <- function(r.km) {
    g <- terra::ifel(G[[sprintf("gap_%gkm", r.km)]] == 1, 1, NA); names(g) <- "gap"; pct <- gap$gap_functional_pct[gap$radius_km == r.km]
    ring <- sf::st_cast(sf::st_boundary(sf::st_buffer(sf::st_geometry(Func), r.km * 1000)), "LINESTRING")
    L.Shade(Wide.Bb) + tmap::tm_shape(crop_bb(g, Wide.Bb)) + tmap::tm_raster(col.scale = tmap::tm_scale_categorical(values = "#fd8d3c", labels = "beyond reach of a working pan"), col.legend = tmap::tm_legend(title = "")) +
      tmap::tm_shape(ring) + tmap::tm_lines(col = "#1f78b4", lwd = 1.5, lty = "dashed") + L.Aoi() + L.Pans() +
      tmap::tm_title(sprintf("Water gap at %g km from the functional pan: %.0f%% of Amaka (%s ha)", r.km, pct, format(round(gap$gap_functional_ha[gap$radius_km == r.km]), big.mark = ",")), size = 0.8)
  }
  save_tm(list(mk(5), mk(10)), "fig_water_gap.png", 2, 11, 5.5)
})
build("fig_pan_maps", {
  Pc <- sf::st_read(file.path(Dir.Out.Rasters, "pan_catchments.shp"), quiet = TRUE); St <- tab("pan_catchment_static.csv"); Dg <- tab("pan_diagnostics.csv")
  for (id in Pans$pan_id) {
    poly <- Pc[Pc$pan_id == id, ]; s <- St[St$pan_id == id, ]; p <- Pans[Pans$pan_id == id, ]; bb <- bb_of(poly, 600)
    Str <- suppressWarnings(sf::st_crop(Streams, bb))
    ttl <- sprintf("%s (%s)\n%s: %.0f ha, mean slope %.1f deg, RPI %.2f, TWI %.1f%s", p$name, gsub("_", " ", sub("existing_", "", p$status)), ifelse(s$stream_fed, "Catchment", "500 m context disc, no catchment"),
                   s$catch_ha, s$slope_mean_deg, s$rpi_now, s$twi_mean, ifelse(isTRUE(s$catch_touches_dem_edge), " (truncated at the DEM edge)", ""))
    m <- L.Shade(bb) + tmap::tm_shape(poly) + tmap::tm_polygons(fill = "#4292c6", fill_alpha = 0.22, col = "#08306b", lwd = 2) +
      (if (nrow(Str)) tmap::tm_shape(Str) + tmap::tm_lines(col = "#2171b5", lwd = 1.2) else NULL) +
      tmap::tm_shape(sf::st_geometry(Aoi.Sf)) + tmap::tm_borders(col = "black", lwd = 2, lty = "dashed") +
      tmap::tm_shape(p) + tmap::tm_symbols(fill = Pan.Cols[[p$state]], col = "white", shape = 21, size = 1) + tmap::tm_text("label", size = 0.8, ymod = -1.1) + tmap::tm_title(ttl, size = 0.7)
    tmap::tmap_save(m, file.path(Dir.Out.Figures, sprintf("fig_pan_%s.png", id)), width = 6, height = 5.5, dpi = 150)
  }
})
build("fig_hypotheses", {
  h <- tab("pan_hypothesis_scores.csv"); pn <- Pans[, c("pan_id", "name", "status")] |> sf::st_drop_geometry()
  d <- dplyr::left_join(h, pn, by = "pan_id") |> dplyr::filter(status == "existing_degraded") |>
    dplyr::arrange(pan_id) |> dplyr::mutate(hypothesis = factor(sub("^H[0-9] ", "", hypothesis), rev(unique(sub("^H[0-9] ", "", hypothesis)))), pan = factor(paste0(name, " (", sub("AMK_", "", pan_id), ")"), unique(paste0(name, " (", sub("AMK_", "", pan_id), ")"))))
  p <- ggplot2::ggplot(d, ggplot2::aes(score, hypothesis)) + ggplot2::geom_col(fill = "#e31a1c", alpha = 0.75, width = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = control_score, shape = "functional dam, same test"), size = 3, colour = "#1f78b4") +
    ggplot2::facet_wrap(~ pan, ncol = 1) + ggplot2::scale_shape_manual(values = c("functional dam, same test" = 18), name = NULL) + ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Evidence score (0.5 = typical of the ranch)", y = NULL, title = "Why might the degraded pans fail? Screening evidence from satellite data",
                  subtitle = wrap("A hypothesis the working dam also scores high on does not explain the failure. The HQ rock-catchment reservoir is not stream-fed and is not scored.", 95)) + theme_amaka
  save_gg(p, "fig_hypotheses.png", 8, 5.5)
})

# Stress and zones 
log_phase("PHASE 4: stress and zones")
build("fig_stress", {
  Lr <- rr("stress_LR.tif"); Dry <- rr("stress_DRY.tif"); Amp <- rr("stress_amplification_DRY_minus_LR.tif")
  lim <- c(min(q_aoi(Lr, .01), q_aoi(Dry, .01)), max(q_aoi(Lr, .99), q_aoi(Dry, .99))); la <- max(abs(c(q_aoi(Amp, .01), q_aoi(Amp, .99)))) * c(-1, 1)
  save_tm(list(map_r(Lr, "Stress index, long rains", "scico.lajolla", lim, mask.aoi = TRUE, legend = "0-1"), map_r(Dry, "Stress index, dry season", "scico.lajolla", lim, mask.aoi = TRUE, legend = "0-1"),
               map_r(Amp, "Dry minus long rains", "scico.vik", la, midpoint = 0, mask.aoi = TRUE, legend = "")), "fig_stress.png", 3, 14, 5.5)
})
build("fig_zones", {
  Z <- rr("recommendation_zones.tif"); Labels <- c("A regeneration feasibility", "B catchment interventions", "C pan protection", "D-wet comfort", "D-dry comfort", "E waterless buffer")
  Zcol <- c("#e6ab02", "#1f78b4", "#e31a1c", "#b2df8a", "#33a02c", "#6a3d9a")
  m <- L.Shade() + tmap::tm_shape(crop_bb(Z, Aoi.Bb)) + tmap::tm_raster(col.scale = tmap::tm_scale_categorical(values = Zcol, levels = Labels), col.legend = tmap::tm_legend(title = "Zone"), col_alpha = 0.85) +
    L.Aoi() + L.Pans() + tmap::tm_title("Recommendation zones (overlaps resolved B > A > C > E > D-dry > D-wet)", size = 0.8)
  tmap::tmap_save(m, file.path(Dir.Out.Figures, "fig_zones.png"), width = 7.5, height = 6.5, dpi = 150)
})
build("fig_zone_areas", {
  d <- tab("recommendation_zone_summary.csv") |> dplyr::mutate(zone = factor(zone, rev(zone)))
  p <- ggplot2::ggplot(d, ggplot2::aes(overlapping_pct_of_aoi, zone)) + ggplot2::geom_col(fill = "grey60", width = 0.6) + ggplot2::geom_col(ggplot2::aes(x = exclusive_pct_of_aoi), fill = "#2171b5", width = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%s ha", format(round(overlapping_ha), big.mark = ","))), hjust = -0.1, size = 3) + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
    ggplot2::labs(x = "% of the Amaka AOI", y = NULL, title = "Area of each recommendation zone", subtitle = wrap("Grey = zone as defined (zones overlap); blue = what is left after overlaps are resolved", 100)) + theme_amaka
  save_gg(p, "fig_zone_areas.png", 8, 3.8)
})

log_phase("SUMMARY")
print(as.data.frame(Fig.Log), row.names = FALSE)
if (any(Fig.Log$status != "ok")) warning(sum(Fig.Log$status != "ok"), " figure(s) failed: see the table above", call. = FALSE)
log_phase("DONE")
