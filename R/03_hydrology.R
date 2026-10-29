# 03_hydrology.R  -  WhiteboxTools terrain and flow routing
# PURPOSE
#   Condition the buffered DEM and derive the hydrological layers every later stage
#   uses: breach, D8 / D-infinity flow, accumulation, streams, Strahler order, specific
#   contributing area, slope, TWI, TPI, hillshade (DESCRIPTION.md 6.3). Also CALIBRATES the
#   stream threshold (9.7) instead of accepting the starting value of 1000 cells.
#
# INPUTS
#   data/dem/amaka_dem_utm37s.tif (02)      data/gee_exports/s2_ndvi_mean_DRY_*.tif (calibration)
#   data/water_pans/*.csv (dam locations used as calibration evidence)
#
# OUTPUTS
#   outputs/rasters/hydrology/_buffered/   full-envelope WBT rasters: dem_breached, d8_pointer,
#       dinf_pointer, flow_accum_d8, sca_dinf, slope_deg, twi, tpi_500m, rtp_500m, hillshade,
#       streams, strahler, streams_strahler.shp
#   outputs/rasters/*.tif   the same layers masked to the Amaka AOI for reporting
#       (twi, slope_deg, sca, flow_accum, streams, strahler, tpi_500m, rtp_500m, hillshade)
#   outputs/tables/03_stream_threshold_calibration.csv, 03_hydrology_summary.csv
#
# RUNTIME NOTES
#   * About a minute. Re-run after changing Cfg.Stream.Threshold or the DEM.
#   * Breach, not fill (9.6): wbt_breach_depressions_least_cost() first, with fill = TRUE only
#     for the pits that cannot be breached within `dist`, so flow stays connected.
#   * Slope, TWI and TPI come from the ORIGINAL DEM (the breached DEM is only for routing).
#   * The terrain is a very flat plain with one rocky inselberg, and GLO-30 is a surface model
#     (canopy included): flow paths on the plain are only as good as the DEM and are checked
#     against imagery in the calibration below.
#   * Flow routing is done by WhiteboxTools only (4.3); terra is used for masking and distances.

if (!exists("Setup.Loaded")) source("R/00_setup.R")
whitebox::wbt_verbose(FALSE)


# PHASE 1: Inputs and helpers
log_phase("PHASE 1: inputs")

Dem.Path <- file.path(Dir.DEM, "amaka_dem_utm37s.tif")
if (!file.exists(Dem.Path)) stop("Run R/02_dem_prep.R first: ", Dem.Path, " is missing.", call. = FALSE)

AOI.Sf  <- read_aoi()
AOI.Vec <- terra::vect(AOI.Sf)
DEM     <- terra::rast(Dem.Path)
Cell.m  <- terra::res(DEM)[1]
Dem.In  <- normalizePath(Dem.Path, winslash = "/")

wb <- function(name) normalizePath(file.path(Dir.Out.Hydro.Buf, name), winslash = "/", mustWork = FALSE)
need <- function(f) { if (!file.exists(f)) stop("WhiteboxTools did not write ", f, call. = FALSE); invisible(f) }
rd <- function(name) {                                   # read a WBT raster and make sure it carries the metric CRS
  r <- terra::rast(wb(name)); if (is.na(terra::crs(r, describe = TRUE)$code)) terra::crs(r) <- paste0("EPSG:", Cfg.CRS.Metric); r
}
Win.Cells <- 2 * floor((Cfg.TPI.Window.m / Cell.m) / 2) + 1          # odd number of cells covering ~500 m
ext1 <- function(r, xy) { v <- terra::extract(r, xy); v[[ncol(v)]] }  # values of a single-layer raster at points (with or without an ID column)
message(sprintf("  DEM %d x %d cells at %g m | TPI window %d cells", terra::ncol(DEM), terra::nrow(DEM), Cell.m, Win.Cells))


# PHASE 2: Conditioning and flow routing 
# dist = 100 cells (3 km): the farthest a breach path may run to reach a lower cell.
log_phase("PHASE 2: breach depressions, D8 / D-infinity flow")

whitebox::wbt_breach_depressions_least_cost(dem = Dem.In, output = wb("dem_breached.tif"), dist = 100, fill = TRUE)
need(wb("dem_breached.tif"))
whitebox::wbt_d8_pointer(dem = wb("dem_breached.tif"), output = wb("d8_pointer.tif"));      need(wb("d8_pointer.tif"))
whitebox::wbt_d_inf_pointer(dem = wb("dem_breached.tif"), output = wb("dinf_pointer.tif")); need(wb("dinf_pointer.tif"))
whitebox::wbt_d8_flow_accumulation(input = wb("dem_breached.tif"), output = wb("flow_accum_d8.tif"), out_type = "cells")
need(wb("flow_accum_d8.tif"))
whitebox::wbt_d_inf_flow_accumulation(input = wb("dem_breached.tif"), output = wb("sca_dinf.tif"),
                                      out_type = "Specific Contributing Area")
need(wb("sca_dinf.tif"))

Breach.Diff <- terra::rast(wb("dem_breached.tif")) - DEM
message(sprintf("  breaching/filling changed %.1f%% of cells (max lowering %.1f m, max raising %.1f m)",
                100 * terra::global(abs(Breach.Diff) > 1e-6, "mean")[[1]],
                -terra::global(Breach.Diff, "min", na.rm = TRUE)[[1]], terra::global(Breach.Diff, "max", na.rm = TRUE)[[1]]))


# PHASE 3: Terrain attributes (from the ORIGINAL DEM) 
log_phase("PHASE 3: slope, TWI, TPI, hillshade")

whitebox::wbt_slope(dem = Dem.In, output = wb("slope_deg.tif"), units = "degrees"); need(wb("slope_deg.tif"))
whitebox::wbt_wetness_index(sca = wb("sca_dinf.tif"), slope = wb("slope_deg.tif"), output = wb("twi.tif")); need(wb("twi.tif"))
whitebox::wbt_diff_from_mean_elev(dem = Dem.In, output = wb("tpi_500m.tif"), filterx = Win.Cells, filtery = Win.Cells)
need(wb("tpi_500m.tif"))                                  # TPI = z - mean(z) in the window, metres (>0 ridge, <0 depression)
whitebox::wbt_relative_topographic_position(dem = Dem.In, output = wb("rtp_500m.tif"), filterx = Win.Cells, filtery = Win.Cells)
need(wb("rtp_500m.tif"))
whitebox::wbt_hillshade(dem = Dem.In, output = wb("hillshade.tif"), azimuth = 315, altitude = 30); need(wb("hillshade.tif"))

Slope.R <- rd("slope_deg.tif"); Twi.R <- rd("twi.tif"); Acc.R <- rd("flow_accum_d8.tif")
message(sprintf("  slope: median %.2f deg, 95th %.2f deg | TWI: median %.1f, range %.1f to %.1f",
                global_q(Slope.R, 0.50),
                global_q(Slope.R, 0.95),
                global_q(Twi.R, 0.50),
                terra::global(Twi.R, "min", na.rm = TRUE)[[1]], terra::global(Twi.R, "max", na.rm = TRUE)[[1]]))


# PHASE 4: Stream threshold calibration (DESCRIPTION.md 9.7) 
# The threshold decides what counts as a stream. Three pieces of evidence for each candidate:
#   (a) the dams: they are built on drainage lines, so the distance from each located pan to the
#       nearest stream cell should be small;
#   (b) Sentinel-2 dry-season NDVI: real channels carry denser woody vegetation, so NDVI in a corridor
#       along the streams should exceed NDVI in a reference ring beside them (plain only, slope < 3 deg,
#       so the inselberg does not dominate);
#   (c) drainage density: a threshold that makes the whole plain a stream network is not credible.
log_phase("PHASE 4: stream threshold calibration")

Ndvi.Dry <- terra::app(terra::rast(gee_file("s2_ndvi_mean_DRY_2019_2024")), "mean", na.rm = TRUE)   # 10 m, mean of 2019-2024
Pans.Df  <- read_pans(); Pans.Df <- Pans.Df[Pans.Df$located, ]
Pans.Sf  <- sf::st_transform(sf::st_as_sf(Pans.Df, coords = c("start_lon", "start_lat"), crs = Cfg.CRS.Geo, remove = FALSE), Cfg.CRS.Metric)
Env.Km2  <- as.numeric(terra::expanse(terra::as.polygons(terra::ext(DEM), crs = terra::crs(DEM)), unit = "km"))
Plain    <- Slope.R < 3

set.seed(42)
Calib <- lapply(Cfg.Stream.Calib.Thresholds, function(thr) {
  f <- wb(sprintf("streams_T%d.tif", thr))
  whitebox::wbt_extract_streams(flow_accum = wb("flow_accum_d8.tif"), output = f, threshold = thr, zero_background = TRUE)
  s  <- terra::rast(f); if (is.na(terra::crs(s, describe = TRUE)$code)) terra::crs(s) <- paste0("EPSG:", Cfg.CRS.Metric)
  sm <- terra::ifel(s > 0, 1, NA)
  d  <- terra::distance(sm)                                           # metres to the nearest stream cell
  n.str <- terra::global(s > 0, "sum", na.rm = TRUE)[[1]]
  cor.c <- terra::ifel(d <= 1.5 * Cell.m & Plain, 1, NA); ring.c <- terra::ifel(d >= 3 * Cell.m & d <= 6 * Cell.m & Plain, 1, NA)
  pick  <- function(m, n = 15000) { xy <- terra::xyFromCell(m, which(!is.na(terra::values(m)))); if (nrow(xy) > n) xy <- xy[sample(nrow(xy), n), , drop = FALSE]; xy }
  v.cor <- ext1(Ndvi.Dry, pick(cor.c)); v.ring <- ext1(Ndvi.Dry, pick(ring.c))
  contrast <- mean(v.cor, na.rm = TRUE) - mean(v.ring, na.rm = TRUE)
  pooled   <- sqrt((var(v.cor, na.rm = TRUE) + var(v.ring, na.rm = TRUE)) / 2)
  dpan <- terra::extract(d, terra::vect(Pans.Sf))[, 2]
  tibble::tibble(threshold_cells = thr, upslope_km2 = thr * Cell.m^2 / 1e6, stream_km = n.str * Cell.m / 1000,
                 drainage_km_per_km2 = n.str * Cell.m / 1000 / Env.Km2,
                 ndvi_corridor = mean(v.cor, na.rm = TRUE), ndvi_ring = mean(v.ring, na.rm = TRUE),
                 ndvi_contrast = contrast, contrast_d = contrast / pooled,
                 pans_within_100m = sum(dpan <= 100), pans_n = length(dpan), pan_median_dist_m = median(dpan),
                 pan_max_dist_m = max(dpan))
})
Calib <- dplyr::bind_rows(Calib)
readr::write_csv(Calib, file.path(Dir.Out.Tables, "03_stream_threshold_calibration.csv"))
print(as.data.frame(Calib |> dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

# Rule: prefer the threshold that best separates real channels from their surroundings in the imagery
# (largest standardised NDVI contrast d); ties within 10% go to the one that puts more dams on a stream.
Best <- Calib |> dplyr::filter(contrast_d >= 0.9 * max(contrast_d, na.rm = TRUE)) |>
  dplyr::arrange(dplyr::desc(pans_within_100m), pan_median_dist_m) |> dplyr::slice(1)
message(sprintf("  best-supported threshold by the rule: %d cells (%.2f km2): contrast d = %.2f, %d of %d dams within 100 m",
                Best$threshold_cells, Best$upslope_km2, Best$contrast_d, Best$pans_within_100m, Best$pans_n))
message("  Cfg.Stream.Threshold (00_setup.R) = ", Cfg.Stream.Threshold, " cells is used for the products below.")


# PHASE 5: Streams, Strahler order, vector log_phase("PHASE 5: streams at the configured threshold")

whitebox::wbt_extract_streams(flow_accum = wb("flow_accum_d8.tif"), output = wb("streams.tif"),
                              threshold = Cfg.Stream.Threshold, zero_background = TRUE);                    need(wb("streams.tif"))
whitebox::wbt_strahler_stream_order(d8_pntr = wb("d8_pointer.tif"), streams = wb("streams.tif"),
                                    output = wb("strahler.tif"), zero_background = TRUE);                  need(wb("strahler.tif"))
whitebox::wbt_raster_streams_to_vector(streams = wb("strahler.tif"), d8_pntr = wb("d8_pointer.tif"),
                                       output = wb("streams_strahler.shp"));                              need(wb("streams_strahler.shp"))

Streams.Sf <- sf::st_read(wb("streams_strahler.shp"), quiet = TRUE)
sf::st_crs(Streams.Sf) <- Cfg.CRS.Metric
names(Streams.Sf)[names(Streams.Sf) %in% c("STRM_VAL", "VALUE")] <- "order"
if (!"order" %in% names(Streams.Sf)) stop("No Strahler order attribute in the stream vector: ", paste(names(Streams.Sf), collapse = ", "), call. = FALSE)
Streams.Sf$length_m <- as.numeric(sf::st_length(Streams.Sf))
Streams.Aoi <- suppressWarnings(sf::st_intersection(Streams.Sf, sf::st_union(sf::st_geometry(AOI.Sf))))
Streams.Aoi$length_m <- as.numeric(sf::st_length(Streams.Aoi))
write_shp(Streams.Sf[, c("order", "length_m")],  file.path(Dir.Out.Hydro.Buf, "streams_strahler_utm.shp"))
write_shp(Streams.Aoi[, c("order", "length_m")], file.path(Dir.Out.Rasters, "streams_aoi.shp"))

Len.Tbl <- dplyr::full_join(
  Streams.Sf |> sf::st_drop_geometry() |> dplyr::group_by(order) |> dplyr::summarise(envelope_km = sum(length_m) / 1000, .groups = "drop"),
  Streams.Aoi |> sf::st_drop_geometry() |> dplyr::group_by(order) |> dplyr::summarise(aoi_km = sum(length_m) / 1000, .groups = "drop"),
  by = "order") |> dplyr::arrange(order) |> dplyr::mutate(dplyr::across(where(is.numeric), ~ tidyr::replace_na(round(.x, 2), 0)))
print(as.data.frame(Len.Tbl), row.names = FALSE)


# PHASE 6: Mask to the AOI for reporting; keep the buffered versions
log_phase("PHASE 6: AOI-masked copies")

mask_aoi <- function(r) terra::mask(terra::crop(r, AOI.Vec), AOI.Vec)
Report.Map <- c(twi = "twi.tif", slope_deg = "slope_deg.tif", sca = "sca_dinf.tif", flow_accum = "flow_accum_d8.tif",
                streams = "streams.tif", strahler = "strahler.tif", tpi_500m = "tpi_500m.tif", rtp_500m = "rtp_500m.tif",
                hillshade = "hillshade.tif")
for (nm in names(Report.Map)) {
  r <- rd(Report.Map[[nm]])
  terra::writeRaster(mask_aoi(r), file.path(Dir.Out.Rasters, paste0(nm, ".tif")), overwrite = TRUE)
}
message("  wrote ", length(Report.Map), " AOI-masked rasters to outputs/rasters/")


# PHASE 7: Summary
log_phase("PHASE 7: summary table")

Aoi.R <- function(name) terra::mask(terra::crop(rd(name), AOI.Vec), AOI.Vec)
gstat <- function(r) c(mean = terra::global(r, "mean", na.rm = TRUE)[[1]], p05 = global_q(r, 0.05),
                       p50 = global_q(r, 0.50), p95 = global_q(r, 0.95))
Hydro.Summary <- dplyr::bind_rows(
  tibble::tibble(item = "dem_source", value = Cfg.DEM.Source),
  tibble::tibble(item = "dem_smoothed", value = as.character(Cfg.DEM.Smooth)),
  tibble::tibble(item = "stream_threshold_cells_used", value = as.character(Cfg.Stream.Threshold)),
  tibble::tibble(item = "stream_threshold_cells_best_by_rule", value = as.character(Best$threshold_cells)),
  tibble::tibble(item = "tpi_window_cells", value = as.character(Win.Cells)),
  tibble::tibble(item = c("aoi_slope_deg_mean", "aoi_slope_deg_p05", "aoi_slope_deg_p50", "aoi_slope_deg_p95"), value = as.character(round(gstat(Aoi.R("slope_deg.tif")), 3))),
  tibble::tibble(item = c("aoi_twi_mean", "aoi_twi_p05", "aoi_twi_p50", "aoi_twi_p95"), value = as.character(round(gstat(Aoi.R("twi.tif")), 2))),
  tibble::tibble(item = c("aoi_tpi500_mean_m", "aoi_tpi500_p05_m", "aoi_tpi500_p50_m", "aoi_tpi500_p95_m"), value = as.character(round(gstat(Aoi.R("tpi_500m.tif")), 2))),
  tibble::tibble(item = "aoi_stream_km_all_orders", value = as.character(round(sum(Len.Tbl$aoi_km), 2))),
  tibble::tibble(item = "aoi_drainage_density_km_per_km2", value = as.character(round(sum(Len.Tbl$aoi_km) / (sum(as.numeric(sf::st_area(AOI.Sf))) / 1e6), 3))),
  tibble::tibble(item = "max_strahler_order", value = as.character(max(Len.Tbl$order)))
)
readr::write_csv(Hydro.Summary, file.path(Dir.Out.Tables, "03_hydrology_summary.csv"))
print(as.data.frame(Hydro.Summary), row.names = FALSE)

# PHASE 7b: Cell-level terrain table (keyed on the RPI-native analysis grid)
# Fine rasters are aggregated into each ~460 m cell with area-weighted means (no resampling, DESCRIPTION.md 9.2);
# stream_frac = share of the cell covered by stream cells, acc_max = largest D8 upslope area in the cell.
log_phase("PHASE 7b: cell-level terrain table")

Grid.List <- get_analysis_grid(); Cells <- Grid.List$cells
Cell.Poly <- cell_polygons(Grid.List$grid)
Terrain.Mean <- c(DEM, rd("twi.tif"), rd("slope_deg.tif"), rd("tpi_500m.tif"), rd("rtp_500m.tif"), (rd("streams.tif") > 0) * 1)
names(Terrain.Mean) <- c("elev_mean", "twi_mean", "slope_mean", "tpi_mean", "rtp_mean", "stream_frac")
Tm <- zonal_stat(Terrain.Mean, Cell.Poly, "mean")
Acc.Max <- zonal_stat(Acc.R, Cell.Poly, "max");  names(Acc.Max) <- "acc_max_cells"
Ord.Max <- zonal_stat(rd("strahler.tif"), Cell.Poly, "max"); names(Ord.Max) <- "strahler_max"
Cells.Terrain <- dplyr::bind_cols(tibble::tibble(cell_id = Cell.Poly$cell_id), Tm, Acc.Max, Ord.Max)
Cells.Terrain$strahler_max[is.na(Cells.Terrain$strahler_max)] <- 0
readr::write_csv(Cells.Terrain, file.path(Dir.Out.Tables, "cells_terrain.csv"))
Chk <- dplyr::left_join(Cells, Cells.Terrain, by = "cell_id") |> dplyr::filter(in_aoi)
message(sprintf("  %d cells; inside the AOI: elevation %.0f-%.0f m, TWI mean %.1f, slope mean %.2f deg, %.0f%% of cells carry a stream",
                nrow(Cells.Terrain), min(Chk$elev_mean), max(Chk$elev_mean), mean(Chk$twi_mean), mean(Chk$slope_mean),
                100 * mean(Chk$stream_frac > 0)))


# PHASE 8: Check figures
log_phase("PHASE 8: check figures")

Hs.R <- rd("hillshade.tif"); names(Hs.R) <- "shade"
Streams.Plot <- Streams.Sf; Streams.Plot$order <- factor(Streams.Plot$order)
Map.Streams <-
  tmap::tm_shape(Hs.R) + tmap::tm_raster(col.scale = tmap::tm_scale_continuous(values = "brewer.greys"), col.legend = tmap::tm_legend_hide()) +
  tmap::tm_shape(Streams.Plot) + tmap::tm_lines(col = "order", col.scale = tmap::tm_scale_categorical(values = c("#9ecae1", "#4292c6", "#2171b5", "#08306b")),
                                                col.legend = tmap::tm_legend(title = "Strahler order"), lwd = 1.3) +
  tmap::tm_shape(sf::st_geometry(AOI.Sf)) + tmap::tm_borders(col = "black", lwd = 2.5) +
  tmap::tm_shape(Pans.Sf) + tmap::tm_symbols(fill = "orange", col = "black", size = 0.7) + tmap::tm_text("pan_id", size = 0.7, ymod = -0.9) +
  tmap::tm_title(sprintf("Streams at %d cells (%.1f km2 upslope), Amaka AOI and pans", Cfg.Stream.Threshold, Cfg.Stream.Threshold * Cell.m^2 / 1e6), size = 0.8)
tmap::tmap_save(Map.Streams, file.path(Dir.Out.Figures, "03_streams_check_map.png"), width = 8, height = 8, dpi = 170)

Calib.Plot <- Calib |>
  tidyr::pivot_longer(c(ndvi_contrast, pans_within_100m, drainage_km_per_km2), names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric = factor(metric, c("pans_within_100m", "ndvi_contrast", "drainage_km_per_km2"),
                                c("Dams within 100 m of a stream (of 4)", "Dry-season NDVI: corridor minus ring", "Drainage density (km per km2)")))
Fig.Calib <- ggplot2::ggplot(Calib.Plot, ggplot2::aes(threshold_cells, value)) +
  ggplot2::geom_line(colour = "grey40") + ggplot2::geom_point(size = 2.5, colour = "#08306b") +
  ggplot2::geom_vline(xintercept = Cfg.Stream.Threshold, linetype = "dashed", colour = "firebrick") +
  ggplot2::scale_x_log10(breaks = Cfg.Stream.Calib.Thresholds) + ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  ggplot2::labs(x = "Stream threshold (D8 cells, log scale)", y = NULL, title = "Stream-threshold calibration (DESCRIPTION.md 9.7)",
                subtitle = "Dashed line = threshold used") + ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(file.path(Dir.Out.Figures, "03_stream_calibration.png"), Fig.Calib, width = 6, height = 7, dpi = 170)
message("  wrote 03_streams_check_map.png and 03_stream_calibration.png")

log_phase("DONE")
