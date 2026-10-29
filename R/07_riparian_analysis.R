# 07_riparian_analysis.R  -  Riparian corridor condition
# PURPOSE
#   Buffer streams of Strahler order >= 2 by Cfg.Riparian.m (125 m) and report their vegetation condition
#   separately from the rest of the ranch: riparian versus interior contrast for every vegetation layer,
#   expressed against the ranch-wide median (DESCRIPTION.md 6.5). Rivers are the dominant baseline attractor for
#   large mammals (Masiaine 2021), so riparian condition is reported on its own. Also a distance-to-stream
#   gradient, which shows WHERE along a corridor the condition changes.
#
# INPUTS
#   outputs/rasters/hydrology/_buffered/streams_strahler_utm.shp (03)
#   outputs/rasters/s2/*.tif (06), data/gee_exports/{worldcover,hansen}_*.tif, outputs/rasters/{rpi_now,pml_*}.tif
#
# OUTPUTS
#   outputs/tables/riparian_sets.csv         area of each riparian set inside the AOI
#   outputs/tables/riparian_condition.csv    metric x set: mean, ranch-wide median, difference, interior contrast
#   outputs/tables/riparian_gradient.csv     condition by distance class from order >= 2 streams
#   outputs/rasters/riparian_mask.tif, riparian_order2plus.shp
#
# RUNTIME NOTES
#   * Under a minute. Only the 10-30 m layers (Sentinel-2, WorldCover, TWI, Hansen) can resolve a 250 m strip;
#     RPI (~463 m) and PML (500 m) are included as the brief asks but each strip pixel is mostly its
#     surroundings (flagged `resolves_strip = FALSE` in the table).
#   * Contrast = mean in the riparian set minus the ranch-wide MEDIAN (brief) and minus the interior mean.
#     std_diff = (riparian - interior) / SD of the layer across the ranch: an effect size, not a p-value
#     (neighbouring pixels are not independent).
#   * Within the AOI the extracted network has no order >= 4 stream (known from 03), so the strongest
#     corridors are order 2 and 3.

if (!exists("Setup.Loaded")) source("R/00_setup.R")

AOI.Sf  <- read_aoi()
AOI.U   <- sf::st_union(sf::st_geometry(AOI.Sf))
Streams <- sf::st_read(file.path(Dir.Out.Hydro.Buf, "streams_strahler_utm.shp"), quiet = TRUE)
Max.Ord <- max(Streams$order)


# PHASE 1: Riparian polygons
log_phase("PHASE 1: riparian polygons")

buf_set <- function(orders) {
  s <- Streams[Streams$order %in% orders, ]
  if (!nrow(s)) return(NULL)
  g <- suppressWarnings(sf::st_intersection(sf::st_union(sf::st_buffer(sf::st_geometry(s), Cfg.Riparian.m)), AOI.U))
  if (!length(g) || all(sf::st_is_empty(g))) NULL else g
}
Sets <- list(order2 = buf_set(2), order3 = buf_set(3), order4plus = if (Max.Ord >= 4) buf_set(4:Max.Ord) else NULL, order2plus = buf_set(2:Max.Ord))
Sets <- Sets[!vapply(Sets, is.null, logical(1))]
Sets$interior <- sf::st_difference(AOI.U, Sets$order2plus)
Sets$aoi      <- AOI.U
Sets.Sf <- sf::st_sf(set = names(Sets), geometry = do.call(c, lapply(Sets, function(g) sf::st_geometry(sf::st_sfc(sf::st_union(g), crs = Cfg.CRS.Metric)))))
Sets.Tbl <- tibble::tibble(set = Sets.Sf$set, area_ha = round(as.numeric(sf::st_area(Sets.Sf)) / 1e4, 1))
Sets.Tbl$pct_of_aoi <- round(100 * Sets.Tbl$area_ha / Sets.Tbl$area_ha[Sets.Tbl$set == "aoi"], 1)
Len <- Streams |> sf::st_drop_geometry() |> dplyr::group_by(order) |> dplyr::summarise(length_km_envelope = sum(length_m) / 1000, .groups = "drop")
Aoi.Len <- suppressWarnings(sf::st_intersection(Streams, AOI.U)); Aoi.Len$len <- as.numeric(sf::st_length(Aoi.Len))
Sets.Tbl$stream_km_in_aoi <- vapply(Sets.Tbl$set, function(s) switch(s, order2 = sum(Aoi.Len$len[Aoi.Len$order == 2]), order3 = sum(Aoi.Len$len[Aoi.Len$order == 3]),
  order4plus = sum(Aoi.Len$len[Aoi.Len$order >= 4]), order2plus = sum(Aoi.Len$len[Aoi.Len$order >= 2]), NA_real_) / 1000, 0)
readr::write_csv(Sets.Tbl, file.path(Dir.Out.Tables, "riparian_sets.csv"))
print(as.data.frame(Sets.Tbl |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
write_shp(sf::st_sf(set = "order2plus", geometry = sf::st_sfc(sf::st_union(Sets$order2plus), crs = Cfg.CRS.Metric)), file.path(Dir.Out.Rasters, "riparian_order2plus.shp"))


# PHASE 2: Riparian versus interior condition, per layer 
log_phase("PHASE 2: riparian vs interior condition")

Dir.S2 <- file.path(Dir.Out.Rasters, "s2")
Hansen <- terra::rast(gee_file("hansen_gfc_v1_11_envelope")); Wc <- terra::rast(gee_file("worldcover_2021_envelope"))
rr <- function(...) terra::rast(file.path(Dir.Out.Rasters, ...))
Metrics <- list(
  ndvi_dry      = list(r = terra::rast(file.path(Dir.S2, "s2_ndvi_DRY_mean.tif")), res = "10 m"),
  ndvi_lr       = list(r = terra::rast(file.path(Dir.S2, "s2_ndvi_LR_mean.tif")),  res = "10 m"),
  bsi_dry       = list(r = terra::rast(file.path(Dir.S2, "s2_bsi_DRY_mean.tif")),  res = "10 m"),
  ndwi_dry      = list(r = terra::rast(file.path(Dir.S2, "s2_ndwi_DRY_mean.tif")), res = "10 m"),
  ndwi_lr       = list(r = terra::rast(file.path(Dir.S2, "s2_ndwi_LR_mean.tif")),  res = "10 m"),
  wc_trees      = list(r = (Wc == 10) * 1, res = "10 m"),
  wc_shrub      = list(r = (Wc == 20) * 1, res = "10 m"),
  wc_grass      = list(r = (Wc == 30) * 1, res = "10 m"),
  twi           = list(r = terra::rast(file.path(Dir.Out.Hydro.Buf, "twi.tif")), res = "30 m"),
  treecover2000 = list(r = Hansen[["treecover2000"]], res = "30 m"),
  tree_loss     = list(r = Hansen[["loss"]], res = "30 m"),
  rpi_now       = list(r = rr("rpi_now.tif"), res = "463 m"),
  rpi_sen_slope = list(r = rr("rpi_sen_slope.tif"), res = "463 m"),
  esfrac_DRY    = list(r = rr("pml_esfrac_DRY.tif"), res = "500 m"),
  esfrac_LR     = list(r = rr("pml_esfrac_LR.tif"), res = "500 m"),
  gpp_DRY       = list(r = rr("pml_gpp_DRY.tif"), res = "500 m"),
  gpp_LR        = list(r = rr("pml_gpp_LR.tif"), res = "500 m"),
  wue_ANNUAL    = list(r = rr("pml_wue_ANNUAL.tif"), res = "500 m"))

aoi_masked <- function(r) { v <- terra::vect(sf::st_transform(AOI.Sf, terra::crs(r))); terra::mask(terra::crop(r, v), v) }
Cond <- dplyr::bind_rows(lapply(names(Metrics), function(m) {
  r <- Metrics[[m]]$r
  z <- zonal_stat(r, Sets.Sf[, "set"], "mean")[[1]]
  names(z) <- NULL
  ram <- aoi_masked(r)
  tibble::tibble(metric = m, resolution = Metrics[[m]]$res, resolves_strip = Metrics[[m]]$res %in% c("10 m", "30 m"), set = Sets.Sf$set, mean = z,
                 ranch_median = global_q(ram, 0.5), ranch_mean = unname(terra::global(ram, "mean", na.rm = TRUE)[[1]]),
                 ranch_sd = unname(terra::global(ram, "sd", na.rm = TRUE)[[1]]))
}))
Cond <- Cond |> dplyr::group_by(metric) |>
  dplyr::mutate(interior_mean = mean[set == "interior"],
                diff_from_median = mean - ranch_median, diff_from_interior = mean - interior_mean,
                std_diff_vs_interior = (mean - interior_mean) / ranch_sd) |> dplyr::ungroup()
readr::write_csv(Cond, file.path(Dir.Out.Tables, "riparian_condition.csv"))
Show <- Cond |> dplyr::filter(set %in% c("order2plus", "order2", "order3")) |>
  dplyr::select(metric, resolution, set, mean, interior_mean, diff_from_interior, std_diff_vs_interior) |>
  dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3)))
print(as.data.frame(Show), row.names = FALSE)


# PHASE 3: Gradient with distance from streams of order >= 2
log_phase("PHASE 3: distance gradient")

Net2 <- sf::st_union(sf::st_geometry(Streams[Streams$order >= 2, ]))
Pts  <- sf::st_sample(AOI.U, size = 12000, type = "regular")
Pts  <- Pts[!sf::st_is_empty(Pts)]
Dist <- as.numeric(sf::st_distance(Pts, Net2))
Bins <- c(0, 60, 125, 250, 500, 1000, Inf); Bin.Lab <- c("0-60 m", "60-125 m", "125-250 m", "250-500 m", "500-1000 m", ">1000 m")
Grad.Vars <- c("ndvi_dry", "ndvi_lr", "bsi_dry", "ndwi_dry", "wc_trees", "wc_shrub", "wc_grass", "twi")
Vals <- vapply(Grad.Vars, function(m) { r <- Metrics[[m]]$r; terra::extract(r, terra::vect(sf::st_transform(sf::st_sf(geometry = Pts), terra::crs(r))))[[2]] }, numeric(length(Pts)))
Grad <- tibble::as_tibble(Vals) |> dplyr::mutate(dist_class = cut(Dist, Bins, labels = Bin.Lab, include.lowest = TRUE)) |>
  dplyr::group_by(dist_class) |> dplyr::summarise(n_points = dplyr::n(), dplyr::across(dplyr::all_of(Grad.Vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
Grad.Cor <- tibble::tibble(metric = Grad.Vars, spearman_with_distance = vapply(Grad.Vars, function(m) suppressWarnings(cor(Dist, Vals[, m], method = "spearman", use = "complete.obs")), 0))
readr::write_csv(Grad, file.path(Dir.Out.Tables, "riparian_gradient.csv"))
readr::write_csv(Grad.Cor, file.path(Dir.Out.Tables, "riparian_gradient_correlation.csv"))
print(as.data.frame(Grad |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)
print(as.data.frame(Grad.Cor |> dplyr::mutate(spearman_with_distance = round(spearman_with_distance, 3))), row.names = FALSE)


# PHASE 4: Mask raster
Dem <- terra::rast(file.path(Dir.DEM, "amaka_dem_utm37s.tif"))
Mask <- terra::rasterize(terra::vect(sf::st_sf(geometry = sf::st_sfc(sf::st_union(Sets$order2plus), crs = Cfg.CRS.Metric))), terra::crop(Dem, terra::vect(AOI.Sf)), field = 1, background = 0)
Mask <- terra::mask(Mask, terra::vect(AOI.Sf)); names(Mask) <- "riparian_order2plus"
terra::writeRaster(Mask, file.path(Dir.Out.Rasters, "riparian_mask.tif"), overwrite = TRUE, datatype = "INT1U")

K <- Cond[Cond$set == "order2plus" & Cond$resolves_strip, ]
message(sprintf("  riparian (order >= 2, %g m) covers %.1f%% of the AOI. Versus the interior: dry-season NDVI %+.3f, BSI %+.3f, tree fraction %+.4f, TWI %+.2f",
                Cfg.Riparian.m, Sets.Tbl$pct_of_aoi[Sets.Tbl$set == "order2plus"],
                K$diff_from_interior[K$metric == "ndvi_dry"], K$diff_from_interior[K$metric == "bsi_dry"],
                K$diff_from_interior[K$metric == "wc_trees"], K$diff_from_interior[K$metric == "twi"]))

log_phase("DONE")
