# 09_fire_history.R  -  Fire history and fire x woody-loss separation
# PURPOSE
#   Burn count and most recent burn year on the 500 m MCD64A1 grid (2001-2024), the year-by-year fire record (burned
#   area, MODIS active-fire detections) next to Hansen tree-cover loss, and the cross-tab of burned-and-lost versus
#   unburned-and-lost woody loss that separates fire-associated loss from the rest (DESCRIPTION.md 6.7). Jones 2025: fire and
#   megaherbivores act synergistically on tree mortality, so the two signatures have to be kept apart.
#
# INPUTS
#   data/gee_exports/{burn_count,most_recent_burn_year,burned_by_year,firms_detections_by_year}_2001_2024.tif   (gee/10)
#   data/gee_exports/hansen_gfc_v1_11_envelope.tif                                                              (gee/08)
#   outputs/rasters/pan_catchments.shp (08; optional - the per-pan rows are skipped without it)
#
# OUTPUTS
#   outputs/rasters/burn_count_2001_2024.tif, most_recent_burn_year.tif (NA = never burned), firms_detections_total.tif
#   outputs/tables/fire_by_year.csv            burned area, active-fire detections and Hansen loss per year (envelope and AOI)
#   outputs/tables/fire_woody_loss_crosstab.csv burned-and-lost vs unburned-and-lost, envelope / AOI / each pan catchment
#   outputs/tables/pan_fire_metrics.csv        fire metrics per pan catchment
#   outputs/tables/cells_fire.csv              per analysis-cell fire and Hansen woody-cover metrics (feeds 10 and 11)
#   outputs/figures/09_fire_by_year.png, 09_fire_map.png
#
# RUNTIME NOTES
#   * A few seconds. Burned area is a 500 m product and Hansen loss a 30 m one: each 30 m pixel is assigned the burn history
#     of the 500 m cell containing it (point-sampling a coarse layer at fine pixels, DESCRIPTION.md 9.2). Fire anywhere in a cell
#     therefore tags every loss pixel in it, which overstates "burned-and-lost" and cannot see small fires.
#   * "burned-and-lost" = a loss pixel whose cell burned in the loss year or up to Cfg.Fire.Loss.Window.Years before it.
#     Only mature-tree pixels (Hansen treecover2000 >= Cfg.Tree.Mature.Pct) are counted as woody loss.
#   * The unburned remainder is "not fire", NOT "elephant": charcoal burning and clearing, drought mortality and mapping noise
#     look the same from orbit. Cause can only be attributed with ground evidence.
#   * To tell fire signal from coincidence, the table gives (i) the share of mature-tree pixel-years exposed to fire = the share
#     of "burned-and-lost" expected by chance, and (ii) the annual loss rate per exposed vs unexposed pixel-year with an exact
#     Poisson rate ratio (exposure_*). The spatial_* columns (loss in cells that ever burned vs never burned) IGNORE timing and
#     are shown only for contrast: they count losses that happened before the fire.
#   * Hansen v1.11 loss in this envelope stops in 2014 and MCD64A1 shows no burns in 2009-2020: a short overlap, so every
#     conclusion here is about 2001-2014.

if (!exists("Setup.Loaded")) source("R/00_setup.R")

# PHASE 1: Inputs
log_phase("PHASE 1: inputs")
Aoi.Sf  <- read_aoi(); Env.Sf <- read_envelope()
Burn.Count <- terra::rast(gee_file("burn_count_2001_2024"))                 # 500 m, times burned 2001-2024
Burn.Last  <- terra::rast(gee_file("most_recent_burn_year_2001_2024"))      # 0 = never burned
Burn.Year  <- terra::rast(gee_file("burned_by_year_2001_2024"))             # 24 layers burned_<year>, 0 / 1
Firms      <- terra::rast(gee_file("firms_detections_by_year_2001_2024"))   # 1 km, MODIS active-fire detections per year
Hansen     <- terra::rast(gee_file("hansen_gfc_v1_11_envelope"))            # 30 m: treecover2000, loss, lossyear, gain, datamask
Yrs <- layer_year(Burn.Year)
stopifnot(identical(Yrs, layer_year(Firms)), identical(Yrs, as.integer(Cfg.Years.Fire)))
Cell500.ha <- prod(terra::res(Burn.Year)) / 1e4; Px30.ha <- prod(terra::res(Hansen)) / 1e4

# the two burn products must agree with each other
BY <- terra::values(Burn.Year)
qc.count <- all(rowSums(BY, na.rm = TRUE) == terra::values(Burn.Count)[, 1])
qc.last  <- all(apply(BY, 1, function(v) if (any(v == 1, na.rm = TRUE)) max(Yrs[which(v == 1)]) else 0) == terra::values(Burn.Last)[, 1])
message(sprintf("  burned_by_year sums to burn_count: %s; its last burned year equals most_recent_burn_year: %s", qc.count, qc.last))
message(sprintf("  years in the exports: %d-%d; Hansen loss years in the envelope: %s",
                min(Yrs), max(Yrs), paste(range(2000 + terra::values(Hansen$lossyear)[terra::values(Hansen$loss)[, 1] %in% 1, 1]), collapse = "-")))


# PHASE 2: Burn rasters (500 m)
log_phase("PHASE 2: burn rasters")
terra::writeRaster(Burn.Count, file.path(Dir.Out.Rasters, "burn_count_2001_2024.tif"), overwrite = TRUE, datatype = "INT1U")
Last.Na <- terra::ifel(Burn.Last > 0, Burn.Last, NA); names(Last.Na) <- "most_recent_burn_year"
terra::writeRaster(Last.Na, file.path(Dir.Out.Rasters, "most_recent_burn_year.tif"), overwrite = TRUE, datatype = "INT2U")
Firms.Total <- sum(Firms); names(Firms.Total) <- "firms_detections_total"
terra::writeRaster(Firms.Total, file.path(Dir.Out.Rasters, "firms_detections_total.tif"), overwrite = TRUE, datatype = "INT2U")
message(sprintf("  envelope: %d of %d cells burned at least once (%.0f%%), up to %d times",
                sum(terra::values(Burn.Count)[, 1] > 0), terra::ncell(Burn.Count), 100 * mean(terra::values(Burn.Count)[, 1] > 0), max(terra::values(Burn.Count)[, 1])))


# PHASE 3: Year-by-year record
log_phase("PHASE 3: year-by-year fire and Hansen loss")
Tc <- terra::values(Hansen$treecover2000)[, 1]; Loss1 <- terra::values(Hansen$loss)[, 1] %in% 1; Ly <- terra::values(Hansen$lossyear)[, 1]; Ly[is.na(Ly)] <- 0
mask30 <- function(poly) !is.na(terra::values(terra::rasterize(terra::vect(sf::st_transform(poly, terra::crs(Hansen))), Hansen, field = 1))[, 1])
Scopes <- list(envelope = mask30(Env.Sf), AOI = mask30(Aoi.Sf))
Polys  <- list(envelope = Env.Sf, AOI = Aoi.Sf)

year_table <- function(label) {
  p <- sf::st_transform(Polys[[label]], terra::crs(Burn.Year)); m <- Scopes[[label]]
  b <- unlist(exactextractr::exact_extract(Burn.Year, p, fun = "sum", progress = FALSE)[1, ]) * Cell500.ha       # coverage-weighted ha
  f <- unlist(exactextractr::exact_extract(Firms, p, fun = "sum", progress = FALSE)[1, ])                         # coverage-weighted detections
  tibble::tibble(scope = label, year = Yrs, burned_ha = unname(b), firms_detections = unname(f),
                 hansen_loss_ha = vapply(Yrs, function(y) sum(m & Loss1 & Ly == y - 2000) * Px30.ha, 0),
                 hansen_loss_canopy_ha = vapply(Yrs, function(y) sum(Tc[m & Loss1 & Ly == y - 2000] / 100) * Px30.ha, 0))
}
Fire.By.Year <- dplyr::bind_rows(lapply(names(Scopes), year_table))
readr::write_csv(Fire.By.Year, file.path(Dir.Out.Tables, "fire_by_year.csv"))
Tot <- Fire.By.Year |> dplyr::group_by(scope) |> dplyr::summarise(burn_years = sum(burned_ha > 0), burned_ha = sum(burned_ha), firms_years = sum(firms_detections > 0.5), firms = sum(firms_detections),
                                                                 loss_ha = sum(hansen_loss_ha), loss_canopy_ha = sum(hansen_loss_canopy_ha), .groups = "drop")
print(as.data.frame(Tot |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)
Aoi.Yr <- Fire.By.Year[Fire.By.Year$scope == "AOI", ]
message(sprintf("  AOI: burned in %s; active-fire detections in %s; Hansen loss in %s",
                paste(Aoi.Yr$year[Aoi.Yr$burned_ha > 0], collapse = ", "), paste(Aoi.Yr$year[Aoi.Yr$firms_detections > 0.5], collapse = ", "), paste(Aoi.Yr$year[Aoi.Yr$hansen_loss_ha > 0], collapse = ", ")))


# PHASE 4: Burned-and-lost versus unburned-and-lost
log_phase("PHASE 4: fire x woody-loss cross-tab")
Xy <- terra::xyFromCell(Hansen, seq_len(terra::ncell(Hansen)))
Cb <- terra::cellFromXY(Burn.Count, Xy)                                      # the 500 m cell of every 30 m pixel
Ever <- terra::values(Burn.Count)[Cb, 1] > 0; Ever[is.na(Ever)] <- FALSE
Win <- rep(FALSE, length(Cb)); Li <- which(Loss1 & !is.na(Cb))               # burned in the loss year or W years before it (loss pixels only)
for (k in 0:Cfg.Fire.Loss.Window.Years) {
  col <- (2000 + Ly[Li] - k) - 2000; ok <- col >= 1 & col <= ncol(BY)
  Win[Li[ok]] <- Win[Li[ok]] | (BY[cbind(Cb[Li[ok]], col[ok])] %in% 1)
}

# Timing-aware test: every mature-tree pixel contributes one pixel-year for each year it is still standing (2001 to its loss year,
# or to 2023 = the end of Hansen v1.11). A pixel-year is EXPOSED when its cell burned in that year or up to W years before. The loss
# rate per exposed pixel-year is compared with the unexposed rate (exact Poisson comparison). Unlike a plain "ever burned" split this
# cannot be fooled by losses that happened before the fire.
exposure_test <- function(risk) {
  idx <- which(risk & !is.na(Cb)); cb <- Cb[idx]; ly <- ifelse(Loss1[idx], Ly[idx], Inf)
  py <- c(0, 0); ev <- c(0, 0)
  for (t in seq_len(2023 - 2000)) {                                         # t = year - 2000 = the column of burned_by_year
    live <- ly >= t                                                         # standing at the start of year t (the loss year itself counts)
    expo <- rep(FALSE, length(idx))
    for (k in 0:Cfg.Fire.Loss.Window.Years) if (t - k >= 1) expo <- expo | (BY[cbind(cb, t - k)] %in% 1)
    py <- py + c(sum(live & expo), sum(live & !expo)); ev <- ev + c(sum(live & expo & ly == t), sum(live & !expo & ly == t))
  }
  pt <- if (all(py > 0) && sum(ev) > 0) stats::poisson.test(ev, py) else NULL
  list(py = py, ev = ev, rr = if (is.null(pt)) NA_real_ else unname(pt$estimate), lo = if (is.null(pt)) NA_real_ else pt$conf.int[1],
       hi = if (is.null(pt)) NA_real_ else pt$conf.int[2], p = if (is.null(pt)) NA_real_ else pt$p.value)
}

crosstab <- function(msk, label) {
  risk <- msk & Tc >= Cfg.Tree.Mature.Pct                                   # mature-tree pixels: the population that can lose woody cover
  lost <- risk & Loss1
  a <- sum(lost & Win); b <- sum(lost & !Win)                               # burned-and-lost / unburned-and-lost
  n11 <- sum(risk & Ever & Loss1); n10 <- sum(risk & Ever & !Loss1); n01 <- sum(risk & !Ever & Loss1); n00 <- sum(risk & !Ever & !Loss1)
  r1 <- if (n11 + n10 > 0) n11 / (n11 + n10) else NA_real_; r0 <- if (n01 + n00 > 0) n01 / (n01 + n00) else NA_real_
  p  <- if ((n11 + n10) > 0 && (n01 + n00) > 0 && (n11 + n01) > 0) stats::fisher.test(matrix(c(n11, n10, n01, n00), 2))$p.value else NA_real_
  ex <- exposure_test(risk)
  stopifnot(ex$ev[1] == a)                                                  # the exposure test and the cross-tab use the same definition of "burned-and-lost"
  tibble::tibble(scope = label,
                 mature_tree_px = sum(risk), mature_tree_ha = sum(risk) * Px30.ha,
                 loss_px_all = sum(msk & Loss1), loss_ha_all = sum(msk & Loss1) * Px30.ha, loss_canopy_ha_all = sum(Tc[msk & Loss1] / 100) * Px30.ha,
                 woody_loss_px = sum(lost), woody_loss_ha = sum(lost) * Px30.ha,
                 burned_and_lost_px = a, unburned_and_lost_px = b, burned_and_lost_pct = ifelse(a + b > 0, 100 * a / (a + b), NA_real_),
                 exposed_pixel_year_share_pct = 100 * ex$py[1] / sum(ex$py),            # share of burned-and-lost expected by chance
                 annual_loss_rate_exposed_pct = ifelse(ex$py[1] > 0, 100 * ex$ev[1] / ex$py[1], NA_real_),
                 annual_loss_rate_unexposed_pct = ifelse(ex$py[2] > 0, 100 * ex$ev[2] / ex$py[2], NA_real_),
                 exposure_rate_ratio = ex$rr, exposure_rr_ci_low = ex$lo, exposure_rr_ci_high = ex$hi, exposure_p = ex$p,
                 spatial_share_in_ever_burned_cells_pct = 100 * mean(Ever[risk]),    # timing IGNORED: cumulative loss in cells that ever burned vs never burned
                 spatial_loss_in_ever_burned_pct = 100 * r1, spatial_loss_in_never_burned_pct = 100 * r0,
                 spatial_risk_ratio = ifelse(is.finite(r0) && r0 > 0, r1 / r0, NA_real_), spatial_fisher_p = p)
}
Cross <- dplyr::bind_rows(lapply(names(Scopes), function(n) crosstab(Scopes[[n]], n)))

# per-pan rows (stream-fed pans: their catchment; others: the context disc)
Pc.f <- file.path(Dir.Out.Rasters, "pan_catchments.shp"); Pc <- NULL
if (file.exists(Pc.f)) {
  Pc <- sf::st_read(Pc.f, quiet = TRUE)
  Cross <- dplyr::bind_rows(Cross, dplyr::bind_rows(lapply(seq_len(nrow(Pc)), function(i) crosstab(mask30(Pc[i, ]), paste0("pan ", Pc$pan_id[i])))))
} else message("  pan_catchments.shp not found: run R/08_water_pans_analysis.R for the per-pan rows")
readr::write_csv(Cross, file.path(Dir.Out.Tables, "fire_woody_loss_crosstab.csv"))
print(as.data.frame(Cross |> dplyr::select(scope, mature_tree_ha, woody_loss_px, burned_and_lost_px, unburned_and_lost_px, burned_and_lost_pct, exposed_pixel_year_share_pct,
                                           annual_loss_rate_exposed_pct, annual_loss_rate_unexposed_pct, exposure_rate_ratio, exposure_rr_ci_low, exposure_rr_ci_high, exposure_p) |>
                      dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)


# PHASE 5: Per-cell table (feeds 10 and 11)
log_phase("PHASE 5: cell table")
Grid <- get_analysis_grid(); Cells <- Grid$cells; Cell.Poly <- cell_polygons(Grid$grid)
Burn.Samp <- sample_cells(c(Burn.Count, Last.Na), Cells)                    # 500 m layers, point-sampled at the cell centroids
names(Burn.Samp) <- c("burn_count", "last_burn_year")
H.Layers <- c(Hansen$treecover2000, Hansen$treecover2000 >= Cfg.Tree.Mature.Pct, Hansen$loss, Hansen$loss * Hansen$treecover2000 / 100)
names(H.Layers) <- c("tc2000_mean", "mature_tree_share", "hansen_loss_share", "hansen_loss_canopy_share")
Hs <- zonal_stat(H.Layers, Cell.Poly, "mean")                               # 30 m layers, area-weighted over the cell polygons
Cells.Fire <- dplyr::bind_cols(tibble::tibble(cell_id = Cell.Poly$cell_id), Burn.Samp, sample_cells(Firms.Total, Cells) |> stats::setNames("firms_detections_total"), Hs)
Cells.Fire$last_burn_year[is.na(Cells.Fire$last_burn_year)] <- NA
Cells.Fire$years_since_burn <- max(Yrs) - Cells.Fire$last_burn_year
Cells.Fire$hansen_loss_ha <- Hs$hansen_loss_share * Cells$area_ha[match(Cells.Fire$cell_id, Cells$cell_id)]
Cells.Fire$hansen_loss_canopy_ha <- Hs$hansen_loss_canopy_share * Cells$area_ha[match(Cells.Fire$cell_id, Cells$cell_id)]
readr::write_csv(Cells.Fire, file.path(Dir.Out.Tables, "cells_fire.csv"))
In.Aoi <- Cells$in_aoi[match(Cells.Fire$cell_id, Cells$cell_id)]
message(sprintf("  %d cells; inside the AOI: %d burned at least once, %d never; mean Hansen tree cover 2000 %.1f%%, %.1f%% of the area holds mature trees",
                nrow(Cells.Fire), sum(Cells.Fire$burn_count[In.Aoi] > 0, na.rm = TRUE), sum(Cells.Fire$burn_count[In.Aoi] == 0, na.rm = TRUE),
                mean(Cells.Fire$tc2000_mean[In.Aoi], na.rm = TRUE), 100 * mean(Cells.Fire$mature_tree_share[In.Aoi], na.rm = TRUE)))


# PHASE 6: Per-pan fire metrics
log_phase("PHASE 6: per-pan fire metrics")
if (!is.null(Pc)) {
  Ever.R <- Burn.Count > 0
  Pan.Fire <- dplyr::bind_rows(lapply(seq_len(nrow(Pc)), function(i) {
    p <- sf::st_transform(Pc[i, ], terra::crs(Burn.Count))
    tibble::tibble(pan_id = Pc$pan_id[i], area_basis = Pc$method[i],
                   burned_share_ever_pct = 100 * exactextractr::exact_extract(Ever.R, p, fun = "mean", progress = FALSE),
                   mean_burn_count = exactextractr::exact_extract(Burn.Count, p, fun = "mean", progress = FALSE),
                   last_burn_year = exactextractr::exact_extract(Last.Na, p, fun = "max", progress = FALSE),
                   firms_detections_per_100ha = 100 * exactextractr::exact_extract(Firms.Total, p, fun = "mean", progress = FALSE) / (prod(terra::res(Firms.Total)) / 1e4))
  })) |> dplyr::left_join(Cross |> dplyr::filter(grepl("^pan ", scope)) |> dplyr::mutate(pan_id = sub("^pan ", "", scope)) |>
                            dplyr::select(pan_id, mature_tree_ha, loss_ha_all, woody_loss_px, burned_and_lost_px, unburned_and_lost_px, burned_and_lost_pct), by = "pan_id")
  readr::write_csv(Pan.Fire, file.path(Dir.Out.Tables, "pan_fire_metrics.csv"))
  print(as.data.frame(Pan.Fire |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
} else message("  skipped (no pan_catchments.shp)")

# PHASE 7: Check figures
log_phase("PHASE 7: check figures")
Yr.Long <- Fire.By.Year |>
  tidyr::pivot_longer(c(burned_ha, firms_detections, hansen_loss_ha), names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric = factor(metric, c("burned_ha", "firms_detections", "hansen_loss_ha"),
                                c("Burned area, MCD64A1 (ha)", "Active-fire detections, MODIS (count)", "Hansen tree-cover loss (ha)")),
                scope = factor(scope, c("AOI", "envelope"), c("Amaka AOI (4,782 ha)", "7 km envelope")))
Fig.Fire <- ggplot2::ggplot(Yr.Long, ggplot2::aes(year, value, fill = metric)) +
  ggplot2::geom_col(show.legend = FALSE) + ggplot2::facet_grid(metric ~ scope, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c("#d95f02", "#b2182b", "#1b7837")) + ggplot2::scale_x_continuous(breaks = seq(2001, 2024, 3)) +
  ggplot2::labs(x = NULL, y = NULL, title = "Fire and woody loss on and around Amaka, 2001-2024",
                subtitle = "Hansen loss ends in 2014; no burned area is mapped inside the envelope after 2008, though active-fire detections continue (2020, 2021, 2024)") +
  ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))
ggplot2::ggsave(file.path(Dir.Out.Figures, "09_fire_by_year.png"), Fig.Fire, width = 9, height = 6.5, dpi = 170)

Hs.R <- terra::rast(file.path(Dir.Out.Hydro.Buf, "hillshade.tif")); names(Hs.R) <- "shade"
Burn.Plot <- terra::ifel(Burn.Count > 0, Burn.Count, NA); names(Burn.Plot) <- "burn_count"
Loss.Pts <- sf::st_as_sf(data.frame(terra::xyFromCell(Hansen, which(Loss1 & Tc >= Cfg.Tree.Mature.Pct))), coords = c("x", "y"), crs = terra::crs(Hansen))
Map.Fire <-
  tmap::tm_shape(Hs.R) + tmap::tm_raster(col.scale = tmap::tm_scale_continuous(values = "brewer.greys"), col.legend = tmap::tm_legend_hide()) +
  tmap::tm_shape(Burn.Plot) + tmap::tm_raster(col.scale = tmap::tm_scale_categorical(values = c("#fdae61", "#d73027")), col_alpha = 0.7, col.legend = tmap::tm_legend(title = "Times burned")) +
  tmap::tm_shape(Loss.Pts) + tmap::tm_dots(fill = "#1b7837", size = 0.25, col = NULL) +
  tmap::tm_shape(sf::st_geometry(Aoi.Sf)) + tmap::tm_borders(col = "black", lwd = 2.5) +
  tmap::tm_title("Burn count 2001-2024 (500 m) and mature-tree loss pixels (green dots, 30 m)", size = 0.8)
tmap::tmap_save(Map.Fire, file.path(Dir.Out.Figures, "09_fire_map.png"), width = 8, height = 8, dpi = 170)
message("  wrote 09_fire_by_year.png and 09_fire_map.png")

log_phase("DONE")
