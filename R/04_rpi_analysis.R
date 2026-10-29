# 04_rpi_analysis.R  -  Resilience / performance (RPI) analysis inside Amaka
# PURPOSE
#   Zonal statistics, spatial pattern and pixel-wise Theil-Sen trends of RPI (and of actual and
#   potential GPP) inside Amaka; flag "degrading" pixels; stratify by TWI quintile; cross-check
#   against the independent MODIS vegetation indices (DESCRIPTION.md 3.2, 6.4).
#
# INPUTS
#   Path.RPI, Path.RPI.GPP.Act, Path.RPI.GPP.Pot, Path.RPI.TempPerf   (native grid, 2000-2022)
#   outputs/tables/cells_terrain.csv (03);  data/gee_exports/modis_{ndvi,evi}_mean_ANNUAL_*.tif
#
# OUTPUTS
#   outputs/tables/rpi_zonal_annual.csv        AOI mean / SD / p10 / p50 / p90 per year (RPI, GPP actual, GPP potential)
#   outputs/tables/rpi_trend_summary.csv       AOI-level Theil-Sen + Mann-Kendall, share of degrading cells
#   outputs/tables/rpi_by_twi_quintile.csv     RPI now / trend by TWI quintile (AOI cells)
#   outputs/tables/crosscheck_modis_vs_rpi.csv RPI / GPP versus MODIS NDVI / EVI (time series and spatial)
#   outputs/tables/cells_rpi.csv               one row per analysis cell (keyed on cell_id)
#   outputs/rasters/rpi_now.tif, rpi_sen_slope.tif, rpi_sen_p.tif, rpi_degrading.tif, gpp_sen_slope.tif
#
# RUNTIME NOTES
#   * Seconds. RPI v2 has 23 layers, 2000-2022 (not 2000-2023 as in the brief).
#   * Everything stays on the RPI native grid (the target grid of DESCRIPTION.md 9.2); no resampling.
#   * "Degrading" = Theil-Sen slope < 0 with Mann-Kendall p < Cfg.Degrade.P (0.1), per cell.
#   * RPI_now (used by the stress index) = mean of Cfg.RPI.Now.Years (2020-2022) to damp single-year noise.
#   * Land-cover stratification is done in 06 (it needs the 10 m land cover); TWI quintiles are done here.
#   * RPI is a performance ratio relative to modelled potential GPP: it says how far a pixel is from what
#     its climate and soil could support, not what caused the gap (9.9: a structural index).

if (!exists("Setup.Loaded")) source("R/00_setup.R")


# PHASE 1: Inputs on the analysis grid
log_phase("PHASE 1: inputs")

Grid.List <- get_analysis_grid(); G <- Grid.List$grid; Cells <- Grid.List$cells
AOI.Sf    <- read_aoi()
Terrain   <- readr::read_csv(file.path(Dir.Out.Tables, "cells_terrain.csv"), show_col_types = FALSE)

crop_native <- function(path) terra::crop(terra::rast(path), terra::ext(G))      # same grid as G by construction
Rpi     <- crop_native(Path.RPI)
Gpp.Act <- crop_native(Path.RPI.GPP.Act)
Gpp.Pot <- crop_native(Path.RPI.GPP.Pot)
Perf    <- crop_native(Path.RPI.TempPerf)
stopifnot("RPI crop is not on the analysis grid" = terra::compareGeom(Rpi[[1]], G, stopOnError = FALSE))
Years <- layer_year(Rpi)
message(sprintf("  %d RPI layers, %d-%d, %d x %d cells; %d cells inside the AOI", terra::nlyr(Rpi), min(Years), max(Years),
                terra::ncol(G), terra::nrow(G), sum(Cells$in_aoi)))


# PHASE 2: Zonal statistics inside Amaka, per year 
log_phase("PHASE 2: AOI statistics per year")

zonal_multi <- function(r, label) {
  ex <- exactextractr::exact_extract(r, sf::st_transform(AOI.Sf, terra::crs(r)), fun = c("mean", "stdev", "quantile"),
                                     quantiles = c(0.1, 0.5, 0.9), progress = FALSE)
  long <- tidyr::pivot_longer(ex, tidyselect::everything(), names_to = c("stat", "layer"),
                              names_pattern = "^(mean|stdev|q10|q50|q90)\\.(.*)$")
  long$year <- layer_year(long$layer)
  out <- tidyr::pivot_wider(long[, c("year", "stat", "value")], names_from = stat, values_from = value)
  dplyr::transmute(out, variable = label, year, mean, sd = stdev, p10 = q10, p50 = q50, p90 = q90)
}
Zonal.Annual <- dplyr::bind_rows(zonal_multi(Rpi, "rpi"), zonal_multi(Gpp.Act, "gpp_actual"), zonal_multi(Gpp.Pot, "gpp_potential")) |>
  dplyr::arrange(variable, year)
stopifnot("pivot_wider produced list columns (DESCRIPTION.md 9.5)" = !any(vapply(Zonal.Annual, is.list, logical(1))))
readr::write_csv(Zonal.Annual, file.path(Dir.Out.Tables, "rpi_zonal_annual.csv"))
print(as.data.frame(Zonal.Annual |> dplyr::filter(variable == "rpi", year %in% c(2000, 2005, 2010, 2015, 2020, 2022)) |>
                      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)


# PHASE 3: Pixel-wise Theil-Sen trends
# trend::sens.slope() per cell through terra::app(); slope = change per year in the layer's own units.
log_phase("PHASE 3: Theil-Sen trends per cell")

sen_fun <- function(v) {
  ok <- !is.na(v)
  if (sum(ok) < 10) return(c(NA_real_, NA_real_))
  s <- tryCatch(trend::sens.slope(v[ok]), error = function(e) NULL)
  if (is.null(s)) return(c(NA_real_, NA_real_))
  c(unname(s$estimates), s$p.value)
}
Sen.Rpi <- terra::app(Rpi, sen_fun);     names(Sen.Rpi) <- c("rpi_sen_slope", "rpi_sen_p")
Sen.Gpp <- terra::app(Gpp.Act, sen_fun); names(Sen.Gpp) <- c("gpp_sen_slope", "gpp_sen_p")
Degrading <- (Sen.Rpi[["rpi_sen_slope"]] < 0) & (Sen.Rpi[["rpi_sen_p"]] < Cfg.Degrade.P); names(Degrading) <- "rpi_degrading"

Rpi.Now    <- terra::mean(Rpi[[Years %in% Cfg.RPI.Now.Years]], na.rm = TRUE);  names(Rpi.Now) <- "rpi_now"
Rpi.First  <- terra::mean(Rpi[[Years %in% min(Years):(min(Years) + 2)]], na.rm = TRUE)
Rpi.Mean   <- terra::mean(Rpi, na.rm = TRUE);                                   names(Rpi.Mean) <- "rpi_mean"
Gpp.Now    <- terra::mean(Gpp.Act[[layer_year(Gpp.Act) %in% Cfg.RPI.Now.Years]], na.rm = TRUE); names(Gpp.Now) <- "gpp_now"
Gpp.Pot.Now <- terra::mean(Gpp.Pot[[layer_year(Gpp.Pot) %in% Cfg.RPI.Now.Years]], na.rm = TRUE); names(Gpp.Pot.Now) <- "gpp_pot_now"
Rpi.Change <- Rpi.Now - Rpi.First;                                              names(Rpi.Change) <- "rpi_change_now_vs_first3"

All.R <- c(Rpi.Now, Rpi.Mean, Rpi.Change, Sen.Rpi, Degrading, Gpp.Now, Gpp.Pot.Now, Sen.Gpp, Perf)
Cells.Rpi <- dplyr::bind_cols(tibble::tibble(cell_id = seq_len(terra::ncell(G))),
                              tibble::as_tibble(terra::values(All.R, dataframe = TRUE)))
Cells.Rpi$rpi_degrading <- as.integer(Cells.Rpi$rpi_degrading)

# Quintiles are defined on the AOI cells so they describe Amaka itself.
Joined <- dplyr::left_join(Cells, Terrain, by = "cell_id") |> dplyr::left_join(Cells.Rpi, by = "cell_id")
quint <- function(x, ref) cut(x, quantile(ref, 0:5 / 5, na.rm = TRUE), include.lowest = TRUE, labels = 1:5)
Cells.Rpi$rpi_now_q <- as.integer(as.character(quint(Joined$rpi_now, Joined$rpi_now[Joined$in_aoi])))
Cells.Rpi$twi_q     <- as.integer(as.character(quint(Joined$twi_mean, Joined$twi_mean[Joined$in_aoi])))
Cells.Rpi$rpi_now_q[!Joined$in_aoi] <- NA; Cells.Rpi$twi_q[!Joined$in_aoi] <- NA
readr::write_csv(Cells.Rpi, file.path(Dir.Out.Tables, "cells_rpi.csv"))

In.Aoi <- Joined$in_aoi
message(sprintf("  AOI cells: RPI_now mean %.3f (range %.3f-%.3f) | Sen slope mean %+.4f per year | %.0f%% degrading (p < %.2f), %.0f%% improving",
                mean(Joined$rpi_now[In.Aoi], na.rm = TRUE), min(Joined$rpi_now[In.Aoi], na.rm = TRUE), max(Joined$rpi_now[In.Aoi], na.rm = TRUE),
                mean(Joined$rpi_sen_slope[In.Aoi], na.rm = TRUE), 100 * mean(Joined$rpi_degrading[In.Aoi] == 1, na.rm = TRUE), Cfg.Degrade.P,
                100 * mean(Joined$rpi_sen_slope[In.Aoi] > 0 & Joined$rpi_sen_p[In.Aoi] < Cfg.Degrade.P, na.rm = TRUE)))


# PHASE 4: AOI-level trend summary
log_phase("PHASE 4: AOI-level trends")

series_trend <- function(df, var) {
  d <- df[df$variable == var, ]; d <- d[order(d$year), ]
  s <- trend::sens.slope(d$mean); m <- trend::mk.test(d$mean)
  first3 <- mean(head(d$mean, 3)); last3 <- mean(tail(d$mean, 3))
  tibble::tibble(variable = var, years = paste(range(d$year), collapse = "-"), mean_all = mean(d$mean),
                 mean_first3 = first3, mean_last3 = last3, change_pct = 100 * (last3 / first3 - 1),
                 sen_slope_per_year = unname(s$estimates), sen_slope_pct_per_year = 100 * unname(s$estimates) / mean(d$mean),
                 mk_p = m$p.value, mk_tau = unname(m$estimates["tau"]))
}
Trend.Summary <- dplyr::bind_rows(lapply(c("rpi", "gpp_actual", "gpp_potential"), series_trend, df = Zonal.Annual))
w <- Joined$aoi_frac[In.Aoi]
Trend.Summary$aoi_cells_degrading_pct <- c(100 * weighted.mean(Joined$rpi_degrading[In.Aoi] == 1, w, na.rm = TRUE), NA, NA)
Trend.Summary$aoi_cells_improving_pct <- c(100 * weighted.mean(Joined$rpi_sen_slope[In.Aoi] > 0 & Joined$rpi_sen_p[In.Aoi] < Cfg.Degrade.P, w, na.rm = TRUE), NA, NA)
readr::write_csv(Trend.Summary, file.path(Dir.Out.Tables, "rpi_trend_summary.csv"))
print(as.data.frame(Trend.Summary |> dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

# Reliability: how well does the model of potential GPP fit here (the RPI denominator)?
message(sprintf("  potential-GPP model fit inside the AOI: rsq mean %.2f, rmse mean %.3f", mean(Joined$rsq[In.Aoi], na.rm = TRUE), mean(Joined$rmse[In.Aoi], na.rm = TRUE)))


# PHASE 5: Stratify by TWI quintile
log_phase("PHASE 5: RPI by TWI quintile")

By.Twi <- Joined[In.Aoi, ] |> dplyr::mutate(twi_q = Cells.Rpi$twi_q[In.Aoi]) |>
  dplyr::group_by(twi_quintile = twi_q) |>
  dplyr::summarise(n_cells = dplyr::n(), twi_mean = mean(twi_mean), rpi_now = mean(rpi_now, na.rm = TRUE),
                   rpi_sen_slope = mean(rpi_sen_slope, na.rm = TRUE), degrading_pct = 100 * mean(rpi_degrading == 1, na.rm = TRUE),
                   gpp_now = mean(gpp_now, na.rm = TRUE), .groups = "drop")
readr::write_csv(By.Twi, file.path(Dir.Out.Tables, "rpi_by_twi_quintile.csv"))
print(as.data.frame(By.Twi |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)


# PHASE 6: Cross-check against MODIS NDVI / EVI 
# An independent greenness record (250 m, 2001-2024). If RPI and GPP track it, the signal is real.
log_phase("PHASE 6: cross-check against MODIS")

modis_aoi <- function(prefix, label) {
  z <- zonal_stat(gee_stack(prefix, "ANNUAL"), AOI.Sf)
  tibble::tibble(year = layer_year(names(z)), value = as.numeric(z[1, ]), variable = label)
}
Mod <- dplyr::bind_rows(modis_aoi("modis_ndvi_mean", "modis_ndvi"), modis_aoi("modis_evi_mean", "modis_evi"))
Ref <- Zonal.Annual |> dplyr::filter(variable %in% c("rpi", "gpp_actual")) |> dplyr::transmute(variable, year, value = mean)
cc <- function(a, b) {
  d <- dplyr::inner_join(dplyr::filter(Ref, variable == a), dplyr::filter(Mod, variable == b), by = "year", suffix = c("_a", "_b"))
  tibble::tibble(level = "AOI mean, per year", a = a, b = b, n = nrow(d), pearson = cor(d$value_a, d$value_b), spearman = cor(d$value_a, d$value_b, method = "spearman"))
}
Cross.Ts <- dplyr::bind_rows(cc("rpi", "modis_ndvi"), cc("rpi", "modis_evi"), cc("gpp_actual", "modis_ndvi"), cc("gpp_actual", "modis_evi"))

Ndvi.St   <- gee_stack("modis_ndvi_mean", "ANNUAL")
Ndvi.Cell <- rowMeans(as.matrix(sample_cells(Ndvi.St[[layer_year(Ndvi.St) %in% Years]], Cells)), na.rm = TRUE)   # centroid sampling (9.2)
Gpp.Cell  <- rowMeans(terra::values(Gpp.Act), na.rm = TRUE)
Rpi.Cell  <- rowMeans(terra::values(Rpi), na.rm = TRUE)
cor2 <- function(a, b, m = "pearson") suppressWarnings(cor(a, b, use = "complete.obs", method = m))
cell_cc <- function(level.txt, a.name, vals, idx) {          # (argument names must not clash with the tibble columns below)
  tibble::tibble(level = level.txt, a = a.name, b = "modis_ndvi", n = sum(idx & is.finite(vals) & is.finite(Ndvi.Cell)),
                 pearson = cor2(vals[idx], Ndvi.Cell[idx]), spearman = cor2(vals[idx], Ndvi.Cell[idx], "spearman"))
}
Spatial <- dplyr::bind_rows(
  cell_cc("cells in AOI, multi-year mean", "gpp_actual", Gpp.Cell, In.Aoi),
  cell_cc("cells in AOI, multi-year mean", "rpi", Rpi.Cell, In.Aoi),
  cell_cc("all envelope cells, multi-year mean", "gpp_actual", Gpp.Cell, Cells$in_env),
  cell_cc("all envelope cells, multi-year mean", "rpi", Rpi.Cell, Cells$in_env))
Cross <- dplyr::bind_rows(Cross.Ts, Spatial)
readr::write_csv(Cross, file.path(Dir.Out.Tables, "crosscheck_modis_vs_rpi.csv"))
print(as.data.frame(Cross |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)


# PHASE 7: Rasters 
log_phase("PHASE 7: rasters (RPI native grid)")

terra::writeRaster(Rpi.Now,   file.path(Dir.Out.Rasters, "rpi_now.tif"),       overwrite = TRUE)
terra::writeRaster(Sen.Rpi[["rpi_sen_slope"]], file.path(Dir.Out.Rasters, "rpi_sen_slope.tif"), overwrite = TRUE)
terra::writeRaster(Sen.Rpi[["rpi_sen_p"]],     file.path(Dir.Out.Rasters, "rpi_sen_p.tif"),     overwrite = TRUE)
terra::writeRaster(Degrading * 1, file.path(Dir.Out.Rasters, "rpi_degrading.tif"), overwrite = TRUE, datatype = "INT1U")
terra::writeRaster(Sen.Gpp[["gpp_sen_slope"]], file.path(Dir.Out.Rasters, "gpp_sen_slope.tif"), overwrite = TRUE)
message("  wrote 5 rasters to outputs/rasters/")

log_phase("DONE")
