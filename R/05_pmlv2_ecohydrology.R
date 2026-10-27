# ==============================================================================
# 05_pmlv2_ecohydrology.R  -  PML_V2 ecohydrology: GPP, Ec, Es, Ei, WUE, plus climate context
# ==============================================================================
# PURPOSE
#   Seasonal (LR / SR / DRY) and annual GPP, Ec, Es, Ei from the gee/09 export (calendar-year
#   seasons, mean daily rates); ET = Ec + Es + Ei; WUE = GPP / ET; Ec/ET and Es/ET (CLAUDE.md 6.4).
#   Es/ET is the bare-soil share of evapotranspiration and feeds the stress index. Also checks the
#   units against the local annual stack (9.4) and builds the climate context: CHIRPS rainfall
#   anomalies (z-score and SPI), SMAP soil moisture, MODIS land-surface temperature and NDVI / EVI.
#
# INPUTS
#   data/gee_exports/pml_{gpp,ec,es,ei}_meanrate_<SEASON>_2001_2020.tif   (gee/09)
#   data/gee_exports/chirps_*, smap_*, modis_lst_day_*, modis_{ndvi,evi}_*  (gee/03-07)
#   Path.PML.Annual (local, annual only: SUM of 8-day rates over Sep-Aug years)
#
# OUTPUTS
#   outputs/tables/pml_seasonal_zonal.csv        AOI per season x year x variable (rates)
#   outputs/tables/pml_summary.csv               per season: multi-year means, totals, trends
#   outputs/tables/crosscheck_pml_local_vs_gee.csv   unit / reducer check against the local stack
#   outputs/tables/climate_context_series.csv, climate_context_summary.csv, climate_vs_pml_correlations.csv
#   outputs/tables/cells_pml.csv                 one row per analysis cell (keyed on cell_id)
#   outputs/rasters/pml_<gpp|et|esfrac|ecfrac|wue>_<SEASON>.tif   multi-year means, PML native grid
#
# RUNTIME NOTES
#   * Seconds. Rates are mm d-1 (Ec, Es, Ei, ET) and gC m-2 d-1 (GPP). Season totals = rate x days
#     (LR 92, SR 92, DRY 181, ANNUAL 365); ratios do not depend on that factor.
#   * Fractions of AOI-level series are RATIOS OF AREA-WEIGHTED MEANS; cell-level fractions are ratios of
#     multi-year cell means. Cells are sampled at RPI centroids (9.2).
#   * Es/ET is highest in the RAINY seasons here (wet soil evaporates), and lower in the dry season,
#     when it is limited by soil water. The index uses the dry-season value as specified, but read it as
#     the bare-soil share of what little water the system moves, not as dryness itself (9.9).
#   * Night-time LST was exported blank (QC mask too strict, known_issues G7), so only day LST is used.
# ==============================================================================

if (!exists("Setup.Loaded")) source("R/00_setup.R")

Grid.List <- get_analysis_grid(); Cells <- Grid.List$cells
AOI.Sf <- read_aoi()
Days <- c(ANNUAL = 365, LR = 92, SR = 92, DRY = 181)
Vars <- c("gpp", "ec", "es", "ei")


# ---- PHASE 1: Load the PML seasonal stacks ---------------------------------------
log_phase("PHASE 1: PML seasonal stacks")

Pml <- lapply(setNames(Seasons.All, Seasons.All), function(s) {
  st <- lapply(setNames(Vars, Vars), function(v) gee_stack(sprintf("pml_%s_meanrate", v), s))
  st$et <- st$ec + st$es + st$ei
  for (k in names(st)) names(st[[k]]) <- sprintf("%s_%s_%d", k, s, layer_year(st[[k]]))
  st
})
Yrs <- layer_year(Pml$ANNUAL$gpp)
message(sprintf("  %d seasons x %d variables x %d years (%d-%d), grid %d x %d at %g m",
                length(Pml), length(Pml$ANNUAL), length(Yrs), min(Yrs), max(Yrs),
                terra::ncol(Pml$ANNUAL$gpp), terra::nrow(Pml$ANNUAL$gpp), terra::res(Pml$ANNUAL$gpp)[1]))


# ---- PHASE 2: AOI series per season, year and variable -----------------------------
log_phase("PHASE 2: AOI zonal series")

Zonal <- dplyr::bind_rows(lapply(Seasons.All, function(s) dplyr::bind_rows(lapply(c(Vars, "et"), function(v) {
  z <- zonal_stat(Pml[[s]][[v]], AOI.Sf)
  tibble::tibble(season = s, year = layer_year(names(z)), variable = v, value = as.numeric(z[1, ]))
}))))
Wide <- tidyr::pivot_wider(Zonal, names_from = variable, values_from = value)
stopifnot("list columns after pivot_wider (9.5)" = !any(vapply(Wide, is.list, logical(1))))
Wide <- Wide |> dplyr::mutate(wue = gpp / et, ec_frac = ec / et, es_frac = es / et, ei_frac = ei / et)
Pml.Zonal <- tidyr::pivot_longer(Wide, c(gpp, ec, es, ei, et, wue, ec_frac, es_frac, ei_frac), names_to = "variable", values_to = "value") |>
  dplyr::arrange(season, variable, year)
readr::write_csv(Pml.Zonal, file.path(Dir.Out.Tables, "pml_seasonal_zonal.csv"))

trend_of <- function(x) { x <- x[is.finite(x)]
  tibble::tibble(sen_slope_per_year = unname(trend::sens.slope(x)$estimates), mk_p = trend::mk.test(x)$p.value) }
Pml.Summary <- Wide |> dplyr::group_by(season) |>
  dplyr::summarise(dplyr::across(c(gpp, ec, es, ei, et, wue, ec_frac, es_frac), list(mean = ~ mean(.x)), .names = "{.col}_{.fn}"),
                   es_frac_recent5 = mean(es_frac[year >= max(year) - 4]), gpp_recent5 = mean(gpp[year >= max(year) - 4]), .groups = "drop") |>
  dplyr::mutate(days = Days[season], et_mm_per_season = et_mean * days, gpp_gC_per_season = gpp_mean * days)
Tr <- Wide |> dplyr::group_by(season) |>
  dplyr::summarise(es_frac_sen = trend_of(es_frac)$sen_slope_per_year, es_frac_mk_p = trend_of(es_frac)$mk_p,
                   gpp_sen = trend_of(gpp)$sen_slope_per_year, gpp_mk_p = trend_of(gpp)$mk_p,
                   et_sen = trend_of(et)$sen_slope_per_year, et_mk_p = trend_of(et)$mk_p, .groups = "drop")
Pml.Summary <- dplyr::left_join(Pml.Summary, Tr, by = "season")
readr::write_csv(Pml.Summary, file.path(Dir.Out.Tables, "pml_summary.csv"))
print(as.data.frame(Pml.Summary |> dplyr::select(season, gpp_mean, et_mean, wue_mean, ec_frac_mean, es_frac_mean, es_frac_recent5, es_frac_sen, es_frac_mk_p, et_mm_per_season) |>
                      dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)


# ---- PHASE 3: Unit and reducer check against the local annual stack (9.4) ---------------
# The local stack is a SUM of ~45 8-day composites of daily rates over a Sep-Aug year; x 8 gives mm or gC m-2
# per year. The GEE export is the MEAN daily rate of calendar-year composites; x 365 gives the same quantity.
# If the two annual totals agree the scale factor (PML_SCALE = 1) and the reducer are right.
log_phase("PHASE 3: check against the local annual stack")

Loc <- zonal_stat(terra::rast(Path.PML.Annual), AOI.Sf)
Loc.Long <- tibble::tibble(layer = names(Loc), value = as.numeric(Loc[1, ])) |>
  dplyr::mutate(variable = tolower(sub("^[0-9]+_([A-Za-z]+)_[0-9]{4}$", "\\1", layer)), year = as.integer(sub(".*_", "", layer))) |>
  dplyr::filter(year %in% Cfg.Years.PML.Local)
Loc.Tot <- Loc.Long |> dplyr::group_by(variable) |> dplyr::summarise(local_total = mean(value) * Cfg.PML.Composite.Days, .groups = "drop")
Gee.Tot <- Wide |> dplyr::filter(season == "ANNUAL", year %in% Cfg.Years.PML.Local) |> dplyr::summarise(dplyr::across(c(gpp, ec, es, ei), mean)) |>
  tidyr::pivot_longer(tidyselect::everything(), names_to = "variable", values_to = "gee_rate") |> dplyr::mutate(gee_total = gee_rate * 365.25)
Series.Cor <- dplyr::inner_join(Loc.Long, Wide |> dplyr::filter(season == "ANNUAL") |> dplyr::select(year, gpp, ec, es, ei) |>
  tidyr::pivot_longer(c(gpp, ec, es, ei), names_to = "variable", values_to = "gee_rate"), by = c("variable", "year")) |>
  dplyr::group_by(variable) |> dplyr::summarise(year_to_year_pearson = cor(value, gee_rate), .groups = "drop")
Check <- Loc.Tot |> dplyr::inner_join(Gee.Tot, by = "variable") |> dplyr::inner_join(Series.Cor, by = "variable") |>
  dplyr::mutate(ratio_gee_over_local = gee_total / local_total) |>
  dplyr::select(variable, local_total, gee_total, ratio_gee_over_local, year_to_year_pearson)
readr::write_csv(Check, file.path(Dir.Out.Tables, "crosscheck_pml_local_vs_gee.csv"))
print(as.data.frame(Check |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
message(sprintf("  GEE / local annual totals, ratio range %.3f-%.3f (about 1 means scale and reducer are right)",
                min(Check$ratio_gee_over_local), max(Check$ratio_gee_over_local)))


# ---- PHASE 4: Multi-year means, rasters and the cell table ----------------------------------
log_phase("PHASE 4: multi-year means and cell table")

Mean.R <- lapply(Seasons.All, function(s) {
  m <- lapply(Pml[[s]], terra::mean, na.rm = TRUE)                       # mean over the 20 years, per variable
  list(gpp = m$gpp, ec = m$ec, es = m$es, ei = m$ei, et = m$et, wue = m$gpp / m$et, ecfrac = m$ec / m$et, esfrac = m$es / m$et)
}) |> setNames(Seasons.All)

for (s in Seasons.All) for (v in c("gpp", "et", "esfrac", "ecfrac", "wue")) {
  r <- Mean.R[[s]][[v]]; names(r) <- paste0(v, "_", s)
  terra::writeRaster(r, file.path(Dir.Out.Rasters, sprintf("pml_%s_%s.tif", v, s)), overwrite = TRUE)
}
Cells.Pml <- tibble::tibble(cell_id = Cells$cell_id)
for (s in Seasons.All) for (v in names(Mean.R[[s]])) {
  r <- Mean.R[[s]][[v]]; names(r) <- paste0(v, "_", s)
  Cells.Pml <- dplyr::bind_cols(Cells.Pml, sample_cells(r, Cells))
}
# the recent five years of dry-season Es/ET, for information
Recent <- layer_year(Pml$DRY$es) >= max(Yrs) - 4
Cells.Pml$esfrac_DRY_recent5 <- sample_cells(terra::mean(Pml$DRY$es[[Recent]]) / terra::mean(Pml$DRY$et[[Recent]]), Cells)[[1]]
readr::write_csv(Cells.Pml, file.path(Dir.Out.Tables, "cells_pml.csv"))
In.Aoi <- Cells$in_aoi
message(sprintf("  AOI cells: Es/ET  ANNUAL %.2f | LR %.2f | SR %.2f | DRY %.2f (cell range in DRY %.2f-%.2f)",
                mean(Cells.Pml$esfrac_ANNUAL[In.Aoi]), mean(Cells.Pml$esfrac_LR[In.Aoi]), mean(Cells.Pml$esfrac_SR[In.Aoi]),
                mean(Cells.Pml$esfrac_DRY[In.Aoi]), min(Cells.Pml$esfrac_DRY[In.Aoi]), max(Cells.Pml$esfrac_DRY[In.Aoi])))


# ---- PHASE 5: Climate context --------------------------------------------------------------------
log_phase("PHASE 5: climate context")

series_aoi <- function(prefix, dataset, variable) dplyr::bind_rows(lapply(Seasons.All, function(s) {
  z <- zonal_stat(gee_stack(prefix, s), AOI.Sf)
  tibble::tibble(dataset = dataset, variable = variable, season = s, year = layer_year(names(z)), value = as.numeric(z[1, ]))
}))
Climate <- dplyr::bind_rows(
  series_aoi("chirps_precip_sum",     "CHIRPS",       "rainfall_mm"),
  series_aoi("smap_rootzone_mean",    "SMAP L4",      "soil_moisture_rootzone"),
  series_aoi("smap_surface_mean",     "SMAP L4",      "soil_moisture_surface"),
  series_aoi("modis_lst_day_mean",    "MODIS MOD11A2", "lst_day_C"),
  series_aoi("modis_ndvi_mean",       "MODIS MOD13Q1", "ndvi"),
  series_aoi("modis_evi_mean",        "MODIS MOD13Q1", "evi"))

# Rainfall anomalies: z-score against the 2000-2024 climatology, and SPI from a gamma fit per season.
spi_gamma <- function(x) {
  ok <- is.finite(x) & x > 0
  fit <- if (sum(ok) >= 10) tryCatch(MASS::fitdistr(x[ok], "gamma"), error = function(e) NULL) else NULL
  if (is.null(fit)) return(rep(NA_real_, length(x)))
  qnorm(pgamma(x, shape = fit$estimate[["shape"]], rate = fit$estimate[["rate"]]))
}
Climate <- Climate |> dplyr::group_by(dataset, variable, season) |>
  dplyr::mutate(z = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE),
                spi = if (variable[1] == "rainfall_mm") spi_gamma(value) else NA_real_,
                drought = variable == "rainfall_mm" & !is.na(spi) & spi <= -1,
                wet = variable == "rainfall_mm" & !is.na(spi) & spi >= 1) |> dplyr::ungroup()
readr::write_csv(Climate, file.path(Dir.Out.Tables, "climate_context_series.csv"))

Clim.Summary <- Climate |> dplyr::group_by(dataset, variable, season) |>
  dplyr::summarise(n_years = sum(is.finite(value)), years = paste(range(year), collapse = "-"), mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
                   sen_slope_per_year = unname(trend::sens.slope(value[is.finite(value)])$estimates), mk_p = trend::mk.test(value[is.finite(value)])$p.value,
                   n_drought_years = sum(drought, na.rm = TRUE), n_wet_years = sum(wet, na.rm = TRUE), .groups = "drop")
readr::write_csv(Clim.Summary, file.path(Dir.Out.Tables, "climate_context_summary.csv"))
print(as.data.frame(Clim.Summary |> dplyr::filter(variable %in% c("rainfall_mm", "lst_day_C"), season %in% c("ANNUAL", "LR", "SR", "DRY")) |>
                      dplyr::select(variable, season, years, mean, sd, sen_slope_per_year, mk_p, n_drought_years, n_wet_years) |>
                      dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
Dry.Years <- Climate |> dplyr::filter(variable == "rainfall_mm", season == "ANNUAL", drought) |> dplyr::pull(year)
message("  annual drought years (SPI <= -1): ", ifelse(length(Dry.Years), paste(Dry.Years, collapse = ", "), "none"))

# Does rainfall explain the PML signals? (AOI means, same calendar-year seasons, 2001-2020)
Rain <- Climate |> dplyr::filter(variable == "rainfall_mm") |> dplyr::select(season, year, rain = value)
Cor.Pml <- Wide |> dplyr::inner_join(Rain, by = c("season", "year")) |> dplyr::group_by(season) |>
  dplyr::summarise(n = dplyr::n(), gpp_vs_rain = cor(gpp, rain), et_vs_rain = cor(et, rain), es_frac_vs_rain = cor(es_frac, rain),
                   gpp_vs_rain_spearman = cor(gpp, rain, method = "spearman"), .groups = "drop")
readr::write_csv(Cor.Pml, file.path(Dir.Out.Tables, "climate_vs_pml_correlations.csv"))
print(as.data.frame(Cor.Pml |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)

log_phase("DONE")
