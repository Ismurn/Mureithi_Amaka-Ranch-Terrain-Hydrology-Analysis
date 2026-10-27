# ==============================================================================
# 06_landcover_studypixels.R  -  Land cover, shrub-to-tree ratio and bare-and-wet flag
# ==============================================================================
# PURPOSE
#   Within-ranch land-cover structure for Amaka. Decision of 2026-10-26: this comes from the 10 m ESA
#   WorldCover map (class fractions) plus Sentinel-2 (bare-soil index, NDVI, NDWI), NOT from the 500 m
#   study pixels, which are ~97% shrubland over Amaka (known_issues D5). The study pixels are kept as the
#   Kenya-wide baseline. This script also TESTS that decision: how much within-ranch variation does each
#   candidate layer actually carry? Derives per analysis cell:
#     - tree / shrub / grass / bare / other fractions (10 m classes aggregated, no resampling),
#     - the shrub-to-tree ratio (McCleery 2025, Jones 2025),
#     - seasonal Sentinel-2 BSI / NDVI / NDWI (means of 2019-2024),
#     - a bare-and-wet flag for Zone A (high BSI AND high TWI, relative to the ranch),
#   and a cross-tab of land cover x TWI quintile x RPI quintile.
#
# INPUTS
#   data/gee_exports/worldcover_2021_envelope.tif, s2_{bsi,ndvi,ndwi}_mean_<SEASON>_2019_2024.tif
#   Path.StudyPixels (baseline comparison); outputs/tables/cells_terrain.csv (03), cells_rpi.csv (04)
#
# OUTPUTS
#   outputs/tables/cells_landcover.csv           one row per analysis cell (keyed on cell_id)
#   outputs/tables/landcover_fractions.csv       Kenya baseline vs Amaka at 500 m vs Amaka at 10 m
#   outputs/tables/landcover_signal_check.csv    within-ranch variation of every candidate layer
#   outputs/tables/crosstab_lc_twi_rpi.csv       dominant class x TWI quintile x RPI quintile (AOI cells)
#   outputs/tables/s2_aoi_series.csv             Sentinel-2 AOI means per season and year
#   outputs/rasters/lc_{trees,shrub,grass,bare}_rpigrid.tif, shrub_tree_ratio.tif, bsi_<SEASON>.tif, bare_and_wet.tif
#
# RUNTIME NOTES
#   * A minute or two (twelve 136 MB Sentinel-2 stacks). Class fractions use exactextractr "frac" on the
#     native 10 m grid; nothing is resampled (CLAUDE.md 9.2).
#   * "High" is RELATIVE to the ranch (top quartile of AOI cells, Cfg.Relative.Q): absolute thresholds select
#     nothing in a landscape that is almost entirely shrubland.
#   * The shrub-to-tree ratio is undefined where tree cover is ~0; it is reported only where trees >= 1% of a
#     cell and paired with an explicit class ("no trees", "shrub-dominated", ...), never divided by zero.
#   * WorldCover v200 (2021) and the study pixels' v100 (2020) are different classifications.
# ==============================================================================

if (!exists("Setup.Loaded")) source("R/00_setup.R")

Grid.List <- get_analysis_grid(); G <- Grid.List$grid; Cells <- Grid.List$cells
AOI.Sf    <- read_aoi()
Cell.Poly <- cell_polygons(G)
Terrain   <- readr::read_csv(file.path(Dir.Out.Tables, "cells_terrain.csv"), show_col_types = FALSE)
Rpi.Tbl   <- readr::read_csv(file.path(Dir.Out.Tables, "cells_rpi.csv"), show_col_types = FALSE)
In.Aoi    <- Cells$in_aoi


# ---- PHASE 1: WorldCover class fractions per cell ---------------------------------
log_phase("PHASE 1: WorldCover 10 m class fractions")

Wc <- terra::rast(gee_file("worldcover_2021_envelope"))
Cls.Names <- c(`10` = "trees", `20` = "shrub", `30` = "grass", `40` = "crop", `50` = "built", `60` = "bare",
               `70` = "snow", `80` = "water", `90` = "wetland", `95` = "mangrove", `100` = "moss")
Fr <- exactextractr::exact_extract(Wc, Cell.Poly, fun = "frac", progress = FALSE)
for (k in names(Cls.Names)) if (!paste0("frac_", k) %in% names(Fr)) Fr[[paste0("frac_", k)]] <- 0
Fr[is.na(Fr)] <- 0
Frac <- tibble::tibble(cell_id = Cell.Poly$cell_id)
for (k in names(Cls.Names)) Frac[[paste0("wc_", Cls.Names[[k]])]] <- Fr[[paste0("frac_", k)]]
Frac$wc_other <- pmax(0, 1 - Frac$wc_trees - Frac$wc_shrub - Frac$wc_grass - Frac$wc_bare)
Frac$wc_dominant <- Cls.Names[as.character(c(10, 20, 30, 40, 50, 60, 80, 90)[max.col(as.matrix(Frac[, paste0("wc_", Cls.Names[c("10", "20", "30", "40", "50", "60", "80", "90")])]), ties.method = "first")])]
Frac$wc_woody <- Frac$wc_trees + Frac$wc_shrub                       # tree + shrub cover (Zone D, siting)

Class.Tbl <- exactextractr::exact_extract(Wc, sf::st_transform(AOI.Sf, terra::crs(Wc)), fun = "frac", progress = FALSE)
message("  WorldCover composition of the AOI (area-weighted): ",
        paste(sprintf("%s %.1f%%", Cls.Names[sub("frac_", "", names(Class.Tbl))], 100 * unlist(Class.Tbl[1, ])), collapse = " | "))

# Shrub-to-tree ratio with explicit handling of "no trees".
Frac$shrub_tree_ratio <- ifelse(Frac$wc_trees >= 0.01, Frac$wc_shrub / Frac$wc_trees, NA_real_)
Frac$shrub_tree_class <- dplyr::case_when(
  Frac$wc_trees < 0.01                    ~ "no trees (<1%)",
  Frac$shrub_tree_ratio > Cfg.Sterile.Ratio ~ "shrub-dominated (ratio > 5)",
  Frac$shrub_tree_ratio >= 1              ~ "mixed (1-5)",
  TRUE                                    ~ "tree-dominated (<1)")


# ---- PHASE 2: Sentinel-2 BSI / NDVI / NDWI per cell and season ------------------------
log_phase("PHASE 2: Sentinel-2 seasonal indices")

S2 <- tibble::tibble(cell_id = Cell.Poly$cell_id)
S2.Series <- list()
Dir.S2 <- file.path(Dir.Out.Rasters, "s2"); dir.create(Dir.S2, showWarnings = FALSE)      # multi-year 10 m means, read by 07 and 09
for (idx in c("bsi", "ndvi", "ndwi")) for (s in Seasons.All) {
  st <- gee_stack(sprintf("s2_%s_mean", idx), s)
  m  <- terra::app(st, "mean", na.rm = TRUE); names(m) <- sprintf("%s_%s", idx, s)
  terra::writeRaster(m, file.path(Dir.S2, sprintf("s2_%s_%s_mean.tif", idx, s)), overwrite = TRUE)
  S2[[sprintf("%s_%s", idx, s)]] <- zonal_stat(m, Cell.Poly)[[1]]
  z <- zonal_stat(st, AOI.Sf)
  S2.Series[[length(S2.Series) + 1]] <- tibble::tibble(dataset = "Sentinel-2", variable = idx, season = s, year = layer_year(names(z)), value = as.numeric(z[1, ]))
  message(sprintf("  %s %-6s done", idx, s))
}
S2.Series <- dplyr::bind_rows(S2.Series)
readr::write_csv(S2.Series, file.path(Dir.Out.Tables, "s2_aoi_series.csv"))


# ---- PHASE 3: Does 10 m carry more within-ranch signal than 500 m? (known_issues D5) ----
# Candidate layers for "bare ground" and "woody cover", compared on the AOI cells only. A layer with almost
# no spread cannot support a min-max-scaled stress term or a relative Zone A threshold.
log_phase("PHASE 3: within-ranch signal check")

Sp <- terra::rast(Path.StudyPixels)                                    # Kenya-wide, 500 m, WGS84
Sp.Cell <- sample_cells(Sp[[c("trees", "shrubland", "grassland", "bare")]], Cells, geo = TRUE)
names(Sp.Cell) <- paste0("sp_", c("trees", "shrub", "grass", "bare"))
Cand <- dplyr::bind_cols(Frac[, c("wc_trees", "wc_shrub", "wc_grass", "wc_bare", "wc_woody")], Sp.Cell, S2[, c("bsi_DRY", "bsi_LR", "ndvi_DRY", "ndvi_LR")])[In.Aoi, ]
Cand.Meta <- tibble::tibble(layer = names(Cand), source = dplyr::case_when(grepl("^wc_", names(Cand)) ~ "WorldCover 10 m, per-cell fraction",
                                                                    grepl("^sp_", names(Cand)) ~ "study pixels 500 m, at cell centre",
                                                                    TRUE ~ "Sentinel-2 10 m, per-cell mean"))
Signal <- dplyr::bind_cols(Cand.Meta, tibble::tibble(
  mean = vapply(Cand, mean, 0, na.rm = TRUE), sd = vapply(Cand, sd, 0, na.rm = TRUE),
  min = vapply(Cand, min, 0, na.rm = TRUE), p10 = vapply(Cand, quantile, 0, 0.1, na.rm = TRUE),
  p90 = vapply(Cand, quantile, 0, 0.9, na.rm = TRUE), max = vapply(Cand, max, 0, na.rm = TRUE)))
Signal$range <- Signal$max - Signal$min
Signal$cv <- Signal$sd / abs(Signal$mean)
readr::write_csv(Signal, file.path(Dir.Out.Tables, "landcover_signal_check.csv"))
print(as.data.frame(Signal |> dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

Bare.Signal <- Signal[Signal$layer %in% c("wc_bare", "sp_bare", "bsi_DRY"), ]
message(sprintf("  bare ground inside the AOI: study pixels range %.3f, WorldCover range %.3f, Sentinel-2 BSI (dry) range %.3f",
                Signal$range[Signal$layer == "sp_bare"], Signal$range[Signal$layer == "wc_bare"], Signal$range[Signal$layer == "bsi_DRY"]))
message(sprintf("  tree cover inside the AOI: study pixels SD %.3f vs WorldCover SD %.3f",
                Signal$sd[Signal$layer == "sp_trees"], Signal$sd[Signal$layer == "wc_trees"]))


# ---- PHASE 4: Bare-and-wet flag (Zone A input) ----------------------------------------------------
# High = top quartile (Cfg.Relative.Q) of the AOI cells; the thresholds come from the AOI and are applied to every
# cell, so pan catchments that extend beyond the AOI are classified on the same footing.
log_phase("PHASE 4: bare-and-wet flag")

Joined <- dplyr::left_join(Cells, Terrain, by = "cell_id") |> dplyr::left_join(Rpi.Tbl, by = "cell_id") |>
  dplyr::left_join(Frac, by = "cell_id") |> dplyr::left_join(S2, by = "cell_id")
q.aoi <- function(x) unname(quantile(x[In.Aoi], Cfg.Relative.Q, na.rm = TRUE))
Thr <- c(twi = q.aoi(Joined$twi_mean), bsi_dry = q.aoi(Joined$bsi_DRY))
Joined$bare_high <- Joined$bsi_DRY >= Thr[["bsi_dry"]]
Joined$twi_high  <- Joined$twi_mean >= Thr[["twi"]]
Joined$bare_and_wet <- Joined$bare_high & Joined$twi_high
message(sprintf("  thresholds from AOI cells: TWI >= %.2f, dry-season BSI >= %.3f | bare-and-wet cells: %d inside the AOI (%.0f%% of AOI cells), %d in the envelope",
                Thr[["twi"]], Thr[["bsi_dry"]], sum(Joined$bare_and_wet[In.Aoi], na.rm = TRUE), 100 * mean(Joined$bare_and_wet[In.Aoi], na.rm = TRUE),
                sum(Joined$bare_and_wet, na.rm = TRUE)))
message(sprintf("  correlation of dry-season BSI with TWI inside the AOI: %.2f", cor(Joined$bsi_DRY[In.Aoi], Joined$twi_mean[In.Aoi], use = "complete.obs")))


# ---- PHASE 5: Cross-tab, baseline table, rasters -----------------------------------------------------------
log_phase("PHASE 5: cross-tab, baseline comparison, rasters")

Cross <- Joined[In.Aoi, ] |>
  dplyr::group_by(wc_dominant, twi_q, rpi_now_q) |>
  dplyr::summarise(n_cells = dplyr::n(), mean_bsi_dry = mean(bsi_DRY, na.rm = TRUE), mean_wc_trees = mean(wc_trees), mean_wc_shrub = mean(wc_shrub),
                   mean_rpi_now = mean(rpi_now, na.rm = TRUE), n_bare_and_wet = sum(bare_and_wet, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(twi_q, rpi_now_q)
readr::write_csv(Cross, file.path(Dir.Out.Tables, "crosstab_lc_twi_rpi.csv"))
message(sprintf("  cross-tab: %d TWI x RPI quintile combinations; dominant WorldCover class(es) in the AOI cells: %s",
                nrow(Cross), paste(unique(Cross$wc_dominant), collapse = ", ")))

# Kenya-wide baseline versus Amaka at 500 m versus Amaka at 10 m
kenya <- terra::global(Sp[[c("trees", "shrubland", "grassland", "bare", "allRangelands")]], "mean", na.rm = TRUE)[[1]]
amaka.sp <- unlist(zonal_stat(Sp[[c("trees", "shrubland", "grassland", "bare", "allRangelands")]], AOI.Sf)[1, ])
wcf <- unlist(Class.Tbl[1, ]); names(wcf) <- Cls.Names[sub("frac_", "", names(wcf))]
getf <- function(n) if (n %in% names(wcf)) unname(wcf[[n]]) else 0
Env.Wc <- unlist(exactextractr::exact_extract(Wc, sf::st_transform(read_envelope(), terra::crs(Wc)), fun = "frac", progress = FALSE)[1, ])
names(Env.Wc) <- Cls.Names[sub("frac_", "", names(Env.Wc))]
getf.env <- function(n) if (n %in% names(Env.Wc)) unname(Env.Wc[[n]]) else 0
Baseline <- tibble::tibble(
  source = c("Study pixels, Kenya-wide mean (500 m, WorldCover v100 2020)", "Study pixels, Amaka AOI (500 m)", "WorldCover v200 2021, Amaka AOI (10 m)", "WorldCover v200 2021, 7 km envelope (10 m)"),
  trees = c(kenya[1], amaka.sp[1], getf("trees"), getf.env("trees")), shrubland = c(kenya[2], amaka.sp[2], getf("shrub"), getf.env("shrub")),
  grassland = c(kenya[3], amaka.sp[3], getf("grass"), getf.env("grass")), bare = c(kenya[4], amaka.sp[4], getf("bare"), getf.env("bare")))
readr::write_csv(Baseline, file.path(Dir.Out.Tables, "landcover_fractions.csv"))
print(as.data.frame(Baseline |> dplyr::mutate(source = substr(source, 1, 52), dplyr::across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)

Out <- Frac |> dplyr::left_join(S2, by = "cell_id") |> dplyr::left_join(Joined[, c("cell_id", "bare_high", "twi_high", "bare_and_wet")], by = "cell_id")
Out$bare_and_wet <- as.integer(Out$bare_and_wet); Out$bare_high <- as.integer(Out$bare_high); Out$twi_high <- as.integer(Out$twi_high)
readr::write_csv(Out, file.path(Dir.Out.Tables, "cells_landcover.csv"))

to_grid <- function(v, nm) { r <- terra::rast(G); terra::values(r) <- v; names(r) <- nm; r }
for (v in c("trees", "shrub", "grass", "bare")) terra::writeRaster(to_grid(Frac[[paste0("wc_", v)]], paste0("lc_", v)),
                                                              file.path(Dir.Out.Rasters, sprintf("lc_%s_rpigrid.tif", v)), overwrite = TRUE)
terra::writeRaster(to_grid(Frac$shrub_tree_ratio, "shrub_tree_ratio"), file.path(Dir.Out.Rasters, "shrub_tree_ratio.tif"), overwrite = TRUE)
for (s in Seasons.All) terra::writeRaster(to_grid(S2[[paste0("bsi_", s)]], paste0("bsi_", s)), file.path(Dir.Out.Rasters, sprintf("bsi_%s.tif", s)), overwrite = TRUE)
terra::writeRaster(to_grid(Out$bare_and_wet, "bare_and_wet"), file.path(Dir.Out.Rasters, "bare_and_wet.tif"), overwrite = TRUE, datatype = "INT1U")
message("  wrote cells_landcover.csv and 10 rasters")

tc <- table(Frac$shrub_tree_class[In.Aoi]); message("  shrub:tree classes inside the AOI: ", paste(sprintf("%s = %d", names(tc), tc), collapse = " | "))

log_phase("DONE")
