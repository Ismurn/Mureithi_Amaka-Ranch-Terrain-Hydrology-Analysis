# Amaka Ranch - Terrain, Hydrology & Vegetation Analysis

This repository contains the R and Google Earth Engine (GEE) code for a terrain, hydrology and vegetation analysis of
**Amaka Ranch and Conservancy** (Taita-Taveta County, Kenya; catalogued as `MUSAMITI` in the TTWCA ranch shapefile),
within the Tsavo Conservation Area and the Kasigau Corridor REDD+ project area.

Use of this code is licensed under [GNU GPL v3](LICENSE).

Contact: Isaac Nduta Mureithi (African Wildlife Foundation) - *IMureithi@awf.org*.

Citation: *to be added once a report or paper exists.*

---

## Project overview

The analysis serves three operational questions, with large herbivores (elephants above all) as the management lens
and water plus cover as the two ranch-scale levers:

1. **Diagnose** why the existing water structures perform as they do (three that do not hold water for long, one working under heavy dry-season pressure).
2. **Monitor** the desilting dam now being worked on, by fixing a defensible pre-intervention baseline for a 6 / 12 / 24 month re-run.
3. **Site** new water pans, ranked with explicit trade-offs, for AWF / TTWCA funding proposals.

Every time-varying dataset is summarised for the **annual**, **long-rains (Mar-May)**, **short-rains (Oct-Dec)** and
**dry-season (Jun-Sep + Jan-Feb)** windows, so dry-season dynamics are not averaged away.

**Status (2026-10-27): version 1 has been run end to end** on the real data. The Word report is written to
`outputs/reports/Amaka_Hydro_Veg_Report_v1.docx` (not versioned). Read its first page before using any result: the analysis is a
satellite screening, the boundary is the coarse TTWCA polygon (20% below the listed size, two of the four located pans fall outside it), there is no ground
truth, and the desilting dam has no coordinates yet.

The full design brief, decisions and pitfalls are in [`DESCRIPTION.md`](DESCRIPTION.md). Differences found between the brief and
the data on disk, the decisions taken and the findings are logged in [`docs/known_issues.md`](docs/known_issues.md); the
method is in [`docs/methodology.md`](docs/methodology.md).

---

## Repository contents

### R scripts (`R/`) - run in numerical order

*Every script has been written and run on the real data. Parameters live in `00_setup.R`.*

| Script | Description |
|---|---|
| *00_setup.R* | Packages, local paths, run configuration (`Cfg.*`), season helpers, shared helpers (`read_pans()`, `get_analysis_grid()`, `zonal_stat()`, `scale_ref()`...), WhiteboxTools. Sourced by every other script. |
| *01_boundary_prep.R* | Verifies the shapefile CRS, builds the AOI (the TTWCA Musamiti polygon as delivered, or a supplied boundary file; `Cfg.AOI.Mode`), builds the 7 km catchment envelope, writes context layers and the GEE boundary. |
| *02_dem_prep.R* | Ingests the GEE DEM, reprojects, clips to the envelope, QA, compares with the phone altitudes at the field points. |
| *03_hydrology.R* | WhiteboxTools breach, D8 / D-infinity flow routing, streams (threshold calibrated), TWI, slope, TPI, Strahler order, cell terrain table. |
| *04_rpi_analysis.R* | RPI zonal statistics, Theil-Sen / Mann-Kendall trends, degrading cells, MODIS cross-check. |
| *05_pmlv2_ecohydrology.R* | Seasonal GPP, ET components, Es/ET, WUE; rainfall, drought years (SPI), SMAP, LST, MODIS context. |
| *06_landcover_studypixels.R* | Within-ranch land cover from 10 m WorldCover + Sentinel-2; shrub-to-tree check; bare-and-wet flag. The 500 m study pixels are the Kenya-wide baseline only. |
| *07_riparian_analysis.R* | Riparian corridor condition and the distance-to-stream gradient. |
| *08_water_pans_analysis.R* | Pan catchments, accessibility sheds, water gap, refugia value, failure-mode hypotheses, demand on the working dam, desilting baseline (needs coordinates). |
| *09_fire_history.R* | Burn count and last burn, year-by-year fire, fire against woody loss (timing-aware exposure test). |
| *10_integrated_stress.R* | Composite stress index (annual, DRY, LR) and recommendation zones A-E. |
| *11_pan_siting.R* | Suitability surface, candidate sites by non-maximum suppression, weight sensitivity, rehabilitation comparison. |
| *12_maps_figures.R* | Report figures (tmap + ggplot2). |

### GEE scripts (`gee/`) - run in the Earth Engine Code Editor

*Isaac ran 01-10 and downloaded the exports into* `data/gee_exports/` *(2026-10-27). He ran them as shared, with one saved change (the band naming in gee/10, already in this repository).
The night-LST export from 04 is empty (`docs/known_issues.md` G7); the script is patched but untested.*

| Script | Exports |
|---|---|
| *01_amaka_boundary_export.js* | The Amaka boundary as an EE asset (all other scripts read it). **Run first.** |
| *02_dem_amaka_export.js* | DEM over the buffered envelope (GLO-30, FABDEM or SRTM). |
| *03_chirps_rainfall.js* | Seasonal rainfall totals. |
| *04_modis_lst.js* | Seasonal day / night land-surface temperature. |
| *05_modis_vi.js* | Seasonal NDVI and EVI (250 m). |
| *06_sentinel2_vi_ndwi.js* | Seasonal NDVI, NDWI and bare-soil index (10 m). |
| *07_smap_soilmoisture.js* | Seasonal root-zone and surface soil moisture. |
| *08_hansen_woodycover.js* | Canopy cover, tree-cover loss, and 10 m WorldCover (the land-cover source). |
| *09_pmlv2_amaka_export.js* | Seasonal GPP, Ec, Es, Ei (the local PML stack is annual only). |
| *10_fire_history.js* | Burn count, most recent burn year, burned-by-year, FIRMS counts. |

### Report

`analysis/Amaka_Analysis.Rmd` reads every number from `outputs/` at knit time and knits to `outputs/reports/Amaka_Hydro_Veg_Report_v1.docx`
(knit it from RStudio, or run it with `RSTUDIO_PANDOC` pointing at RStudio's pandoc folder; `docs/known_issues.md` E6).

### Data

Nothing under `data/` or `outputs/` is versioned (except `data/water_pans/water_pans.csv`, `data/water_pans/amaka_field_points.csv` and `.gitkeep` placeholders).
Local inputs live on OneDrive and are declared in `R/00_setup.R`; they are catalogued in
[`docs/data_sources.md`](docs/data_sources.md).

| File | Description |
|---|---|
| *data/water_pans/water_pans.csv* | The five known pans (three degraded, one working under pressure, one being desilted without coordinates). Candidate sites are appended only on Isaac's confirmation (`Cfg.Siting.Append`). Schema in `DESCRIPTION.md` 3.3 plus `field_id`. |
| *data/water_pans/amaka_field_points.csv* | Points from the field visit (GPS, elevation, notes, proposed pan status) and the dam centres read from imagery. |

---

## Getting started

Clone the repository, open `research-repository-template.Rproj` in RStudio, and check the paths in Phase 3 of
`R/00_setup.R` (or set the `AMAKA_DATASETS_DIR` environment variable to the `Datasets_and_Shapefiles` folder).
Required: R >= 4.4 and the packages listed in `00_setup.R` (installed on first run by `pacman`).

## Running the analysis

1. `source("R/01_boundary_prep.R")` - writes `data/boundary/`, the CRS check outputs and the GEE boundary block.
2. **CRS check** (brief 9.1): *done* - confirmed as EPSG:32737 by overlay on imagery (`docs/known_issues.md` C1). To re-check on another machine or with a new boundary, overlay `outputs/figures/crs_check/*.kml` on imagery at zoom 14 or higher; if the file tag is wrong, set `Cfg.Ranches.CRS.Assumed` in `00_setup.R` and re-run step 1.
3. In the Code Editor: create the asset folder, run and export `gee/01_amaka_boundary_export.js` (instructions in its header), then `gee/02_dem_amaka_export.js` and the seasonal scripts. Save your edited copies back to `gee/`.
4. Download the Drive exports (`Amaka_GEE_Exports/`) into `data/gee_exports/`.
5. Run `R/02` to `R/12` in order (each sources `00_setup.R`; none takes more than a few minutes).
6. Add coordinates for the desilting dam (AMK_P05) to `data/water_pans/water_pans.csv` and re-run `R/08` to write its pre-intervention baseline.
7. Knit `analysis/Amaka_Analysis.Rmd`.

If the boundary changes (`Cfg.AOI.Mode <- "file"`), re-run from step 1: the GEE exports cover a 7 km envelope, so they only need repeating if the new boundary leaves it.
