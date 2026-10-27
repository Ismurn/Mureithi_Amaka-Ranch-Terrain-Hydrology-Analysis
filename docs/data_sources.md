# Data sources

Every dataset the pipeline reads or will read. Local paths are declared once, in `R/00_setup.R`
(Phase 3); GEE assets are declared at the top of each `gee/*.js` script.

## 1. Local inputs (OneDrive, never copied into the repository)

Root: `Datasets_and_Shapefiles/` (override with the `AMAKA_DATASETS_DIR` environment variable).

| Object in `00_setup.R` | File | Content | Grid / CRS |
|---|---|---|---|
| `Path.Ranches` | `Tsavo_Conservation_Areas/RanchesTTWCA.shp` | 28 TTWCA ranch polygons. **Amaka = `R_NAME == "MUSAMITI"`.** | tagged EPSG:32737, confirmed by overlay on imagery (`known_issues.md` C1) |
| `Path.SettledAreas` | `Tsavo_Conservation_Areas/SettledAreas.shp` | 9 settled-area polygons (one has an empty name) | tagged EPSG:32737, same source as the ranches (C1) |
| `Path.TsavoEast` / `Path.TsavoWest` | `WDPA_Amboseli_Tsavo/National Parks&Reserves/Tsavo{East,West}-WDPA_..._shp/..._shp_0/..._shp-polygons.shp` | WDPA national park polygons (Apr 2025) | EPSG:4326 |
| `Path.Kenya.Admin1` | `Country_Borders/kenya_regional_admin1_boundary_wgs84.shp` | Kenya regional boundaries | EPSG:4326 |
| `Path.RPI` | `RPI Datasets (GuyLomax)/rpi_ea_v2.tif` | Relative Performance Index, layers `rpi.2000` ... `rpi.2022` (23) | EPSG:4326, 0.004167 deg (~460 m), East Africa |
| `Path.RPI.GPP.Act` | `.../actual_gpp_ea_v2.tif` | Actual GPP, `GPP.2000` ... `GPP.2022` | as above |
| `Path.RPI.GPP.Pot` | `.../potential_gpp_ea_v2.tif` | Potential GPP, `potential_gpp_predicted.2000` ... | as above |
| `Path.RPI.TempPerf` | `.../rpi_ea_temporal_performance.tif` | Model skill: `mae`, `rmse`, `rsq` | as above |
| `Path.PML.Annual` | `PMLv2_and_IGBPpixels/PML_V2_Kenya_Annual_Stack.tif` | **Annual** GPP, Ec, Es, Ei, 2001-2020, 80 layers named `<i>_<var>_<year>`. Each = **sum** of ~45 8-day composites of daily rates (`CAS/IGSNRR/PML/V2_v017`) over a **Sep-Aug** year; x 8 for mm / gC m-2 per year; 2020 looks partial. See `known_issues.md` D3. | EPSG:4326, 0.004492 deg (~500 m), Kenya |
| `Path.PML.Mean` | `PMLv2_and_IGBPpixels/PML_V2_Kenya_Drive_Export.tif` | 4 layers: GPP, Ec, Es, Ei = mean of the annual sums above (includes partial 2020) | as above |
| `Path.StudyPixels` | `PMLv2_and_IGBPpixels/studyPixels_Kenya_Export.tif` | `igbp` (MODIS MCD12Q1 mode 2001-2023), `trees`, `shrubland`, `grassland`, `bare`, `allRangelands` (ESA WorldCover **v100**, 2020, fractions at 500 m). Kept only where IGBP is 6-10 and `allRangelands` >= 0.9, within three Kenyan ecoregions. In Amaka it is ~97% shrubland (`known_issues.md` D5). | EPSG:4326, 0.004492 deg |
| `Path.Pans` | `data/water_pans/water_pans.csv` (in the repo) | Existing / candidate pans, schema in `DESCRIPTION.md` 3.3 plus `field_id` (known_issues F17). Five rows: AMK_P01 Mickles 2, AMK_P02 HQ reservoir, AMK_P03 Mickles 1 (degraded), AMK_P04 (functional under pressure), AMK_P05 the desilting dam (no coordinates). | WGS84 lon / lat |
| (field file) | `data/water_pans/amaka_field_points.csv` (in the repo) | Isaac's iPhone-GPS visit: 7 located points (borehole, 3 dams, rock-catchment reservoir, 2 scout posts) and 2 unsurveyed features (desilting dam, planned EL tank), with the original DMS strings, elevation, `position_type` (close for HQ and scout posts, standoff for dams and borehole) and `uncertainty_m_est` (50 / 150 m, estimates), image-derived dam centres `lat_img` / `lon_img` for PT02, PT04 (reservoir), PT05, PT07 (Esri imagery of 2021-10-23, not surveyed; originals untouched), and *proposed* pan type / status. See `known_issues.md` section 6. | WGS84 lon / lat |

## 2. Google Earth Engine exports

Account `ee-imnduta1`; scripts saved under `users/imnduta1/Ismurn/Amaka/`; exports go to the Drive folder
`Amaka_GEE_Exports/` and are downloaded into `data/gee_exports/`. Everything is exported over the **buffered
envelope** (`Cfg.Buffer.m` + 1 km pad), EPSG:32737. **Isaac ran scripts 01-10 as shared (one saved change: the band naming in `gee/10`) and downloaded the 63 GeoTIFFs (2026-10-27)**; the
one defect is the night-LST export, which is all NoData (known_issues G7; `gee/04` is patched but untested).

Seasons (same calendar year): **ANNUAL** Jan-Dec, **LR** Mar-May, **SR** Oct-Dec, **DRY** Jan-Feb + Jun-Sep.

Seasonal files are named `<prefix>_<SEASON>_<y0>_<y1>.tif`, one band per year named `<label>_<SEASON>_<year>`.

| Script | Asset | Variable(s) -> output prefix (band label) | Reducer | Export scale | Years |
|---|---|---|---|---|---|
| `01_amaka_boundary_export.js` | (literal coordinates from R) | boundary -> EE asset `Amaka_Ranch` | - | - | - |
| `02_dem_amaka_export.js` | `COPERNICUS/DEM/GLO30` (or FABDEM / SRTM) | `amaka_dem_<glo30\|fabdem\|srtm>_envelope` (`DEM`) | mosaic, bilinear | 30 m | static |
| `03_chirps_rainfall.js` | `UCSB-CHG/CHIRPS/DAILY` | `chirps_precip_sum` (`precip_mm`) | **sum** (mm) | 5566 m native | 2000-2024 |
| `04_modis_lst.js` | `MODIS/061/MOD11A2` | `modis_lst_day_mean` (`lst_day_C`), `modis_lst_night_mean` (`lst_night_C`, **empty in the first export**, G7) | mean, day QC = good, night QC <= other quality (patched, untested) | 1000 m native | 2001-2024 |
| `05_modis_vi.js` | `MODIS/061/MOD13Q1` | `modis_ndvi_mean` (`ndvi`), `modis_evi_mean` (`evi`) | mean, SummaryQA <= 1 | 250 m native | 2001-2024 |
| `06_sentinel2_vi_ndwi.js` | `COPERNICUS/S2_SR_HARMONIZED` | `s2_ndvi_mean` (`ndvi`), `s2_ndwi_mean` (`ndwi`), `s2_bsi_mean` (`bsi`, bare-soil index: continuous bare proxy for Zone A) | mean, SCL mask, cloud < 60% | 10 m | 2019-2024 |
| `07_smap_soilmoisture.js` | `NASA/SMAP/SPL4SMGP/008` | `smap_rootzone_mean` (`sm_rootzone`), `smap_surface_mean` (`sm_surface`) | mean | 9000 m native | 2016-2024 |
| `08_hansen_woodycover.js` | `UMD/hansen/global_forest_change_2023_v1_11`; `ESA/WorldCover/v200` (**the within-ranch land-cover source**); optional `ESA/WorldCover/v100` | `hansen_gfc_v1_11_envelope` (5 bands), `worldcover_2021_envelope`, optional `worldcover_2020_v100_envelope` | static | 30 m / 10 m | static |
| `09_pmlv2_amaka_export.js` | `CAS/IGSNRR/PML/V2_v017` (same as the user's local PML scripts) | `pml_{gpp,ec,es,ei}_meanrate` | **mean of daily rates** (totals = rate x days, in R); calendar-year seasons | 500 m native | 2001-2020 (clamped to full years) |
| `10_fire_history.js` | `MODIS/061/MCD64A1`; `FIRMS` | `burn_count`, `most_recent_burn_year`, `burned_by_year`, `firms_detections_by_year` | max within year, sum across years | 500 m / 1000 m | 2001-2024 |

Sentinel-2 surface reflectance starts globally in late 2018; MODIS LST / VI start Feb 2000 and SMAP L4 Mar 2015, so
the first complete year is used. The scripts clamp the requested years to full-coverage years automatically.

Grid handling in R follows brief 9.2: every raster stays on its native grid; analysis tables come from point-sampling
at RPI pixel centroids, and any resampling is done once, onto a documented target grid.
