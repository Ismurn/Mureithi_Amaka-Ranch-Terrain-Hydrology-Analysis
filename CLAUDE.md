# Amaka Ranch — Terrain, Hydrology & Vegetation Analysis

**Repo:** `Ismurn/Mureithi_Amaka-Ranch-Terrain-Hydrology-Analysis`
**Lead:** Isaac Mureithi (AWF / OPALS, Univ. of Exeter)
**Version:** v2 (26 October 2026) — supersedes v1 of the same date
**Status at handover:** Brief locked. Three existing pans (statuses known) and one desilting dam in progress; coordinates to be supplied. Code to be built from scratch.

**What changed from v1:** restructured around three operational use cases (diagnose / monitor / site); new seasonality layer across all time-varying datasets; new scripts for riparian analysis, fire history, and pan siting; expanded water-pans CSV schema; expanded recommendation zones (now A–E, with explicit trade-offs); new §13 "Literature grounding" citing the four papers that justify each analytical choice; expanded known-pitfalls section.

---

## 0. Purpose of this document

This is the handover brief for **Claude Code** to build the full R + Google Earth Engine (GEE) analytical pipeline for Amaka Ranch and Conservancy. The brief fixes the design decisions already made, lists the datasets and tools to use, specifies the repo structure, codifies the house conventions, flags known pitfalls from the two closest prior pipelines (Tsavo West Rhino Habitat Intelligence; Kitobo Forest Conservation Plan), and grounds every analytical choice in the peer-reviewed literature (§13). Everything downstream should refer back here.

---

## 1. Project objectives

The analysis serves **three distinct operational use cases** at Amaka Ranch, each tied to a different stakeholder question:

1. **Diagnose existing water-pan performance.** Three pans currently exist on Amaka: two are degraded and do not hold water for long, one holds water but is under immense pressure in the dry season. The pipeline must characterise *why* each existing pan is in its current condition (contributing catchment area, upslope vegetation stress, sediment supply signals, accessibility-shed demand) so management interventions target the right failure mode.

2. **Monitor an intervention in progress.** A desilting dam is currently being worked on. The pipeline must establish a defensible pre-intervention baseline for its catchment (vegetation condition, land-cover fractions, Es/ET ratio, soil-moisture context) so a re-run at 6/12/24 months post-intervention produces before/after evidence suitable for donor reporting and adaptive management.

3. **Site new water pans for future fundraising.** The pipeline must identify and rank candidate locations for new pans, with explicit multi-criteria rationale (water-gap filling, large-herbivore comfort, HWC avoidance, waterless-buffer preservation, woody-sacrifice-zone avoidance), to feed directly into fundraising proposals written by AWF/TTWCA.

Underlying all three: characterising vegetation condition, hydrological structure, and large-herbivore dynamics across the ranch, with elephants as the dominant management lens because elephant movement between Tsavos and the ranches is the dominant large-mammal driver in this landscape, and water + cover are the two ranch-scale levers that matter most.

Amaka sits within the Kasigau Corridor REDD+ project area and the broader TTWCA space; the outputs are designed to be operationally useful to Amaka management, TTWCA, and AWF Conservation Planning simultaneously.

---

## 2. Study area

**Amaka Ranch and Conservancy** — community-managed rangeland within the Tsavo Conservation Area, Taita-Taveta County, Kenya. Governance context: under the Kasigau Corridor REDD+ project area. Neighbouring units: other TTWCA ranches, Tsavo East and West NPs, and the Tsavo Settled Areas.

Isolate the AOI from `RanchesTTWCA.shp` using `R_NAME == "AMAKA RANCH"`, write it to `data/boundary/Amaka_Ranch.shp`, and treat that file as the canonical AOI. Verify CRS before anything else (see §9.1).

For hydrology, do **not** clip the DEM to Amaka alone — clip to a **buffered catchment envelope** of 5–10 km around the ranch so upslope contributing area is retained.

**Season definitions (used throughout the pipeline):**

- **Long rains (LR):** March–May
- **Short rains (SR):** October–December
- **Dry season (DRY):** June–September + January–February

All time-varying datasets below produce three seasonal aggregates in addition to annual, so that dry-season dynamics (the operational pressure point per Valls-Fox 2018 and Masiaine 2021) are not washed out by annual means.

---

## 3. Data inventory

### 3.1 Local inputs Isaac already holds

Paths below are placeholders — Isaac will edit `R/00_setup.R` with the real local paths at run time. Treat these as "sitting on OneDrive".

| Dataset | Format | Current location (indicative) |
|---|---|---|
| Kenya regional admin boundary | `.shp` (WGS84) | `TCA_GAE-Shapefiles/Country_Borders/` |
| Tsavo East / West NPs (WDPA) | `.shp` | `TCA_GAE-Shapefiles/WDPA_Amboseli_Tsavo/National Parks&Reserves/` |
| Tsavo ranches (`RanchesTTWCA.shp`) | `.shp` | `TCA_GAE-Shapefiles/Tsavo_Conservation_Areas/` |
| Tsavo settled areas | `.shp` | `TCA_GAE-Shapefiles/Tsavo_Conservation_Areas/` |
| **RPI v2** (actual GPP, potential GPP, RPI, temporal performance) | 4 × `.tif` | OneDrive — path to be set in `00_setup.R` |
| **PML_V2 Kenya** (GPP, Ec, Es, Ei) | multiband `.tif` | OneDrive — path to be set in `00_setup.R` |
| **Study pixels Kenya** (IGBP class + tree/shrub/grass/bare/rangeland fractions) | multiband `.tif` | OneDrive — path to be set in `00_setup.R` |

All of these must be loadable through `R/00_setup.R` so no other script hard-codes a Windows path.

### 3.2 To be acquired via GEE

Export to Drive and ingest from `data/gee_exports/`. Each time-varying dataset exports **four aggregates** (annual + LR + SR + DRY) unless noted.

| Dataset | GEE asset | Role |
|---|---|---|
| **DEM** — Copernicus GLO-30 (primary) | `COPERNICUS/DEM/GLO30` | Hydrology base. Fallback: SRTM (`USGS/SRTMGL1_003`) or FABDEM (`projects/sat-io/open-datasets/FABDEM`). Preferred because it is bare-earth-corrected relative to SRTM in vegetated terrain. (Not time-varying — single export.) |
| **CHIRPS daily rainfall** | `UCSB-CHG/CHIRPS/DAILY` | Climate context, drought flagging (SPI / rainfall anomaly). 2000–2024 to match RPI window. Seasonal sums. |
| **MODIS LST (MOD11A2, 8-day, 1 km)** | `MODIS/061/MOD11A2` | Surface thermal stress, day/night. Seasonal means. |
| **MODIS NDVI/EVI (MOD13Q1, 250 m)** | `MODIS/061/MOD13Q1` | Vegetation phenology, cross-check for RPI/PML. Seasonal means. |
| **Sentinel-2 NDVI/NDWI (10 m, 2015–)** | `COPERNICUS/S2_SR_HARMONIZED` | Fine-scale greenness and water presence; current state of pans and riparian strips. Seasonal means. |
| **SMAP L4 soil moisture (9 km)** | `NASA/SMAP/SPL4SMGP/008` | Root-zone moisture proxy. Seasonal means. |
| **ERA5-Land** | `ECMWF/ERA5_LAND/DAILY_AGGR` | For reference ET (PET) and air T, if needed. Seasonal means. |
| **Hansen Global Forest Change** | `UMD/hansen/global_forest_change_2023_v1_11` | Woody cover baseline + loss, for elephant-impact signature and sterile-landscape diagnostic. (Not time-varying in the seasonal sense — annual loss layer.) |
| **ESA WorldCover 2021 (10 m)** | `ESA/WorldCover/v200/2021` | Cross-check against study pixels. (Single epoch.) |
| **MODIS burned area (MCD64A1, 500 m)** | `MODIS/061/MCD64A1` | Fire history for burn-count-per-pixel and most-recent-burn-year; separates fire-driven from elephant-driven woody loss (Jones 2025). |
| **NASA FIRMS active fire** | `FIRMS` | Finer-dated active fire points, 2001–present. |

All GEE scripts live under the account `ee-imnduta1`, saved in the Code Editor at `users/imnduta1/Ismurn/Amaka/`.

### 3.3 To be supplied by Isaac — `data/water_pans/water_pans.csv`

Expanded schema (v2):

| Column | Required | Values / notes |
|---|---|---|
| `pan_id` | yes | Short unique string (e.g. `AMK_P01`) |
| `name` | yes | Human-readable name used by ranch management |
| `longitude` | yes | WGS84 decimal degrees |
| `latitude` | yes | WGS84 decimal degrees |
| `type` | yes | `natural_pan` / `earth_dam` / `borehole` / `other` |
| `status` | **yes** | One of: `existing_functional_under_pressure` / `existing_degraded` / `existing_desilting` / `candidate_new` / `candidate_baseline_pre_intervention` |
| `season_functional` | optional | `perennial` / `seasonal_lr` / `seasonal_sr` / `rarely` / `never_recently` |
| `notes` | optional | Free text — age, construction type, known issues, management history |

The pipeline must tolerate the file being empty on first run and re-ingest it in a dedicated stage (§6.6). Candidate rows marked `candidate_new` are generated by the siting script (§6.9), not supplied by Isaac — Isaac supplies the three/four currently known existing pans with their status flags.

---

## 4. Tools & environment

### 4.1 R stack (Isaac's house conventions — do not deviate)

- Package management: `pacman::p_load(...)` only. No `library()` chains.
- Core spatial: `sf`, `terra`, `tidyverse`, `tmap`.
- Zonal statistics: `exactextractr::exact_extract()`.
- Trend tests: `trend::sens.slope()` (Theil-Sen).
- Visualisation: `scico`, `viridis`, `patchwork`.
- Hydrology: `whitebox` R package (WhiteboxTools). Call `whitebox::install_whitebox()` once in setup; wrap `wbt_init()` in `00_setup.R`.
- Naming: dot notation for objects (`Amaka.Boundary`, `RPI.Clip`), snake_case for file names.
- Code style: **heavily commented, clear phase structure**. Each script opens with a block comment stating purpose, inputs, outputs, and runtime notes.

### 4.2 GEE

- Account: `ee-imnduta1` (legacy `users/imnduta1`).
- Scripts saved under `users/imnduta1/Ismurn/Amaka/`.
- Export targets: Google Drive folder `Amaka_GEE_Exports/` → ingested from `data/gee_exports/`.
- Default export CRS: `EPSG:32737` (WGS84 / UTM Zone 37S) at 30 m scale for raster outputs. 10 m for Sentinel-2.
- **Seasonality handled server-side.** Each time-varying GEE script defines LR / SR / DRY date filters at the top, loops over years, and exports one annual plus three seasonal composites per year. Reducer choice (sum vs mean) is documented in each script's header.

### 4.3 Hydrology

- Primary engine: **WhiteboxTools** via the `whitebox` R package. Do not fall back to `terra::terrain()` for flow routing — WhiteboxTools' breach/fill routines are what we rely on.
- Buffered catchment: 5–10 km around Amaka boundary (parameterise, default 7 km).

---

## 5. Repository structure (target)

```
Mureithi_Amaka-Ranch-Terrain-Hydrology-Analysis/
├── README.md                         # Project overview + how to run
├── CLAUDE.md                         # This brief (v2)
├── .gitignore                        # Excludes data/, outputs/, *.tif, *.Rproj.user
├── R/
│   ├── 00_setup.R                    # pacman::p_load, local paths, wbt_init, season date helpers
│   ├── 01_boundary_prep.R            # Extract Amaka, verify CRS, write Amaka_Ranch.shp, build buffered envelope
│   ├── 02_dem_prep.R                 # Ingest Copernicus DEM, clip to buffered catchment
│   ├── 03_hydrology.R                # WhiteboxTools: breach, D8, flow accum, streams, TWI, sub-catchments
│   ├── 04_rpi_analysis.R             # RPI zonal + spatial pattern + Theil-Sen trends inside Amaka
│   ├── 05_pmlv2_ecohydrology.R       # GPP, Ec, Es, Ei; seasonal ratios (Ec/ET, Es/ET); WUE
│   ├── 06_landcover_studypixels.R    # IGBP + tree/shrub/grass/bare fractions; shrub-to-tree ratio
│   ├── 07_riparian_analysis.R        # (NEW) 100–150 m buffers on streams of order ≥ 2; riparian condition
│   ├── 08_water_pans_analysis.R      # Pan pour-points → sub-catchments → per-pan diagnostics (expanded)
│   ├── 09_fire_history.R             # (NEW) Burn-count-per-pixel, most-recent-burn-year, fire × woody loss
│   ├── 10_integrated_stress.R        # Composite stress index; overlay with hydrology; recommendation zones A–E
│   ├── 11_pan_siting.R               # (NEW) Multi-criteria suitability + ranked candidate pan sites
│   └── 12_maps_figures.R             # tmap-based figure factory for the report
├── gee/
│   ├── 01_amaka_boundary_export.js
│   ├── 02_dem_amaka_export.js
│   ├── 03_chirps_rainfall.js         # Seasonal
│   ├── 04_modis_lst.js               # Seasonal
│   ├── 05_modis_vi.js                # Seasonal
│   ├── 06_sentinel2_vi_ndwi.js       # Seasonal
│   ├── 07_smap_soilmoisture.js       # Seasonal
│   ├── 08_hansen_woodycover.js
│   ├── 09_pmlv2_amaka_export.js      # Optional re-export at matching temporal reducer, seasonal
│   └── 10_fire_history.js            # (NEW) MCD64A1 burned area + FIRMS points
├── data/
│   ├── boundary/                     # Amaka_Ranch.shp + context (NPs, settled areas)
│   ├── rpi/                          # (gitignored — points to OneDrive)
│   ├── pml_v2/                       # (gitignored)
│   ├── study_pixels/                 # (gitignored)
│   ├── dem/                          # Copernicus DEM clipped to buffered catchment
│   ├── gee_exports/                  # Seasonal composites (gitignored, except small ones)
│   └── water_pans/
│       └── water_pans.csv            # Populated by Isaac with 3 existing + 1 desilting; candidates added by §6.9
├── outputs/
│   ├── rasters/                      # Derived rasters (TWI, stress index, suitability, etc.)
│   ├── figures/                      # PNG/PDF
│   ├── tables/                       # CSV / XLSX (incl. pan diagnostics, candidate sites)
│   └── reports/                      # Knitted RMarkdown
├── docs/
│   ├── methodology.md                # Narrative method, with equations and citations
│   ├── data_sources.md               # Catalogue of every dataset used
│   ├── literature_grounding.md       # (NEW) Paper-by-paper mapping to analytical choices (mirrors §13)
│   └── known_issues.md               # Running log of bugs hit and fixes applied
└── analysis/
    └── Amaka_Analysis.Rmd            # Narrative notebook that sources R/ scripts
```

Keep raw and derived data out of git — only code, docs, and small CSVs (like `water_pans.csv`) live under version control.

---

## 6. Analytical pipeline (ordered)

Each sub-section corresponds to one script in `R/`. Scripts must be independently runnable given that `R/00_setup.R` has been sourced.

### 6.1 Boundary preparation (`01_boundary_prep.R`)

- Load `RanchesTTWCA.shp`.
- Filter `R_NAME == "AMAKA RANCH"` (keep exact string case-insensitive via `toupper()`).
- Verify CRS with `st_crs()`. The Tsavo shapefiles have historically come in both WGS84 and Arc 1960 / UTM 37S (EPSG:21037). **Confirm numerically** by comparing a known distance against Google Earth; see §9.1.
- Reproject to `EPSG:4326` for mapping and to `EPSG:32737` (WGS84/UTM 37S) for metric operations. Save both.
- Build the **buffered catchment envelope**: `st_buffer(Amaka, dist = 7000)` in metric CRS, then `st_bbox()` → the DEM clipping extent.
- Also load TsavoE, TsavoW, settled areas, and Tsavo ranches for context layers in maps.

### 6.2 DEM preparation (`02_dem_prep.R`)

- Ingest Copernicus GLO-30 clipped to the buffered catchment (from `gee/02_dem_amaka_export.js`).
- Reproject to `EPSG:32737` if not already.
- Mild smoothing optional: `wbt_gaussian_filter()` with `sigma = 1.0` — document whether used.
- Export: `data/dem/amaka_dem_utm37s.tif`.

### 6.3 Hydrology (`03_hydrology.R`)

Order of operations in WhiteboxTools (all paths are temp files inside `outputs/rasters/hydrology/`):

1. **Breach depressions** — `wbt_breach_depressions_least_cost()`. Prefer breach over fill in arid landscapes so real pans are not artificially flooded.
2. **D8 flow pointer** — `wbt_d8_pointer()`.
3. **D-infinity flow pointer** — `wbt_d_inf_pointer()` (for TWI).
4. **D8 flow accumulation** — `wbt_d8_flow_accumulation()`.
5. **Stream extraction** — `wbt_extract_streams()` with threshold calibrated against Sentinel-2 imagery and known seasonal rivers. Start with `threshold = 1000` cells and iterate.
6. **Strahler stream order** — `wbt_strahler_stream_order()`.
7. **Streams as vector** — `wbt_raster_streams_to_vector()`.
8. **Specific contributing area** — `wbt_d_inf_flow_accumulation(out_type = "SpecificContributingArea")`.
9. **Slope** — `wbt_slope()` in degrees.
10. **TWI** — `wbt_wetness_index(sca, slope)`.
11. **TPI (relative topographic position)** — `wbt_relative_topographic_position()` at a 500 m window.
12. **Hillshade** for base maps — `wbt_hillshade()`.

After hydrology, mask all derived rasters to the Amaka boundary for reporting (keep buffered versions in a `_buffered/` subfolder for diagnostics).

### 6.4 Vegetation condition (`04_rpi_analysis.R`, `05_pmlv2_ecohydrology.R`, `06_landcover_studypixels.R`)

**RPI (`04`):**
- Crop and mask the three RPI rasters and the temporal performance raster to Amaka.
- Zonal stats inside Amaka: mean, SD, 10th/50th/90th percentile — annual.
- Pixel-wise Theil-Sen slope on annual RPI (or on the actual-GPP stack if the per-year layers are available): wrap `trend::sens.slope()` via `terra::app()`.
- Flag "degrading" pixels: significant negative Sen slope at p < 0.1.
- Stratify by land-cover class (from study pixels) and by TWI quintile.

**PML_V2 (`05`):**
- Load the PML_V2 Kenya stack and crop to Amaka.
- Derive **three seasonal sums** (LR, SR, DRY) for GPP, Ec, Es, Ei, plus annual. Confirm the time step of the stack first — units matter (daily vs 8-day vs monthly).
- Compute ET = Ec + Es + Ei for each season.
- Derive ratios per season:
  - `WUE_season = GPP_season / ET_season` (water-use efficiency)
  - `Ec_frac_season = Ec_season / ET_season` (plant-driven water flux share — low in degraded land)
  - `Es_frac_season = Es_season / ET_season` (bare-soil evaporation share — **high in degraded land; this is a key stress signal, strongest in DRY**)
- Zonal + maps per season.

**Study pixels / land cover (`06`):**
- Crop to Amaka. Report fractions tree / shrub / grass / bare / rangeland.
- **Derive shrub-to-tree ratio** as a new indicator layer (McCleery 2025, Jones 2025 — shrub-dominance signals active megaherbivore impact history or sterile-landscape conditions).
- Cross-tab IGBP class × TWI quintile × RPI quintile.
- Flag "bare-and-wet" pixels (high TWI + high bare %) — prime regen-feasibility candidates for Zone A.
- **Decision (2026-10-26):** within-ranch fractions, the shrub-to-tree ratio and the bare-and-wet flag are derived from 10 m WorldCover + Sentinel-2, not from the 500 m study pixels (which are ~97% shrubland over Amaka). The study pixels remain the Kenya-wide baseline. See §10 and `docs/known_issues.md` D5.

### 6.5 Riparian buffer analysis (`07_riparian_analysis.R`) — NEW

- Buffer all extracted streams of Strahler order ≥ 2 by 100–150 m (parameterise, default 125 m) → riparian polygons.
- Compute per-polygon zonal stats on all vegetation layers (RPI, PML seasonal, land cover, shrub-to-tree ratio).
- Compute riparian-vs-interior condition contrast: for each metric, express the riparian value as a difference from the ranch-wide median.
- Rationale (Masiaine 2021): rivers are the dominant baseline-condition visitation attractor for large mammals; riparian corridor condition is likely the single most informative vegetation signal on the ranch and must be reported separately.

### 6.6 Water-pans analysis (`08_water_pans_analysis.R`) — EXPANDED

**This stage expects `data/water_pans/water_pans.csv` to exist and have at least one row. If empty, the script writes a `NO_PANS_YET` marker to `outputs/tables/` and exits cleanly.**

For **each pan, regardless of status**:

1. Snap pan point to the nearest stream pixel (`wbt_jenson_snap_pour_points()` or `wbt_snap_pour_points()`). Use a snap distance of 100 m — document choice.
2. Delineate contributing catchment with `wbt_watershed()`.
3. For that catchment, compute per season (LR / SR / DRY / annual):
   - Area (ha), mean slope, mean TWI, max flow accumulation.
   - Mean RPI, Sen slope on RPI (is it degrading?).
   - Mean ET, Es/ET (bare-soil stress), WUE.
   - Mean soil moisture (SMAP), rainfall sum (CHIRPS).
   - Land-cover fractions and shrub-to-tree ratio.
   - Distance from pan to Amaka boundary, nearest Settled Area, and nearest neighbour pan.
4. Compute **accessibility shed at four radii** (3, 5, 10, 15 km) — the area within elephant daily range of this pan. Union across all functional pans gives the ranch-wide "watered area"; the complement is the ranch's **water-gap** (input to §6.9).
5. Compute **refugia value** per pan = f(distance from ranch boundary, distance from Settled Areas, distance from main through-roads). Interior pans distant from boundaries function as refugia during pastoral stress events (Masiaine 2021).

Then **status-specific diagnostics**:

- **For `existing_degraded` pans (the two that do not hold water):** report the hypotheses in rank order based on EO evidence:
  - *Catchment too small to fill:* contributing area below some threshold proportional to CHIRPS annual sum.
  - *Siltation from upslope overgrazing:* high bare fraction + high Es/ET + declining RPI in the upslope catchment.
  - *Flashy runoff from degraded vegetation:* low RPI + low TWI + steep slopes in the catchment.
  - *Terrain geometry wrong for retention:* mean slope too high or too low; TPI positive (ridge) rather than negative (depression).
- **For `existing_functional_under_pressure` pan (the one that holds):** quantify the demand it is absorbing — how many other functional pans sit within its 5 km, 10 km, 15 km shed? If it is the only functional pan within elephant range, the "immense dry-season pressure" is explained quantitatively (this one pan is doing the work of several).
- **For `existing_desilting` pan:** record the full per-season baseline table. Write it to `outputs/tables/desilting_baseline_PRE_INTERVENTION_<date>.csv` and flag it clearly for a 6/12/24 month re-run.

Outputs: `outputs/tables/pan_catchment_metrics.csv`, `outputs/tables/pan_diagnostics.csv`, `outputs/rasters/pan_catchments.tif`, `outputs/rasters/water_gap.tif`, map per pan + a landscape overview.

### 6.7 Fire history (`09_fire_history.R`) — NEW

- Ingest MCD64A1 burned-area and FIRMS active-fire points from `gee/10_fire_history.js`.
- Produce two rasters at 500 m: `burn_count_2001_2024.tif` (number of times burned) and `most_recent_burn_year.tif`.
- Overlay with Hansen tree-cover loss: cross-tab "burned-and-lost" vs "unburned-and-lost" to separate fire-driven from likely elephant-driven woody loss (Jones 2025: fire-megaherbivore synergy amplifies tree mortality; the two signatures must be distinguished).
- Report per pan catchment (feeds into §6.8).

### 6.8 Integrated stress & recommendations (`10_integrated_stress.R`) — EXPANDED

**Composite stress index** per pixel (and summarised per pan catchment, per season) — all terms min-max scaled to [0, 1]:

```
StressIndex =  0.30 * (1 - RPI_now)
             + 0.20 * is_degrading_RPI          # 0/1 from Sen slope
             + 0.20 * Es_frac_DRY               # dry-season bare-soil ET share
             + 0.15 * bare_fraction             # from study pixels
             + 0.15 * (1 - normalised_TWI)      # water-poor terrain
```

Compute **three stress surfaces**: annual, dry-season, long-rains. The dry-season surface is the one that drives management recommendations; the long-rains surface serves as the "best-case baseline" against which to measure stress amplification.

Weights are a starting point — expose them as parameters. Record the version used in each run inside the output filename (e.g. `stress_DRY_v1_030_020_020_015_015.tif`).

**Recommendation zones (overlay layer, A–E):**

- **Zone A — Regeneration feasibility.** High TWI × high bare % × non-degrading RPI (water present, substrate stressed but trajectory recoverable). **Reframed from v1:** the zone identifies where regeneration is physically *possible*; actualising regeneration requires ecological interventions outside EO scope (seed addition, protection from livestock, grass biomass reduction). Jones 2025 shows recruitment is bottlenecked by cryptic herbivores, grass competition, and abiotic factors regardless of large-herbivore pressure.
- **Zone B — Water-catchment interventions** (bunds, sand dams, gully plugs): low-order streams with high upslope stress + proximity to a pan.
- **Zone C — Pan protection / enclosure**: pan catchments in the top quartile of StressIndex and with a negative RPI Sen slope. **Trade-off pairing (new):** every Zone C recommendation carries a counter-metric reporting the expected effect on small-mammal regulation and tree recruitment, per McCleery 2025 (megaherbivore exclusion releases cryptic herbivores that suppress recruitment). The output never prescribes exclusion as unambiguously positive.
- **Zone D — Large-herbivore comfort** (renamed from "Elephant comfort" per McCleery 2025 — shrub cover and water drive mesoherbivore activity as well, elephant is the dominant case but not the only one). Split into two layers:
  - **D-wet:** dispersed, pan-independent large-herbivore use across the ranch during wet seasons.
  - **D-dry:** concentrated-around-pans pressure during the dry season — the operationally urgent layer, driving HWC and out-of-ranch movement risk. Criteria: distance-to-functional-pan < 5 km AND shrub + tree cover > 25% AND TWI above median.
- **Zone E — Waterless retention / heterogeneity buffer** (NEW, from Valls-Fox 2018 and Jones 2025). Areas intentionally far from any pan, with intact woody cover and acceptable RPI despite isolation. These function as reserves of mature tree structure, seed sources, and habitat for low-density sensitive grazers; placing a new pan here would eliminate the exact buffer the ecology needs. **Zone E is an exclusion criterion for the siting step (§6.9).**

**Per Zone C recommendation, report a "woody sacrifice zone" counter-metric:** expected Hansen loss within the pan's elephant accessibility shed, so pan placement decisions are made with the woody-cover cost visible.

### 6.9 Pan siting (`11_pan_siting.R`) — NEW

**Multi-criteria suitability analysis producing fundraising-ready candidate pan sites.**

Suitability surface = weighted sum of (all layers min-max scaled to [0, 1]):

```
Suitability =  0.25 * is_in_water_gap                 # outside existing accessibility sheds
             + 0.20 * terrain_feasibility             # TWI, flow accumulation, slope in viable range
             + 0.15 * is_in_Zone_D_dry                # would retain large herbivores inside Amaka
             + 0.10 * refugia_value                   # interior position, far from boundaries
             + 0.10 * normalised_shrub_tree_cover     # cover available for retained herbivores
             - 0.20 * is_in_Zone_E                    # EXCLUSION: do not place in waterless buffer
             - 0.10 * woody_sacrifice_risk            # avoid compounding elephant impact on mature trees
             - 0.10 * proximity_to_settled_areas      # avoid drawing elephants into HWC zones
```

Weights are v1 defaults — expose as function arguments.

From the suitability surface, extract the **top N candidate sites** (default N = 8) using a non-maximum-suppression algorithm with a minimum inter-site spacing of 3 km. For each candidate, record:

- Coordinates (longitude/latitude, WGS84)
- Contributing catchment area (ha) and estimated annual water yield (from CHIRPS × runoff coefficient)
- Water-gap area that this site would newly cover (ha)
- Expected large-herbivore comfort gain (overlap with Zone D-dry)
- Woody sacrifice zone (expected mature-tree-area within elephant accessibility shed)
- Distance to nearest Settled Area and ranch boundary
- Short textual rationale ("Fills a 1,200 ha water gap on the north-east interior; low HWC exposure; moderate cover; does not conflict with Zone E")

Outputs:
- `outputs/rasters/pan_suitability.tif` — continuous suitability surface
- `outputs/tables/candidate_pan_sites.csv` — ranked top N with full attribute table
- `outputs/figures/pan_siting_map.png` — the fundraising-ready map (existing pans with status + candidate sites ranked)
- Append candidate rows (`status = "candidate_new"`) to `data/water_pans/water_pans.csv` on an Isaac confirm.

### 6.10 Figures (`12_maps_figures.R`)

Figure factory using `tmap` + `patchwork`. Named outputs with stable file names so RMarkdown can embed them. Default palettes: `scico::scico(palette = "roma")` for RPI, `"davos"` for ET components, `"lajolla"` for stress, `"batlow"` for TWI, `"vik"` for suitability. Hillshade underlay on all terrain maps. Standard panel: dry-season vs wet-season side by side for every seasonal layer.

---

## 7. GEE scripts (`gee/`)

Each script:

- Opens with a header block: author, date, purpose, inputs, outputs, export settings, **season date filters, and reducer choice with justification**.
- Imports the Amaka boundary as an EE asset (uploaded once from `01_amaka_boundary_export.js`).
- Exports to Drive folder `Amaka_GEE_Exports/`.
- Uses `scale = 30` by default, `scale = 10` only for Sentinel-2.
- Uses `crs = "EPSG:32737"`.
- Clips to the buffered catchment envelope, not the ranch boundary (so downstream R code can decide).
- Seasonal scripts loop over year × season and export one composite per combination.

Priority order for Claude Code to implement:

1. `01_amaka_boundary_export.js` — upload Amaka as an EE FeatureCollection asset for all other scripts to import.
2. `02_dem_amaka_export.js` — needed to start hydrology.
3. `03_chirps_rainfall.js`, `05_modis_vi.js`, `06_sentinel2_vi_ndwi.js` — climate and vegetation context (seasonal).
4. `10_fire_history.js` — needed before §6.7.
5. `08_hansen_woodycover.js` — needed before §6.8.
6. Everything else when §6 needs it.

---

## 8. Deliverables

- Full repo pushed to `Ismurn/Mureithi_Amaka-Ranch-Terrain-Hydrology-Analysis` by Isaac (not by Claude Code).
- `analysis/Amaka_Analysis.Rmd` knitting to a Word document inside `outputs/reports/Amaka_Hydro_Veg_Report_v1.docx`.
- Key figures as standalone PNGs in `outputs/figures/`.
- Pan-catchment metrics table in `outputs/tables/pan_catchment_metrics.csv` and `.xlsx`.
- **Pan diagnostics table** `outputs/tables/pan_diagnostics.csv` — per-pan failure-mode or pressure-mode diagnosis.
- **Desilting pre-intervention baseline** `outputs/tables/desilting_baseline_PRE_INTERVENTION_<date>.csv` — flagged for 6/12/24 month re-run.
- **Candidate pan sites table** `outputs/tables/candidate_pan_sites.csv` and the siting map `outputs/figures/pan_siting_map.png` — the fundraising-ready outputs.
- `outputs/rasters/stress_DRY.tif`, `stress_LR.tif`, `stress_annual.tif`, `recommendation_zones.tif` (categorical: A/B/C/D-wet/D-dry/E), `pan_suitability.tif`, `water_gap.tif`.

---

## 9. Known pitfalls (carry-overs from Tsavo West and Kitobo work + new entries from the literature grounding)

### 9.1 CRS verification is non-optional

Tsavo shapefiles historically arrive in **Arc 1960 / UTM 37S (EPSG:21037)** despite WDPA files being WGS84. The Kitobo work hit this exact bug — survey points assumed to be WGS84 were actually Arc 1960, giving a ~250 m offset until caught.

In `01_boundary_prep.R`:

```r
# Sanity check: compare a known pair of corner points' metric distance
# with Google Earth's measurement for the same two points. A metric offset
# > ~50 m on a ~1 km baseline means the CRS tag is wrong — do NOT just
# st_set_crs(), reproject deliberately.
```

### 9.2 Raster grid alignment

RPI, PML_V2, study pixels, Copernicus DEM, MODIS, CHIRPS, Sentinel-2, SMAP all have different native grids. **Do not blindly `terra::resample()` everything to a common grid** — reprojection artefacts corrupt zonal statistics.

Preferred approach (from the Tsavo West rhino pipeline):

- Keep each raster in its native CRS/grid.
- Use RPI pixel centroids as sample points (`terra::xyFromCell()`).
- Point-sample all other rasters at those centroids with `terra::extract()`.
- Build the analysis data frame from the sampled points.

For continuous maps and zonal summaries, resample **once** to a defined target grid (RPI's native grid is the natural choice), and document that choice.

### 9.3 `gridDist()` returning all zeros

If you compute distance-to-water using `terra::gridDist()`, make sure the input is a **proper binary raster** (`1` for water, `NA` for non-water — not `0`). The Tsavo rhino work lost half a day to this. Fix pattern:

```r
water.binary <- ifel(streams > 0, 1, NA)
dist.to.water <- gridDist(water.binary, target = 1)
```

### 9.4 PML_V2 temporal reducer

If re-exporting PML_V2 from GEE for Amaka, confirm the temporal reducer matches the units. The Tsavo rhino pipeline silently used the wrong reducer once and had to re-export. Document in the GEE script header which reducer was used and why.

### 9.5 `pivot_wider()` list columns

When tidying zonal outputs with multiple metrics per polygon, `pivot_wider()` can return list columns if names are not unique. Verify with `glimpse()` and unnest explicitly.

### 9.6 Breach vs fill

Arid landscapes with real pans need `wbt_breach_depressions_least_cost()`, not `wbt_fill_depressions()`. Filling will erase the pans themselves.

### 9.7 Stream threshold calibration

The flow-accumulation threshold for stream extraction is landscape-specific. Calibrate by overlaying extracted streams on Sentinel-2 and against known seasonal rivers (Tsavo River, Voi River tributaries) before accepting the output. Iterate on `threshold` in the WBT call.

### 9.8 EO is daytime-biased — misses elephant temporal niche shifts (NEW, Valls-Fox 2018)

Elephant bulls accommodate cattle presence primarily through *temporal* avoidance (nocturnal shifts) rather than spatial avoidance. Satellite overpasses are daytime-biased, so the stress index captures spatial responses but is systematically blind to temporal ones. Interpretation of Zone D for elephant bulls specifically should carry this caveat.

### 9.9 Stress index is a structural index — misses cryptic non-structural pathways (NEW, McCleery 2025)

McCleery 2025 demonstrated that non-structural pathways (behavioural facilitation, game trails, moribund grass removal, forb composition shifts, seed/insect availability, predator-prey behavioural dynamics) are more than twice as strong as vegetation-structure pathways in explaining mesoherbivore facilitation by megaherbivores. Our stress index is a vegetation-structure index. Any interpretation attributing causality to vegetation structure alone will underestimate megaherbivore influence. Flag in the methodology write-up when discussing recommendation zones.

### 9.10 No ground-truth exclosures on Amaka or TTWCA (NEW)

No calibration opportunity for the EO stress signal is available on-ranch. The pipeline proceeds on EO alone; findings should be presented with appropriate uncertainty, and recommendations should prefer actions that are robust to a reasonable range of EO signal error.

### 9.11 Human shield effect — distance-to-settlements is not monotonic for wildlife (NEW, Masiaine 2021)

Distance-to-Settled-Areas is used throughout the pipeline as a baseline-condition metric (settled areas repel wildlife under normal conditions). Under acute pastoralist-stress events the sign can reverse — wildlife may move *toward* human habitations to reduce exposure to armed poachers (Masiaine 2021). The layer's interpretation is condition-dependent, not absolute.

### 9.12 Mega- vs large-herbivore attribution care (NEW, Jones 2025 + McCleery 2025)

Our framing is elephant-centric because elephants dominate megaherbivore biomass in Tsavo. But Jones 2025 shows some species (S. birrea) are *facilitated* by large herbivores via cryptic-herbivore suppression. Avoid attributing all woody dynamics to elephants alone in interpretation; flag species-specific nuance when discussing recruitment and woody-cover change.

---

## 10. Decisions deferred to Isaac

These are flagged for Isaac to confirm before or during the Claude Code run. Defaults are listed; Claude Code should use the default unless Isaac has edited this section.

| Decision | Default | Notes |
|---|---|---|
| DEM source | Copernicus GLO-30 | Switch to FABDEM if Copernicus looks noisy in vegetated patches. |
| Catchment buffer radius | 7 km | 5 km if DEM processing is slow; 10 km if Amaka is small and sits high in its catchment. |
| Stream extraction threshold | start at 1000 cells, iterate | Calibrate against Sentinel-2. |
| Stress index weights | 0.30 / 0.20 / 0.20 / 0.15 / 0.15 | Reasonable starting point; expose as function arguments. |
| Pan snap distance | 100 m | May need to grow for pans sited on broad floodplains. |
| RPI time dimension | annual, 2000–2022 | Corrected 2026-10-26: `rpi_ea_v2.tif` has 23 layers (`rpi.2000`–`rpi.2022`), not 2000–2023. |
| **Season date filters** | LR = Mar–May, SR = Oct–Dec, DRY = Jun–Sep + Jan–Feb | Standard southern Kenya bimodal rainfall pattern. Confirm against Amaka management's own working seasons. |
| **Riparian buffer width** | 125 m | 100 m if riparian strips are narrow; 150 m where seasonal floodplains are wider. |
| **Pan accessibility shed radii** | 3, 5, 10, 15 km | Report all four; 10–15 km is the dry-season elephant daily-range reference. |
| **Pan suitability weights** (§6.9) | as listed in §6.9 | v1 defaults; sweep when first candidate maps are reviewed. |
| **Minimum inter-candidate spacing** | 3 km | Prevents candidate clustering; adjust if ranch geometry is tight. |
| **Number of candidate pan sites to return** | 8 | Report top 8; AWF/TTWCA selects for proposals. |
| **Shrub-to-tree ratio threshold for "sterile landscape" flag** | ratio > 5 | Shrub cover more than five times tree cover in a pan catchment with high Hansen loss. Adjust after first run. |
| **Amaka AOI** | **The TTWCA Musamiti polygon exactly as delivered (4,782 ha)** — decided 2026-10-27 (a hull-based AOI was tried and rejected) | The polygon is coarse: ~20% below the 14,815 acres TTWCA lists for Amaka, and 3 of the 7 field points (PT05–PT07) fall in the neighbouring Kambanga polygon; accepted, pans outside it are still analysed and reported as outside. Replace with Amaka's real boundary if supplied (`Cfg.AOI.Mode = "file"`). Overrides §2's "`R_NAME == "AMAKA RANCH"`" filter (the ranch is `MUSAMITI` in the file). See `docs/known_issues.md` F8, F9. |
| **Land-cover source for within-ranch structure** (bare fraction, shrub-to-tree ratio, Zone A) | **10 m ESA WorldCover v200 + Sentinel-2 (BSI, NDVI)** — decided 2026-10-26 | The 500 m study pixels are ~97% shrubland over Amaka (no pixel > 5% bare) and stay as the Kenya-wide baseline only. "High bare" / "high TWI" thresholds are relative to the ranch. See `docs/known_issues.md` D5. |

---

## 11. Roles

- **Isaac Mureithi** — scientific lead, boundary/data/pan-coordinate provider, repo owner (pushes to GitHub manually), decision-maker on §10, interface with Amaka management and AWF Conservation Planning.
- **Claude Code** — implements the pipeline in R and GEE according to this brief; writes all scripts, maintains `docs/known_issues.md` and `docs/literature_grounding.md`, and keeps the README current.
- **Claude (chat)** — authored this brief; available for mid-stream clarifications and interpretation of results.

---

## 12. First actions for Claude Code

1. Create the directory skeleton in §5. Stub every R and GEE file with its header block.
2. Write `R/00_setup.R` with `pacman::p_load(sf, terra, tidyverse, tmap, exactextractr, whitebox, trend, scico, viridis, patchwork, rmarkdown, janitor, lubridate)`, a config block for local paths (RPI, PML_V2, study pixels, OneDrive root), **season date helpers** (functions `season_of_date()`, `is_dry()`, `is_lr()`, `is_sr()`), and `whitebox::install_whitebox()` guarded by `if (!wbt_init(check_version = FALSE))`.
3. Implement `R/01_boundary_prep.R` end-to-end, including the CRS sanity check in §9.1.
4. Implement `gee/01_amaka_boundary_export.js` and `gee/02_dem_amaka_export.js` so the DEM is in hand before hydrology.
5. Pause and report back with the knitted intermediate outputs (CRS verification figure, DEM clip map) before moving to §6.3.
6. From then on, work through §6.3 → §6.9 in order, keeping `docs/known_issues.md` and `docs/literature_grounding.md` current as issues are hit and literature cross-references are implemented.
7. §6.9 (pan siting) is the last analytical script and must be run after all stress and recommendation zones are produced.

---

## 13. Literature grounding (NEW in v2)

Every substantive analytical choice in the brief is grounded in one or more of four recent peer-reviewed papers. The table below maps each paper to the §6 sections it informs.

| Paper | Journal, Year | Core finding(s) carried into the brief | Sections it informs |
|---|---|---|---|
| **Valls-Fox et al.** — "Water and cattle shape habitat selection by wild herbivores at the edge of a protected area" | *Animal Conservation*, 2018 | Dry season drives aggregation around water; buffalo-cattle overlap increases fivefold in dry season; elephant bulls use temporal avoidance (nocturnal), buffalo use spatial avoidance; "virtual fence" effect of cattle herding; waterless buffers as a management lever for landscape heterogeneity | §2 (season definitions); §6.4 (seasonal PML ratios); §6.6 (seasonal pan metrics, accessibility shed, refugia value); §6.8 (D-wet / D-dry split; Zone E); §6.9 (siting weights); §9.8 (temporal blind spot) |
| **Jones et al.** — "Mega- and large-herbivores influence survival but not recruitment rates of African savanna trees" | *Biological Conservation*, 2025 | Megaherbivores (elephants) limit tree *survival*; structural defences (spines) deter large herbivores but not elephants; recruitment is bottlenecked by cryptic herbivores, grass competition, and abiotic factors regardless of large-herbivore pressure; realised λ rarely exceeds 1 even with exclusion; water heterogeneity buffers tree impact | §6.4 (shrub-to-tree ratio as sterile-landscape precursor); §6.7 (fire × woody loss separation); §6.8 (Zone A reframed as "feasibility"; Zone C trade-off pairing; woody sacrifice zone counter-metric); §9.12 (attribution care) |
| **McCleery et al.** — "Megaherbivores facilitate large grazing herbivores and suppress small mammals through vegetation structure and cryptic pathways" | *Biological Conservation*, 2025 | Non-structural pathways >2× stronger than structural pathways; megaherbivores suppress small mammals (which eat seeds/seedlings) through cryptic routes; mesoherbivore activity responds positively to shrub cover and negatively to tall tree canopy; mesoherbivores cannot replace megaherbivores | §6.4 (shrub-to-tree ratio as derived indicator); §6.8 (Zone D renamed "large herbivore comfort"; Zone C trade-off pairing closing the Jones loop); §9.9 (structural-index limitation) |
| **Masiaine et al.** — "Landscape-level changes to large mammal space use in response to a pastoralist incursion" | *Ecological Indicators*, 2021 | Rivers are the dominant baseline visitation attractor; cattle competitively exclude grazers at water; elephants held river access until overwhelmed by 40,000 head; interior dams functioned as refugia; "human shield" effect reversed distance-to-settlement sign under acute stress | §6.5 (riparian buffer analysis — new script); §6.6 (refugia value metric); §9.11 (human shield / non-monotonic distance-to-settlements) |

**No calibration ground-truth is available** (no exclosures or controlled grazing plots on Amaka or neighbouring TTWCA ranches — §9.10). The literature grounding is therefore the primary defence of analytical choices in the eventual write-up.
