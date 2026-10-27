# Methodology

Narrative method for the Amaka Ranch analysis, one section per stage. Section numbers in brackets refer to `CLAUDE.md`;
`known_issues.md` records what the data did to the plan. Every stage has been run on the real data (2026-10-27).

## 1. Study area and area of interest [2, 6.1]

Amaka Ranch and Conservancy is taken from `RanchesTTWCA.shp` as the feature(s) with `R_NAME` in
`{"MUSAMITI", "AMAKA RANCH"}` (the ranch is catalogued under its former name). Names are compared after upper-casing,
trimming and collapsing repeated spaces. The feature is dissolved to a single valid MULTIPOLYGON and written to
`data/boundary/Amaka_Ranch.shp` (EPSG:32737, canonical) and `Amaka_Ranch_wgs84.shp` (EPSG:4326, mapping and GEE).

**The AOI is the TTWCA polygon exactly as delivered.** It is coarse: four corners, 4,782 ha against the 14,815 acres
(about 5,995 ha) TTWCA lists for Amaka, and 3 of the 7 points Isaac recorded on the ranch (PT05, PT06, PT07) fall in
the neighbouring Kambanga polygon (`docs/known_issues.md` F2, F8). A hull of the field points was tried and rejected
(F9); there is no better boundary available, so the polygon is used as is. Pans outside it are still analysed and are
reported as outside the AOI (F15). `Cfg.AOI.Mode` is `"musamiti"` (this polygon) or `"file"` (a boundary supplied by Amaka
management or TTWCA, the intended end state); every output carries the mode.

Seasons are defined on calendar months of the same year: long rains (LR) March-May, short rains (SR) October-December,
dry season (DRY) January-February plus June-September; ANNUAL is the calendar year.

## 2. Coordinate reference system check [9.1]

Metric operations use WGS 84 / UTM 37S (EPSG:32737). The ranch and settled-area shapefiles are tagged EPSG:32737, but
Tsavo files have previously been Arc 1960 / UTM 37S (EPSG:21037). Two checks are run on the coordinate *numbers*,
read under each hypothesis (H1 = file tag, H2 = Arc 1960):

1. **Numeric edge test.** Vertices of all ranch polygons lying within 1.5 km of the Tsavo East / West boundary
   (WDPA, WGS84) under either hypothesis form one common vertex set. For each hypothesis the distance of every vertex
   to the park boundary is computed; the median distance and the share within 100 m are compared. Where a ranch edge
   is also a park edge the correct hypothesis puts vertices on the park line. This is indicative only (the two sources
   are digitised independently); here it was inconclusive (C2).
2. **Visual overlay.** Amaka is written as a KML under H1 and under H2, with its four extreme vertices tabulated as
   latitude / longitude. The hypothesis whose outline lies on real features in Google Earth is correct. Position is
   tested rather than distance, because a datum shift is nearly a pure translation. Isaac's overlay confirmed H1 (C1).

The assumption actually used is `Cfg.Ranches.CRS.Assumed`; the coordinates are **relabelled** (never silently
`st_set_crs`-ed) in `read_tsavo_layer()` and then reprojected.

## 3. Catchment envelope and DEM [4.3, 6.1, 6.2]

Hydrology needs the upslope contributing area, so the DEM is never clipped to the ranch. The Amaka polygon is buffered by
`Cfg.Buffer.m` (default 7 km) in EPSG:32737 and the **bounding box of the buffer** is the clip extent. GEE exports use the
same buffer plus a 1 km pad; `02_dem_prep.R` then crops to the exact envelope and checks that the DEM covers it.

DEM source is Copernicus GLO-30 (FABDEM or SRTM selectable). GLO-30 is a surface model; this matters for slope and flow
routing in woody patches. It is compared with the phone altitudes at the field points (RMSE 3.8 m, `02_dem_vs_field_elevations.csv`).

## 4. Land-cover source [6.4, decision 2026-10-26]

The 500 m "study pixels" (MODIS IGBP rangeland classes AND ESA WorldCover v100 rangeland fraction >= 0.9) are the
Kenya-wide baseline, but over Amaka they are 96.6% shrubland, 2.9% grassland, 0.3% trees and 0.2% bare, with no pixel above
5% bare. They cannot support within-ranch metrics. Within-ranch structure therefore comes from the 10 m WorldCover v200 map
(class fractions computed onto the RPI native grid and onto pan catchments with `exact_extract`, never by resampling) and
from Sentinel-2 (a bare-soil index and NDVI as continuous proxies). The 10 m map shows the same picture (94.0% shrubland,
0% bare), so BSI is the bare term and "high" is always relative to the ranch.

## 5. Analysis grids and resolution [9.2]

Two grids are used. The **analysis grid** is the RPI native grid (0.004167 degrees, ~460 m) over the envelope: 49 x 45 =
2,205 cells, each with a `cell_id` that keys every cell-level table (`outputs/tables/cells_*.csv`). A cell belongs to the AOI
(`in_aoi`) when at least half of it lies inside (222 cells). Fine rasters (30 m, 10 m) are averaged over each cell polygon
with area-weighted zonal means; coarse rasters are point-sampled at the cell centroid. No raster is resampled to
make another fit. The **30 m DEM grid** carries the hydrology, the stress surface, the recommendation zones and the siting;
there the coarse layers (RPI, PML) are *replicated* with nearest-neighbour (no interpolation) and the fine layers keep their
detail, and cell tables are area-weighted means of the 30 m surfaces.

## 6. Hydrology [6.3, 9.6, 9.7]

`wbt_breach_depressions_least_cost` (distance 100 cells, remaining sinks filled) conditions the DEM, because filling would erase
real pans. D8 pointer and flow accumulation (in cells) and D-infinity specific catchment area follow. Slope is in degrees; the
topographic wetness index is TWI = ln(SCA / tan(slope)); TPI is the elevation minus the mean of a 500 m window (17 x 17
cells; `wbt_diff_from_mean_elev`), and RTP the related relative topographic position. Slope, TWI and TPI come from the original DEM; the
breached DEM is used only for routing. Streams are cells with accumulation >= `Cfg.Stream.Threshold`
(2,000 cells = 1.8 km2), chosen from a calibration over 250-5,000 cells that scores (a) how many of the four dams lie within 100 m
of a stream, (b) the dry-season NDVI contrast between 100 m corridors and a ring beyond them, and (c) drainage density
(`03_stream_threshold_calibration.csv`; known_issues H1). Strahler order and a vector network are derived from the same raster.

## 7. Vegetation condition: RPI [6.4]

For each year, the area-weighted mean, SD and 10th / 50th / 90th percentile of RPI over the AOI cells. The trend of the AOI
mean is the Theil-Sen slope (`trend::sens.slope`) with a Mann-Kendall test; per-cell Sen slopes and p-values mark **degrading**
cells (negative slope, p < `Cfg.Degrade.P` = 0.1) and improving cells. `RPI_now` is the mean of 2020-2022. RPI and actual GPP are
cross-checked against MODIS NDVI and EVI, per year and across cells, because RPI is a ratio of actual to modelled potential GPP and
need not track greenness.

## 8. Seasonal ecohydrology and climate context [6.4, 2]

PML-V2 comes from `gee/09` as the **mean of daily rates** over each calendar-year season; totals are rate x days (LR 92, SR 92,
DRY 181, ANNUAL 365). ET = Ec + Es + Ei; Es/ET is the soil-evaporation share (the "bare-soil" ratio of the brief) and
WUE = GPP / ET. The local 80-layer annual stack (sum of 8-day rates over Sep-Aug) is used only to check scale and reducer against the
GEE totals (ratios 0.954-1.014). Climate context series for each season: CHIRPS rainfall with the standardised precipitation index
(gamma fit, `MASS::fitdistr`; drought SPI <= -1, wet SPI >= 1), MODIS day LST and NDVI / EVI, SMAP soil moisture; Sen slopes, and
Pearson / Spearman correlations of PML GPP, ET and Es/ET with rainfall.

## 9. Land cover and riparian corridors [6.4, 6.5]

Land cover: WorldCover class fractions (`exact_extract`, `"frac"`) per cell and per catchment; the shrub-to-tree ratio is reported only
where trees are at least 1% of the area (218 of 222 AOI cells are not); Sentinel-2 BSI, NDVI and NDWI as season means.
**Bare-and-wet** cells are top-quartile TWI and top-quartile dry-season BSI (`Cfg.Relative.Q`). Riparian: streams of Strahler
order >= 2 buffered by `Cfg.Riparian.m` (125 m) and clipped to the AOI; for each metric the mean in each riparian set minus the ranch
median and minus the interior mean, and the difference in ranch standard deviations (an effect size, not a p-value, since neighbouring
pixels are not independent). Only 10-30 m layers can resolve a 250 m strip; RPI and PML are included and flagged. A distance
gradient uses 12,000 regular points and distance classes from the nearest order >= 2 stream.

## 10. Water pans [6.6]

`read_pans()` joins `water_pans.csv` to the field file on `field_id` and takes the **image-read centre** as the start where there is one,
else the recorded coordinates. **Stream-fed pans** (`Cfg.Pan.StreamFed.Types`: earth dams, natural pans) are snapped to the stream network
(Jenson method, `Cfg.Snap.Img.m` = 60 m for image-read centres, else the stand-off uncertainty; widened in steps to 300 m if no stream is
reached) and their catchment delineated with `wbt_watershed` on the D8 pointer. Other types get a `Cfg.Pan.Context.m` (500 m) context disc
(known_issues F12). Per catchment and season: area, slope, TWI, flow accumulation, RPI and Sen slope, GPP, ET, Es/ET, WUE, SMAP soil
moisture, CHIRPS rainfall, day LST, land-cover fractions, BSI, and an annual yield = area x rainfall x `Cfg.Runoff.Coef` (assumed 0.10).
Distances to the boundary (negative outside), Settled Areas, Tsavo NP and the nearest pan; **sheds** at 3 / 5 / 10 / 15 km; the **water gap** =
AOI cells farther than r km from every functional pan (primary radius `Cfg.WaterGap.Radius.km` = 5); **refugia value** = mean of
min(d / 3 km, 1) for the boundary and min(d / 15 km, 1) for settled areas (no road layer).

**Failure-mode hypotheses** for degraded stream-fed pans: H1 catchment too small = 1 - area / area of the working dam; H2 siltation =
mean of CDFs for dry-season BSI, dry-season Es/ET and RPI decline; H3 flashy runoff = mean of CDFs for low RPI, low TWI and steep slope;
H4 terrain geometry = mean of the slope-viability penalty (`Cfg.Terrain.Slope.Deg` trapezoid) and the CDF of TPI at the site. Each CDF is
the normal CDF of (catchment mean - ranch mean) / SD of ranch means over windows the size of the catchment, so 0.5 is typical. The working
dam is scored by the same tests as a control. The working dam's **demand** is the number of other functional and non-functional pans within
5 / 10 / 15 km and its shed share of the AOI. The desilting dam's baseline is written only when it has coordinates.

## 11. Fire history [6.7]

MCD64A1 burn count and most recent burn year (500 m), burned-by-year (24 layers), and FIRMS detections (1 km). Each 30 m Hansen pixel takes
the burn history of its 500 m cell. A mature-tree pixel (treecover2000 >= `Cfg.Tree.Mature.Pct` = 10) lost in year L is **burned-and-lost**
when its cell burned in L or in the `Cfg.Fire.Loss.Window.Years` (1) year before. Fire association is tested with every mature-tree pixel-year
(2001 to its loss year, or 2023): exposed when the cell burned in that year or the year before; the loss rate per exposed against unexposed
pixel-year is compared with an exact Poisson test. The share of exposed pixel-years is the share of burned-and-lost expected by chance.

## 12. Stress index and recommendation zones [6.8]

StressIndex = 0.30 (1 - RPI_now) + 0.20 is_degrading + 0.20 Es_frac + 0.15 bare + 0.15 (1 - TWI*), with weights in `Cfg.Stress.Weights`
(version recorded in the file name, e.g. `stress_DRY_v1_030_020_020_015_015.tif`), for ANNUAL, DRY and LR. Es_frac and the bare term (Sentinel-2 BSI)
are season-specific. Each term is min-max scaled between the 2nd and 98th percentile of the AOI (`Cfg.Stress.Clip`), clamped to [0, 1]
(`scale_ref()`); a flat layer is set to 0.5 with a warning. Stress amplification = DRY - LR.

Zones (all thresholds relative to the AOI): **A** top-quartile TWI and top-quartile dry BSI and not degrading; **B** Strahler order 1-2 reaches
with top-quartile mean upslope stress (D8 mass flux of the DRY stress divided by flow accumulation; checked: loading 1 reproduces accumulation
exactly) within `Cfg.Zone.B.Pan.km` of a pan or inside a pan catchment; **C** pan catchments whose mean DRY stress is in the top quartile of
catchment-sized windows of the ranch and whose RPI Sen slope is negative; **D-wet** shrub + tree cover > `Cfg.Zone.Cover.Min` and TWI above the
median; **D-dry** D-wet within `Cfg.Zone.D.Dist.km` of a functional pan; **E** farther than `Cfg.Zone.E.Dist.km` from every functional pan, with
that cover and RPI_now at or above the `Cfg.Zone.E.Rpi.Q` quantile. Output: a categorical raster (overlaps resolved B > A > C > E > D-dry > D-wet) and a
bitmask raster (A=1 B=2 C=4 D-wet=8 D-dry=16 E=32). Each pan carries a woody-sacrifice counter-metric (mature-tree area and Hansen loss in its
`Cfg.Woody.Shed.km` shed) and the McCleery-2025 trade-off note.

## 13. Pan siting [6.9]

Suitability = 0.25 water_gap + 0.20 terrain + 0.15 zone_d_dry + 0.10 refugia + 0.10 cover - 0.20 zone_E - 0.10 woody_risk - 0.10 settled
(`Cfg.Suit.Weights`, an argument of `Run.Pan.Siting()`), each term on [0, 1] on the 30 m grid: **water_gap** = mean over `Cfg.Siting.Gap.Radii.km` of
"farther than r from every functional pan"; **terrain** = geometric mean of the slope trapezoid (5 x 5 smoothed), the flow-accumulation trapezoid
(`Cfg.Terrain.Acc.Cells`) and scaled TWI; **zone_d_dry** = share of the 5 km neighbourhood that is D-wet habitat; **refugia** as above; **cover** = shrub + tree fraction;
**woody_risk** = scaled mature-tree share within the shed; **settled** = 1 - min(d / 15 km, 1). Hard rules: inside the AOI and >= `Cfg.Siting.Edge.m` from its
boundary, outside Zone E, terrain >= `Cfg.Siting.Min.Terrain`. Candidates are picked by **non-maximum suppression** (spacing `Cfg.Siting.MinSpace.m`) on the
suitability smoothed over `Cfg.Siting.Smooth.m`. If fewer than N sites fit, the run is repeated at half the spacing as a tagged alternative. For each candidate:
coordinates, snapped catchment and yield, water-gap area newly covered (5 and 10 km), comfort habitat newly within 5 km of water, mature trees and Hansen loss in the
5 km shed, distances, existing pans upslope, and a stability score (share of 16 runs with each weight halved or raised by half that pick a site within half the spacing).
Existing non-functional pans get the same gap-gain numbers (`rehabilitation_gap_gain.csv`). Candidate rows are appended to `water_pans.csv` only when `Cfg.Siting.Append` is TRUE.

## 14. Figures and report [6.10]

`12_maps_figures.R` builds the report figures with tmap (hillshade underlay) and ggplot2, with stable file names and scico palettes (roma RPI, davos ET, lajolla stress,
batlow TWI, vik diverging). `analysis/Amaka_Analysis.Rmd` reads every number from `outputs/` at knit time and renders `outputs/reports/Amaka_Hydro_Veg_Report_v1.docx`.

## 15. Interpretation limits to carry into the write-up [9.8-9.12]

The stress index is a vegetation-structure index and daytime-biased; it cannot see non-structural or nocturnal pathways. No on-ranch exclosure data exist to calibrate it.
Distance to settlements is not monotonic for wildlife. Woody dynamics should not be attributed to elephants alone; the unburned remainder of woody loss is "not fire".
Hypothesis scores are screening evidence, not diagnoses. Yields rest on an assumed runoff coefficient. Details and the numbers behind each statement are in `known_issues.md` section 7.
