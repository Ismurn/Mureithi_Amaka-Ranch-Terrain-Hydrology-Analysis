# 00_setup.R  -  Amaka Ranch Terrain, Hydrology & Vegetation Analysis
# PURPOSE
#   Single entry point every other script sources. Loads packages, resolves the
#   project root, holds ALL local file paths (so no other script hard-codes a
#   Windows path), holds the run configuration (DESCRIPTION.md section 10 defaults),
#   defines the season helpers (section 2) and initialises WhiteboxTools.
#
# INPUTS
#   None read here. Paths to the OneDrive datasets are only *declared* (Phase 3);
#   existence is reported, never enforced, so the setup never blocks a script
#   that does not need a given dataset.
#
# OUTPUTS
#   Objects in the global environment (dot notation):
#     Dir.*      project folders          Path.*   local input files
#     Cfg.*      run configuration        Seasons.Def + season helper functions
#   Creates the outputs/ sub-folders if they are missing.
#
# RUNTIME NOTES
#   * Run from the RStudio project (or any folder below it): the project root is
#     found by walking up until R/00_setup.R is seen.
#   * To point at the data on another machine without editing this file, set the
#     environment variable AMAKA_DATASETS_DIR (e.g. in ~/.Renviron).
#   * Other scripts start with:  if (!exists("Setup.Loaded")) source("R/00_setup.R")
#   * First run installs missing packages (pacman) and the WhiteboxTools binary.


# PHASE 1: Packages 
# House convention: pacman::p_load only, no library() chains.
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  sf, terra, tidyverse, tmap, exactextractr, whitebox, trend,
  scico, viridis, patchwork, rmarkdown, janitor, lubridate
)

# terra::extract() is masked by tidyr::extract() once tidyverse is attached.
# Scripts must call terra::extract() explicitly (DESCRIPTION.md section 9.2 workflow).
tmap::tmap_mode("plot")


# PHASE 2: Project root and folders
find_project_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(d, "R", "00_setup.R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) {
      stop("Project root not found: no R/00_setup.R at or above '", start, "'.",
           " Open the .Rproj or setwd() into the repository first.", call. = FALSE)
    }
    d <- parent
  }
}

Dir.Root          <- find_project_root()
Dir.R             <- file.path(Dir.Root, "R")
Dir.GEE           <- file.path(Dir.Root, "gee")
Dir.Docs          <- file.path(Dir.Root, "docs")
Dir.Analysis      <- file.path(Dir.Root, "analysis")

Dir.Data          <- file.path(Dir.Root, "data")
Dir.Boundary      <- file.path(Dir.Data, "boundary")
Dir.Boundary.Ctx  <- file.path(Dir.Boundary, "context")
Dir.RPI           <- file.path(Dir.Data, "rpi")
Dir.PML           <- file.path(Dir.Data, "pml_v2")
Dir.StudyPixels   <- file.path(Dir.Data, "study_pixels")
Dir.DEM           <- file.path(Dir.Data, "dem")
Dir.GEE.Exports   <- file.path(Dir.Data, "gee_exports")
Dir.Pans          <- file.path(Dir.Data, "water_pans")

Dir.Out           <- file.path(Dir.Root, "outputs")
Dir.Out.Rasters   <- file.path(Dir.Out, "rasters")
Dir.Out.Hydro     <- file.path(Dir.Out.Rasters, "hydrology")   # WBT working files (03)
Dir.Out.Hydro.Buf <- file.path(Dir.Out.Hydro, "_buffered")     # full-envelope hydrology rasters (03)
Dir.Out.Figures   <- file.path(Dir.Out, "figures")
Dir.Out.Tables    <- file.path(Dir.Out, "tables")
Dir.Out.Reports   <- file.path(Dir.Out, "reports")

for (d in c(Dir.Boundary, Dir.Boundary.Ctx, Dir.DEM, Dir.GEE.Exports, Dir.Pans,
            Dir.Out.Rasters, Dir.Out.Hydro, Dir.Out.Hydro.Buf, Dir.Out.Figures, Dir.Out.Tables,
            Dir.Out.Reports)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# PHASE 3: Local input paths  (EDIT HERE ONLY) 
# Everything below sits on OneDrive; nothing is copied into the repository.
Dir.Datasets <- Sys.getenv(
  "AMAKA_DATASETS_DIR",
  unset = "C:/Users/IsaacM/OneDrive - African Wildlife Foundation/R WORK/Datasets_and_Shapefiles"
)

# -- Vectors --
Path.Ranches      <- file.path(Dir.Datasets, "Tsavo_Conservation_Areas", "RanchesTTWCA.shp")
Path.SettledAreas <- file.path(Dir.Datasets, "Tsavo_Conservation_Areas", "SettledAreas.shp")
Path.TsavoEast    <- file.path(Dir.Datasets, "WDPA_Amboseli_Tsavo", "National Parks&Reserves",
                               "TsavoEast-WDPA_WDOECM_Apr2025_Public_752_shp",
                               "WDPA_WDOECM_Apr2025_Public_752_shp_0",
                               "WDPA_WDOECM_Apr2025_Public_752_shp-polygons.shp")
Path.TsavoWest    <- file.path(Dir.Datasets, "WDPA_Amboseli_Tsavo", "National Parks&Reserves",
                               "TsavoWest-WDPA_WDOECM_Apr2025_Public_19564_shp",
                               "WDPA_WDOECM_Apr2025_Public_19564_shp_0",
                               "WDPA_WDOECM_Apr2025_Public_19564_shp-polygons.shp")
Path.Kenya.Admin1 <- file.path(Dir.Datasets, "Country_Borders", "kenya_regional_admin1_boundary_wgs84.shp")

# RPI v2 (Guy Lomax; East Africa, 0.004167 deg, 23 annual layers 2000-2022)
Dir.RPI.Src       <- file.path(Dir.Datasets, "RPI Datasets (GuyLomax)")
Path.RPI          <- file.path(Dir.RPI.Src, "rpi_ea_v2.tif")                    # rpi.2000 ... rpi.2022
Path.RPI.GPP.Act  <- file.path(Dir.RPI.Src, "actual_gpp_ea_v2.tif")             # GPP.2000 ... GPP.2022
Path.RPI.GPP.Pot  <- file.path(Dir.RPI.Src, "potential_gpp_ea_v2.tif")          # potential_gpp_predicted.2000 ...
Path.RPI.TempPerf <- file.path(Dir.RPI.Src, "rpi_ea_temporal_performance.tif")  # mae, rmse, rsq

# PML_V2 and study pixels (Kenya, 0.004492 deg ~ 500 m, WGS84)
Dir.PML.Src       <- file.path(Dir.Datasets, "PMLv2_and_IGBPpixels")
# 80 layers <i>_<GPP|Ec|Es|Ei>_<2001..2020>. Each is a SUM of ~45 8-day composites of daily RATES
# (CAS/IGSNRR/PML/V2_v017) over a Sep-Aug hydrological year, so x 8 gives mm or gC m-2 per year.
# The 2020 layer looks like a partial year (known_issues D3).
Path.PML.Annual   <- file.path(Dir.PML.Src, "PML_V2_Kenya_Annual_Stack.tif")
Path.PML.Mean     <- file.path(Dir.PML.Src, "PML_V2_Kenya_Drive_Export.tif")    # 4 layers: GPP, Ec, Es, Ei
Path.StudyPixels  <- file.path(Dir.PML.Src, "studyPixels_Kenya_Export.tif")     # igbp, trees, shrubland, grassland, bare, allRangelands
                                                                                # Kenya-wide BASELINE only: ~97% shrubland over Amaka (known_issues D5)

# 10 m land cover for within-ranch structure: exported by gee/08, downloaded to data/gee_exports/
Path.WorldCover     <- file.path(Dir.GEE.Exports, "worldcover_2021_envelope.tif")      # ESA WorldCover v200 (2021), band Map
Path.WorldCover.V100 <- file.path(Dir.GEE.Exports, "worldcover_2020_v100_envelope.tif") # optional: the version behind the study pixels

# Project-side files
Path.Pans         <- file.path(Dir.Pans, "water_pans.csv")
Path.FieldPoints  <- file.path(Dir.Pans, "amaka_field_points.csv")              # iPhone-GPS visit, 7 located points
Path.Amaka.UTM    <- file.path(Dir.Boundary, "Amaka_Ranch.shp")                 # canonical AOI (EPSG:32737), written by 01
Path.Amaka.WGS84  <- file.path(Dir.Boundary, "Amaka_Ranch_wgs84.shp")           # mapping / GEE upload (EPSG:4326)


# PHASE 4: Run configuration  (DESCRIPTION.md section 10 defaults)
# Coordinate reference systems 
Cfg.CRS.Geo     <- 4326L    # WGS84, mapping
Cfg.CRS.Metric  <- 32737L   # WGS84 / UTM 37S, all metric operations
Cfg.CRS.Arc1960 <- 21037L   # Arc 1960 / UTM 37S, the historic mis-tag (section 9.1)

# What the numbers inside RanchesTTWCA.shp / SettledAreas.shp actually ARE.
# The files carry an EPSG:32737 tag, but section 9.1 warns Tsavo shapefiles have
# arrived as Arc 1960 / UTM 37S before (Kitobo: ~250 m offset). Keep 32737 until
# the Google Earth overlay from 01_boundary_prep.R shows otherwise; if it shows
# an offset set this to 21037L, re-run 01, and flip Cfg.CRS.Confirmed to TRUE
# either way once you have looked at the overlay.
Cfg.Ranches.CRS.Assumed <- 32737L
Cfg.CRS.Confirmed       <- TRUE

# Area of interest
# Amaka Ranch and Conservancy is catalogued under its old name in the TTWCA file.
Cfg.Amaka.Field <- "R_NAME"
Cfg.Amaka.Names <- c("MUSAMITI", "AMAKA RANCH")   # first = what the file holds today
Cfg.Buffer.m    <- 7000                           # catchment envelope; 5000-10000 allowed

# Which boundary counts as "Amaka" (decision 2026-10-27, Isaac; known_issues F8, F9) -----
#   "musamiti" (current): the TTWCA polygon exactly as delivered. It is coarse (4 corners, 4,782 ha
#                against the 14,815 acres TTWCA lists for Amaka) and field points PT05-PT07 fall in
#                the neighbouring Kambanga polygon. Accepted for now. (A hull of the field points was
#                tried and rejected.)
#   "file"     : a boundary supplied by Amaka management / TTWCA; set Cfg.AOI.File (shp, gpkg, kml).
#                The intended end state: change the mode, re-run 01.
Cfg.AOI.Mode          <- "musamiti"
Cfg.AOI.File          <- NA_character_
Cfg.AOI.Listed.acres  <- 14815                    # TTWCA member list (ttwcakenya.com/our-members), for sanity checks

# DEM (02) 
Cfg.DEM.Source  <- "glo30"   # file = data/gee_exports/amaka_dem_<source>_envelope.tif ("glo30" | "fabdem" | "srtm")
Cfg.DEM.Smooth  <- FALSE     # optional wbt_gaussian_filter(); record the choice, it changes TWI
Cfg.DEM.Sigma   <- 1.0

# Land cover (06, 10) 
# Decision 2026-10-26: within-ranch structure (bare fraction, shrub:tree, Zone A) comes from 10 m
# ESA WorldCover plus Sentinel-2 (BSI, NDVI); the 500 m study pixels are the Kenya-wide baseline only.
Cfg.LC.Source <- "worldcover10m"
WorldCover.Classes <- c(trees = 10L, shrubland = 20L, grassland = 30L, cropland = 40L, builtup = 50L,
                        bare = 60L, water = 80L, herbaceous_wetland = 90L, mangroves = 95L)

# Time windows 
Cfg.Years.RPI       <- 2000:2022   # rpi_ea_v2.tif has 23 layers (brief said 2000-2023)
Cfg.Years.PML.Local <- 2001:2019   # local annual stack, full Sep-Aug years only (2020 looks partial)
Cfg.PML.Composite.Days <- 8        # local annual layers x this = mm / gC m-2 per year
Cfg.Years.GEE       <- 2000:2024   # default climate / vegetation exports

# Hydrology (03, 07, 08) 
Cfg.Stream.Threshold <- 2000       # D8 cells (1.8 km2). Calibrated in 03 (known_issues H1): the sparsest network that still
                                   # puts all four dams on a stream and keeps most of the riparian NDVI contrast (brief start value: 1000)
Cfg.Riparian.m       <- 125        # 100-150 m
Cfg.Snap.m           <- 100        # pan pour-point snap distance
Cfg.TPI.Window.m     <- 500
Cfg.Shed.Radii.km    <- c(3, 5, 10, 15)

# Composite stress index (10): order = RPI deficit, degrading, Es/ET DRY, bare, TWI deficit
Cfg.Stress.Weights <- c(rpi = 0.30, degrading = 0.20, es_dry = 0.20, bare = 0.15, twi = 0.15)

# Pan siting (11)
Cfg.Suit.Weights <- c(water_gap = 0.25, terrain = 0.20, zone_d_dry = 0.15, refugia = 0.10,
                      cover = 0.10, zone_e = -0.20, woody_risk = -0.10, settled = -0.10)
Cfg.Siting.N          <- 8
Cfg.Siting.MinSpace.m <- 3000

# Flags 
Cfg.Sterile.Ratio <- 5             # shrub:tree ratio threshold for "sterile landscape"
Cfg.Degrade.P     <- 0.1           # Sen-slope significance for "degrading" pixels

# Analysis parameters added when 02-12 were implemented (2026-10-27) 
# Years actually present in the GEE exports (calendar-year seasons, see docs/data_sources.md)
Cfg.Years.PML.Seasonal <- 2001:2020   # gee/09 (the asset runs to Dec 2020)
Cfg.Years.S2           <- 2019:2024   # gee/06
Cfg.Years.SMAP         <- 2016:2024   # gee/07
Cfg.Years.MODIS        <- 2001:2024   # gee/04, gee/05
Cfg.Years.CHIRPS       <- 2000:2024   # gee/03
Cfg.Years.Fire         <- 2001:2024   # gee/10
Cfg.RPI.Now.Years      <- 2020:2022   # "RPI_now" in the stress index = mean of the last three RPI years
Seasons.All            <- c("ANNUAL", "LR", "SR", "DRY")

# "High" is defined RELATIVE to the ranch (top quartile of cells inside the AOI): inside Amaka nearly
# every cell is shrubland with almost no bare ground, so absolute cut-offs select nothing (known_issues D5).
Cfg.Relative.Q         <- 0.75

# Hydrology / pans
Cfg.Stream.Calib.Thresholds <- c(250, 500, 1000, 2000, 5000)   # D8 cells tried in the 9.7 calibration (03)
Cfg.Snap.Img.m         <- 60          # snap distance when a pan starts from an image-read centre (+-5 m imagery, +-15 m reading)
Cfg.WaterGap.Radius.km <- 5           # radius of the primary water-gap layer (all four radii are written)
Cfg.Runoff.Coef        <- 0.10        # annual runoff coefficient for water-yield estimates (an ASSUMPTION; no gauge data)
Cfg.Tree.Mature.Pct    <- 10          # Hansen treecover2000 at or above this = "mature tree" cell (canopy > 5 m)
Cfg.Woody.Shed.km      <- 5           # shed radius for the woody-sacrifice metrics
Cfg.Refugia.Ref.km     <- c(boundary = 3, settled = 15, roads = 10)   # distance at which each refugia term saturates
Cfg.Pan.StreamFed.Types <- c("earth_dam", "natural_pan")   # pan types that get a snapped pour point and a delineated catchment (08).
                                      # Other types (rock-catchment reservoir, borehole, ...) are NOT stream-fed: snapping them to the
                                      # nearest channel would hand them a river basin (the HQ reservoir would get 13,000 ha), so they
                                      # get a context disc instead (known_issues F12)
Cfg.Pan.Context.m      <- 500         # radius of that context disc
Cfg.Fire.Loss.Window.Years <- 1       # a Hansen loss pixel counts as "burned-and-lost" when its 500 m cell burned in the loss year or
                                      # up to this many years before it (post-fire mortality is delayed; Hansen dates loss by growing season)

# Recommendation zones (10) and pan siting (11). All "high / low" cut-offs are relative to the ranch (Cfg.Relative.Q).
Cfg.Stress.Clip        <- c(0.02, 0.98)   # stress terms are min-max scaled between these AOI percentiles (a few extreme pixels would otherwise
                                          # squash everything else); values outside are clamped to 0 / 1
Cfg.Zone.Cover.Min     <- 0.25        # Zones D and E: shrub + tree cover fraction above this (DESCRIPTION.md 6.8)
Cfg.Zone.D.Dist.km     <- 5           # Zone D-dry: within this distance of a functional pan
Cfg.Zone.E.Dist.km     <- 10          # Zone E: farther than this from every functional pan (dry-season elephant reference is 10-15 km)
Cfg.Zone.E.Rpi.Q       <- 0.25        # Zone E "acceptable RPI": RPI_now at or above this quantile of the ranch
Cfg.Zone.B.Orders      <- 1:2         # Zone B: low-order streams
Cfg.Zone.B.Pan.km      <- 3           # Zone B: within this distance of a pan, or inside a pan catchment
Cfg.Siting.Gap.Radii.km <- c(3, 5, 10, 15)   # radii averaged in the siting water-gap term (11); a single value gives the brief's literal binary
Cfg.Siting.Min.Terrain <- 0.5         # a candidate needs at least this terrain feasibility (hard eligibility rule, on top of the weighted sum)
Cfg.Siting.Edge.m      <- 300         # ... and must lie at least this far inside the AOI boundary (the polygon is coarse: known_issues F8)
Cfg.Siting.Smooth.m    <- 300         # suitability is smoothed over this radius before candidates are picked (single-pixel peaks are noise)

# Pan siting terrain feasibility (11): a trapezoid per variable = 0 outside the outer pair,
# 1 between the inner pair. Farm dams in dry savanna typically sit on 0.3-20 km2 of upslope area on gentle ground.
Cfg.Terrain.Slope.Deg  <- c(0.2, 0.6, 3.0, 6.0)               # degrees
Cfg.Terrain.Acc.Cells  <- c(300, 2000, 20000, 150000)         # D8 contributing area, 30 m cells (0.27, 1.8, 18, 135 km2)
Cfg.Siting.Append      <- FALSE       # append candidate rows to water_pans.csv only after Isaac confirms


# PHASE 5: Season helpers (section 2) 
# LR  = Mar-May,  SR = Oct-Dec,  DRY = Jun-Sep + Jan-Feb of the SAME calendar year.
# "ANNUAL" is the calendar year. The three seasons cover 12 months with no gap.
Seasons.Def <- tibble::tibble(
  season = c("LR", "SR", "DRY"),
  label  = c("Long rains", "Short rains", "Dry season"),
  months = list(3:5, 10:12, c(1:2, 6:9))
)

# Accepts Date / POSIXct / ISO character, or an integer month 1-12.
.month_of <- function(x) if (is.numeric(x)) as.integer(x) else lubridate::month(x)

season_of_date <- function(x) {
  m <- .month_of(x)
  dplyr::case_when(
    m %in% 3:5          ~ "LR",
    m %in% 10:12        ~ "SR",
    m %in% c(1:2, 6:9)  ~ "DRY",
    TRUE                ~ NA_character_
  )
}
is_lr  <- function(x) season_of_date(x) == "LR"
is_sr  <- function(x) season_of_date(x) == "SR"
is_dry <- function(x) season_of_date(x) == "DRY"

# Months belonging to a season code; "ANNUAL" returns 1:12.
season_months <- function(season) {
  if (identical(season, "ANNUAL")) return(1:12)
  Seasons.Def$months[[match(season, Seasons.Def$season)]]
}


# PHASE 6: Small shared helpers 
# Announce a phase in the console so long scripts are easy to follow.
log_phase <- function(txt) message("\n== ", txt, " ", strrep("=", max(3, 70 - nchar(txt))))

# Read one of the Tsavo-Conservation-Area shapefiles (ranches, settled areas) under
# the CRS assumption in Cfg.Ranches.CRS.Assumed. This is a deliberate RELABEL of
# the coordinate numbers, announced in the console - never a silent st_set_crs().
read_tsavo_layer <- function(path, crs.assumed = Cfg.Ranches.CRS.Assumed) {
  x   <- sf::st_read(path, quiet = TRUE)
  tag <- sf::st_crs(x)$epsg
  if (is.na(tag) || tag != crs.assumed) {
    message("  [CRS] ", basename(path), ": file tag ", ifelse(is.na(tag), "none", paste0("EPSG:", tag)),
            " relabelled as EPSG:", crs.assumed, " (Cfg.Ranches.CRS.Assumed)")
    x <- suppressWarnings(sf::st_set_crs(x, crs.assumed))
  }
  sf::st_make_valid(x)
}

# Write a shapefile, replacing any previous copy (all sidecar files).
write_shp <- function(x, path) {
  sf::st_write(x, path, delete_layer = TRUE, quiet = TRUE)
  invisible(path)
}

# Build the Amaka AOI (EPSG:32737) from the dissolved TTWCA Musamiti polygon, per Cfg.AOI.Mode
# (see Phase 4). Returns list(geom, mode, src, label). A supplied file with no CRS is refused:
# a CRS is never guessed (DESCRIPTION.md 9.1).
build_aoi <- function(ttwca.geom, mode = Cfg.AOI.Mode, file = Cfg.AOI.File) {
  clean <- function(g) {
    g <- sf::st_make_valid(sf::st_union(sf::st_zm(g, drop = TRUE, what = "ZM")))
    if (any(sf::st_geometry_type(g) == "GEOMETRYCOLLECTION")) g <- sf::st_collection_extract(g, "POLYGON")
    sf::st_cast(g, "MULTIPOLYGON")
  }
  geom <- clean(ttwca.geom)

  if (identical(mode, "musamiti")) {
    src <- "TTWCA Musamiti polygon as delivered"

  } else if (identical(mode, "file")) {
    if (is.na(file) || !path_exists(file)) stop("Cfg.AOI.File not found: ", file, call. = FALSE)
    x <- sf::st_read(file, quiet = TRUE)
    if (is.na(sf::st_crs(x))) stop("The supplied boundary has no CRS; set it deliberately, never guess: ", file, call. = FALSE)
    if (!all(sf::st_geometry_type(x) %in% c("POLYGON", "MULTIPOLYGON"))) stop("The supplied boundary must be polygons: ", file, call. = FALSE)
    geom <- clean(sf::st_transform(sf::st_geometry(x), Cfg.CRS.Metric))
    src  <- paste0("supplied boundary: ", basename(file))

  } else {
    stop("Cfg.AOI.Mode must be 'musamiti' or 'file', not '", mode, "'", call. = FALSE)
  }

  list(geom = geom, mode = mode, src = src, label = "Amaka Ranch and Conservancy")
}

# file.exists() reports FALSE for paths over 260 characters on Windows even when the
# file is fine and GDAL can read it (the Tsavo West WDPA path is 265). Fall back to
# the \\?\ long-path form before declaring a file missing.
path_exists <- function(p) {
  ok <- file.exists(p)
  long <- !ok & nchar(p) >= 250 & .Platform$OS.type == "windows"
  if (any(long)) {
    ok[long] <- file.exists(paste0("\\\\?\\", gsub("/", "\\\\", normalizePath(p[long], winslash = "/", mustWork = FALSE))))
  }
  ok
}

# Soft check of declared inputs: lists what is missing, never stops the run.
check_inputs <- function() {
  p <- c(Ranches = Path.Ranches, SettledAreas = Path.SettledAreas, TsavoEast = Path.TsavoEast,
         TsavoWest = Path.TsavoWest, KenyaAdmin1 = Path.Kenya.Admin1, RPI = Path.RPI,
         RPI.GPP.Act = Path.RPI.GPP.Act, RPI.GPP.Pot = Path.RPI.GPP.Pot,
         RPI.TempPerf = Path.RPI.TempPerf, PML.Annual = Path.PML.Annual,
         PML.Mean = Path.PML.Mean, StudyPixels = Path.StudyPixels)
  tibble::tibble(input = names(p), exists = path_exists(p), path = unname(p))
}


# PHASE 6b: Shared helpers for the analysis scripts (02-12) 
# -- Boundary / envelope 
read_aoi      <- function() sf::st_read(Path.Amaka.UTM, quiet = TRUE)
read_envelope <- function() sf::st_read(file.path(Dir.Boundary, sprintf("Amaka_Envelope_%gkm.shp", Cfg.Buffer.m / 1000)), quiet = TRUE)

# GEE exports (data/gee_exports/) 
# One seasonal stack: <prefix>_<SEASON>_<y0>_<y1>.tif, band names as exported (e.g. es_DRY_2001).
gee_stack <- function(prefix, season = NULL) {
  pat <- if (is.null(season)) sprintf("^%s_[0-9]{4}_[0-9]{4}\\.tif$", prefix)
         else sprintf("^%s_%s_[0-9]{4}_[0-9]{4}\\.tif$", prefix, season)
  f <- list.files(Dir.GEE.Exports, pattern = pat, full.names = TRUE)
  if (length(f) != 1) stop("Expected one GEE export matching ", pat, " in data/gee_exports, found ", length(f), call. = FALSE)
  terra::rast(f)
}
gee_file <- function(name) {
  f <- file.path(Dir.GEE.Exports, paste0(name, ".tif"))
  if (!file.exists(f)) stop("GEE export not found: ", f, call. = FALSE)
  f
}
# Year of each layer from its name: es_DRY_2001 -> 2001, rpi.2000 -> 2000, burned_2001 -> 2001
layer_year <- function(r) as.integer(sub(".*[._]([0-9]{4})$", "\\1", if (is.character(r)) r else names(r)))

# Scaling 
# Min-max scale to [0, 1] (DESCRIPTION.md 6.8 / 6.9). A (nearly) constant layer is NOT stretched, because
# the stretch would only amplify noise: fractions already on [0, 1] are kept as they are
# (already01 = TRUE); anything else becomes 0.5 with a warning. min.range = smallest span worth stretching.
minmax01 <- function(x, already01 = FALSE, min.range = 0, label = "") {
  rng  <- if (inherits(x, "SpatRaster")) unlist(terra::global(x, c("min", "max"), na.rm = TRUE)) else range(x, na.rm = TRUE)
  span <- unname(rng[2] - rng[1])
  if (!is.finite(span) || span <= 0 || span < min.range) {
    if (already01) return(if (inherits(x, "SpatRaster")) terra::clamp(x, 0, 1) else pmin(pmax(x, 0), 1))
    warning("minmax01(", label, "): range ", signif(span, 3), " is too small to stretch; set to 0.5", call. = FALSE)
    return(if (inherits(x, "SpatRaster")) x * 0 + 0.5 else ifelse(is.na(x), NA_real_, 0.5))
  }
  (x - unname(rng[1])) / span
}

# Trapezoid membership: 0 outside p[1]..p[4], 1 between p[2]..p[3], linear ramps between.
trapezoid <- function(x, p) {
  stopifnot(length(p) == 4, all(diff(p) > 0))
  if (inherits(x, "SpatRaster")) {
    up <- terra::clamp((x - p[1]) / (p[2] - p[1]), 0, 1); dn <- terra::clamp((p[4] - x) / (p[4] - p[3]), 0, 1)
    return(min(up, dn))                                   # terra's Summary method: cell-wise minimum of the two rasters
  }
  pmin(pmax((x - p[1]) / (p[2] - p[1]), 0), 1, pmax((p[4] - x) / (p[4] - p[3]), 0))
}

# Quantile of a raster layer (terra::global() needs a NAMED function, not an anonymous one).
global_q <- function(r, p) unname(terra::global(r, quantile, probs = p, na.rm = TRUE)[[1]])

# Min-max scaling RELATIVE to the ranch (10, 11): between two quantiles of the values inside `ref` (the AOI mask raster), clamped to
# [0, 1], so that a few extreme pixels do not squash the rest. A (nearly) flat layer is not stretched: it becomes 0.5 with a warning.
scale_limits <- function(x, ref, clip = Cfg.Stress.Clip) c(global_q(terra::mask(x, ref), clip[1]), global_q(terra::mask(x, ref), clip[2]))
scale_ref <- function(x, ref, clip = Cfg.Stress.Clip, min.range = 0.01, label = "") {
  lim <- scale_limits(x, ref, clip)
  if (!is.finite(diff(lim)) || diff(lim) < min.range) {
    warning("scale_ref(", label, "): range ", signif(diff(lim), 3), " is too small to stretch; set to 0.5", call. = FALSE)
    return(x * 0 + 0.5)
  }
  terra::clamp((x - lim[1]) / diff(lim), 0, 1)
}

# Zonal statistics 
# Area-weighted statistic of every layer of `r` inside each polygon (polygons are reprojected to the
# raster CRS first). Returns one row per polygon, one column per layer.
zonal_stat <- function(r, polys, fun = "mean") {
  polys <- sf::st_transform(polys, terra::crs(r))
  out <- exactextractr::exact_extract(r, polys, fun = fun, progress = FALSE)
  if (is.null(dim(out))) out <- data.frame(value = out)
  if (terra::nlyr(r) == 1) names(out) <- names(r) else names(out) <- sub(paste0("^", fun, "\\."), "", names(out))
  tibble::as_tibble(out)
}
# Point-sample every layer of `r` at the analysis-cell centroids (DESCRIPTION.md 9.2: native grids, no resampling).
sample_cells <- function(r, cells, geo = FALSE) {
  xy <- if (geo) cbind(cells$lon, cells$lat) else cbind(cells$x, cells$y)
  v  <- terra::extract(r, xy)
  v  <- v[, setdiff(names(v), "ID"), drop = FALSE]
  tibble::as_tibble(v)
}

# Analysis grid: the RPI native grid over the catchment envelope (DESCRIPTION.md 9.2)
# Every cell-level table in the pipeline is keyed on cell_id from this grid. Cached on first use.
cell_polygons <- function(g) {
  p <- sf::st_as_sf(terra::as.polygons(g, dissolve = FALSE))
  sf::st_transform(p, Cfg.CRS.Metric)
}
get_analysis_grid <- function(refresh = FALSE) {
  f.r <- file.path(Dir.Out.Rasters, "analysis_grid_rpi.tif")
  f.t <- file.path(Dir.Out.Tables, "analysis_cells.csv")
  if (!refresh && file.exists(f.r) && file.exists(f.t)) {
    return(list(grid = terra::rast(f.r), cells = readr::read_csv(f.t, show_col_types = FALSE)))
  }
  env <- read_envelope(); aoi <- sf::st_union(sf::st_geometry(read_aoi()))
  env.ll <- sf::st_transform(env, Cfg.CRS.Geo)
  g <- terra::crop(terra::rast(Path.RPI)[[1]], terra::ext(terra::vect(env.ll)), snap = "out")
  terra::values(g) <- seq_len(terra::ncell(g)); names(g) <- "cell_id"
  poly <- cell_polygons(g)
  ctr  <- sf::st_centroid(sf::st_geometry(poly))
  ctr.ll <- sf::st_coordinates(sf::st_transform(ctr, Cfg.CRS.Geo))
  ctr.xy <- sf::st_coordinates(ctr)
  frac <- function(zone) {
    hit <- lengths(sf::st_intersects(poly, zone)) > 0
    out <- numeric(nrow(poly))
    if (any(hit)) out[hit] <- as.numeric(sf::st_area(sf::st_intersection(sf::st_geometry(poly)[hit], zone))) / as.numeric(sf::st_area(poly))[hit]
    out
  }
  cells <- tibble::tibble(cell_id = poly$cell_id, lon = ctr.ll[, 1], lat = ctr.ll[, 2], x = ctr.xy[, 1], y = ctr.xy[, 2],
                          area_ha = round(as.numeric(sf::st_area(poly)) / 1e4, 3),
                          aoi_frac = round(frac(aoi), 4), env_frac = round(frac(sf::st_union(sf::st_geometry(env))), 4))
  cells$in_aoi <- cells$aoi_frac >= 0.5      # a cell "belongs to Amaka" when at least half of it lies inside the AOI
  cells$in_env <- cells$env_frac >= 0.5
  terra::writeRaster(g, f.r, overwrite = TRUE, datatype = "INT4S")
  readr::write_csv(cells, f.t)
  list(grid = g, cells = cells)
}

# Pans 
# water_pans.csv joined with the image-read centres and stand-off uncertainty of the field file.
# Start position = image-read centre where there is one, else the recorded coordinates.
read_pans <- function() {
  p <- readr::read_csv(Path.Pans, show_col_types = FALSE)
  fp <- if (file.exists(Path.FieldPoints)) readr::read_csv(Path.FieldPoints, show_col_types = FALSE) else NULL
  keep <- c("lat_img", "lon_img", "img_offset_m", "img_confidence", "uncertainty_m_est", "position_type")
  if (!is.null(fp) && "field_id" %in% names(p)) {
    j <- fp[, intersect(c("id", keep), names(fp))]
    p <- dplyr::left_join(p, j, by = c("field_id" = "id"))
  }
  for (k in keep) if (!k %in% names(p)) p[[k]] <- if (k %in% c("img_confidence", "position_type")) NA_character_ else NA_real_
  p$start_lon <- dplyr::coalesce(p$lon_img, p$longitude)
  p$start_lat <- dplyr::coalesce(p$lat_img, p$latitude)
  p$start_source <- dplyr::case_when(!is.na(p$lon_img) ~ "image-read centre", !is.na(p$longitude) ~ "recorded coordinates", TRUE ~ "none")
  p$snap_m   <- ifelse(!is.na(p$lon_img), Cfg.Snap.Img.m, pmax(Cfg.Snap.m, dplyr::coalesce(p$uncertainty_m_est, Cfg.Snap.m)))
  p$located  <- !is.na(p$start_lon)
  p$functional <- p$status %in% "existing_functional_under_pressure"
  p
}

# Tables
# Minimal .xlsx writer (no extra package: only 'zip', which officer already brings in).
write_xlsx_min <- function(df, path, sheet = "data") {
  if (!requireNamespace("zip", quietly = TRUE)) { warning("package 'zip' not available; not written: ", path, call. = FALSE); return(invisible(FALSE)) }
  df  <- as.data.frame(df, stringsAsFactors = FALSE)
  esc <- function(s) { s <- gsub("&", "&amp;", s, fixed = TRUE); s <- gsub("<", "&lt;", s, fixed = TRUE)
                       s <- gsub(">", "&gt;", s, fixed = TRUE); gsub("[\u0001-\u0008\u000b\u000c\u000e-\u001f]", "", s) }
  colref <- function(i) { s <- ""; while (i > 0) { r <- (i - 1) %% 26; s <- paste0(LETTERS[r + 1], s); i <- (i - 1) %/% 26 }; s }
  cell <- function(v, ref) {
    if (length(v) == 0 || is.na(v)) return("")
    if (is.numeric(v)) return(sprintf('<c r="%s"><v>%s</v></c>', ref, format(v, digits = 15, scientific = FALSE, trim = TRUE)))
    sprintf('<c r="%s" t="inlineStr"><is><t xml:space="preserve">%s</t></is></c>', ref, esc(as.character(v)))
  }
  hdr  <- paste0(vapply(seq_along(df), function(j) cell(names(df)[j], paste0(colref(j), 1)), ""), collapse = "")
  rows <- c(sprintf('<row r="1">%s</row>', hdr),
            vapply(seq_len(nrow(df)), function(i)
              sprintf('<row r="%d">%s</row>', i + 1, paste0(vapply(seq_along(df), function(j) cell(df[[j]][i], paste0(colref(j), i + 1)), ""), collapse = "")), ""))
  tmp <- tempfile("xlsx"); dir.create(tmp); dir.create(file.path(tmp, "_rels")); dir.create(file.path(tmp, "xl", "_rels"), recursive = TRUE); dir.create(file.path(tmp, "xl", "worksheets"))
  w <- function(rel, txt) writeLines(txt, file.path(tmp, rel), useBytes = TRUE)
  xml <- '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
  w("[Content_Types].xml", paste0(xml, '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>'))
  w("_rels/.rels", paste0(xml, '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>'))
  w("xl/workbook.xml", paste0(xml, sprintf('<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="%s" sheetId="1" r:id="rId1"/></sheets></workbook>', esc(substr(sheet, 1, 31)))))
  w("xl/_rels/workbook.xml.rels", paste0(xml, '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>'))
  w("xl/worksheets/sheet1.xml", paste0(xml, '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>', paste0(rows, collapse = ""), '</sheetData></worksheet>'))
  if (file.exists(path)) file.remove(path)
  zip::zipr(normalizePath(path, winslash = "/", mustWork = FALSE), files = c("[Content_Types].xml", "_rels", "xl"), root = tmp)
  unlink(tmp, recursive = TRUE)
  invisible(TRUE)
}


# PHASE 7: WhiteboxTools
# Install the binary once, then point the R wrapper at it. wbt_init() reports
# whether a working binary was found.
.wbt.ok <- isTRUE(tryCatch(whitebox::wbt_init(check_version = FALSE), error = function(e) FALSE))
if (!.wbt.ok) {
  message("WhiteboxTools binary not found - installing (one-off).")
  whitebox::install_whitebox()
  .wbt.ok <- isTRUE(tryCatch(whitebox::wbt_init(check_version = FALSE), error = function(e) FALSE))
}
if (!.wbt.ok) warning("WhiteboxTools could not be initialised; 03_hydrology.R will not run.", call. = FALSE)
rm(.wbt.ok)


# PHASE 8: Report 
Inputs.Status <- check_inputs()
if (!all(Inputs.Status$exists)) {
  message("00_setup.R: missing inputs (edit Phase 3 or set AMAKA_DATASETS_DIR):")
  print(Inputs.Status[!Inputs.Status$exists, c("input", "path")], n = Inf)
}
if (!Cfg.CRS.Confirmed) {
  message("00_setup.R: Tsavo-shapefile CRS is UNCONFIRMED (assumed EPSG:", Cfg.Ranches.CRS.Assumed,
          "). Run 01_boundary_prep.R and check the Google Earth overlay (DESCRIPTION.md 9.1).")
}

if (identical(Cfg.AOI.Mode, "musamiti")) {
  message("00_setup.R: AOI = the TTWCA Musamiti polygon (coarse; ~20% below the size TTWCA lists, and ",
          "3 of the 7 field points fall outside it): docs/known_issues.md F8.")
}

Setup.Loaded <- TRUE
message("00_setup.R loaded. Project root: ", Dir.Root)
