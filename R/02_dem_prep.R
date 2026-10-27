# ==============================================================================
# 02_dem_prep.R  -  Ingest the GEE DEM and clip it to the catchment envelope
# ==============================================================================
# PURPOSE
#   Take the DEM exported by gee/02_dem_amaka_export.js, put it on EPSG:32737,
#   clip it to the 7 km catchment envelope (NOT to Amaka alone, so upslope
#   contributing area is kept), run basic quality checks, optionally smooth it,
#   and write the DEM that 03_hydrology.R consumes.
#
# INPUTS
#   data/gee_exports/amaka_dem_<Cfg.DEM.Source>_envelope.tif   (from GEE)
#   data/boundary/Amaka_Envelope_<N>km.shp, Amaka_Ranch.shp    (from 01)
#
# OUTPUTS
#   data/dem/amaka_dem_utm37s.tif             DEM for hydrology, EPSG:32737
#   outputs/tables/02_dem_prep_log.csv        provenance and QA of this run
#   outputs/figures/02_dem_clip_map.png       DEM + hillshade + AOI + envelope
#
# RUNTIME NOTES
#   * Seconds. Re-run whenever the DEM source, buffer or smoothing option changes.
#   * Smoothing (Cfg.DEM.Smooth) is OFF by default. If you switch it on the choice
#     is written to the log, because it changes TWI and stream extraction.
#   * GLO-30 is a surface model (it includes tree canopy); FABDEM is the forest-
#     and building-corrected derivative. Cfg.DEM.Source switches between them.
# ==============================================================================

if (!exists("Setup.Loaded")) source("R/00_setup.R")


# ---- PHASE 1: Locate the GEE export ----------------------------------------
log_phase("PHASE 1: locate DEM export")

DEM.Src.Path <- file.path(Dir.GEE.Exports, sprintf("amaka_dem_%s_envelope.tif", Cfg.DEM.Source))
if (!file.exists(DEM.Src.Path)) {
  found <- list.files(Dir.GEE.Exports, pattern = "dem.*\\.tif$", ignore.case = TRUE)
  stop("DEM export not found: ", DEM.Src.Path,
       "\nRun gee/02_dem_amaka_export.js, then download the GeoTIFF from Drive folder 'Amaka_GEE_Exports' ",
       "into data/gee_exports/.",
       ifelse(length(found), paste0("\nDEM-like files already there: ", paste(found, collapse = ", ")), ""),
       call. = FALSE)
}

Envelope.Path <- file.path(Dir.Boundary, sprintf("Amaka_Envelope_%gkm.shp", Cfg.Buffer.m / 1000))
Amaka.Envelope <- sf::st_read(Envelope.Path, quiet = TRUE)
Amaka.UTM      <- sf::st_read(Path.Amaka.UTM, quiet = TRUE)
Amaka.Buffer   <- sf::st_read(file.path(Dir.Boundary, sprintf("Amaka_Buffer_%gkm.shp", Cfg.Buffer.m / 1000)), quiet = TRUE)


# ---- PHASE 2: Load, reproject, clip ----------------------------------------
log_phase("PHASE 2: load, reproject, clip")

DEM.Raw <- terra::rast(DEM.Src.Path)[[1]]
names(DEM.Raw) <- "elev"
Src.EPSG <- suppressWarnings(as.integer(terra::crs(DEM.Raw, describe = TRUE)$code))
message(sprintf("  source: EPSG:%s, %g x %g %s, %d x %d cells",
                Src.EPSG, terra::res(DEM.Raw)[1], terra::res(DEM.Raw)[2],
                ifelse(isTRUE(terra::is.lonlat(DEM.Raw)), "deg", "m"), terra::ncol(DEM.Raw), terra::nrow(DEM.Raw)))

if (is.na(Src.EPSG) || Src.EPSG != Cfg.CRS.Metric) {
  message("  reprojecting to EPSG:", Cfg.CRS.Metric, " (bilinear, 30 m)")
  DEM.UTM <- terra::project(DEM.Raw, paste0("EPSG:", Cfg.CRS.Metric), method = "bilinear", res = 30)
} else {
  DEM.UTM <- DEM.Raw
}

# The DEM must cover the whole envelope or the upslope area is silently truncated.
Env.Ext <- terra::ext(terra::vect(Amaka.Envelope))
Dem.Ext <- terra::ext(DEM.UTM)
Covers  <- Dem.Ext$xmin <= Env.Ext$xmin + 30 & Dem.Ext$xmax >= Env.Ext$xmax - 30 &
           Dem.Ext$ymin <= Env.Ext$ymin + 30 & Dem.Ext$ymax >= Env.Ext$ymax - 30
if (!Covers) {
  warning("The DEM does not fully cover the ", Cfg.Buffer.m / 1000, " km envelope; ",
          "re-export with a larger EXPORT_PAD_M in gee/02.", call. = FALSE)
}
DEM.UTM <- terra::crop(DEM.UTM, terra::vect(Amaka.Envelope), snap = "out")


# ---- PHASE 3: Quality checks -------------------------------------------------
log_phase("PHASE 3: quality checks")

Elev.Stats <- terra::global(DEM.UTM, c("min", "max", "mean", "sd"), na.rm = TRUE)
Na.Frac    <- terra::global(is.na(DEM.UTM), "mean")[[1]]
message(sprintf("  elevation %.0f-%.0f m (mean %.0f, sd %.0f) | NA cells: %.2f%%",
                Elev.Stats$min, Elev.Stats$max, Elev.Stats$mean, Elev.Stats$sd, 100 * Na.Frac))
if (Na.Frac > 0.01) stop("More than 1% of the DEM is NA - check the GEE export.", call. = FALSE)
if (Elev.Stats$min < -20 || Elev.Stats$max > 2500) {
  warning("Elevation range looks implausible for Taita-Taveta lowlands - check units / source.", call. = FALSE)
}


# ---- PHASE 3b: Check against the elevations read off the phone in the field ------
# A free validation (known_issues F1). The phone's GPS altitude is good to ~10-30 m and GLO-30 is a surface
# model (canopy included), so only gross disagreement means anything; the table is for the record.
log_phase("PHASE 3b: DEM versus field elevations")

if (file.exists(Path.FieldPoints)) {
  Field.Elev <- readr::read_csv(Path.FieldPoints, show_col_types = FALSE)
  Field.Elev <- Field.Elev[!is.na(Field.Elev$latitude) & !is.na(Field.Elev$elevation_m), ]
  Pts.Rec <- sf::st_transform(sf::st_as_sf(Field.Elev, coords = c("longitude", "latitude"), crs = Cfg.CRS.Geo), Cfg.CRS.Metric)
  Field.Elev$dem_at_recorded_m <- terra::extract(DEM.UTM, terra::vect(Pts.Rec))[, 2]
  has.img <- !is.na(Field.Elev$lat_img)
  Field.Elev$dem_at_image_m <- NA_real_
  if (any(has.img)) {
    Pts.Img <- sf::st_transform(sf::st_as_sf(Field.Elev[has.img, ], coords = c("lon_img", "lat_img"), crs = Cfg.CRS.Geo), Cfg.CRS.Metric)
    Field.Elev$dem_at_image_m[has.img] <- terra::extract(DEM.UTM, terra::vect(Pts.Img))[, 2]
  }
  Field.Elev$diff_recorded_m <- round(Field.Elev$dem_at_recorded_m - Field.Elev$elevation_m, 1)
  Field.Elev$diff_image_m    <- round(Field.Elev$dem_at_image_m - Field.Elev$elevation_m, 1)
  Elev.Check <- Field.Elev[, c("id", "name", "position_type", "elevation_m", "dem_at_recorded_m", "diff_recorded_m",
                               "dem_at_image_m", "diff_image_m")]
  readr::write_csv(Elev.Check, file.path(Dir.Out.Tables, "02_dem_vs_field_elevations.csv"))
  print(as.data.frame(Elev.Check |> dplyr::mutate(name = substr(name, 1, 24), dplyr::across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)
  message(sprintf("  DEM minus phone altitude at the recorded points: mean %+.1f m, RMSE %.1f m (n = %d)",
                  mean(Field.Elev$diff_recorded_m), sqrt(mean(Field.Elev$diff_recorded_m^2)), nrow(Field.Elev)))
} else {
  message("  no field points file; skipped")
}


# ---- PHASE 4: Optional smoothing -------------------------------------------
log_phase("PHASE 4: optional smoothing")

Dem.Out.Path <- file.path(Dir.DEM, "amaka_dem_utm37s.tif")
if (isTRUE(Cfg.DEM.Smooth)) {
  Tmp.In  <- file.path(Dir.Out.Hydro, "dem_prep_unsmoothed.tif")
  Tmp.Out <- file.path(Dir.Out.Hydro, "dem_prep_smoothed.tif")
  terra::writeRaster(DEM.UTM, Tmp.In, overwrite = TRUE)
  whitebox::wbt_gaussian_filter(input = Tmp.In, output = Tmp.Out, sigma = Cfg.DEM.Sigma)
  DEM.UTM <- terra::rast(Tmp.Out)
  names(DEM.UTM) <- "elev"
  message("  Gaussian filter applied, sigma = ", Cfg.DEM.Sigma)
} else {
  message("  no smoothing (Cfg.DEM.Smooth = FALSE)")
}


# ---- PHASE 5: Write ------------------------------------------------------------
log_phase("PHASE 5: write DEM and log")

terra::writeRaster(DEM.UTM, Dem.Out.Path, overwrite = TRUE, datatype = "FLT4S",
                   gdal = c("COMPRESS=LZW"))

Dem.Log <- tibble::tibble(
  run_date    = format(Sys.Date()),
  source      = Cfg.DEM.Source,
  source_file = basename(DEM.Src.Path),
  source_epsg = Src.EPSG,
  reprojected = is.na(Src.EPSG) || Src.EPSG != Cfg.CRS.Metric,
  smoothed    = isTRUE(Cfg.DEM.Smooth),
  sigma       = ifelse(isTRUE(Cfg.DEM.Smooth), Cfg.DEM.Sigma, NA_real_),
  buffer_m    = Cfg.Buffer.m,
  res_m       = terra::res(DEM.UTM)[1],
  ncol        = terra::ncol(DEM.UTM), nrow = terra::nrow(DEM.UTM),
  elev_min    = Elev.Stats$min, elev_max = Elev.Stats$max, elev_mean = round(Elev.Stats$mean, 1),
  na_pct      = round(100 * Na.Frac, 3),
  covers_envelope = Covers
)
readr::write_csv(Dem.Log, file.path(Dir.Out.Tables, "02_dem_prep_log.csv"))
message("  wrote ", Dem.Out.Path)


# ---- PHASE 6: DEM clip map ----------------------------------------------------
log_phase("PHASE 6: DEM clip map")

# Hillshade for the base map only; flow routing stays with WhiteboxTools (4.3).
Slope.Rad  <- terra::terrain(DEM.UTM, "slope",  unit = "radians")
Aspect.Rad <- terra::terrain(DEM.UTM, "aspect", unit = "radians")
Hill.Map   <- terra::shade(Slope.Rad, Aspect.Rad, angle = 40, direction = 315)
names(Hill.Map) <- "shade"

Map.DEM <-
  tmap::tm_shape(Hill.Map) +
  tmap::tm_raster(col.scale = tmap::tm_scale_continuous(values = "brewer.greys"),
                  col.legend = tmap::tm_legend_hide()) +
  tmap::tm_shape(DEM.UTM) +
  tmap::tm_raster(col = "elev", col_alpha = 0.55,
                  col.scale  = tmap::tm_scale_continuous(values = "scico.batlow"),
                  col.legend = tmap::tm_legend(title = "Elevation (m)")) +
  tmap::tm_shape(Amaka.Buffer) + tmap::tm_borders(col = "grey20", lwd = 1, lty = "dashed") +
  tmap::tm_shape(Amaka.Envelope) + tmap::tm_borders(col = "blue", lwd = 1.5) +
  tmap::tm_shape(Amaka.UTM) + tmap::tm_borders(col = "black", lwd = 2.5) +
  tmap::tm_title(sprintf("Amaka DEM clip (%s, %g km envelope)", toupper(Cfg.DEM.Source), Cfg.Buffer.m / 1000),
                 size = 0.9) +
  tmap::tm_credits("Black: Amaka | dashed: buffer | blue: DEM clip extent", position = c("left", "bottom"), size = 0.6)

Fig.Path <- file.path(Dir.Out.Figures, "02_dem_clip_map.png")
tmap::tmap_save(Map.DEM, Fig.Path, width = 7, height = 7, dpi = 200)
message("  wrote ", basename(Fig.Path))

log_phase("DONE")
