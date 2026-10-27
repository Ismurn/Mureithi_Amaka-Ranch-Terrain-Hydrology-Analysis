# ==============================================================================
# 01_boundary_prep.R  -  Extract Amaka, verify CRS, build the catchment envelope
# ==============================================================================
# PURPOSE
#   (1) Verify the coordinate system of the Tsavo ranch / settled-area shapefiles
#       BEFORE anything is built on them (CLAUDE.md 9.1).
#   (2) Isolate Amaka Ranch and Conservancy from RanchesTTWCA.shp. The ranch is
#       catalogued there under its old name, R_NAME == "MUSAMITI".
#   (3) Build the AOI per Cfg.AOI.Mode: "musamiti" (current) is the TTWCA polygon exactly as
#       delivered; "file" is a boundary supplied by Amaka management. The TTWCA polygon (4,782 ha)
#       is ~20% below the 14,815 acres TTWCA lists for Amaka and 3 of the 7 field points fall
#       outside it (docs/known_issues.md F8); this is accepted for now. Write the AOI (EPSG:32737
#       for metric work, EPSG:4326 for maps and GEE) and the 7 km buffered catchment envelope.
#   (4) Write the context layers used in every map (Tsavo East / West NPs,
#       settled areas, TTWCA ranches).
#   (5) Refresh the boundary block inside gee/01_amaka_boundary_export.js so the
#       GEE boundary asset is built from exactly this geometry.
#
# INPUTS  (paths declared in 00_setup.R)
#   Path.Ranches, Path.SettledAreas, Path.TsavoEast, Path.TsavoWest
#
# OUTPUTS
#   data/boundary/Amaka_Ranch.shp              canonical AOI, EPSG:32737
#   data/boundary/Amaka_Ranch_wgs84.shp        EPSG:4326 (maps, GEE upload)
#   data/boundary/Amaka_Buffer_<N>km.shp       buffered AOI, EPSG:32737
#   data/boundary/Amaka_Envelope_<N>km.shp     bbox of the buffer (DEM clip extent)
#   data/boundary/Amaka_Envelope_<N>km_wgs84.shp
#   data/boundary/context/*.shp                NPs, settled areas, ranches (EPSG:32737)
#   outputs/figures/crs_check/*.kml            Amaka under each CRS hypothesis
#   outputs/figures/01_crs_verification.png    overview + zoom, both hypotheses
#   outputs/figures/01_aoi.png                 the AOI, the TTWCA polygon and the field points
#   data/boundary/context/Musamiti_TTWCA_polygon.shp   the TTWCA polygon as delivered
#   outputs/tables/01_crs_np_edge_test.csv     numeric edge test (Phase 2)
#   outputs/tables/01_crs_corner_check.csv     corner coordinates for Google Earth
#   outputs/tables/01_boundary_summary.csv
#   gee/01_amaka_boundary_export.js            boundary block refreshed in place
#
# RUNTIME NOTES
#   * Seconds. Re-run whenever Cfg.Ranches.CRS.Assumed or Cfg.Buffer.m changes.
#   * The CRS is NOT silently "fixed". Phase 2 gives numeric evidence; the
#     decision is made by looking at the KML overlays in Google Earth, then
#     setting Cfg.Ranches.CRS.Assumed / Cfg.CRS.Confirmed in 00_setup.R.
# ==============================================================================

if (!exists("Setup.Loaded")) source("R/00_setup.R")


# ---- PHASE 1: Load source layers --------------------------------------------
log_phase("PHASE 1: load source layers")

Ranches.Raw <- sf::st_read(Path.Ranches,      quiet = TRUE)   # as delivered (file tag untouched)
Settled.Raw <- sf::st_read(Path.SettledAreas, quiet = TRUE)
TsavoE.Raw  <- sf::st_read(Path.TsavoEast,    quiet = TRUE)   # WDPA, WGS84
TsavoW.Raw  <- sf::st_read(Path.TsavoWest,    quiet = TRUE)

CRS.Tags <- tibble::tibble(
  layer    = c("RanchesTTWCA", "SettledAreas", "TsavoEast (WDPA)", "TsavoWest (WDPA)"),
  epsg_tag = c(sf::st_crs(Ranches.Raw)$epsg, sf::st_crs(Settled.Raw)$epsg,
               sf::st_crs(TsavoE.Raw)$epsg,  sf::st_crs(TsavoW.Raw)$epsg),
  n        = c(nrow(Ranches.Raw), nrow(Settled.Raw), nrow(TsavoE.Raw), nrow(TsavoW.Raw))
)
print(CRS.Tags)

# The WDPA parks are the independent reference below, so they must be WGS84.
stopifnot("WDPA reference layers are not tagged EPSG:4326" =
            all(CRS.Tags$epsg_tag[3:4] == Cfg.CRS.Geo))

# What the delivered ranch file claims to be (H1) versus the known failure (H2).
Tag.Ranches <- sf::st_crs(Ranches.Raw)$epsg
if (is.na(Tag.Ranches)) Tag.Ranches <- Cfg.CRS.Metric

# Read the numbers of `geom` as if they were in `epsg.assumed`, return them in UTM 37S.
as_utm37 <- function(geom, epsg.assumed) {
  g <- suppressWarnings(sf::st_set_crs(sf::st_geometry(geom), epsg.assumed))
  sf::st_transform(g, Cfg.CRS.Metric)
}


# ---- PHASE 2: CRS verification, part 1 - numeric edge test ------------------
# Evidence, not proof. Where a ranch edge is also a park edge (many TTWCA ranches
# abut Tsavo East / West), the ranch vertices should sit ON the WDPA park boundary
# if the CRS assumption is right. Compare the two hypotheses on the SAME vertices.
# Caveat: the two sources are digitised independently, so absolute distances are
# noisy; only the DIFFERENCE between H1 and H2 carries information.
log_phase("PHASE 2: CRS verification (numeric edge test)")

TsavoNP.UTM <- c(sf::st_geometry(sf::st_transform(TsavoE.Raw, Cfg.CRS.Metric)),
                 sf::st_geometry(sf::st_transform(TsavoW.Raw, Cfg.CRS.Metric))) |>
  sf::st_make_valid() |>
  sf::st_union()
NP.Edge <- sf::st_boundary(TsavoNP.UTM)

Ranch.Vertices <- sf::st_cast(sf::st_cast(sf::st_geometry(Ranches.Raw), "MULTIPOLYGON"),
                              "POINT", warn = FALSE)
Vtx.H1 <- as_utm37(Ranch.Vertices, Tag.Ranches)           # H1: coordinates are what the file tag says
Vtx.H2 <- as_utm37(Ranch.Vertices, Cfg.CRS.Arc1960)       # H2: coordinates are really Arc 1960 / UTM 37S

Dist.H1 <- as.numeric(sf::st_distance(Vtx.H1, NP.Edge))
Dist.H2 <- as.numeric(sf::st_distance(Vtx.H2, NP.Edge))
Near    <- Dist.H1 < 1500 | Dist.H2 < 1500                # same vertex set for both hypotheses

CRS.Edge.Test <- tibble::tibble(
  hypothesis      = c(paste0("H1: file tag EPSG:", Tag.Ranches), "H2: Arc 1960 / UTM 37S (EPSG:21037)"),
  n_vertices      = sum(Near),
  median_dist_m   = round(c(median(Dist.H1[Near]), median(Dist.H2[Near]))),
  pct_within_100m = round(100 * c(mean(Dist.H1[Near] < 100), mean(Dist.H2[Near] < 100)), 1)
)
print(CRS.Edge.Test)

# Rule 1: share of vertices within 100 m of the park edge.
Rule.Share <- dplyr::case_when(
  CRS.Edge.Test$pct_within_100m[2] > 2 * CRS.Edge.Test$pct_within_100m[1] &
    CRS.Edge.Test$pct_within_100m[2] - CRS.Edge.Test$pct_within_100m[1] > 5 ~ "H2",
  CRS.Edge.Test$pct_within_100m[1] > 2 * CRS.Edge.Test$pct_within_100m[2] &
    CRS.Edge.Test$pct_within_100m[1] - CRS.Edge.Test$pct_within_100m[2] > 5 ~ "H1",
  TRUE ~ "none"
)
# Rule 2 (stricter): vertices ON the park line (< 60 m) under one hypothesis but far from it
# (> 250 m) under the other. These are the only places where a ~300 m datum shift is
# unambiguous. A lean is called only when BOTH rules agree.
Disc.H1 <- sum(Near & Dist.H1 < 60 & Dist.H2 > 250)
Disc.H2 <- sum(Near & Dist.H2 < 60 & Dist.H1 > 250)
Rule.Disc <- dplyr::case_when(
  Disc.H2 >= 3 & Disc.H2 >= 2 * Disc.H1 ~ "H2",
  Disc.H1 >= 3 & Disc.H1 >= 2 * Disc.H2 ~ "H1",
  TRUE ~ "none"
)
CRS.Edge.Lean <- dplyr::case_when(
  Rule.Share == "H2" & Rule.Disc == "H2" ~ "leans H2 (Arc 1960)",
  Rule.Share == "H1" & Rule.Disc == "H1" ~ "leans H1 (file tag)",
  TRUE ~ "inconclusive"
)
message(sprintf("  Edge test %s | rule 1 (share within 100 m): %s | rule 2 (vertices fitting only one hypothesis): H1 %d vs H2 %d -> %s",
                CRS.Edge.Lean, Rule.Share, Disc.H1, Disc.H2, Rule.Disc))
message("  Indicative only: ranch and park edges are digitised independently. A visual check on imagery at zoom >= 14 decides.")
readr::write_csv(dplyr::mutate(CRS.Edge.Test, n_fit_only_this_hypothesis = c(Disc.H1, Disc.H2),
                               rule_share = Rule.Share, rule_discriminating = Rule.Disc, verdict = CRS.Edge.Lean),
                 file.path(Dir.Out.Tables, "01_crs_np_edge_test.csv"))


# ---- PHASE 3: Isolate Amaka (Musamiti) --------------------------------------
log_phase("PHASE 3: isolate Amaka Ranch")

Ranches.UTM <- sf::st_transform(read_tsavo_layer(Path.Ranches), Cfg.CRS.Metric)

# Normalise names: upper case, trimmed, internal runs of spaces collapsed
# (the file holds e.g. "MUGENO  RANCH" with a double space).
norm_name <- function(x) gsub("\\s+", " ", toupper(trimws(as.character(x))))

Amaka.Match <- norm_name(Ranches.UTM[[Cfg.Amaka.Field]]) %in% norm_name(Cfg.Amaka.Names)
if (!any(Amaka.Match)) {
  Name.Pool <- sort(unique(norm_name(Ranches.UTM[[Cfg.Amaka.Field]])))
  Close     <- agrep(Cfg.Amaka.Names[1], Name.Pool, max.distance = 0.3, value = TRUE)
  stop("No feature with ", Cfg.Amaka.Field, " in {", paste(Cfg.Amaka.Names, collapse = ", "),
       "}.\nClosest names: ", ifelse(length(Close), paste(Close, collapse = ", "), "none"),
       "\nAll names: ", paste(Name.Pool, collapse = " | "),
       "\nFix Cfg.Amaka.Names in R/00_setup.R.", call. = FALSE)
}

Amaka.Src <- sf::st_zm(Ranches.UTM[Amaka.Match, ], drop = TRUE, what = "ZM")
message("  Matched ", nrow(Amaka.Src), " feature(s): ",
        paste(unique(as.character(Amaka.Src[[Cfg.Amaka.Field]])), collapse = "; "))

# Dissolve the TTWCA polygon to one valid MULTIPOLYGON.
Amaka.Geom <- sf::st_make_valid(sf::st_union(sf::st_geometry(Amaka.Src)))
if (any(sf::st_geometry_type(Amaka.Geom) == "GEOMETRYCOLLECTION")) {
  Amaka.Geom <- sf::st_collection_extract(Amaka.Geom, "POLYGON")
}
Amaka.Geom <- sf::st_cast(Amaka.Geom, "MULTIPOLYGON")
Amaka.Geom.TTWCA <- Amaka.Geom        # the polygon exactly as TTWCA draws it (kept for reference and maps)

Musamiti.Area.ha <- as.numeric(sf::st_area(Amaka.Geom.TTWCA)) / 1e4
Src.Area.ha      <- suppressWarnings(max(as.numeric(unlist(
  sf::st_drop_geometry(Amaka.Src)[intersect(c("HECTARES", "AREA_HA"), names(Amaka.Src))])), na.rm = TRUE))
message(sprintf("  TTWCA Musamiti polygon: %.0f ha from geometry | %s ha in the source table",
                Musamiti.Area.ha, ifelse(is.finite(Src.Area.ha), sprintf("%.0f", Src.Area.ha), "n/a")))

# The AOI proper, per Cfg.AOI.Mode (see 00_setup.R): musamiti / file.
AOI           <- build_aoi(Amaka.Geom.TTWCA)
Amaka.Geom    <- AOI$geom
Amaka.Area.ha <- as.numeric(sf::st_area(Amaka.Geom)) / 1e4
Listed.ha     <- Cfg.AOI.Listed.acres * 0.40468564
message(sprintf("  AOI mode '%s' (%s): %.0f ha = %+.0f%% against the %.0f ha (%s acres) TTWCA lists for Amaka",
                AOI$mode, AOI$src, Amaka.Area.ha, 100 * (Amaka.Area.ha / Listed.ha - 1), Listed.ha,
                format(Cfg.AOI.Listed.acres, big.mark = ",")))
if (AOI$mode == "musamiti") {
  message("  Note: the TTWCA polygon is coarse and smaller than the size TTWCA lists (docs/known_issues.md F8).")
}

Amaka.UTM <- sf::st_sf(
  ranch_id  = "AMAKA",
  name      = AOI$label,
  src_rname = paste(unique(as.character(Amaka.Src[[Cfg.Amaka.Field]])), collapse = ";"),
  area_ha   = round(Amaka.Area.ha, 1),
  aoi_mode  = AOI$mode,
  aoi_src   = AOI$src,
  crs_epsg  = Cfg.Ranches.CRS.Assumed,
  crs_conf  = ifelse(Cfg.CRS.Confirmed, "yes", "no"),
  geometry  = Amaka.Geom
)
Amaka.WGS84  <- sf::st_transform(Amaka.UTM, Cfg.CRS.Geo)
Musamiti.UTM <- sf::st_sf(r_name = "MUSAMITI", geometry = Amaka.Geom.TTWCA)


# ---- PHASE 3b: CRS verification, part 2 - Google Earth overlay --------------
# The decisive check. A datum mix-up is a (near) pure translation, so the distance
# between two corners is unchanged and cannot reveal it; position can. Open both
# KML files in Google Earth: the outline that sits on the real boundary features
# (fence line, track, change in vegetation) is the correct interpretation.
log_phase("PHASE 3b: Google Earth overlay for the CRS decision")

Amaka.Raw.Geom <- sf::st_geometry(Ranches.Raw)[Amaka.Match] |> sf::st_union()
Amaka.H1.WGS   <- sf::st_transform(as_utm37(Amaka.Raw.Geom, Tag.Ranches),      Cfg.CRS.Geo)
Amaka.H2.WGS   <- sf::st_transform(as_utm37(Amaka.Raw.Geom, Cfg.CRS.Arc1960), Cfg.CRS.Geo)

Dir.CRS.Check <- file.path(Dir.Out.Figures, "crs_check")
dir.create(Dir.CRS.Check, showWarnings = FALSE, recursive = TRUE)
write_kml <- function(geom, label, file) {
  tryCatch({
    sf::st_write(sf::st_sf(name = label, geometry = geom), file, driver = "KML",
                 delete_dsn = file.exists(file), quiet = TRUE)
    message("  wrote ", basename(file))
  }, error = function(e) message("  KML not written (", basename(file), "): ", conditionMessage(e)))
}
write_kml(Amaka.H1.WGS, paste0("Amaka if coordinates are EPSG:", Tag.Ranches, " (file tag)"),
          file.path(Dir.CRS.Check, paste0("Amaka_H1_EPSG", Tag.Ranches, ".kml")))
write_kml(Amaka.H2.WGS, "Amaka if coordinates are Arc 1960 / UTM 37S (EPSG:21037)",
          file.path(Dir.CRS.Check, "Amaka_H2_Arc1960_EPSG21037.kml"))

# Four extreme vertices (N, S, E, W) under both hypotheses, lat/lon for pasting
# straight into the Google Earth search box.
Amaka.Pts.Raw <- sf::st_cast(sf::st_cast(Amaka.Raw.Geom, "MULTIPOLYGON"), "POINT", warn = FALSE)
Amaka.XY      <- sf::st_coordinates(Amaka.Pts.Raw)[, 1:2]
Corner.Idx    <- c(North = which.max(Amaka.XY[, 2]), South = which.min(Amaka.XY[, 2]),
                   East  = which.max(Amaka.XY[, 1]), West  = which.min(Amaka.XY[, 1]))
Corner.H1 <- sf::st_transform(as_utm37(Amaka.Pts.Raw[Corner.Idx], Tag.Ranches),     Cfg.CRS.Geo)
Corner.H2 <- sf::st_transform(as_utm37(Amaka.Pts.Raw[Corner.Idx], Cfg.CRS.Arc1960), Cfg.CRS.Geo)

CRS.Corner.Check <- tibble::tibble(
  vertex   = names(Corner.Idx),
  H1_lat   = round(sf::st_coordinates(Corner.H1)[, 2], 6),
  H1_lon   = round(sf::st_coordinates(Corner.H1)[, 1], 6),
  H2_lat   = round(sf::st_coordinates(Corner.H2)[, 2], 6),
  H2_lon   = round(sf::st_coordinates(Corner.H2)[, 1], 6),
  shift_m  = round(as.numeric(sf::st_distance(Corner.H1, Corner.H2, by_element = TRUE)))
)
print(CRS.Corner.Check)
readr::write_csv(CRS.Corner.Check, file.path(Dir.Out.Tables, "01_crs_corner_check.csv"))
message("  H1 vs H2 differ by ~", round(mean(CRS.Corner.Check$shift_m)),
        " m at Amaka. Anything above ~50 m means the tag must be settled before use (9.1).")


# ---- PHASE 4: Buffered catchment envelope -----------------------------------
# Hydrology needs the upslope contributing area, so the DEM is clipped to the
# bounding box of a buffer around the ranch, never to the ranch itself.
log_phase("PHASE 4: buffered catchment envelope")

Buffer.km     <- Cfg.Buffer.m / 1000
Amaka.Buffer  <- sf::st_buffer(Amaka.UTM, dist = Cfg.Buffer.m)
Envelope.BBox <- sf::st_bbox(Amaka.Buffer)
Amaka.Envelope <- sf::st_sf(
  buffer_m = Cfg.Buffer.m,
  geometry = sf::st_as_sfc(Envelope.BBox)
)
Amaka.Envelope.WGS84 <- sf::st_transform(Amaka.Envelope, Cfg.CRS.Geo)

message(sprintf("  Envelope (EPSG:%d): x %.0f-%.0f, y %.0f-%.0f  (%.1f x %.1f km)",
                Cfg.CRS.Metric, Envelope.BBox["xmin"], Envelope.BBox["xmax"],
                Envelope.BBox["ymin"], Envelope.BBox["ymax"],
                (Envelope.BBox["xmax"] - Envelope.BBox["xmin"]) / 1000,
                (Envelope.BBox["ymax"] - Envelope.BBox["ymin"]) / 1000))


# ---- PHASE 5: Write boundary + context layers -------------------------------
log_phase("PHASE 5: write boundary and context layers")

write_shp(Amaka.UTM,   Path.Amaka.UTM)
write_shp(Amaka.WGS84, Path.Amaka.WGS84)
write_shp(Amaka.Buffer,   file.path(Dir.Boundary, sprintf("Amaka_Buffer_%gkm.shp",   Buffer.km)))
write_shp(Amaka.Envelope, file.path(Dir.Boundary, sprintf("Amaka_Envelope_%gkm.shp", Buffer.km)))
write_shp(Amaka.Envelope.WGS84, file.path(Dir.Boundary, sprintf("Amaka_Envelope_%gkm_wgs84.shp", Buffer.km)))

# Context layers, all EPSG:32737, trimmed to a few attributes.
trim_ctx <- function(x, keep) {
  x <- sf::st_zm(x, drop = TRUE, what = "ZM")
  x <- x[, intersect(keep, names(x))]
  sf::st_cast(sf::st_make_valid(x), "MULTIPOLYGON")
}
Ctx.Ranches <- trim_ctx(Ranches.UTM, c("R_NAME", "OWNERSHIP"))
Ctx.Settled <- trim_ctx(sf::st_transform(read_tsavo_layer(Path.SettledAreas), Cfg.CRS.Metric), "R_NAME")
Ctx.TsavoE  <- trim_ctx(sf::st_transform(TsavoE.Raw, Cfg.CRS.Metric), c("NAME", "WDPAID", "DESIG_ENG", "IUCN_CAT"))
Ctx.TsavoW  <- trim_ctx(sf::st_transform(TsavoW.Raw, Cfg.CRS.Metric), c("NAME", "WDPAID", "DESIG_ENG", "IUCN_CAT"))

write_shp(Ctx.Ranches, file.path(Dir.Boundary.Ctx, "Tsavo_Ranches.shp"))
write_shp(Ctx.Settled, file.path(Dir.Boundary.Ctx, "Tsavo_Settled_Areas.shp"))
write_shp(Ctx.TsavoE,  file.path(Dir.Boundary.Ctx, "Tsavo_East_NP.shp"))
write_shp(Ctx.TsavoW,  file.path(Dir.Boundary.Ctx, "Tsavo_West_NP.shp"))
# The TTWCA polygon as delivered, kept apart from the AOI so it is never lost.
write_shp(Musamiti.UTM, file.path(Dir.Boundary.Ctx, "Musamiti_TTWCA_polygon.shp"))

# Summary table
Boundary.Summary <- tibble::tibble(
  ranch_id         = "AMAKA",
  aoi_mode         = AOI$mode,
  aoi_source       = AOI$src,
  source_name      = Amaka.UTM$src_rname,
  area_ha_geometry = round(Amaka.Area.ha, 1),
  area_ha_musamiti_polygon = round(Musamiti.Area.ha, 1),
  area_ha_listed_ttwca     = round(Listed.ha, 1),
  aoi_vs_listed_pct        = round(100 * (Amaka.Area.ha / Listed.ha - 1), 1),
  area_ha_source   = ifelse(is.finite(Src.Area.ha), Src.Area.ha, NA_real_),
  n_vertices       = nrow(sf::st_coordinates(Amaka.Geom)),
  buffer_m         = Cfg.Buffer.m,
  envelope_xmin    = unname(Envelope.BBox["xmin"]), envelope_xmax = unname(Envelope.BBox["xmax"]),
  envelope_ymin    = unname(Envelope.BBox["ymin"]), envelope_ymax = unname(Envelope.BBox["ymax"]),
  crs_assumed      = Cfg.Ranches.CRS.Assumed,
  crs_confirmed    = Cfg.CRS.Confirmed,
  h1_h2_shift_m    = round(mean(CRS.Corner.Check$shift_m)),
  edge_test        = CRS.Edge.Lean
)
readr::write_csv(Boundary.Summary, file.path(Dir.Out.Tables, "01_boundary_summary.csv"))
print(t(Boundary.Summary))


# ---- PHASE 6: Refresh the GEE boundary block --------------------------------
# gee/01_amaka_boundary_export.js carries the AOI as literal coordinates between two
# marker comments, so the boundary asset in GEE is built from this exact geometry
# (no manual shapefile upload). Everything outside the markers is left untouched.
log_phase("PHASE 6: refresh gee/01_amaka_boundary_export.js")

ring_js  <- function(m) paste0("[", paste(sprintf("[%.6f,%.6f]", m[, 1], m[, 2]), collapse = ","), "]")
poly_js  <- function(p) paste0("    [", paste(vapply(p, ring_js, character(1)), collapse = ",\n     "), "]")
coords_js <- function(multipoly) paste0("[\n", paste(vapply(multipoly, poly_js, character(1)), collapse = ",\n"), "\n  ]")

update_gee_boundary <- function(js.path, geom.wgs84, area.ha, aoi) {
  if (!file.exists(js.path)) { warning("GEE script not found: ", js.path); return(invisible(FALSE)) }
  txt   <- readLines(js.path, warn = FALSE)
  begin <- grep("AMAKA_BOUNDARY_BEGIN", txt, fixed = TRUE)
  end   <- grep("AMAKA_BOUNDARY_END",   txt, fixed = TRUE)
  if (length(begin) != 1 || length(end) != 1 || end <= begin) {
    warning("Marker comments missing in ", js.path); return(invisible(FALSE))
  }
  g <- sf::st_geometry(sf::st_cast(geom.wgs84, "MULTIPOLYGON"))[[1]]
  block <- c(
    sprintf("// >>> AMAKA_BOUNDARY_BEGIN  (written by R/01_boundary_prep.R on %s - do not edit by hand)", format(Sys.Date())),
    sprintf("var AMAKA_CRS_STATUS = '%s';", ifelse(Cfg.CRS.Confirmed, "CONFIRMED", "UNCONFIRMED")),
    sprintf("var AMAKA_AOI_MODE   = '%s';   // MUSAMITI = the TTWCA polygon as delivered; FILE = a supplied boundary", toupper(aoi$mode)),
    sprintf("var AMAKA_AOI_SRC    = '%s';", aoi$src),
    sprintf("var AMAKA_CRS_EPSG   = %d;      // CRS the source shapefile coordinates were taken to be", Cfg.Ranches.CRS.Assumed),
    sprintf("var AMAKA_R_AREA_HA  = %.1f;   // planar area in EPSG:32737, for the sanity check in GEE", area.ha),
    sprintf("var AMAKA_BUFFER_M_R = %d;     // buffer used in R (Cfg.Buffer.m)", as.integer(Cfg.Buffer.m)),
    paste0("var AMAKA_COORDS = ", coords_js(g), ";"),
    "// <<< AMAKA_BOUNDARY_END"
  )
  writeLines(c(txt[seq_len(begin - 1)], block, txt[-seq_len(end)]), js.path, useBytes = TRUE)
  message(sprintf("  wrote %d vertices into %s", nrow(sf::st_coordinates(g)), basename(js.path)))
  invisible(TRUE)
}
update_gee_boundary(file.path(Dir.GEE, "01_amaka_boundary_export.js"), Amaka.WGS84, Amaka.Area.ha, AOI)


# ---- PHASE 7: CRS verification figure ---------------------------------------
log_phase("PHASE 7: CRS verification figure")

Hyp.Lines <- sf::st_sf(
  hypothesis = factor(c(paste0("H1: file tag (EPSG:", Tag.Ranches, ")"), "H2: Arc 1960 (EPSG:21037)")),
  geometry   = c(sf::st_transform(Amaka.H1.WGS, Cfg.CRS.Metric), sf::st_transform(Amaka.H2.WGS, Cfg.CRS.Metric))
)
Ctx.NP <- rbind(Ctx.TsavoE, Ctx.TsavoW)[, 0]

Map.Overview <-
  tmap::tm_shape(Ctx.Ranches) + tmap::tm_polygons(fill = "grey94", col = "grey65", lwd = 0.4) +
  tmap::tm_shape(Ctx.Settled) + tmap::tm_polygons(fill = "khaki", fill_alpha = 0.6, col = "grey60", lwd = 0.3) +
  tmap::tm_shape(Ctx.NP)      + tmap::tm_polygons(fill = "darkseagreen3", fill_alpha = 0.6, col = "darkgreen", lwd = 0.5) +
  tmap::tm_shape(Amaka.UTM)   + tmap::tm_polygons(fill = "orange", fill_alpha = 0.6, col = "black", lwd = 1) +
  tmap::tm_shape(Musamiti.UTM) + tmap::tm_borders(col = "blue", lwd = 1.2, lty = "dashed") +
  tmap::tm_title("A. Amaka AOI (orange) in the TTWCA landscape; TTWCA Musamiti polygon dashed blue", size = 0.8)

Map.Zoom <-
  tmap::tm_shape(Ctx.Ranches, bbox = Envelope.BBox) + tmap::tm_polygons(fill = "grey94", col = "grey65", lwd = 0.4) +
  tmap::tm_shape(Ctx.Settled) + tmap::tm_polygons(fill = "khaki", fill_alpha = 0.6, col = "grey60", lwd = 0.3) +
  tmap::tm_shape(Ctx.NP)      + tmap::tm_polygons(fill = "darkseagreen3", fill_alpha = 0.6, col = "darkgreen", lwd = 0.5) +
  tmap::tm_shape(Hyp.Lines)   + tmap::tm_borders(col = "hypothesis",
                                                col.scale  = tmap::tm_scale_categorical(values = c("black", "red")),
                                                col.legend = tmap::tm_legend(title = "Boundary if the coordinates are:"),
                                                lwd = 2) +
  tmap::tm_title("B. Both CRS hypotheses (envelope extent)", size = 0.9)

Fig.Path <- file.path(Dir.Out.Figures, "01_crs_verification.png")
tmap::tmap_save(tmap::tmap_arrange(Map.Overview, Map.Zoom, ncol = 2), Fig.Path,
                width = 12, height = 6, dpi = 200)
message("  wrote ", basename(Fig.Path))


# ---- PHASE 7b: AOI figure ----------------------------------------------------
# The picture to look at when judging the AOI: the AOI, the TTWCA Musamiti polygon, the
# neighbouring Kambanga polygon, and the field points (how many fall outside the AOI).
log_phase("PHASE 7b: AOI figure")

Kambanga.UTM <- Ctx.Ranches[norm_name(Ctx.Ranches$R_NAME) == "KAMBANGA RANCH", ]
Field.UTM <- NULL
if (file.exists(Path.FieldPoints)) {
  Field.Pts <- readr::read_csv(Path.FieldPoints, show_col_types = FALSE)
  Field.Pts <- Field.Pts[!is.na(Field.Pts$latitude), ]
  Field.UTM <- sf::st_transform(sf::st_as_sf(Field.Pts, coords = c("longitude", "latitude"), crs = Cfg.CRS.Geo),
                                Cfg.CRS.Metric)
}
# Frame the AOI and, when there are field points, every one of them (some lie outside the AOI).
Fig.Frame <- sf::st_geometry(sf::st_buffer(Amaka.UTM, 2500))
if (!is.null(Field.UTM)) Fig.Frame <- c(Fig.Frame, sf::st_buffer(sf::st_union(sf::st_geometry(Field.UTM)), 1500))
Map.AOI <-
  tmap::tm_shape(Ctx.Ranches, bbox = sf::st_bbox(Fig.Frame)) +
  tmap::tm_polygons(fill = "grey94", col = "grey65", lwd = 0.4) +
  tmap::tm_shape(Kambanga.UTM) + tmap::tm_polygons(fill = "khaki", fill_alpha = 0.5, col = "grey40", lwd = 1) +
  tmap::tm_shape(Amaka.UTM)    + tmap::tm_polygons(fill = "orange", fill_alpha = 0.5, col = "black", lwd = 2.2) +
  tmap::tm_shape(Musamiti.UTM) + tmap::tm_borders(col = "blue", lwd = 1.5, lty = "dashed")
if (!is.null(Field.UTM)) {
  Map.AOI <- Map.AOI +
    tmap::tm_shape(Field.UTM) + tmap::tm_symbols(fill = "white", col = "black", size = 0.6) +
    tmap::tm_text("id", size = 0.7, ymod = -0.9)
  Field.Note <- sprintf(" | field points outside the AOI: %d of %d",
                        sum(lengths(sf::st_intersects(Field.UTM, Amaka.UTM)) == 0), nrow(Field.UTM))
} else {
  Field.Note <- ""
}
Map.AOI <- Map.AOI + tmap::tm_title(
  sprintf("Amaka AOI: %.0f ha (%s). Black = AOI, blue dashed = TTWCA Musamiti polygon, yellow = Kambanga (TTWCA)%s",
          Amaka.Area.ha, AOI$src, Field.Note), size = 0.65)
Fig.AOI.Path <- file.path(Dir.Out.Figures, "01_aoi.png")
tmap::tmap_save(Map.AOI, Fig.AOI.Path, width = 9, height = 7, dpi = 200)
message("  wrote ", basename(Fig.AOI.Path))


# ---- PHASE 8: Hand-over message ----------------------------------------------
log_phase("DONE")
message("Amaka AOI (", AOI$mode, "): ", sprintf("%.0f ha", Amaka.Area.ha), " | envelope ", Buffer.km, " km | CRS ",
        ifelse(Cfg.CRS.Confirmed, "CONFIRMED", "UNCONFIRMED"), " (assumed EPSG:", Cfg.Ranches.CRS.Assumed, ")")
if (AOI$mode == "musamiti") {
  message("If Amaka's real boundary arrives: set Cfg.AOI.Mode = 'file' and Cfg.AOI.File in R/00_setup.R, ",
          "re-run this script, then re-paste gee/01 and re-run its asset export.")
}
if (!Cfg.CRS.Confirmed) {
  message("NEXT: open both KMLs in ", Dir.CRS.Check, " in Google Earth, decide, then set ",
          "Cfg.Ranches.CRS.Assumed / Cfg.CRS.Confirmed in R/00_setup.R and re-run this script.")
}
