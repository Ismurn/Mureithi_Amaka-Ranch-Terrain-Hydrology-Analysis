/**
 * 01_amaka_boundary_export.js
 * ----------------------------------------------------------------------------
 * Project : Amaka Ranch - Terrain, Hydrology & Vegetation Analysis
 * Author  : Isaac Mureithi (AWF / OPALS, Univ. of Exeter)
 * Written : 2026-10-26      GEE account: ee-imnduta1
 * Save as : users/imnduta1/Ismurn/Amaka/01_amaka_boundary_export
 *
 * PURPOSE
 *   Create the Amaka Ranch boundary as an Earth Engine asset. Every other script
 *   in this folder reads that asset, so this one is run ONCE (and again only if
 *   the boundary or the CRS decision changes).
 *
 *   NOTE The boundary is TTWCA's Musamiti polygon exactly as delivered (AMAKA_AOI_MODE below).
 *   It is coarse: 4 corners, about 20% smaller than the 14,815 acres TTWCA lists for Amaka, and
 *   3 of the 7 field points lie outside it (docs/known_issues.md F8). If Amaka's real boundary is
 *   supplied, set Cfg.AOI.Mode = 'file' in R/00_setup.R, re-run R/01_boundary_prep.R, paste the
 *   refreshed script here and re-run the asset export (overwrite the asset: every other
 *   script reads it).
 *
 * HOW IT WORKS
 *   The boundary is NOT uploaded by hand. R/01_boundary_prep.R writes the polygon
 *   (WGS84 lon/lat) straight into the block between the AMAKA_BOUNDARY markers
 *   below, so the asset always comes from the same geometry as the R analysis.
 *   Do not edit between the markers; re-run the R script instead.
 *
 * BEFORE YOU RUN
 *   1. Assets tab -> NEW -> Folder, create:  users/imnduta1/Ismurn/Amaka
 *      (Export cannot create the folder itself.)
 *   2. If your assets live under the newer path, change AMAKA_ASSET below, e.g.
 *      'projects/ee-imnduta1/assets/Amaka/Amaka_Ranch'. Use the SAME string in
 *      every other script in this folder.
 *   3. Click Run, check the printed area against the R value, then open the
 *      Tasks tab and press RUN on 'Amaka_Ranch_boundary_asset'.
 *
 * INPUTS    AMAKA_COORDS (below, from R)
 * OUTPUTS   Asset  AMAKA_ASSET  (FeatureCollection, 1 feature)
 *             properties: ranch_id, name, src_rname, aoi_mode, aoi_src, crs_epsg, crs_status, r_area_ha
 * SEASONS / REDUCERS   not applicable (static geometry)
 * ----------------------------------------------------------------------------
 */

// ---- 1. CONFIG ---------------------------------------------------------------
var AMAKA_ASSET      = 'users/imnduta1/Ismurn/Amaka/Amaka_Ranch';
var EXPORT_PAD_M     = 1000;   // extra margin added to the buffered envelope in the other scripts
var RUN_ASSET_EXPORT = true;   // false = preview only, no export task is created

// ---- 2. BOUNDARY (written by R - do not edit between the markers) --------------
// >>> AMAKA_BOUNDARY_BEGIN  (written by R/01_boundary_prep.R on 2026-10-27 - do not edit by hand)
var AMAKA_CRS_STATUS = 'CONFIRMED';
var AMAKA_AOI_MODE   = 'MUSAMITI';   // MUSAMITI = the TTWCA polygon as delivered; FILE = a supplied boundary
var AMAKA_AOI_SRC    = 'TTWCA Musamiti polygon as delivered';
var AMAKA_CRS_EPSG   = 32737;      // CRS the source shapefile coordinates were taken to be
var AMAKA_R_AREA_HA  = 4781.6;   // planar area in EPSG:32737, for the sanity check in GEE
var AMAKA_BUFFER_M_R = 7000;     // buffer used in R (Cfg.Buffer.m)
var AMAKA_COORDS = [
    [[[38.764003,-3.881540],[38.814991,-3.882460],[38.817026,-3.950824],[38.759332,-3.956707],[38.764003,-3.881540]]]
  ];
// <<< AMAKA_BOUNDARY_END

// ---- 3. BUILD THE FEATURE -------------------------------------------------------
if (AMAKA_COORDS.length === 0) {
  throw new Error('AMAKA_COORDS is empty. Run R/01_boundary_prep.R first and paste the refreshed script.');
}

// Planar (non-geodesic) edges in lon/lat reproduce the R polygon vertex-for-vertex.
var amakaGeom = ee.Geometry.MultiPolygon(AMAKA_COORDS, 'EPSG:4326', false);

var amaka = ee.FeatureCollection([
  ee.Feature(amakaGeom, {
    ranch_id:   'AMAKA',
    name:       'Amaka Ranch and Conservancy',
    src_rname:  'MUSAMITI',
    aoi_mode:   AMAKA_AOI_MODE,
    aoi_src:    AMAKA_AOI_SRC,
    crs_epsg:   AMAKA_CRS_EPSG,
    crs_status: AMAKA_CRS_STATUS,
    r_area_ha:  AMAKA_R_AREA_HA
  })
]);

// ---- 4. SANITY CHECKS --------------------------------------------------------------
// GEE measures area on the ellipsoid, R in UTM 37S; they should agree to well under 1%.
var gee_area_ha = amakaGeom.area(1).divide(1e4);
print('Area in GEE (ha):', gee_area_ha);
print('Area from R (ha):', AMAKA_R_AREA_HA);
print('Difference (%):', gee_area_ha.subtract(AMAKA_R_AREA_HA).divide(AMAKA_R_AREA_HA).multiply(100));
print('CRS status carried from R:', AMAKA_CRS_STATUS, '(assumed EPSG:' + AMAKA_CRS_EPSG + ')');
print('AOI mode carried from R:', AMAKA_AOI_MODE, '-', AMAKA_AOI_SRC);
if (AMAKA_AOI_MODE === 'MUSAMITI') {
  print('Note: the TTWCA polygon is coarse and smaller than the size TTWCA lists for Amaka (known_issues F8).');
}

// The buffered catchment envelope that scripts 02-10 export over.
var envelope = amakaGeom.buffer(AMAKA_BUFFER_M_R + EXPORT_PAD_M, 10).bounds(10);

Map.centerObject(amakaGeom, 11);
Map.setOptions('HYBRID');
Map.addLayer(ee.FeatureCollection([ee.Feature(envelope)]).style({color: '0000ff', fillColor: '00000000', width: 2}),
             {}, 'Export envelope (buffer + pad)');
Map.addLayer(amaka.style({color: 'ff8800', fillColor: 'ff880033', width: 2}), {}, 'Amaka boundary');

// ---- 5. EXPORT THE ASSET ---------------------------------------------------------------
if (RUN_ASSET_EXPORT) {
  Export.table.toAsset({
    collection:  amaka,
    description: 'Amaka_Ranch_boundary_asset',
    assetId:     AMAKA_ASSET
  });
}
