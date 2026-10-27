/**
 * 02_dem_amaka_export.js
 * ----------------------------------------------------------------------------
 * Project : Amaka Ranch - Terrain, Hydrology & Vegetation Analysis
 * Author  : Isaac Mureithi (AWF / OPALS, Univ. of Exeter)
 * Written : 2026-10-26      GEE account: ee-imnduta1
 * Save as : users/imnduta1/Ismurn/Amaka/02_dem_amaka_export
 *
 * PURPOSE
 *   Export the DEM over the BUFFERED catchment envelope (not the ranch alone) so
 *   that upslope contributing area is kept for the WhiteboxTools hydrology
 *   (R/03_hydrology.R). Run 01_amaka_boundary_export.js first.
 *
 * INPUTS
 *   AMAKA_ASSET (from script 01)
 *   DEM_SOURCE selects one of:
 *     'GLO30'  COPERNICUS/DEM/GLO30                  (primary; 30 m; a SURFACE model, i.e. includes canopy)
 *     'FABDEM' projects/sat-io/open-datasets/FABDEM  (GLO-30 with forests / buildings removed;
 *                                                     check its licence before using results in a proposal)
 *     'SRTM'   USGS/SRTMGL1_003                      (fallback; 30 m)
 *
 * OUTPUTS  (Drive folder Amaka_GEE_Exports/)
 *   amaka_dem_<glo30|fabdem|srtm>_envelope.tif   single band 'DEM', metres, float32
 *   -> download into data/gee_exports/ ; R/02_dem_prep.R reads it.
 *
 * EXPORT SETTINGS   EPSG:32737, 30 m, bilinear resampling, region = buffered envelope
 * SEASONS / REDUCERS  not applicable (static surface; tiles are mosaicked, not reduced)
 * ----------------------------------------------------------------------------
 */

// ---- 1. CONFIG ---------------------------------------------------------------
var AMAKA_ASSET   = 'users/imnduta1/Ismurn/Amaka/Amaka_Ranch';   // from 01_amaka_boundary_export.js
var BUFFER_M      = 7000;     // keep equal to Cfg.Buffer.m in R/00_setup.R (5000-10000 allowed)
var EXPORT_PAD_M  = 1000;     // margin so R can clip the exact UTM envelope
var DEM_SOURCE    = 'GLO30';  // 'GLO30' | 'FABDEM' | 'SRTM'
var EXPORT_CRS    = 'EPSG:32737';
var EXPORT_SCALE  = 30;
var DRIVE_FOLDER  = 'Amaka_GEE_Exports';
var RUN_EXPORT    = true;     // false = preview only

// ---- 2. AOI ----------------------------------------------------------------------
var amaka    = ee.FeatureCollection(AMAKA_ASSET);
var envelope = amaka.geometry().buffer(BUFFER_M + EXPORT_PAD_M, 10).bounds(10);

// ---- 3. DEM ------------------------------------------------------------------------
var key = DEM_SOURCE.toLowerCase();   // used in the file name; R reads Cfg.DEM.Source = this string
var dem;

if (DEM_SOURCE === 'GLO30') {
  var tiles30 = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM').filterBounds(envelope);
  // A mosaic has no fixed projection; give it the tiles' native one so the
  // bilinear resampling below happens from the true 1 arc-second grid.
  dem = tiles30.mosaic().setDefaultProjection(tiles30.first().projection());
} else if (DEM_SOURCE === 'FABDEM') {
  var tilesFab = ee.ImageCollection('projects/sat-io/open-datasets/FABDEM').select('b1').filterBounds(envelope);
  dem = tilesFab.mosaic().setDefaultProjection(tilesFab.first().projection());
} else if (DEM_SOURCE === 'SRTM') {
  dem = ee.Image('USGS/SRTMGL1_003').select('elevation');
} else {
  throw new Error('DEM_SOURCE must be GLO30, FABDEM or SRTM');
}
dem = dem.rename('DEM').resample('bilinear').toFloat();

// ---- 4. SANITY CHECKS -----------------------------------------------------------------
// Elevation over the ranch and over the whole envelope. Taita-Taveta lowlands should be
// a few hundred metres; anything wild means the wrong asset or band.
var statReducer = ee.Reducer.min().combine(ee.Reducer.max(), '', true)
                    .combine(ee.Reducer.mean(), '', true);
print('DEM source:', DEM_SOURCE);
print('Elevation over Amaka (m):',
      dem.reduceRegion({reducer: statReducer, geometry: amaka.geometry(), scale: 30, maxPixels: 1e9}));
print('Elevation over envelope (m):',
      dem.reduceRegion({reducer: statReducer, geometry: envelope, scale: 30, maxPixels: 1e9}));

Map.centerObject(amaka, 11);
Map.setOptions('HYBRID');
Map.addLayer(ee.Terrain.hillshade(dem), {min: 0, max: 255}, 'Hillshade', true, 0.6);
Map.addLayer(dem, {min: 300, max: 900, palette: ['#2b6a99', '#9bc57a', '#f2e394', '#c0794a', '#ffffff']},
             'DEM (m)', true, 0.5);
Map.addLayer(ee.FeatureCollection([ee.Feature(envelope)]).style({color: '0000ff', fillColor: '00000000', width: 2}),
             {}, 'Export envelope');
Map.addLayer(amaka.style({color: 'ff8800', fillColor: '00000000', width: 2}), {}, 'Amaka boundary');

// ---- 5. EXPORT ------------------------------------------------------------------------
if (RUN_EXPORT) {
  Export.image.toDrive({
    image:          dem,
    description:    'amaka_dem_' + key + '_envelope',
    folder:         DRIVE_FOLDER,
    fileNamePrefix: 'amaka_dem_' + key + '_envelope',
    region:         envelope,
    crs:            EXPORT_CRS,
    scale:          EXPORT_SCALE,
    maxPixels:      1e9,
    fileFormat:     'GeoTIFF'
  });
}
