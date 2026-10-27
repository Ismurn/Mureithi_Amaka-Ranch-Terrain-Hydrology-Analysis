/**
 * 10_fire_history.js
 * ----------------------------------------------------------------------------
 * Project : Amaka Ranch - Terrain, Hydrology & Vegetation Analysis
 * Author  : Isaac Mureithi (AWF / OPALS, Univ. of Exeter)
 * Written : 2026-10-26      GEE account: ee-imnduta1
 * Save as : users/imnduta1/Ismurn/Amaka/10_fire_history
 *
 * PURPOSE
 *   Fire history, so fire-driven woody loss can be separated from likely elephant-driven
 *   loss (Jones 2025; DESCRIPTION.md 6.7): burn count per pixel, most recent burn year, burned-by-year
 *   stack, and (optionally) MODIS active-fire detection counts.
 *   Run 01_amaka_boundary_export.js first.
 *
 * INPUTS
 *   MODIS/061/MCD64A1  monthly burned area, 500 m, band BurnDate (day of year; > 0 = burned,
 *                      0 = not burned, negative = unmapped / unburnable)
 *   FIRMS              daily active-fire raster, 1 km: confidence (0-100)  [optional]
 * OUTPUTS  (Drive folder Amaka_GEE_Exports/)
 *   burn_count_<y0>_<y1>.tif            number of YEARS with at least one burn (0 = never)
 *   most_recent_burn_year_<y0>_<y1>.tif last year burned, 0 = never burned
 *   burned_by_year_<y0>_<y1>.tif        one 0/1 band per year, band burned_<year>
 *   firms_detections_by_year_<y0>_<y1>.tif  active-fire detections per year, band firms_<year>  [optional]
 *
 * SEASON FILTERS / REDUCERS
 *   Not seasonal: fire is summarised per calendar year. Within a year the monthly burn masks are
 *   combined with MAX (a pixel burned in two months of one year still counts once), then years are
 *   SUMMED for the count and MAX-ed against the year number for the most recent burn.
 *   FIRMS: detections with confidence >= FIRMS_MIN_CONF are SUMMED (a count of detections).
 *
 * EXPORT SETTINGS   EPSG:32737, 500 m for MCD64A1 (native 463 m), 1000 m for FIRMS (native 1 km).
 *
 * NOTE  FIRMS in Earth Engine is a 1 km raster, not a table of fire points. If you need the
 *       original points, download them from NASA FIRMS and add them as a separate layer in R.
 * ----------------------------------------------------------------------------
 */

// ---- 1. CONFIG ---------------------------------------------------------------
var AMAKA_ASSET      = 'users/imnduta1/Ismurn/Amaka/Amaka_Ranch';   // from 01_amaka_boundary_export.js
var BUFFER_M         = 7000;     // keep equal to Cfg.Buffer.m in R/00_setup.R
var EXPORT_PAD_M     = 1000;     // margin so R can clip the exact UTM envelope
var YEAR_START       = 2001;     // MCD64A1 starts Nov 2000; 2001 is the first complete year
var YEAR_END         = 2024;
var EXPORT_CRS       = 'EPSG:32737';
var BURN_SCALE       = 500;
var FIRMS_SCALE      = 1000;
var DRIVE_FOLDER     = 'Amaka_GEE_Exports';
var EXPORT_FIRMS     = true;     // optional active-fire detection counts
var FIRMS_MIN_CONF   = 50;       // detection confidence threshold (0-100)
var RUN_EXPORTS      = true;     // false = preview only

// ---- 2. AOI ---------------------------------------------------------------------
var amaka    = ee.FeatureCollection(AMAKA_ASSET);
var envelope = amaka.geometry().buffer(BUFFER_M + EXPORT_PAD_M, 10).bounds(10);

function seq(a, b) {
  var out = [];
  for (var i = a; i <= b; i++) { out.push(i); }
  return out;
}
var years = seq(YEAR_START, YEAR_END);

// ---- 3. BURNED AREA (MCD64A1) ---------------------------------------------------------
var ba = ee.ImageCollection('MODIS/061/MCD64A1').select('BurnDate');

// Generic band name collection for math operations (sum / max)
var annualGeneric = years.map(function (y) {
  return ba.filterDate(y + '-01-01', (y + 1) + '-01-01')
           .map(function (im) { return im.gt(0); })
           .max()
           .unmask(0)
           .rename('burned')
           .toByte();
});

// Unique band names for multi-band stacking
var annualNamed = years.map(function (y) {
  return ba.filterDate(y + '-01-01', (y + 1) + '-01-01')
           .map(function (im) { return im.gt(0); })
           .max()
           .unmask(0)
           .rename('burned_' + y)
           .toByte();
});

var burnedByYear = ee.Image.cat(annualNamed);

var burnCount = ee.ImageCollection.fromImages(annualGeneric).sum().rename('burn_count').toByte();

var lastBurn = ee.ImageCollection.fromImages(
  annualGeneric.map(function (img, i) { return img.multiply(years[i]).rename('year').toInt16(); })
).max().rename('most_recent_burn_year');

// ---- 4. SANITY CHECKS ----------------------------------------------------------------------
print('Years summarised:', YEAR_START + '-' + YEAR_END);
print('Burn count over Amaka (min / max / mean years burned):',
      burnCount.reduceRegion({
        reducer: ee.Reducer.min().combine(ee.Reducer.max(), '', true).combine(ee.Reducer.mean(), '', true),
        geometry: amaka.geometry(), scale: BURN_SCALE, maxPixels: 1e9}));

Map.centerObject(amaka, 10);
Map.setOptions('HYBRID');
Map.addLayer(burnCount.selfMask(), {min: 1, max: 10, palette: ['#ffffb2', '#fd8d3c', '#bd0026', '#400010']},
             'Burn count (years burned)');
Map.addLayer(lastBurn.selfMask(), {min: YEAR_START, max: YEAR_END, palette: ['#4575b4', '#fee090', '#d73027']},
             'Most recent burn year', false);
Map.addLayer(amaka.style({color: 'ff8800', fillColor: '00000000', width: 2}), {}, 'Amaka boundary');

// ---- 5. EXPORT -------------------------------------------------------------------------------
var tag = YEAR_START + '_' + YEAR_END;
var queued = 0;

function exportBurn(img, name) {
  queued += 1;
  if (!RUN_EXPORTS) { return; }
  Export.image.toDrive({
    image: img, description: name, folder: DRIVE_FOLDER, fileNamePrefix: name,
    region: envelope, crs: EXPORT_CRS, scale: BURN_SCALE, maxPixels: 1e10, fileFormat: 'GeoTIFF'
  });
}
exportBurn(burnCount, 'burn_count_' + tag);
exportBurn(lastBurn, 'most_recent_burn_year_' + tag);
exportBurn(burnedByYear, 'burned_by_year_' + tag);

// ---- 6. ACTIVE FIRE DETECTIONS (FIRMS, optional) ---------------------------------------------
if (EXPORT_FIRMS) {
  var firms = ee.ImageCollection('FIRMS').select('confidence');
  var firmsYears = years.map(function (y) {
    return firms.filterDate(y + '-01-01', (y + 1) + '-01-01')
                .map(function (im) { return im.gte(FIRMS_MIN_CONF); })
                .sum()
                .unmask(0)
                .rename('firms_' + y)
                .toInt16();
  });
  queued += 1;
  if (RUN_EXPORTS) {
    Export.image.toDrive({
      image: ee.Image.cat(firmsYears), description: 'firms_detections_by_year_' + tag, folder: DRIVE_FOLDER,
      fileNamePrefix: 'firms_detections_by_year_' + tag, region: envelope, crs: EXPORT_CRS,
      scale: FIRMS_SCALE, maxPixels: 1e10, fileFormat: 'GeoTIFF'
    });
  }
}
print('Export tasks queued - run each from the Tasks tab:', queued);
