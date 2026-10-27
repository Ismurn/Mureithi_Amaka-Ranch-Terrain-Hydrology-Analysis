/**
 * 05_modis_vi.js
 * ----------------------------------------------------------------------------
 * Project : Amaka Ranch - Terrain, Hydrology & Vegetation Analysis
 * Author  : Isaac Mureithi (AWF / OPALS, Univ. of Exeter)
 * Written : 2026-10-26      GEE account: ee-imnduta1
 * Save as : users/imnduta1/Ismurn/Amaka/05_modis_vi
 *
 * PURPOSE
 *   NDVI and EVI seasonal and annual means: vegetation phenology and an independent
 *   cross-check on RPI / PML. Run 01_amaka_boundary_export.js first.
 *
 * INPUTS
 *   MODIS/061/MOD13Q1 (16-day, 250 m): NDVI, EVI (DN x 0.0001), SummaryQA
 * OUTPUTS  (Drive folder Amaka_GEE_Exports/)
 *   modis_ndvi_mean_<ANNUAL|LR|SR|DRY>_<y0>_<y1>.tif   band ndvi_<SEASON>_<year>
 *   modis_evi_mean_<ANNUAL|LR|SR|DRY>_<y0>_<y1>.tif    band evi_<SEASON>_<year>
 *
 * SEASON FILTERS (DESCRIPTION.md section 2; months of the SAME calendar year)
 *   ANNUAL Jan-Dec | LR Mar-May | SR Oct-Dec | DRY Jan-Feb + Jun-Sep
 *   A 16-day composite is assigned to the season of its START date (error <= 15 days).
 *   MOD13Q1 starts Feb 2000, so the default start year is 2001 (first complete year).
 *
 * REDUCER
 *   MEAN of the 16-day composites (an index, not an amount). Pixels with SummaryQA > 1
 *   (snow / cloudy) are masked; good (0) and marginal (1) data are kept.
 *
 * EXPORT SETTINGS
 *   EPSG:32737, scale 250 m = MODIS native; nearest-neighbour from the sinusoidal grid.
 * ----------------------------------------------------------------------------
 */

// ---- 1. CONFIG (edit here) -----------------------------------------------------
var AMAKA_ASSET  = 'users/imnduta1/Ismurn/Amaka/Amaka_Ranch';   // from 01_amaka_boundary_export.js
var BUFFER_M     = 7000;     // keep equal to Cfg.Buffer.m in R/00_setup.R
var EXPORT_PAD_M = 1000;     // margin so R can clip the exact UTM envelope
var YEAR_START   = 2001;     // clamped to years the dataset fully covers
var YEAR_END     = 2024;
var EXPORT_CRS   = 'EPSG:32737';
var EXPORT_SCALE = 250;     // metres; MOD13Q1 native 250 m
var DRIVE_FOLDER = 'Amaka_GEE_Exports';
var EXPORT_MODE  = 'stack';  // 'stack' = 4 files (one per season, a band per year); 'per_year' = one file per year x season
var RUN_EXPORTS  = true;     // false = preview only, no tasks created

// ---- 2. SEASONS ---------------------------------------------------------------------
var SEASONS = {
  ANNUAL: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
  LR:     [3, 4, 5],
  SR:     [10, 11, 12],
  DRY:    [1, 2, 6, 7, 8, 9]
};
var SEASON_ORDER = ['ANNUAL', 'LR', 'SR', 'DRY'];

// ---- 3. AOI -----------------------------------------------------------------------------
var amaka    = ee.FeatureCollection(AMAKA_ASSET);
var envelope = amaka.geometry().buffer(BUFFER_M + EXPORT_PAD_M, 10).bounds(10);

// ---- 4. HELPERS (same in every seasonal script) -------------------------------------------
var TASKS_QUEUED = 0;

function seq(a, b) {
  var out = [];
  for (var i = a; i <= b; i++) { out.push(i); }
  return out;
}

function monthFilter(months) {
  var fs = months.map(function (m) { return ee.Filter.calendarRange(m, m, 'month'); });
  return fs.length === 1 ? fs[0] : ee.Filter.or.apply(null, fs);
}

// Images of one year and one season (ANNUAL = the whole calendar year).
function seasonYear(coll, year, season) {
  var c = coll.filter(ee.Filter.calendarRange(year, year, 'year'));
  if (season !== 'ANNUAL') { c = c.filter(monthFilter(SEASONS[season])); }
  return c;
}

// Years in [YEAR_START, YEAR_END] that the collection covers completely, so a half
// year can never be exported as if it were a full season.
function effectiveYears(coll) {
  var tmin = new Date(coll.aggregate_min('system:time_start').getInfo());
  var tmax = new Date(coll.aggregate_max('system:time_start').getInfo());
  var first = tmin.getUTCFullYear() + ((tmin.getUTCMonth() > 0 || tmin.getUTCDate() > 20) ? 1 : 0);
  var last  = tmax.getUTCFullYear() - (tmax.getUTCMonth() < 11 ? 1 : 0);
  var y0 = Math.max(YEAR_START, first);
  var y1 = Math.min(YEAR_END, last);
  if (y1 < y0) { throw new Error('No complete years between ' + YEAR_START + ' and ' + YEAR_END); }
  print('Years exported (requested ' + YEAR_START + '-' + YEAR_END + '):', y0 + '-' + y1);
  return seq(y0, y1);
}

function queueExport(img, name) {
  TASKS_QUEUED += 1;
  if (!RUN_EXPORTS) { return; }
  Export.image.toDrive({
    image:          img,
    description:    name,
    folder:         DRIVE_FOLDER,
    fileNamePrefix: name,
    region:         envelope,
    crs:            EXPORT_CRS,
    scale:          EXPORT_SCALE,
    maxPixels:      1e10,
    fileFormat:     'GeoTIFF'
  });
}

// opts: coll, band, reducer(ImageCollection -> Image), prefix, label, years
function exportSeasonal(opts) {
  var years = opts.years;
  SEASON_ORDER.forEach(function (season) {
    var layers = years.map(function (y) {
      return opts.reducer(seasonYear(opts.coll, y, season).select(opts.band))
               .rename(opts.label + '_' + season + '_' + y)
               .toFloat();
    });
    var tag = opts.prefix + '_' + season;
    if (EXPORT_MODE === 'stack') {
      queueExport(ee.Image.cat(layers), tag + '_' + years[0] + '_' + years[years.length - 1]);
    } else {
      layers.forEach(function (img, i) { queueExport(img, tag + '_' + years[i]); });
    }
  });
}

// ---- 5. DATA -------------------------------------------------------------------------------
var vi = ee.ImageCollection('MODIS/061/MOD13Q1').filterBounds(envelope).map(function (im) {
  var good = im.select('SummaryQA').lte(1);
  return im.select(['NDVI', 'EVI']).multiply(0.0001).updateMask(good)
           .rename(['ndvi', 'evi']).copyProperties(im, ['system:time_start']);
});
var years = effectiveYears(vi);

// ---- 6. PREVIEW ----------------------------------------------------------------------------
var lastYear = years[years.length - 1];
Map.centerObject(amaka, 10);
Map.setOptions('HYBRID');
Map.addLayer(seasonYear(vi, lastYear, 'DRY').select('ndvi').mean(),
             {min: 0.1, max: 0.7, palette: ['#a6611a', '#dfc27d', '#80cdc1', '#018571']}, 'DRY NDVI ' + lastYear);
Map.addLayer(amaka.style({color: 'ff8800', fillColor: '00000000', width: 2}), {}, 'Amaka boundary');

// ---- 7. EXPORT ------------------------------------------------------------------------------
exportSeasonal({coll: vi, band: 'ndvi', reducer: function (c) { return c.mean(); },
                prefix: 'modis_ndvi_mean', label: 'ndvi', years: years});
exportSeasonal({coll: vi, band: 'evi', reducer: function (c) { return c.mean(); },
                prefix: 'modis_evi_mean', label: 'evi', years: years});
print('Export tasks queued - run each from the Tasks tab:', TASKS_QUEUED);

