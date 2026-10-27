/**
 * 04_modis_lst.js
 * ----------------------------------------------------------------------------
 * Project : Amaka Ranch - Terrain, Hydrology & Vegetation Analysis
 * Author  : Isaac Mureithi (AWF / OPALS, Univ. of Exeter)
 * Written : 2026-10-26      GEE account: ee-imnduta1
 * Save as : users/imnduta1/Ismurn/Amaka/04_modis_lst
 *
 * PURPOSE
 *   Day and night land-surface temperature (deg C), seasonal and annual means: thermal
 *   stress on the ranch and in each pan catchment. Run 01_amaka_boundary_export.js first.
 *
 * INPUTS
 *   MODIS/061/MOD11A2 (8-day, 1 km): LST_Day_1km, LST_Night_1km (DN x 0.02 = Kelvin), QC_Day, QC_Night
 * OUTPUTS  (Drive folder Amaka_GEE_Exports/)
 *   modis_lst_day_mean_<ANNUAL|LR|SR|DRY>_<y0>_<y1>.tif    band lst_day_C_<SEASON>_<year>
 *   modis_lst_night_mean_<ANNUAL|LR|SR|DRY>_<y0>_<y1>.tif  band lst_night_C_<SEASON>_<year>
 *
 * SEASON FILTERS (CLAUDE.md section 2; months of the SAME calendar year)
 *   ANNUAL Jan-Dec | LR Mar-May | SR Oct-Dec | DRY Jan-Feb + Jun-Sep
 *   An 8-day composite is assigned to the season of its START date (error <= 7 days).
 *   MOD11A2 starts Feb 2000, so the default start year is 2001 (first complete year).
 *
 * REDUCER
 *   MEAN of the 8-day composites. Temperature is a state variable, not an amount, so it
 *   is averaged, never summed. Day pixels are masked unless QC bits 0-1 = 00 (good quality,
 *   QC_MAX_DAY = 0). The first run applied the same rule to night and the night files came back
 *   EMPTY (every pixel NoData; known_issues G7), so night now also keeps 'other quality'
 *   (QC_MAX_NIGHT = 1). UNTESTED: check that the night files are not all NoData before using them.
 *
 * EXPORT SETTINGS
 *   EPSG:32737, scale 1000 m = MODIS native. 1 km cells are nearest-neighbour from the
 *   sinusoidal grid, so positions are good to about half a pixel.
 * ----------------------------------------------------------------------------
 */

// ---- 1. CONFIG (edit here) -----------------------------------------------------
var AMAKA_ASSET  = 'users/imnduta1/Ismurn/Amaka/Amaka_Ranch';   // from 01_amaka_boundary_export.js
var BUFFER_M     = 7000;     // keep equal to Cfg.Buffer.m in R/00_setup.R
var EXPORT_PAD_M = 1000;     // margin so R can clip the exact UTM envelope
var YEAR_START   = 2001;     // clamped to years the dataset fully covers
var YEAR_END     = 2024;
var EXPORT_CRS   = 'EPSG:32737';
var EXPORT_SCALE = 1000;     // metres; MOD11A2 native 1 km
var DRIVE_FOLDER = 'Amaka_GEE_Exports';
var EXPORT_MODE  = 'stack';  // 'stack' = 4 files (one per season, a band per year); 'per_year' = one file per year x season
var RUN_EXPORTS  = true;     // false = preview only, no tasks created
var QC_MAX_DAY   = 0;   // 0 = good quality only, 1 = also 'other quality'
var QC_MAX_NIGHT = 1;   // night with 0 exported as all NoData (known_issues G7); 1 is untested

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
var lst = ee.ImageCollection('MODIS/061/MOD11A2').filterBounds(envelope).map(function (im) {
  var dayOk   = im.select('QC_Day').bitwiseAnd(3).lte(QC_MAX_DAY);
  var nightOk = im.select('QC_Night').bitwiseAnd(3).lte(QC_MAX_NIGHT);
  var day   = im.select('LST_Day_1km').multiply(0.02).subtract(273.15).updateMask(dayOk).rename('lst_day');
  var night = im.select('LST_Night_1km').multiply(0.02).subtract(273.15).updateMask(nightOk).rename('lst_night');
  return day.addBands(night).copyProperties(im, ['system:time_start']);
});
var years = effectiveYears(lst);

// ---- 6. PREVIEW ----------------------------------------------------------------------------
var lastYear = years[years.length - 1];
Map.centerObject(amaka, 9);
Map.setOptions('HYBRID');
Map.addLayer(seasonYear(lst, lastYear, 'DRY').select('lst_day').mean(),
             {min: 25, max: 45, palette: ['#2c7bb6', '#ffffbf', '#d7191c']}, 'DRY day LST ' + lastYear + ' (C)');
Map.addLayer(amaka.style({color: 'ff8800', fillColor: '00000000', width: 2}), {}, 'Amaka boundary');

// ---- 7. EXPORT ------------------------------------------------------------------------------
exportSeasonal({coll: lst, band: 'lst_day', reducer: function (c) { return c.mean(); },
                prefix: 'modis_lst_day_mean', label: 'lst_day_C', years: years});
exportSeasonal({coll: lst, band: 'lst_night', reducer: function (c) { return c.mean(); },
                prefix: 'modis_lst_night_mean', label: 'lst_night_C', years: years});
print('Export tasks queued - run each from the Tasks tab:', TASKS_QUEUED);

