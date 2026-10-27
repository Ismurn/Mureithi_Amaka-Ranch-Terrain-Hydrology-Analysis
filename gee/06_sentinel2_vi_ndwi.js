/**
 * 06_sentinel2_vi_ndwi.js
 * ----------------------------------------------------------------------------
 * Project : Amaka Ranch - Terrain, Hydrology & Vegetation Analysis
 * Author  : Isaac Mureithi (AWF / OPALS, Univ. of Exeter)
 * Written : 2026-10-26      GEE account: ee-imnduta1
 * Save as : users/imnduta1/Ismurn/Amaka/06_sentinel2_vi_ndwi
 *
 * PURPOSE
 *   Fine-scale greenness (NDVI), surface-water signal (NDWI) and a bare-soil index (BSI) at
 *   10 m: current state of pans and riparian strips, the imagery against which the stream
 *   threshold is calibrated (CLAUDE.md 9.7), and - with ESA WorldCover (08) - the source of
 *   within-ranch land-cover structure (bare fraction, Zone A). The 500 m study pixels are
 *   ~97% shrubland over Amaka, so they carry almost no within-ranch signal (known_issues D5).
 *   Run 01_amaka_boundary_export.js first.
 *
 * INPUTS
 *   COPERNICUS/S2_SR_HARMONIZED: B2, B3, B4, B8, B11, SCL (scene classification), CLOUDY_PIXEL_PERCENTAGE
 *   NDVI = (B8 - B4) / (B8 + B4)      NDWI (McFeeters) = (B3 - B8) / (B3 + B8)
 *   BSI  = ((B11 + B4) - (B8 + B2)) / ((B11 + B4) + (B8 + B2))   higher = more exposed soil.
 *   BSI is a continuous proxy, not a cover fraction: read it relative to the ranch (e.g. the
 *   dry-season top quartile), not against an absolute threshold. B11 is 20 m and is resampled.
 * OUTPUTS  (Drive folder Amaka_GEE_Exports/)
 *   s2_ndvi_mean_<ANNUAL|LR|SR|DRY>_<y0>_<y1>.tif   band ndvi_<SEASON>_<year>
 *   s2_ndwi_mean_<ANNUAL|LR|SR|DRY>_<y0>_<y1>.tif   band ndwi_<SEASON>_<year>
 *   s2_bsi_mean_<ANNUAL|LR|SR|DRY>_<y0>_<y1>.tif    band bsi_<SEASON>_<year>
 *
 * SEASON FILTERS (CLAUDE.md section 2; months of the SAME calendar year)
 *   ANNUAL Jan-Dec | LR Mar-May | SR Oct-Dec | DRY Jan-Feb + Jun-Sep
 *   Surface reflectance is only globally available from late 2018, so the default start
 *   year is 2019. Earlier years (2016-2018) are sparse; lower YEAR_START at your own risk.
 *
 * REDUCER
 *   MEAN of per-scene indices after cloud masking. Scenes with CLOUDY_PIXEL_PERCENTAGE >=
 *   MAX_CLOUDY_PCT are dropped; remaining pixels are kept only if SCL is vegetation (4),
 *   not vegetated (5), water (6) or unclassified (7). Use a median instead by changing the
 *   reducer in section 7 if short-lived cloud edges survive in the mean.
 *
 * EXPORT SETTINGS
 *   EPSG:32737, scale 10 m. Three indices x four seasons = 12 stack files, each about 6 bands
 *   x 25 MB; lower EXPORT_SCALE to 20 or 30 if Drive space is tight.
 * ----------------------------------------------------------------------------
 */

// ---- 1. CONFIG (edit here) -----------------------------------------------------
var AMAKA_ASSET  = 'users/imnduta1/Ismurn/Amaka/Amaka_Ranch';   // from 01_amaka_boundary_export.js
var BUFFER_M     = 7000;     // keep equal to Cfg.Buffer.m in R/00_setup.R
var EXPORT_PAD_M = 1000;     // margin so R can clip the exact UTM envelope
var YEAR_START   = 2019;     // clamped to years the dataset fully covers
var YEAR_END     = 2024;
var EXPORT_CRS   = 'EPSG:32737';
var EXPORT_SCALE = 10;     // metres; Sentinel-2 native 10 m
var DRIVE_FOLDER = 'Amaka_GEE_Exports';
var EXPORT_MODE  = 'stack';  // 'stack' = 4 files (one per season, a band per year); 'per_year' = one file per year x season
var RUN_EXPORTS  = true;     // false = preview only, no tasks created
var MAX_CLOUDY_PCT = 60;   // drop scenes cloudier than this (scene-level metadata)

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
function prepS2(im) {
  var scl = im.select('SCL');
  var ok  = scl.eq(4).or(scl.eq(5)).or(scl.eq(6)).or(scl.eq(7));
  var ndvi = im.normalizedDifference(['B8', 'B4']).rename('ndvi');
  var ndwi = im.normalizedDifference(['B3', 'B8']).rename('ndwi');
  var bsi  = im.expression('((B11 + B4) - (B8 + B2)) / ((B11 + B4) + (B8 + B2))', {
    B2: im.select('B2'), B4: im.select('B4'), B8: im.select('B8'), B11: im.select('B11')
  }).rename('bsi');
  return ndvi.addBands(ndwi).addBands(bsi).updateMask(ok).copyProperties(im, ['system:time_start']);
}
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
           .filterBounds(envelope)
           .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', MAX_CLOUDY_PCT))
           .map(prepS2);
var years = effectiveYears(s2);

// ---- 6. PREVIEW ----------------------------------------------------------------------------
var lastYear = years[years.length - 1];
Map.centerObject(amaka, 12);
Map.setOptions('HYBRID');
Map.addLayer(seasonYear(s2, lastYear, 'DRY').select('ndvi').mean(),
             {min: 0.1, max: 0.7, palette: ['#a6611a', '#dfc27d', '#80cdc1', '#018571']}, 'DRY NDVI ' + lastYear);
Map.addLayer(seasonYear(s2, lastYear, 'LR').select('ndwi').mean(),
             {min: -0.6, max: 0.3, palette: ['#8c510a', '#f6e8c3', '#c7eae5', '#01665e']}, 'LR NDWI ' + lastYear, false);
Map.addLayer(seasonYear(s2, lastYear, 'DRY').select('bsi').mean(),
             {min: -0.1, max: 0.2, palette: ['#1a9850', '#ffffbf', '#d73027']}, 'DRY bare-soil index ' + lastYear, false);
Map.addLayer(amaka.style({color: 'ff8800', fillColor: '00000000', width: 2}), {}, 'Amaka boundary');

// ---- 7. EXPORT ------------------------------------------------------------------------------
exportSeasonal({coll: s2, band: 'ndvi', reducer: function (c) { return c.mean(); },
                prefix: 's2_ndvi_mean', label: 'ndvi', years: years});
exportSeasonal({coll: s2, band: 'ndwi', reducer: function (c) { return c.mean(); },
                prefix: 's2_ndwi_mean', label: 'ndwi', years: years});
exportSeasonal({coll: s2, band: 'bsi', reducer: function (c) { return c.mean(); },
                prefix: 's2_bsi_mean', label: 'bsi', years: years});
print('Export tasks queued - run each from the Tasks tab:', TASKS_QUEUED);

