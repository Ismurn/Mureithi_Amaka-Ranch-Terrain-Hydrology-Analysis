/**
 * 09_pmlv2_amaka_export.js
 * ----------------------------------------------------------------------------
 * Project : Amaka Ranch - Terrain, Hydrology & Vegetation Analysis
 * Author  : Isaac Mureithi (AWF / OPALS, Univ. of Exeter)
 * Written : 2026-10-26      GEE account: ee-imnduta1
 * Save as : users/imnduta1/Ismurn/Amaka/09_pmlv2_amaka_export
 *
 * PURPOSE
 *   Seasonal GPP, Ec, Es and Ei. REQUIRED, not optional: the PML_V2 stack held locally
 *   (PML_V2_Kenya_Annual_Stack.tif) is annual only (80 layers = 4 variables x 2001-2020), so the
 *   dry-season Es/ET term of the stress index cannot be derived from it (CLAUDE.md 6.4, 6.8).
 *   Run 01_amaka_boundary_export.js first.
 *
 * INPUTS
 *   CAS/IGSNRR/PML/V2_v017, 8-day, 500 m: GPP (gC m-2 d-1), Ec, Es, Ei (mm d-1).
 *   This is the asset your studyPixels / PML_V2_Kenya scripts use, so the seasonal exports
 *   are on the same product as the local annual stack.
 *   PML_SCALE = 1: the stored values are already physical daily rates. Checked on the local
 *   stack (your Sep-Aug sums over Amaka give GPP ~106 and ET ~80-90 per year, i.e. ~2.4 gC m-2 d-1
 *   and ~1.8 mm d-1 over ~45 composites); section 5 prints the same check from the asset itself.
 *   The asset's end date is not known from here: the local 2020 layer looks like a partial year.
 *   Section 5 prints the last image date; the year clamp below uses only fully covered years.
 * OUTPUTS  (Drive folder Amaka_GEE_Exports/)
 *   pml_gpp_meanrate_<ANNUAL|LR|SR|DRY>_<y0>_<y1>.tif   band gpp_<SEASON>_<year>   (gC m-2 d-1)
 *   pml_ec_meanrate_...   pml_es_meanrate_...   pml_ei_meanrate_...   (mm d-1)
 *
 * SEASON FILTERS (CLAUDE.md section 2; months of the SAME calendar year)
 *   ANNUAL Jan-Dec | LR Mar-May | SR Oct-Dec | DRY Jan-Feb + Jun-Sep
 *   An 8-day composite is assigned to the season of its START date (error <= 7 days).
 *
 * REDUCER
 *   MEAN of the daily RATES carried by the 8-day composites (CLAUDE.md 9.4: the reducer must
 *   match the units). Season totals are rate x days in the season, done in R/05 (LR 92, SR 92,
 *   DRY 181-182, ANNUAL 365-366 days). Ratios such as Es/ET = Es / (Ec + Es + Ei) and WUE =
 *   GPP / ET do not depend on that factor, which is what the stress index uses.
 *   NB your local PML_V2_Kenya_Annual_Stack.tif is a SUM of these daily rates over ~45 composites
 *   per Sep-Aug year: multiply it by 8 to get mm (or gC m-2) per year. It is not directly
 *   comparable to the calendar-year ANNUAL layer exported here.
 *
 * EXPORT SETTINGS
 *   EPSG:32737, scale 500 m = PML_V2 native.
 * ----------------------------------------------------------------------------
 */

// ---- 1. CONFIG (edit here) -----------------------------------------------------
var AMAKA_ASSET  = 'users/imnduta1/Ismurn/Amaka/Amaka_Ranch';   // from 01_amaka_boundary_export.js
var BUFFER_M     = 7000;     // keep equal to Cfg.Buffer.m in R/00_setup.R
var EXPORT_PAD_M = 1000;     // margin so R can clip the exact UTM envelope
var YEAR_START   = 2001;     // clamped to years the dataset fully covers
var YEAR_END     = 2024;
var EXPORT_CRS   = 'EPSG:32737';
var EXPORT_SCALE = 500;     // metres; PML_V2 native 500 m
var DRIVE_FOLDER = 'Amaka_GEE_Exports';
var EXPORT_MODE  = 'stack';  // 'stack' = 4 files (one per season, a band per year); 'per_year' = one file per year x season
var RUN_EXPORTS  = true;     // false = preview only, no tasks created
var PML_ASSET = 'CAS/IGSNRR/PML/V2_v017';   // same asset as your studyPixels / PML_V2_Kenya scripts
var PML_SCALE = 1;           // stored values are already physical rates (see header); do not use 0.01

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
var pmlRaw = ee.ImageCollection(PML_ASSET);

// Verification block: read these in the Console BEFORE exporting.
print('PML bands:', pmlRaw.first().bandNames());
print('PML first image date:', ee.Date(pmlRaw.aggregate_min('system:time_start')));
print('PML last image date:',  ee.Date(pmlRaw.aggregate_max('system:time_start')));
print('PML first image (open it: look for scale / units in the band properties):', pmlRaw.first());

var pml = pmlRaw.filterBounds(envelope).map(function (im) {
  return im.select(['GPP', 'Ec', 'Es', 'Ei']).multiply(PML_SCALE)
           .rename(['gpp', 'ec', 'es', 'ei']).copyProperties(im, ['system:time_start']);
});
var years = effectiveYears(pml);

// Plausibility check: mean daily GPP over Amaka in the long rains of the last year.
// Expect roughly 0.5-5 gC m-2 d-1 (about 2.4 as an annual mean). Hundreds would mean the scale
// factor is wrong; hundredths would mean it was applied twice.
var lastYear = years[years.length - 1];
print('Mean daily GPP over Amaka, LR ' + lastYear + ' (gC m-2 d-1):',
      seasonYear(pml, lastYear, 'LR').select('gpp').mean()
        .reduceRegion({reducer: ee.Reducer.mean(), geometry: amaka.geometry(), scale: 500, maxPixels: 1e9}));

// ---- 6. PREVIEW ----------------------------------------------------------------------------
Map.centerObject(amaka, 10);
Map.setOptions('HYBRID');
Map.addLayer(seasonYear(pml, lastYear, 'DRY').select('es').mean(),
             {min: 0, max: 1.5, palette: ['#ffffcc', '#fd8d3c', '#800026']}, 'DRY Es ' + lastYear + ' (mm/d)');
Map.addLayer(amaka.style({color: 'ff8800', fillColor: '00000000', width: 2}), {}, 'Amaka boundary');

// ---- 7. EXPORT ------------------------------------------------------------------------------
exportSeasonal({coll: pml, band: 'gpp', reducer: function (c) { return c.mean(); },
                prefix: 'pml_gpp_meanrate', label: 'gpp', years: years});
exportSeasonal({coll: pml, band: 'ec', reducer: function (c) { return c.mean(); },
                prefix: 'pml_ec_meanrate', label: 'ec', years: years});
exportSeasonal({coll: pml, band: 'es', reducer: function (c) { return c.mean(); },
                prefix: 'pml_es_meanrate', label: 'es', years: years});
exportSeasonal({coll: pml, band: 'ei', reducer: function (c) { return c.mean(); },
                prefix: 'pml_ei_meanrate', label: 'ei', years: years});
print('Export tasks queued - run each from the Tasks tab:', TASKS_QUEUED);

