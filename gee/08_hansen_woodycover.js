/**
 * 08_hansen_woodycover.js
 * ----------------------------------------------------------------------------
 * Project : Amaka Ranch - Terrain, Hydrology & Vegetation Analysis
 * Author  : Isaac Mureithi (AWF / OPALS, Univ. of Exeter)
 * Written : 2026-10-26      GEE account: ee-imnduta1
 * Save as : users/imnduta1/Ismurn/Amaka/08_hansen_woodycover
 *
 * PURPOSE
 *   (1) Woody-cover baseline and tree-cover loss (elephant / fire impact signature and the
 *       "woody sacrifice zone" counter-metric, CLAUDE.md 6.7-6.8).
 *   (2) 10 m ESA WorldCover land cover: the SOURCE of within-ranch land-cover structure (tree /
 *       shrub / grass / bare fractions, shrub-to-tree ratio, bare-and-wet Zone A). Decision of
 *       2026-10-26: the 500 m study pixels are ~97% shrubland over Amaka and stay as the
 *       Kenya-wide baseline only (known_issues D5).
 *   Run 01_amaka_boundary_export.js first.
 *
 * INPUTS
 *   UMD/hansen/global_forest_change_2023_v1_11 : treecover2000, loss, lossyear, gain, datamask
 *   ESA/WorldCover/v200  (2021)  band Map   - used for the analysis
 *   ESA/WorldCover/v100  (2020)  band Map   - OPTIONAL: the version behind the study pixels, for a
 *                                             like-for-like comparison with the Kenya-wide baseline
 *   WorldCover classes: 10 trees, 20 shrubland, 30 grassland, 40 cropland, 50 built-up,
 *                       60 bare / sparse vegetation, 80 permanent water, 90 herbaceous wetland, 95 mangroves
 * OUTPUTS  (Drive folder Amaka_GEE_Exports/)
 *   hansen_gfc_v1_11_envelope.tif      5 bands, int16, 30 m:
 *       treecover2000 (% canopy cover in 2000) | loss (0/1, 2001-2023) |
 *       lossyear (0 = none, 1-23 = 2001-2023)  | gain (0/1, 2000-2012) | datamask (0 no data, 1 land, 2 water)
 *   worldcover_2021_envelope.tif       band 'Map', class codes above, 10 m
 *   worldcover_2020_v100_envelope.tif  band 'Map', 10 m   (only if EXPORT_WORLDCOVER_V100)
 *
 * SEASON FILTERS / REDUCERS
 *   Not applicable: all are static layers (loss is annual in its own 'lossyear' band, not seasonal).
 *
 * EXPORT SETTINGS   EPSG:32737; 30 m for Hansen (native), 10 m for WorldCover; nearest neighbour
 *                   (all bands are categorical or integer percentages). Fractions at coarser
 *                   grids are computed in R from the 10 m classes, not here.
 *
 * CAUTION  (a) Hansen "tree cover" is canopy > 5 m. In shrub-dominated savanna it under-detects loss of
 *          small woody plants, so treat loss here as the loss of TALL woody cover only
 *          (consistent with the mature-tree framing of the woody sacrifice zone).
 *          (b) WorldCover is a classification, not a measurement: it may also find Amaka almost all
 *          shrubland. Section 4 prints the class counts so you can see this on the first run; if
 *          class 60 (bare) is nearly absent, use the continuous Sentinel-2 bare-soil index from 06.
 *          (c) v100 and v200 are different classifications, not a change-detection pair.
 * ----------------------------------------------------------------------------
 */

// ---- 1. CONFIG ---------------------------------------------------------------
var AMAKA_ASSET            = 'users/imnduta1/Ismurn/Amaka/Amaka_Ranch';   // from 01_amaka_boundary_export.js
var BUFFER_M               = 7000;     // keep equal to Cfg.Buffer.m in R/00_setup.R
var EXPORT_PAD_M           = 1000;     // margin so R can clip the exact UTM envelope
var EXPORT_CRS             = 'EPSG:32737';
var DRIVE_FOLDER           = 'Amaka_GEE_Exports';
var EXPORT_WORLDCOVER_V100 = false;    // also export the 2020 / v100 map (study-pixel version)
var RUN_EXPORTS            = true;     // false = preview only

// ---- 2. AOI ---------------------------------------------------------------------
var amaka    = ee.FeatureCollection(AMAKA_ASSET);
var envelope = amaka.geometry().buffer(BUFFER_M + EXPORT_PAD_M, 10).bounds(10);

// ---- 3. DATA ----------------------------------------------------------------------
var gfc    = ee.Image('UMD/hansen/global_forest_change_2023_v1_11');
var hansen = gfc.select(['treecover2000', 'loss', 'lossyear', 'gain', 'datamask']).toInt16();

var wc21 = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map').toInt16();
var wc20 = ee.ImageCollection('ESA/WorldCover/v100').first().select('Map').toInt16();

// ---- 4. SANITY CHECKS -------------------------------------------------------------------
// Pixel counts per WorldCover class at 10 m (one pixel = 0.01 ha). Compare the shares with the
// 500 m study pixels over Amaka: shrubland 96.6%, grassland 2.9%, trees 0.3%, bare 0.2%.
function classCounts(img, geom, label) {
  var h = img.reduceRegion({
    reducer: ee.Reducer.frequencyHistogram(), geometry: geom, scale: 10, maxPixels: 1e9});
  print(label, h.get('Map'));
}
classCounts(wc21, amaka.geometry(), 'WorldCover v200 (2021) pixel counts over Amaka:');
classCounts(wc21, envelope,         'WorldCover v200 (2021) pixel counts over the envelope:');
classCounts(wc20, amaka.geometry(), 'WorldCover v100 (2020) pixel counts over Amaka:');

print('Mean canopy cover in 2000 over Amaka (%):',
      gfc.select('treecover2000').reduceRegion({
        reducer: ee.Reducer.mean(), geometry: amaka.geometry(), scale: 30, maxPixels: 1e9}));
print('Share of Amaka pixels with recorded loss 2001-2023 (0-1):',
      gfc.select('loss').reduceRegion({
        reducer: ee.Reducer.mean(), geometry: amaka.geometry(), scale: 30, maxPixels: 1e9}));

Map.centerObject(amaka, 12);
Map.setOptions('HYBRID');
// The catalogue image carries its own class palette, which is lost once the band is cast, so the
// map layer is taken from the original image.
Map.addLayer(ee.ImageCollection('ESA/WorldCover/v200').first(), {bands: ['Map']}, 'WorldCover 2021', true, 0.7);
Map.addLayer(gfc.select('treecover2000'), {min: 0, max: 60, palette: ['#ffffe5', '#78c679', '#004529']},
             'Canopy cover 2000 (%)', false, 0.7);
Map.addLayer(gfc.select('lossyear').selfMask(), {min: 1, max: 23, palette: ['#fee08b', '#d73027', '#67001f']},
             'Loss year (1 = 2001 ... 23 = 2023)', false);
Map.addLayer(amaka.style({color: 'ff8800', fillColor: '00000000', width: 2}), {}, 'Amaka boundary');

// ---- 5. EXPORT ------------------------------------------------------------------------------
function exportStatic(img, name, scale) {
  if (!RUN_EXPORTS) { return; }
  Export.image.toDrive({
    image:          img,
    description:    name,
    folder:         DRIVE_FOLDER,
    fileNamePrefix: name,
    region:         envelope,
    crs:            EXPORT_CRS,
    scale:          scale,
    maxPixels:      1e10,
    fileFormat:     'GeoTIFF'
  });
}

exportStatic(hansen, 'hansen_gfc_v1_11_envelope', 30);
exportStatic(wc21,   'worldcover_2021_envelope', 10);
if (EXPORT_WORLDCOVER_V100) {
  exportStatic(wc20, 'worldcover_2020_v100_envelope', 10);
}
