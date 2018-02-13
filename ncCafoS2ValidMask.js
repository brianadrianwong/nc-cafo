/*
This script was a relativley simple way to delineate CAFO rooftops in North Carolina. My primary interest is to
learn of temporal statistics with respect to each CAFO via the Landsat archive, however, there was no good way 
to properly mask out surrounding pixels, hence this pre-processing step.

Here, I've used the higher resolution Sentinel-2 imagery to address this by 1) creating a cloud-free composite,
2) calculate an NDVI composite, 3) create a binary mask through regionally optimized thresholds to isolate low
NDVI pixels, 4) clipped the resulting mask to a 150m search radius per CAFO, and 5) finally export the merged 
image to my GEE Assets.

~Brian Wong, bw161@duke.edu
02-12-2018
*/

//===========================================================================//
//======================= DEFINE VARIABLES & IMPORTS  =======================//
//===========================================================================//

//1. Seven study areas:
var geometry1 = /* color: #d63000 */ee.Geometry.Polygon(
        [[[-81.57348630717024, 36.589045409164065],
          [-81.57348630717024, 36.58684006055065],
          [-81.66961667826399, 35.91127532595572],
          [-81.29058835795149, 35.77546276080151],
          [-80.68084714701399, 35.73534250498251],
          [-80.01068113138899, 36.133415215393455],
          [-80.11230466654524, 36.573606645651786]]]),
    geometry2 = /* color: #98ff00 */ee.Geometry.Polygon(
        [[[-81.67098996927962, 35.907938462340745],
          [-82.59656242951996, 35.92406533270156],
          [-84.03167722513899, 35.15357671429774],
          [-81.75201413920149, 35.18051990168808],
          [-80.59729163282122, 35.09421448715172],
          [-80.29220578959212, 35.380629515528916],
          [-80.68222043802962, 35.73088345023301],
          [-81.28921506693587, 35.77100595359329]]]),
    geometry3 = /* color: #0b4a8b */ee.Geometry.Polygon(
        [[[-80.668487213552, 35.739731655904514],
          [-80.283965729177, 35.381119385351894],
          [-80.28458907210467, 35.38138694860043],
          [-79.85961880534887, 35.31279062775469],
          [-79.44763150066137, 35.24552565471178],
          [-79.04113736003637, 35.66612938874201],
          [-78.65112271159887, 36.53161691511408],
          [-80.09857144206762, 36.57353772829965],
          [-80.0038143619895, 36.134455033540426]]]),
    geometry4 = /* color: #ffc82d */ee.Geometry.Polygon(
        [[[-80.6671142578125, 35.034470859677846],
          [-80.82504272460938, 34.88140130091527],
          [-80.79598980313091, 34.82458455667537],
          [-80.44601440429688, 34.69530850470841],
          [-79.53277587890625, 34.72691746228301],
          [-79.45106506347656, 35.23381857983907],
          [-80.28533935546875, 35.3778302200345]]]),
    geometry5 = /* color: #00ffff */ee.Geometry.Polygon(
        [[[-77.77907747772838, 35.703544744052685],
          [-78.54195649407768, 35.57742729895932],
          [-78.81318675023505, 35.06366643133464],
          [-78.34076474172741, 34.57040294924255],
          [-78.20764597566244, 34.31859481528724],
          [-77.472140904002, 34.60771197483773],
          [-77.14803753928152, 34.86502527920525],
          [-77.07250444724968, 35.663945272798564]]]),
    geometry6 = /* color: #bf04c2 */ee.Geometry.Polygon(
        [[[-78.63327026367188, 36.540444239088316],
          [-79.03564453125, 35.65050805311388],
          [-78.55165477494836, 35.58384777848075],
          [-77.77890999786462, 35.705944843652496],
          [-76.5472412109375, 35.63711590493754],
          [-75.8935546875, 36.544857523189904]]]),
    geometry7 = /* color: #d63000 */ee.Geometry.Polygon(
        [[[-78.55499275960028, 35.564602022223355],
          [-79.05349739827216, 35.63718176156033],
          [-79.43801888264716, 35.24334870647593],
          [-79.53597297764725, 34.700247072855376],
          [-78.98116168200102, 34.215688117781745],
          [-78.21471783124099, 34.32418300823886],
          [-78.34260623886007, 34.568409675096554],
          [-78.81591805256903, 35.06257366596004]]]);

// 2. CAFOS
var cafoPts = ee.FeatureCollection('users/bawong/ejLab/ncPoultryCafoPts')

// 3. North Carolina polygon
var nc_bounds = ee.FeatureCollection('ft:17aT9Ud-YnGiXdXEJUyycH2ocUqreOeKGbzCkUw')
  .filter(ee.Filter.eq('name', 'North Carolina'))

// 4. Temporal Date filters
var startDate = '01-01'
var endDate = '12-31'

//===========================================================================//
//============================ DEFINE FUNCTIONS  ============================//
//===========================================================================//

// 1. NDVI
function ndvi(image){
    return image.addBands(image.normalizedDifference(['B8', 'B4']).rename('ndvi'))
}

// 2. Buffer for CAFOs
function buffer(feature){
  var radius = 150
  var bufferedCafo = feature.buffer(radius)
  return bufferedCafo
}

// 3. Sentinel-2 QA band-based cloud mask
function maskS2clouds(image) {
  var qa = image.select('QA60');
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = ee.Number(2).pow(10).int();
  var cirrusBitMask = ee.Number(2).pow(11).int();
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));
  // Return the masked and scaled data.
  return image.updateMask(mask).divide(10000);
}

//===========================================================================//
//============================== THE ANALYSIS  ==============================//
//===========================================================================//

// 1. Map the function over one year of data and take the median.
var s2_composite = ee.ImageCollection('COPERNICUS/S2')
                  .filterDate('2016-'+startDate, '2017-'+endDate)
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))// Pre-filter to get less cloudy scenes
                  .map(maskS2clouds)
                  .median()
                  .clip(nc_bounds)

// 2. Trick to turn reduced median composite back into an image collection.               
var s2_ic = ee.ImageCollection.fromImages([s2_composite])

// 3. For each geometry, map the NDVI function, clip to AOI and select for optimized threshold
var s2_ndvi_composite_geom1 = s2_ic.map(ndvi).select('ndvi').median().clip(geometry1).lte(0.3)
var s2_ndvi_composite_geom2 = s2_ic.map(ndvi).select('ndvi').median().clip(geometry2).lte(0.27)
var s2_ndvi_composite_geom3 = s2_ic.map(ndvi).select('ndvi').median().clip(geometry3).lte(0.275)
var s2_ndvi_composite_geom4 = s2_ic.map(ndvi).select('ndvi').median().clip(geometry4).lte(0.26)
var s2_ndvi_composite_geom5 = s2_ic.map(ndvi).select('ndvi').median().clip(geometry5).lte(0.225)
var s2_ndvi_composite_geom6 = s2_ic.map(ndvi).select('ndvi').median().clip(geometry6).lte(0.235)
var s2_ndvi_composite_geom7 = s2_ic.map(ndvi).select('ndvi').median().clip(geometry7).lte(0.23)

// 4. Clip each composite to buffered 150m CAFO radius.
var cafoBuffered = cafoPts.map(buffer)
var ndvi_derived_mask_geom1 = s2_ndvi_composite_geom1.clipToCollection(cafoBuffered)
var ndvi_derived_mask_geom2 = s2_ndvi_composite_geom2.clipToCollection(cafoBuffered)
var ndvi_derived_mask_geom3 = s2_ndvi_composite_geom3.clipToCollection(cafoBuffered)
var ndvi_derived_mask_geom4 = s2_ndvi_composite_geom4.clipToCollection(cafoBuffered)
var ndvi_derived_mask_geom5 = s2_ndvi_composite_geom5.clipToCollection(cafoBuffered)
var ndvi_derived_mask_geom6 = s2_ndvi_composite_geom6.clipToCollection(cafoBuffered)
var ndvi_derived_mask_geom7 = s2_ndvi_composite_geom7.clipToCollection(cafoBuffered)

// 5. Join all the separate images into a single image.
var final_valid_mask = ee.ImageCollection.fromImages([ndvi_derived_mask_geom1,
  ndvi_derived_mask_geom2, ndvi_derived_mask_geom3, ndvi_derived_mask_geom4,
  ndvi_derived_mask_geom5, ndvi_derived_mask_geom6, ndvi_derived_mask_geom7]).median()

// 6. Export results to GEE assets.
Export.image.toAsset({
  image: final_valid_mask,
  description: 'cafo_valid_mask_export',
  assetId: 'users/bawong/ejLab/cafoNdviDerivedValidMask_20180212',
  region: nc_bounds,
  scale: 10,
  maxPixels: 1e13
})

// 7. Lastly, let's take a look at our results for visual confirmation.
var imageVisParam = {"opacity":1,"bands":["ndvi"],"palette":["ffffff","ff1901"]};
Map.addLayer(cafoBuffered, {}, 'Buffered CAFOs')
Map.addLayer(final_valid_mask, imageVisParam, 'NDVI-derived Final Valid CAFO Mask', true, 0.5)