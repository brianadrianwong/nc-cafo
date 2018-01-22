// 1. CAFOS
var cafoPts = ee.FeatureCollection("users/bawong/ejLab/ncPoultryCafoPts");

// 2. NC
var nc_bounds = ee.FeatureCollection('ft:17aT9Ud-YnGiXdXEJUyycH2ocUqreOeKGbzCkUw')
  .filter(ee.Filter.eq('name', 'North Carolina'))

// 3. Dates
var startDate = '01-01'
var endDate = '12-31'

// 4. NDVI
var ndvi = function(image) {
    return image.addBands(image.normalizedDifference(['B8', 'B4']).rename('ndvi'))
};

// 5. NDWI
var ndwi = function(image) {
    return image.addBands(image.normalizedDifference(['B3', 'B8']).rename('ndwi'))
};

// Load Sentinel-2 TOA reflectance data.
var s2 = ee.ImageCollection('COPERNICUS/S2')
  //.median()
//print(s2)

// Function to mask clouds using the Sentinel-2 QA band.
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

// Map the function over one year of data and take the median.
var composite = s2.filterDate('2016-'+startDate, '2017-'+endDate)
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                  .map(maskS2clouds)
                  .median()
                  .clip(nc_bounds)
            
                  
var composite2 = ee.ImageCollection.fromImages([composite])
var ndviComp = composite2.map(ndvi).select('ndvi').median().lte(0.2)
var ndwiComp = composite2.map(ndwi).select('ndwi')
print(ndviComp)
var imageVisParam = {"opacity":1,"bands":["ndvi"],"palette":["ffffff","ff1901"]};
Map.addLayer(ndviComp, imageVisParam, 'ndvi', true, 0.5)
Map.addLayer(ndwiComp, {}, 'ndwi', false)

var buffer = function(feature){
  var bufferDist = 150
  var bufferedPt = feature.buffer(bufferDist)
  return bufferedPt
}

var cafoBuffered = cafoPts.map(buffer)
Map.addLayer(cafoBuffered, {}, 'cafoBuffered')

var ndviDerivedFinalMask = ndviComp.clipToCollection(cafoBuffered)
Map.addLayer(ndviDerivedFinalMask, imageVisParam, 'final mask', true, 0.5)

Export.image.toAsset({
  image: ndviDerivedFinalMask,
  description: 'exportTask',
  assetId: 'users/bawong/ejLab/cafoNdviDerivedValidMask',
  region: nc_bounds,
  scale: 10,
  maxPixels: 1e13
})