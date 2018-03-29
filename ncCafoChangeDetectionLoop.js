/* Change Detection Script

Author: Brian Wong, bw161@duke.edu
2018-02-26

GEE Link: https://code.earthengine.google.com/6755a4a3adc4c5ac0f0c5d68a746b4e3
*/

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////---DEFINE VARIABLES---///////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

// 1. CAFO feature collection
var cafoPts = ee.FeatureCollection("users/bawong/ejLab/ncPoultryCafoPts") // point ft

// 2. NC
var nc_bounds = ee.FeatureCollection('ft:17aT9Ud-YnGiXdXEJUyycH2ocUqreOeKGbzCkUw')
  .filter(ee.Filter.eq('name', 'North Carolina'))
  
// 3. buffer variable fed thru script for search radius
var testPoint = ee.Geometry.Point([-80.24823932310002,35.086285156311874]) // 2011 expansion

var buffer_distance = 150
var bufferedPoint = testPoint.buffer(buffer_distance)
Map.addLayer(bufferedPoint, {}, 'test-roi')

// 4. Dates
var startDate = '04-01'
var endDate = '10-30'

// 5. Tasseled Cap component
//var tcComponent = 'brightness'
//var tcComponent = 'greenness'
//var tcComponent = 'wetness'
//var tcComponent = 'ndvi'
//var tcComponent = 'gi'
var tcComponent = 'tc_bg'
//var tcComponent = 'tc_bg_ndvi'
//var tcComponent = 'tc_bgw'

// 6. Define an Array of Tasseled Cap coefficients
var coefficients = ee.Array([
  [0.3037, 0.2793, 0.4743, 0.5585, 0.5082, 0.1863],
  [-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800],
  [0.1509, 0.1973, 0.3279, 0.3406, -0.7112, -0.4572],
  [-0.8242, 0.0849, 0.4392, -0.0580, 0.2012, -0.2768],
  [-0.3280, 0.0549, 0.1075, 0.1855, -0.4357, 0.8085],
  [0.1084, -0.9022, 0.4120, 0.0573, -0.0251, 0.0238]
]);

// 7. visualization variables
var vizParams = {
  bands: ['brightness', 'greenness', 'wetness'],
  min: -0.1, max: [0.5, 0.1, 0.1]};
  
// 8. the valid mask
var validMask = ee.Image('users/bawong/ejLab/cafoNdviDerivedValidMask_20180212')

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////---DEFINE FUNCTIONS---///////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

//1. Cloud Mask for SR products using the pixel_qa band
var maskClouds = function(image){
  var mask = image.select('pixel_qa').bitwiseAnd(2).neq(0);
  var clearImage = image.mask(mask);
  return clearImage
}

//2. Cleaner tool // based on Matt Hancher L5 cleaner function on forums
var outlierValue = 16400
var cleaner = function(outlierValue) {
  return function(image){
      // Compute the minimum mask across all bands, to eliminate known-risky edge pixels.
      var min_mask = image.mask().reduce(ee.Reducer.min());
      // Compute a mask that eliminates fully-saturated (likey non-physical) sensor values.
      var sat_mask = image.reduce(ee.Reducer.max()).lt(outlierValue);
      // Combine these masks and erode neighboring pixels just a bit.
      var new_mask = min_mask.min(sat_mask).focal_min(3);
      return image.updateMask(new_mask);
  };
};

//3. NDVI
var ndvi = function(image) {
    return image.addBands(image.normalizedDifference(['B4', 'B3']).rename('ndvi'))
};

//4. Tasseled Cap
var tc = function(image){
  var scene = image.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B7'])//.copyProperties(['id', 'system:time_start:'])
  var arrayImage1D = scene.toArray(); // Make an Array Image, with a 1-D Array per pixel.
  var arrayImage2D = arrayImage1D.toArray(1); // Make an Array Image with a 2-D Array per pixel, 6x1.
  var componentsImage = ee.Image(coefficients) // Do a matrix multiplication: 6x6 times 6x1.
    .matrixMultiply(arrayImage2D)
    .arrayProject([0])// Get rid of the extra dimensions.
    .arrayFlatten([['brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth']])
    .copyProperties(image/*['id','system:time_start' ,'CLOUD_COVER']*/); // Pulls all 'ordinary' metadata from original image to output
  return image.addBands(componentsImage)
}
 
//5. Tasseled Cap Band Ratio: Brightness to Greenness
var tcRatio = function(image){
  var scene = image.addBands(image.normalizedDifference(['brightness', 'greenness']).rename('tc_bg'))
  return scene
}

//6. Tasseled Cap Band Ratio: Brightness-Greenness to NDVI
var tcBgNdviRatio = function(image){
  var scene = image.addBands(image.normalizedDifference(['tc_bg', 'ndvi']).rename('tc_bg_ndvi'))
  return scene
}

//7. Tasseled Cap Band Ratio: Brightness-Greenness-Wetness
var tcBGWindex = function(image){
  var b = image.select('brightness')
  var g = image.select('greenness')
  var w = image.select('wetness')
  var scene = image.addBands((b.subtract(g).subtract(w)).divide((b.add(g).add(w))).rename('tc_bgw'))
  return scene
}

//8. Update Mask to ignore 0's due to no data spots
var nonZeroPixels = function(image){
  var i = image.select(tcComponent).gte(0.01)
  var nonZeroPixels = image.mask(i)
  return nonZeroPixels
}

var greenIndex = function(image) {
    return image.addBands(image.select('B4').divide(image.select('B2')).rename('gi'))
};

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////---BUILD IMAGE COLLECTIONS---///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

var tcImageCollection2010 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
  .filterBounds(nc_bounds)
  .filterDate('2010-'+startDate, '2010-'+endDate)
  //.map(cleaner(outlierValue))
  .map(maskClouds)
  .map(tc)
  .map(ndvi)
  .map(tcRatio)
  .map(tcBgNdviRatio)
  .map(tcBGWindex)
  .map(greenIndex)
  .map(nonZeroPixels)
  .median()
  .set({"year": 2010})
  .updateMask(validMask)
  
var tcImageCollection2011 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
  .filterBounds(nc_bounds)
  .filterDate('2011-'+startDate, '2011-'+endDate)
  //.map(cleaner(outlierValue))
  .map(maskClouds)
  .map(tc)
  .map(ndvi)
  .map(tcRatio)
  .map(tcBgNdviRatio)
  .map(tcBGWindex)
  .map(greenIndex)
  .map(nonZeroPixels)
  .median()  
  .set({"year": 2011})
  .updateMask(validMask)

var tcImageCollection2012 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
  .filterBounds(nc_bounds)
  .filterDate('2012-'+startDate, '2012-'+endDate)
  //.map(cleaner(outlierValue))
  .map(maskClouds)
  .map(tc)
  .map(ndvi)
  .map(tcRatio)
  .map(tcBgNdviRatio)
  .map(tcBGWindex)
  .map(greenIndex)
  .map(nonZeroPixels)
  .median()
  .set({"year": 2012})
  .updateMask(validMask)

var tcImageCollection2013 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
  .filterBounds(nc_bounds)
  .filterDate('2013-'+startDate, '2013-'+endDate)
  //.map(cleaner(outlierValue))
  .map(maskClouds)
  .map(tc)
  .map(ndvi)
  .map(tcRatio)
  .map(tcBgNdviRatio)
  .map(tcBGWindex)
  .map(greenIndex)
  .map(nonZeroPixels)
  .median()
  .set({"year": 2013})
  .updateMask(validMask)
  
var tcImageCollection2014 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
  .filterBounds(nc_bounds)
  .filterDate('2014-'+startDate, '2014-'+endDate)
  //.map(cleaner(outlierValue))
  .map(maskClouds)
  .map(tc)
  .map(ndvi)
  .map(tcRatio)
  .map(tcBgNdviRatio)
  .map(tcBGWindex)
  .map(greenIndex)
  .map(nonZeroPixels)
  .median()
  .set({"year": 2014})
  .updateMask(validMask)
  
var tcImageCollection2015 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
  .filterBounds(nc_bounds)
  .filterDate('2015-'+startDate, '2015-'+endDate)
  //.map(cleaner(outlierValue))
  .map(maskClouds)
  .map(tc)
  .map(ndvi)
  .map(tcRatio)
  .map(tcBgNdviRatio)
  .map(tcBGWindex)
  .map(greenIndex)
  .map(nonZeroPixels)
  .median()
  .set({"year": 2015})
  .updateMask(validMask)
  
var tcImageCollection2016 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
  .filterBounds(nc_bounds)
  .filterDate('2016-'+startDate, '2016-'+endDate)
  //.map(cleaner(outlierValue))
  .map(maskClouds)
  .map(tc)
  .map(ndvi)
  .map(tcRatio)
  .map(tcBgNdviRatio)
  .map(tcBGWindex)
  .map(greenIndex)
  .map(nonZeroPixels)
  .median()
  .set({"year": 2016})
  .updateMask(validMask)
  
var tcImageCollection2017 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
  .filterBounds(nc_bounds)
  .filterDate('2017-'+startDate, '2017-'+endDate)
  //.map(cleaner(outlierValue))
  .map(maskClouds)
  .map(tc)
  .map(ndvi)
  .map(tcRatio)
  .map(tcBgNdviRatio)
  .map(tcBGWindex)
  .map(greenIndex)
  .map(nonZeroPixels)
  .median()
  .set({"year": 2017})
  .updateMask(validMask)

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////---ANALYSIS FOR SPATIAL STATS---////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

// var annual_img = tcImageCollection2010
// var annual_img = tcImageCollection2011
// var annual_img = tcImageCollection2012
// var annual_img = tcImageCollection2013
// var annual_img = tcImageCollection2014
// var annual_img = tcImageCollection2015
// var annual_img = tcImageCollection2016
var annual_img = tcImageCollection2017

var ic = ee.ImageCollection.fromImages([annual_img])

function buff (f) {
  var radius = 150
  return f.buffer(radius)
}

var cafos_buffered = cafoPts.map(buff)

function addYear(image) {
  var yr = ee.Image.constant(image.get('year'));
  var int = yr.toInt();
  return(image.addBands(int.rename('yr')));
}

var ic_yr_separated = ic.map(addYear);

var get_stats = function(image) {
  var reduced = image.reduceRegions({
    collection: cafos_buffered,
    reducer: ee.Reducer.median(),
    scale:30 
  });
  return(reduced);
};

var out_table = ic_yr_separated.map(get_stats).flatten();

//////////////////////////////////////////////////////////////////////////////////
///////////////////////////////---EXPORT RESULTS---///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

Export.table.toDrive({
  collection: out_table,
  description: "ejLab_cafo_table_exp_yr_",
  fileFormat: "CSV"
});