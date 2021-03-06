///---TEST POINT LONG-LAT---/// 
//random samples generated from GEE:

//var testPoint = ee.Geometry.Point([-81.71978572559999,35.36723263510879]) // 2011
//var testPoint = ee.Geometry.Point([-77.77700561229997,34.87806276711418])
//var testPoint = ee.Geometry.Point([-79.99871537069997,34.844988819914576]) 
//var testPoint = ee.Geometry.Point([-80.59407977900001,36.189869904900284])
//var testPoint = ee.Geometry.Point([-80.50903996450003,35.205065677710564]) // looks like expansion
//var testPoint = ee.Geometry.Point([-77.28304413099998,34.98331078711302]) //
//var testPoint = ee.Geometry.Point([-76.8140278017,36.14907838810069])
//var testPoint = ee.Geometry.Point([-79.15733683740004,36.276812572099395])
//var testPoint = ee.Geometry.Point([-79.58915315899999,35.87388910380345])
//var testPoint = ee.Geometry.Point([-79.11579299810002,34.398446746419765])

//var testPoint = ee.Geometry.Point([-80.33897945610006,34.90696912761388])
//var testPoint = ee.Geometry.Point([-79.74023075419998,35.116498692611486])
//var testPoint = ee.Geometry.Point([-79.51511031030003,35.62567426010602]) //too small - all < 0.7
var testPoint = ee.Geometry.Point([-81.46429249080003,35.47817250360757]) // 2013
//var testPoint = ee.Geometry.Point([-81.01944367800004,35.728322018004995]) // invalid
//var testPoint = ee.Geometry.Point([-78.4363702533,35.00218976811276])
//var testPoint = ee.Geometry.Point([-77.90864644080003,35.32543224210927])
//var testPoint = ee.Geometry.Point([-78.16335814589999,36.038728857101745]) // unclear
//var testPoint = ee.Geometry.Point([-79.7008878771,35.506341322107275]) // unclear unless us 75
//var testPoint = ee.Geometry.Point([-77.54050368220003,35.35845011110884])

//var testPoint = ee.Geometry.Point([-81.05814465689996,35.883279987403334])
//var testPoint = ee.Geometry.Point([-80.24823932310002,35.086285156311874]) // 2011 exansion
//var testPoint = ee.Geometry.Point([-78.35544702440001,35.102260803111704])
//var testPoint = ee.Geometry.Point([-78.8087865579,34.438381587019244])
//var testPoint = ee.Geometry.Point([-79.69099862609995,35.174693247510845])
//var testPoint = ee.Geometry.Point([-81.37385956169997,35.90306899430311]) // to few pixels
//var testPoint = ee.Geometry.Point([-80.2081209281,35.024803868012526])
//var testPoint = ee.Geometry.Point([-76.96584126999997,36.12957949210089])
//var testPoint = ee.Geometry.Point([-77.98550458300004,36.002221632102106])
//var testPoint = ee.Geometry.Point([-77.98177798490003,35.265906776109915])

////////////////////////////  
///---DEFINE VARIABLES---///
////////////////////////////

// 1. CAFOS
var cafoPts = ee.FeatureCollection("users/bawong/ejLab/ncPoultryCafoPts");

// 2. NC
var nc_bounds = ee.FeatureCollection('ft:17aT9Ud-YnGiXdXEJUyycH2ocUqreOeKGbzCkUw')
  .filter(ee.Filter.eq('name', 'North Carolina'))
  
// 3. buffer variable fed thru script for search radius
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
//var validMask = ee.Image("users/bawong/ejLab/cafoValidMask")
var validMask = ee.Image("users/bawong/ejLab/cafoNdviDerivedValidMask")

////////////////////////////
///---DEFINE FUNCTIONS---///
////////////////////////////

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

// 8 update mask
var nonZeroPixels = function(image){
  var i = image.select(tcComponent).gte(0.01)
  var nonZeroPixels = image.mask(i)
  return nonZeroPixels
}

var greenIndex = function(image) {
    return image.addBands(image.select('B4').divide(image.select('B2')).rename('gi'))
};

////////////////////////////////////////////////////////////////////////////////
//-- APPLYING FUNCTIONS TO ANNUAL IMAGE COLLECTIONS FOR ANNUAL TC COMPONENT --//
////////////////////////////////////////////////////////////////////////////////

// var tcImageCollection2007 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
//   .filterBounds(nc_bounds)
//   .filterDate('2007-'+startDate, '2007-'+endDate)
//   //.map(cleaner(outlierValue))
//   .map(maskClouds)
//   .map(tc)
//   .map(ndvi)
//   .map(tcRatio)
//   .map(tcBgNdviRatio)
//   .map(tcBGWindex)
//   .map(greenIndex)
//   .map(nonZeroPixels)
//   .median()
//   .set({"year": 2007})
//   .mask(validMask)

// var tcImageCollection2008 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
//   .filterBounds(nc_bounds)
//   .filterDate('2008-'+startDate, '2008-'+endDate)
//   //.map(cleaner(outlierValue))
//   .map(maskClouds)
//   .map(tc)
//   .map(ndvi)
//   .map(tcRatio)
//   .map(tcBgNdviRatio)
//   .map(tcBGWindex)
//   .map(greenIndex)
//   .map(nonZeroPixels)
//   .median()
//   .set({"year": 2008})
//   .mask(validMask)

// var tcImageCollection2009 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
//   .filterBounds(nc_bounds)
//   .filterDate('2009-'+startDate, '2009-'+endDate)
//   //.map(cleaner(outlierValue))
//   .map(maskClouds)
//   .map(tc)
//   .map(ndvi)
//   .map(tcRatio)
//   .map(tcBgNdviRatio)
//   .map(tcBGWindex)
//   .map(greenIndex)
//   .map(nonZeroPixels)
//   .median()
//   .set({"year": 2009})
//   .mask(validMask)

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
  //.max()
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
  //.max()
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
  //.max()
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
  //.max()
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
  //.max()
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
  //.max()
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
  //.max()
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
  //.max()
  .median()
  .set({"year": 2017})
  .updateMask(validMask)

/////////////////////////////////////////////////////////////////////////////////
///---APPLY SPATIAL REDUCTION TO POINT OF INTEREST FOR TIME SERIES ANALYSIS---///
/////////////////////////////////////////////////////////////////////////////////

var perc = 50 // settled on median essentially
// var bufferedMedianTC_2007 = tcImageCollection2007.select(tcComponent).clip(bufferedPoint).set({"year": 2007})
// var bufferedMedianTC_2008 = tcImageCollection2008.select(tcComponent).clip(bufferedPoint).set({"year": 2008})
// var bufferedMedianTC_2009 = tcImageCollection2009.select(tcComponent).clip(bufferedPoint).set({"year": 2009})
var bufferedMedianTC_2010 = tcImageCollection2010.select(tcComponent).clip(bufferedPoint).set({"year": 2010})
var bufferedMedianTC_2011 = tcImageCollection2011.select(tcComponent).clip(bufferedPoint).set({"year": 2011})
var bufferedMedianTC_2012 = tcImageCollection2012.select(tcComponent).clip(bufferedPoint).set({"year": 2012})
var bufferedMedianTC_2013 = tcImageCollection2013.select(tcComponent).clip(bufferedPoint).set({"year": 2013})
var bufferedMedianTC_2014 = tcImageCollection2014.select(tcComponent).clip(bufferedPoint).set({"year": 2014})
var bufferedMedianTC_2015 = tcImageCollection2015.select(tcComponent).clip(bufferedPoint).set({"year": 2015})
var bufferedMedianTC_2016 = tcImageCollection2016.select(tcComponent).clip(bufferedPoint).set({"year": 2016})
var bufferedMedianTC_2017 = tcImageCollection2017.select(tcComponent).clip(bufferedPoint).set({"year": 2017})

/////////////////////////////////////////////////////////////////////////
//-- REDUCE ALL IMAGES INTO 1 IMAGE COLLECTION FOR PLOTTING PURPOSES --//
/////////////////////////////////////////////////////////////////////////
var collection = ee.ImageCollection.fromImages([bufferedMedianTC_2010, bufferedMedianTC_2011,
bufferedMedianTC_2012, bufferedMedianTC_2013, bufferedMedianTC_2014, 
bufferedMedianTC_2015, bufferedMedianTC_2016, bufferedMedianTC_2017])

///---PLOTTING TIME SERIES--///
var plot = ui.Chart.image.seriesByRegion(collection, bufferedPoint, ee.Reducer.percentile([75]), tcComponent, 30, 'year')
//var plot = ui.Chart.image.seriesByRegion(collection, bufferedPoint, ee.Reducer.median(), tcComponent, 30, 'year')
print(plot)
var plot2 = ui.Chart.image.series(collection, bufferedPoint, ee.Reducer.median(), 30, 'year')
print(plot2)
///---VISUALIZATION---///
//var imageVisParam = {"opacity":1,"bands":["tc_bg"],"min":0,"palette":["ffffff","bababa","03dcff"]};
var imageVisParam = {"opacity":1,"bands":[tcComponent],"min":0,"max":1,"palette":["ffffff","bababa","ff1105"]};

// Map.addLayer(tcImageCollection2007, imageVisParam, '2007', false)
// Map.addLayer(tcImageCollection2008, imageVisParam, '2008', false)
// Map.addLayer(tcImageCollection2009, imageVisParam, '2009', false)
Map.addLayer(tcImageCollection2010, imageVisParam, '2010', true, 0.5)
Map.addLayer(tcImageCollection2011, imageVisParam, '2011', false)
Map.addLayer(tcImageCollection2012, imageVisParam, '2012', false)
Map.addLayer(tcImageCollection2013, imageVisParam, '2013', false)
Map.addLayer(tcImageCollection2014, imageVisParam, '2014', false)
Map.addLayer(tcImageCollection2015, imageVisParam, '2015', false)
Map.addLayer(tcImageCollection2016, imageVisParam, '2016', false)
Map.addLayer(tcImageCollection2017, imageVisParam, '2017', true, 0.5)

Map.addLayer(cafoPts, {}, 'cafo Pts')
Map.centerObject(bufferedPoint, 15)