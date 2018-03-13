/* I used this simple script for reporting out my methods of this projcet. One particular
need was to know how many images were access and/or used. This is a simple GEE method to
do so.

Author: Brian Wong, bw161@duke.edu
2018-02-26

GEE Link: https://code.earthengine.google.com/45225d46e6778ea9311424f103bc62de
*/

// Define AOI:
var nc_bounds = ee.FeatureCollection('ft:17aT9Ud-YnGiXdXEJUyycH2ocUqreOeKGbzCkUw')
  .filter(ee.Filter.eq('name', 'North Carolina'))

// Specify the Input Collections. I used 2 so switched them out to evaluate for both.
// 1.) ESA's Sentinel-1 used for Valid Mask Creation:
//var inputCollection_raw = ee.ImageCollection('COPERNICUS/S2')
// 2.) Landsat-7 SR used for Change Detection analysis:
var inputCollection_raw = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')

// I used two date ranges, the first for Valid Mask Creation & the second for Change Detection.
// var start = '01-01'
// var end = '12-31'
var start = '04-01'
var end = '10-30'

// Filter the annual Input Collection. Separated by year due to non-full year use for Change Detection.
var inputCollection_2010 = inputCollection_raw
  .filterBounds(nc_bounds)
  .filterDate('2010-'+start, '2010-'+end)
var inputCollection_2011 = inputCollection_raw
  .filterBounds(nc_bounds)
  .filterDate('2011-'+start, '2011-'+end)
var inputCollection_2012 = inputCollection_raw
  .filterBounds(nc_bounds)
  .filterDate('2012-'+start, '2012-'+end)
var inputCollection_2013 = inputCollection_raw
  .filterBounds(nc_bounds)
  .filterDate('2013-'+start, '2013-'+end)
var inputCollection_2014 = inputCollection_raw
  .filterBounds(nc_bounds)
  .filterDate('2014-'+start, '2014-'+end)
var inputCollection_2015 = inputCollection_raw
  .filterBounds(nc_bounds)
  .filterDate('2015-'+start, '2015-'+end)
var inputCollection_2016 = inputCollection_raw
  .filterBounds(nc_bounds)
  .filterDate('2016-'+start, '2016-'+end)
var inputCollection_2017 = inputCollection_raw
  .filterBounds(nc_bounds)
  .filterDate('2017-'+start, '2017-'+end)

var ic = ee.ImageCollection(inputCollection_2010.merge(inputCollection_2011)
  .merge(inputCollection_2012).merge(inputCollection_2013).merge(inputCollection_2014)
  .merge(inputCollection_2015).merge(inputCollection_2016).merge(inputCollection_2017))

// Reduce full image collection by the number of images in each.
//print(ic)
var countImg = ic.reduce(ee.Reducer.count())
    
  // View the countImg
Map.addLayer(countImg,{"opacity":1,"bands":["B1_count"],"min":1,"max":127,"palette":["ffc762","ff5454","ff0584"]},'Image Count Scene', true, 0.9);

var ncOutline = ee.Image().toByte().paint(nc_bounds.geometry(),1 ,4);
Map.addLayer(ncOutline, {"opacity":0.9,"bands":["constant"],"max":3,"palette":["2700ff","ffffff"]}, 'NC Study Area')

Map.centerObject(nc_bounds, 7)
Map.setOptions('SATELLITE')

// Check out the stats of the collection:
print('Total Number of Scenes:', ic.size())
print('Total Number of Bytes in Collection:', ic.reduceColumns(ee.Reducer.sum(), ['system:asset_size']))

// Results:
// Landsat-7: 1,497 images = approximately 467GB used
// Sentinel-2: 3,156 images = approximately 2,272GB used