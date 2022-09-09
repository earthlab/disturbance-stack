// Climate Disturbance Stack Analysis: Step 1, data wrangling
// Tyler McIntosh, CU Boulder Earth Lab, 9/9/2022

// This script takes in the TerraClimate dataset, as well as a domain of interest to which results will be clipped to.
// Datasets are wrangled in such a way as to enable future analysis as outlined by
// Hammond et al. 2022, which requires monthly averages for each variable of interest over the entire timeframe
// as well as the minimum/maximum of each variable in a given year, at a given location, in addition to the month in which that max/min occurred.
// The script outputs two image collections: one for monthly averages (one image per variable), and the other for min/max/months (one image per year).
// The image collections are then output as a series named geotiffs to a folder in the user's google drive.

//////////// DATA

//Data imports copied from GEE
var terraclimate = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE"),
    neondomains = ee.FeatureCollection("users/tymc5571/NEON_Domains");



//////////// SETUP

//Google Drive folder for outputs
var foldername = 'GEE_Exports';

//Domain of interest. Here, relevant NEON domains
var domain = neondomains.filter("DomainName == 'Northern Rockies' || DomainName == 'Great Basin' || \
DomainName == 'Pacific Northwest' || DomainName == 'Pacific Southwest' || DomainName == 'Desert Southwest' || \
DomainName == 'Southern Rockies / Colorado Plateau'");

//View terraclimate data structure
print('All TerraClimate', terraclimate);

//Select only variables of interest. Variables shown below, along with negative or positive direction indicating warm/dry
// TMAX = tmmx, VPD = vpd, CWD = def, SOIL M = soil, PPT = pr, PDSI = pdsi
// Positive variable: tmmx, vpd, def
// Negative variables: soil, pr, pdsi
var vars = terraclimate.select('tmmx', 'vpd', 'def', 'soil', 'pr', 'pdsi');


////////////// AVERAGES OVER TIME
//Create list of images, one for each month of the year, with 6 bands each containing the average of a given variable over the entire timeframe 

//Create list of months
var months = ee.List.sequence(1, 12);

//Function to pull monthly data over entire time frame, average, and reset image properties/band names
var getmonthmean = function(month) {
  var means = vars.filter(ee.Filter.calendarRange(month,month,'month')).mean();
  means = means.rename(['tmmx_mean', 'vpd_mean', 'def_mean', 'soil_mean', 'pr_mean', 'pdsi_mean']);
  means = means.set('month', ee.Number(month).toInt());
  return means.clip(domain);
};

//Map function over list of months, print output, test by visualizing a sample variable
var monthmeans = ee.ImageCollection.fromImages(months.map(getmonthmean));
print('Monthly means image collection:', monthmeans);
Map.addLayer(monthmeans.filter(ee.Filter.eq('month',7)).select('tmmx_mean'), {min: -100, max: 500, palette:['blue', 'red']}, 'temp', false);


//////////////// MIN/MAX/MONTHS
//Create imagecollections for each variable of interest, with an image for each year containing two bands:
// 1) the minimum or maximum of the variable that occurred at that location during the given year, and 2) the month of that value's occurrence

//List of years to map over
var years = ee.List.sequence(1958, 2021);

//Function to create a min reducer that retains band names & values for all other bands
function argminReduce(imageCollection) {
  var bandNames = imageCollection.first().bandNames();
  return imageCollection.reduce(
    ee.Reducer.min(bandNames.size())
      .setOutputs(bandNames)
  );
}

//Function to create a max reducer that retains band names & values for all other bands
function argmaxReduce(imageCollection) {
  var bandNames = imageCollection.first().bandNames();
  return imageCollection.reduce(
    ee.Reducer.max(bandNames.size())
      .setOutputs(bandNames)
  );
}

//Function to add the image month as a band called 'month'
function addmonthband(image) {
  var monthVal = image.date().get('month');
  return image.addBands(ee.Image.constant(monthVal).toInt().rename('month'));
}

//Functions to get negative variable minimums and month of those minimums, rename appropriately
var getprecipmin = function(year) {
  var workingyear = vars.select('pr').filter(ee.Filter.calendarRange(year, year, 'year'));
  var minimums = argminReduce(workingyear.map(addmonthband));
  minimums = minimums.rename(['pr_min', 'pr_month']);
  minimums = minimums.set('year', ee.Number(year).toInt()); //set 'year' property here, which is carried over to all others when bands are added in combine function
  return minimums;
};
var getpdsimin = function(year) {
  var workingyear = vars.select('pdsi').filter(ee.Filter.calendarRange(year, year, 'year'));
  var minimums = argminReduce(workingyear.map(addmonthband));
  minimums = minimums.rename(['pdsi_min', 'pdsi_month']);
  return minimums;
};
var getsoilmin = function(year) {
  var workingyear = vars.select('soil').filter(ee.Filter.calendarRange(year, year, 'year'));
  var minimums = argminReduce(workingyear.map(addmonthband));
  minimums = minimums.rename(['soil_min', 'soil_month']);
  return minimums;
};

// Functions to get positive variable maximums and month of those maximums, rename appropriately
var gettempmax = function(year) {
  var workingyear = vars.select('tmmx').filter(ee.Filter.calendarRange(year, year, 'year'));
  var maximums = argmaxReduce(workingyear.map(addmonthband));
  maximums = maximums.rename(['tmmx_max', 'tmmx_month']);
  return maximums;
};
var getvpdmax = function(year) {
  var workingyear = vars.select('vpd').filter(ee.Filter.calendarRange(year, year, 'year'));
  var maximums = argmaxReduce(workingyear.map(addmonthband));
  maximums = maximums.rename(['vpd_max', 'vpd_month']);
  return maximums;
};
var getdeficitmax = function(year) {
  var workingyear = vars.select('def').filter(ee.Filter.calendarRange(year, year, 'year'));
  var maximums = argmaxReduce(workingyear.map(addmonthband));
  maximums = maximums.rename(['def_max', 'def_month']);
  return maximums;
};

//Run functions
var precipmins = years.map(getprecipmin);
var pdsimins = years.map(getpdsimin);
var soilmins = years.map(getsoilmin);
var tempmaxs = years.map(gettempmax);
var vpdmaxs = years.map(getvpdmax);
var deficitmaxs = years.map(getdeficitmax);

////////////// COMBINE
//Combine all into a single image collection

//List of list image numbers to map over
var imageNums = ee.List.sequence(0, 63);

//Function to combine all variable bands into single images & clip to domain of interest
var combine = function(imageNum) {
  var combined = ee.Image(precipmins.get(imageNum));
  combined = combined
    .addBands(ee.Image(pdsimins.get(imageNum)))
    .addBands(soilmins.get(imageNum))
    .addBands(tempmaxs.get(imageNum))
    .addBands(vpdmaxs.get(imageNum))
    .addBands(deficitmaxs.get(imageNum));
  return combined.clip(domain);
}; 

//Run function, print output, visualize sample to check
var allvariables = ee.ImageCollection.fromImages(imageNums.map(combine));
print('All variable relevant min/max & month of occurence:', allvariables);
Map.addLayer(allvariables.filter(ee.Filter.eq('year', 1958)).select('tmmx_month'), {min: 0, max: 12, palette:['blue', 'red']}, 'tmmx-month', false);



////////////// EXPORT
// Export all images of interest using batch download from fitoprincipe repo: https://github.com/fitoprincipe/geetools-code-editor
// To auto-run all tasks directly in GEE console, can follow instructions at bottom of this blog:
// https://benny.istan.to/blog/20220319-batch-task-execution-in-google-earth-engine-code-editor
// which utilizes code from https://github.com/kongdd/gee_monkey

//Export monthly means images
var batch = require('users/fitoprincipe/geetools:batch');
batch.Download.ImageCollection.toDrive(monthmeans, foldername,
  {region: domain,
  name: 'MonthlyMeans_Month{month}',
  scale: 4638.3
  });

//Export all variable images
var batch = require('users/fitoprincipe/geetools:batch');
batch.Download.ImageCollection.toDrive(allvariables, foldername,
  {region: domain,
  name: 'AllVariables_{year}',
  scale: 4638.3
  });

