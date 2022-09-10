#Climate disturbance data analysis
#Tyler McIntosh, 9/9/2022

#This analysis is based off of Hammond et al 2022.
#Data inputs for this script were wrangled in GEE,
#copy of script found in RProject workflow document 02_ClimateDataWrangling.js.

# Analysis per Hammond et al 2022:
# "For each climate variable (TMAX, VPD, CWD, SOIL M, PPT, PDSI, we calculated the anomaly 
# (for example, TMAX of the warmest month in the year of mortality - TMAX average for the 
# 61 values of the same month from 1958–2019). To provide cross-comparison between metrics
# with different scales and across disparate climates, anomalies were standardized into 
# zscores, such that time series had a mean of 0 and a standard deviation of 1 based on 
# using the entire period of record 1958–2019."

# Variable averages for each month across the time period, as well as max/mins/month in 
# which occurred all part of GEE outputs.

#This script uses the following naming conventions wherever possible:
# lowerCamelCase for variables
# period.separated for functions
# underscore_separated for files

############ SETUP WORKSPACE ############

###Standard libraries, unneeded libraries commented out

#Standard libraries
library(tidyverse) #Tidyverse!
library(here) #Relative path best practices
#library(kableExtra) #Table creation
#library(knitr) #For use with R markdown

#Geographic libraries
library(terra) #New raster data package
#library(sf) #New vector data package
#library(tmap) #Thematic mapping
#library(tmaptools) #Supports Tmap
#library(leaflet) #For interactive web mapping
#library(landscapemetrics) #Fragstats alternative

###Clean workspace
rm(list=ls()) #Ensure empty workspace
here() #Check here location

################## LOAD DATASETS #############
#Get file names of severity data
monthFileNames <- list.files(path = here("data", "climate", "GEE_terraclimate_prewrangled"), 
                               pattern ="MonthlyMeans*", full.names = TRUE)
yearVarFileNames <- list.files(path = here("data", "climate", "GEE_terraclimate_prewrangled"), 
                               pattern ="AllVariables*", full.names = TRUE)

#Create SpatRasterDataset of monthly rasters
monthMeans <- vector(mode = "list", length = 1) #Initialize empty list
for (file in monthFileNames) {
  monthData <- rast(file)
  monthMeans <- append(monthMeans, monthData)
  
}
monthMeans <- monthMeans[-1]
names(monthMeans) <- month.name
monthmeans <- sds(monthmeans)  #turn list into spatrasterdataset

#Create SpatRasterDataset of yearly variable rasters
yearVarData <- vector(mode = "list", length = 1) #Initialize empty list
for (file in yearVarFileNames) {
  yearDats <- rast(file)
  yearVarData <- append(yearVarData, yearDats)
  
}
yearVarData <- yearVarData[-1]
names(yearVarData) <- as.character(seq(1958, 2021))
yearVarData <- sds(yearVarData) #turn list into spatrasterdataset

#Compare rasters of both types to ensure same geometries
compareGeom(monthMeans$January, yearVarData$`1958`)

################ CALCULATE ANOMALY ###############

#calculate anomaly for precipitation for one year first, test method. Then apply over entire list
testyear <- yearVarData$`1958`

