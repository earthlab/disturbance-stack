#Climate disturbance data analysis
#Tyler McIntosh, 9/9/2022

#This analysis is based off of Hammond et al 2022.
#Data inputs for this script were wrangled in GEE,
#copy of script found in RProject workflow document 02_ClimateDataWrangling.js.

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

#Create list of monthly raster stacks
monthMeans <- vector(mode = "list", length = 1) #Initialize empty list
for (file in monthFileNames) {
  monthData <- rast(file)
  monthMeans <- append(monthMeans, monthData)
  
}
monthMeans <- monthMeans[-1]
names(monthMeans) <- month.name

#Create list of yearly variable raster stacks
yearVarData <- vector(mode = "list", length = 1) #Initialize empty list
for (file in yearVarFileNames) {
  yearDats <- rast(file)
  yearVarData <- append(yearVarData, yearDats)
  
}
yearVarData <- yearVarData[-1]
names(yearVarData) <- as.character(seq(1958, 2021))

################ PERFORM ANALYSIS ###############
#Analysis per Hammond et al 2022


