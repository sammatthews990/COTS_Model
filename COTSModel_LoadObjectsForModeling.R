##########################
#  Script for modeling COTS outbreaks in the Great Barrier Reef
#    This script sets up the data structures necessary for modeling GOTS outbreaks in the GBR 
#      
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
# 
#  05 April 2015 -- started scripting

##########################



###############################
# Prepare the "XYZ"-type dataframe (each row represents attributes for one population)
###############################

setwd(ENVDATA_DIRECTORY)
PopData <- read.table("ENV_dataGBR.txt",header=T,sep="\t")
head(PopData)

     #  harvest key global variables from the input file
NPOPS <- nrow(PopData)
UNIQUEREEFIDS <<- unique(PopData$REEF_ID)
NREEFS <<- length(UNIQUEREEFIDS)

#### determine reef percent for each cell of interest in the 1km grid

reefpercent <- ReadRaster("reefPercentRaster.asc",projection=projection)
PopData[,'PercentReef'] <- extract(reefpercent,PopData[,c('x','y')])   # add a new column to the envData (XYZ format)
rm(reefpercent) # remove from memory (save space!)


###############################
# MAKE UP DATA INPUTS
###############################

initDensityA <- round(rnorm(NREEFS,10000,1000))    # density per km2 of reef habitat
initDensityS <- round(rnorm(NREEFS,1000,100))










