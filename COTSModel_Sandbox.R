########################
# COTSModel_Sandbox - Rough Working Out etc
########################
#   Authors: Kevin Shoemaker, Sam Matthews

PopData = read.table("PopData.txt", header = T, sep="\t")

setwd(SPATIALDATA_DIRECTORY)
reefpercent <- ReadRaster("reefpercentraster.asc",projection=projection,plot=F)
plot(reefpercent)


