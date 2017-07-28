##############################################
###  GIS functions for the COTS population model
#################################### 



###################
#  ReadRaster
##########
# OBJECTIVE:
#    Generic function for reading in ascii raster layers for COTS population model.
#    uses the "raster" package.
# PARAMS:
#    - rastername: name of the raster layer to be loaded into memory (must be .ASC file)
#        NOTE: this raster layer must be stored in the spatial layers project resource file in Dropbox
#    - projection: cartographic projection in PROJ4 style
#    - plot: logical (T or F) indicating whether the layer should be plotted
# RETURNS:
#    - raster object (see "raster" package for more details)
###################

ReadRaster <- function(rastername,projection=projection,plot=F){
  setwd(SPATIALDATA_DIRECTORY)
  newraster <- readGDAL(rastername,p4s=projection)
  newraster <- raster(newraster)
  if(plot) plot(newraster)
  return(newraster)
}

###################
#  readShapefile
##########
# OBJECTIVE:
#    Generic function for reading in ESRI shapefile layers for COTS population model.
#    uses the "raster" package.
# PARAMS:
#    - shapefilename: name of the shapefile layer to be loaded into memory
#    - projection: cartographic projection in PROJ4 style
#    - plot: logical (T or F) indicating whether the layer should be plotted
# RETURNS:
#    - Spatial points, lines, or polygon object (see "sp" package for more details)
###################
readShapefile <- function(shapefilename,projection=projection,plot=FALSE){
  setwd(SPATIALDATA_DIRECTORY)
  shapefilename=gsub(".shp","",shapefilename)
  datadir <- sprintf("%s\\%s",SPATIALDATA_DIRECTORY,shapefilename)   
  setwd(datadir)
  newshapefile <- readOGR(datadir,layer=shapefilename)
  if(plot) plot(newshapefile)
  return(newshapefile)
}





