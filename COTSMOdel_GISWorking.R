##########################
#  R script: generate key GIS layers for modeling COTS population in the GBR 
#
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
#  
#  08 April 2015 -- started scripting
#
#  26 September 2016 -- redefined GIS layers to comply with GBRMPA MArineBioregions_WGS84 Shapefile
#
#         URGENT: ENV DATA does not cover 2/3 of our reef site
#         URGENT: Reef percent raster does not seem to line up with shapefile in many areas

##########################


###############################
#  PREPARE KEY GIS LAYERS
###############################

################
### create "reefraster" layer: represents percent of each grid cell with reef coverage. 

### UPDATE 28 September 2009 -- Convert to using the MarineBioregions_WGS84 Shapefile as Basis
################

### Create template raster for the GBR (proper extent and resolution)
################

setwd(ENVDATA_DIRECTORY)    # first load the xyz data for each population of interest 
EnvData <- read.table("ENV_dataGBR.txt",header=T,sep="\t")


#NPOPS <- nrow(EnvData)  

minx <- min(EnvData$x)-0.005  # farthest west
miny <- min(EnvData$y)-0.005  # farthest south
maxx <- max(EnvData$x)+0.005  # farthest east
maxy <- max(EnvData$y)+0.005  # farthest north
studyRegion <- extent(minx,maxx,miny,maxy)
cat(sprintf("Western margin is %s degrees, southern margin is %s degrees, eastern margin is %s degrees, and northern margin is %s degrees",minx,miny,maxx,maxy))

plot(studyRegion, main="Study Region (blank)")

# save study region extent (as "extent" object)
setwd(SPATIALDATA_DIRECTORY)
save(studyRegion,file="studyRegion.RData")


template <- raster(ext=studyRegion,resolution=0.01,vals=0, crs = projection)   # create empty raster
#template001 <- raster(ext=studyRegion,resolution=0.001,vals=0, crs = projection)   # create empty raster at 0.001 res
#plot(template)
setwd(SPATIALDATA_DIRECTORY)
writeRaster(template,filename="templateRaster.asc",format="ascii",overwrite=T)   # write to file

##################
## create reefraster from GBRMPA Classification at 0.001 resolution
##################
reefshape <- readShapefile("MarineBioregions_WGS84", projection = projection, plot=T)
reefraster001 <- rasterize(reefshape, template001, field = "REEF_ID", crs = projection)
setwd(SPATIALDATA_DIRECTORY)
writeRaster(reefraster,filename="reefraster001.asc",format="ascii",overwrite=T)   # write to file

summary(reefshape)
sum(reefpercent@data@values>50)
##################
## Create ReefPercent Raster (aggregate reefraster to 0.01 resolution to calculate percent reef in each 0.01 cell)
##################

reefraster01 <- reclassify(reefraster001,rcl=c(NA,NA,0, -Inf,0.5,0, 0.501,Inf,1))  # temporary layer: reclassify to binary
compute_percentReef <- function(t,na.rm=TRUE) sum(t,na.rm)  ## aggregation function
reefpercent <- aggregate(reefraster01, fact=10, fun=compute_percentReef, expand=TRUE, na.rm=TRUE) 
reefpercent <- reefpercent-1    # for some reason, the result ends up between 1 and 101- reformulate for proper percent
plot(reefpercent)
setwd(SPATIALDATA_DIRECTORY)
writeGDAL(reefpercent, fnam="reefpercent.tif", drivername = "GTiff", type = "Float32")
writeRaster(reefpercent,filename="reefPercentRaster.tif",format="GTiff",overwrite=T)   # write to file



##################
## Create ReefBioregion ReefID raster files  
##################

reefbioregion <- rasterize(reefshape, template, field = "BIOREGION")
writeRaster(reefbioregion,filename="reefbioregion.asc",format="ascii",overwrite=T)

reefID <- rasterize(reefshape, template, field = "OBJECTID")
writeRaster(reefID,filename="reefID.asc",format="ascii",overwrite=T)

##################
## create a XYZ grid from reefshape keeping all fields
##################
reefpercent.df <- data.frame(data.frame(x = coordinates(reefpercent)[,1], y = coordinates(reefpercent)[,2], reefpercent = reefpercent@data@values))
reefpercent.df <- subset(reefpercent.df, reefpercent>0)
coords <- reefpercent.df[,1:2]
sp <- SpatialPoints(coords = coords, proj4string = CRS(projection))

#extract vals from shapefile for every pixel that contains reef
vals <- over(sp, reefshape)
reefs <- cbind(coordinates(sp), vals, reefpercent = reefpercent.df[,3])

#create subset of pixels that contain reef but aren't assigned attributes
reef.NA <- subset(reefs, is.na(REEF_ID))
#creat subset that do have attributes
reef.YES <- subset(reefs, !is.na(REEF_ID))
#convert both the spatial points dataframes to determine the geographic distance between them
coordinates(reef.NA) <- ~x+y; coordinates(reef.YES) <- ~x+y
Gdist <- gDistance(reef.NA,reef.YES, byid = T)
min.d <- apply(Gdist, 2, which.min) #creates vector of the minimum distances which we use to index the cords 
min.gd <- apply(Gdist, 2, min) #creates a vecor of the min dist so we can see how far they have been displaced

reef.NA$min.gd <- min.gd
summary(min.gd)
#Double check that the matches are ok
reef.NAdf <- as.data.frame(reef.NA)
reef.YESdf <- as.data.frame(reef.YES)
coord.matches <- cbind(reef.NAdf[,1:2],reef.YESdf[min.d,1:2], min.gd)

#For now I'll just discard reefs that have no info from the associated polygons

##################
#Check how all shapefiles, rasters and environmental data line up
##################
#
# THis section is used to make sure we have all the correct data for every reef cell
# At the moment a lot of reef cells -- especially those with high reef percentage are not represented in our enviro data
# This means we may have to interpolate the dataset further - we only have data for about 1/3
envcoords <- EnvData[,2:3]
envpoints <- SpatialPoints(coords = envcoords, proj4string = CRS(projection))


plotzoom <- function (zoomextent) {
  
  
  plot(reefpercent, xlim=c(zoomextent[1],zoomextent[2]), ylim=c(zoomextent[3],zoomextent[4]))
  #plot(reefraster001, xlim=c(zoomextent[1],zoomextent[2]), ylim=c(zoomextent[3],zoomextent[4]), add = TRUE)
  plot(reefshape, xlim=c(zoomextent[1],zoomextent[2]), ylim=c(zoomextent[3],zoomextent[4]), add = TRUE)
  points(sp, xlim=c(zoomextent[1],zoomextent[2]), ylim=c(zoomextent[3],zoomextent[4]), pch=20, col="blue", cex=0.5)
  points(envpoints, xlim=c(zoomextent[1],zoomextent[2]), ylim=c(zoomextent[3],zoomextent[4]), pch=17, col="red", cex=0.5)
}
windows()
plotzoom(c(151.4, 151.8, -23.03,-22.55))


plot(reefpercent)
plot(reefshape, add=T)
points(envpoints, pch=17, col="blue", cex=0.)


#merge env data with reef data

data <- merge(reef.YESdf, EnvData, by = c("x", "y"))
setwd(DATA_DIRECTORY)
write.table(data, "ENVData_all.txt", row.names = F, sep = "\t")




totalCells <- length(template@data@values)    # 1,335,372 cells
cat(sprintf("total cells in study region raster is %s",totalCells))


################
### Create Reef ID layer for the GBR
################

## rasterize an XYS variable.
reefID <- as.numeric(PopData$REEF_ID)
ndx <- !is.na(reefID)
reefIDraster <- rasterize(PopData[,c('x','y')][ndx,], template, field=reefID[ndx])   # memory limitations

setwd(SPATIALDATA_DIRECTORY)
writeRaster(reefIDraster,filename="reefIDRaster.asc",format="ascii",overwrite=T)   # write to file
plot(reefIDraster)

rm(reefIDraster)


ReadRaster("reefraster2.asc", projection = projection, plot=t)
