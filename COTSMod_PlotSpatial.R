#### Functions

scale <- function(x) {
  x <- (x- min(na.omit(x))) / (max(na.omit(x))- min(na.omit(x)))
  x
}

addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

northarrow <- function(loc,size,bearing=0,cols,cex=1,...) {
  # checking arguments
  if(missing(loc)) stop("loc is missing")
  if(missing(size)) stop("size is missing")
  # default colors are white and black
  if(missing(cols)) cols <- rep(c("white","black"),8)
  # calculating coordinates of polygons
  radii <- rep(size/c(1,4,2,4),4)
  x <- radii[(0:15)+1]*cos((0:15)*pi/8+bearing)+loc[1]
  y <- radii[(0:15)+1]*sin((0:15)*pi/8+bearing)+loc[2]
  # drawing polygons
  for (i in 1:15) {
    x1 <- c(x[i],x[i+1],loc[1])
    y1 <- c(y[i],y[i+1],loc[2])
    polygon(x1,y1,col=cols[i])
  }
  # drawing the last polygon
  polygon(c(x[16],x[1],loc[1]),c(y[16],y[1],loc[2]),col=cols[16])
  # drawing letters
  b <- c("E","N","W","S")
  for (i in 0:3) text((size+par("cxy")[1])*cos(bearing+i*pi/2)+loc[1],
                      (size+par("cxy")[2])*sin(bearing+i*pi/2)+loc[2],b[i+1],
                      cex=cex)
}

scalebar <- function(loc,length,unit="km",division.cex=.8,...) {
  if(missing(loc)) stop("loc is missing")
  if(missing(length)) stop("length is missing")
  x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
  y <- c(0,length/(10*3:1))+loc[2]
  cols <- rep(c("black","white"),2)
  for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
  for (i in 1:5) segments(x[i],y[2],x[i],y[3])
  labels <- (x[c(1,3)]-loc[1])*100
  labels <- append(labels,paste((x[5]-loc[1])*100,unit))
  text(x[c(1,3,5)],y[4],labels=labels,adj=.5,cex=division.cex)
}


# setwd("/export/home/q-z/smatthew/COTS_SDM")



library(PBSmapping)
library(maptools)
library(RColorBrewer)
library(colorRamps)
library(SDMTools)
library(plotfunctions)
# Import and rotate spatial layers
shape <- importShapefile("Data/SDE_OWNER_crcgis_land250/SDE_OWNER_crcgis_land250", readDBF=FALSE)
shape2 <- importShapefile("Data/ManagementAreas/Management_Areas_of_the_GBRMP__Poly_", readDBF=FALSE)

shape.poly <- rgdal::readOGR("Data/SDE_OWNER_crcgis_land250/SDE_OWNER_crcgis_land250.shp", p4s = "+proj=longlat +datum=WGS84")
shape.elide <- elide(shape.poly, rotate=45)
shape45 <- SpatialPolygons2PolySet(shape.elide)

shape2.poly <- rgdal::readOGR("Data/ManagementAreas/Management_Areas_of_the_GBRMP__Poly_.shp")
shape2.poly <- spTransform(shape2.poly, crs(shape.poly))
shape2.elide <- elide(shape2.poly, rotate=45)
shape245 <- SpatialPolygons2PolySet(shape2.elide)

reefs.poly <- rgdal::readOGR("Data/Indicative_Reef_boundary/Indicative_Reef_boundary.shp", p4s = "+proj=longlat +datum=WGS84")
reefs.elide <- elide(reefs.poly, rotate=45)
reefs45 <- SpatialPolygons2PolySet(reefs.elide)

grid <- SpatialPoints(cbind(data.grid$lon, data.grid$lat), proj4string = CRS("+proj=longlat +datum=WGS84"))
grid.elide <- elide(grid, bb=bbox(shape.poly), rotate=45)
grid45 <- as.data.frame(grid.elide)

gridln <- gridlines(shape.poly)
gridln.elide <- elide(gridln, bb=bbox(shape.poly), rotate=45)
gridln45 <- SpatialLines2PolySet(gridln.elide)

# obs <- SpatialPoints(cbind(COTS_Manta$X_COORD, COTS_Manta$Y_COORD), proj4string = CRS("+proj=longlat +datum=WGS84"))
# obs.elide <- elide(obs, bb=bbox(shape.poly), rotate=45)
# obs45 <- as.data.frame(obs.elide)

# pst.grid <- SpatialPoints(cbind(data.WQ$X, data.WQ$Y), proj4string = CRS("+proj=longlat +datum=WGS84"))
# pst.elide45 <- elide(pst.grid, bb=bbox(shape.poly), rotate=45)
# pst45 <- as.data.frame(pst.elide45)


