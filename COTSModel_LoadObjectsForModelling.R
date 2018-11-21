##########################
#  Script for modeling COTS outbreaks in the Great Barrier Reef
#    This script sets up the data structures necessary for modeling GOTS outbreaks in the GBR 
#      
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
# 
#  05 April 2015 -- started scripting

##########################


setwd(DATA_DIRECTORY)
 
PopData <- read.csv("Environmental/Environmental_Data.csv", header = TRUE)[1:npops,] # we will just use a subset for testing
COTS.data <- read.csv("Disturbance/CoTS_data.csv", header = TRUE)[1:npops,]
cotsconn = read.csv("Connectivity/cotsconn_nodes.csv", header = T)
colnames(cotsconn)[2] = "REEF_NAME"
COTS.ConnMat = as.matrix(read.csv("Connectivity/cotsconn.csv", header = T), nrow = npops)
reefs = read.csv("Connectivity/reefs.csv", header = F) %>% dplyr::pull(1)
row.names(COTS.ConnMat) = reefs; colnames(COTS.ConnMat)= reefs
COTS.ConnMat[1:10, 1:10]
setwd(BASE_DIRECTORY)
load(file = "R_Objects/ConnMatFull.Rdata")
#Remove rows and columns nolonger in dataset
setwd(ENVDATA_DIRECTORY)
gone = read.csv("data.grid.csv", header=T) %>% 
  dplyr::anti_join(PopData, by = c("lon", "lat")) %>% dplyr::pull(PIXEL_ID)
ConnMat.full = ConnMat.full[-gone, -gone]

setwd(SPATIALDATA_DIRECTORY)
shp = rgdal::readOGR("Indicative_Reef_boundary/Indicative_Reef_boundary.shp")
zone = rgdal::readOGR("Great_Barrier_Reef_Marine_Park_Zoning/Great_Barrier_Reef_Marine_Park_Zoning.shp")
shp2 = rgdal::readOGR("SDE_OWNER_crcgis_land250/SDE_OWNER_crcgis_land250.shp")
cotsconn = read.csv("../Chapter4_CoTS Connectivity/cotsconn_nodes.csv", header = T)
colnames(cotsconn)[2] = "REEF_NAME"

# Load Connectivity Matrix
setwd(BASE_DIRECTORY)
#load("R_Objects/ProbDistance.Rdata")
#load("R_Objects/ConnMatFull.Rdata")
#ConnMat = ConnMat.full[1:npops, 1:npops]
#rm(ConnMat.full)

# Load Fertilisation vs Density Models Estimates
load("R_Objects/FvDParams.Rdata")

# Load Environmental Data Frame
setwd(ENVDATA_DIRECTORY)
#PercentReef = read.csv("data.grid.csv", header=T)[,c(2,3,5)] %>% inner_join(PopData[,2:3], by = c("lon", "lat"))
data.grid = read.csv("data.grid.csv", header=T)[,c(2,3,5,44:50)] %>% 
  dplyr::inner_join(PopData, by = c("lon", "lat"))
data.grid = data.grid[c(11, 1:3, 12:51,4:10)]
WQ = PopData$Primary + PopData$Secondary + PopData$Tertiary
PopData = data.grid[1:7]

# Gompertz parameters: disturbance & WQ effect sizes (as per MacNeil et al. in prep)
bleaching.mn.sd <- c(-0.19, 0.01)
COTS.mn.sd <- c(-0.54, 0.04)
disease.mn.sd <- c(-0.13, 0.01)
storms.mn.sd <- c(-0.64, 0.01)
unknown.mn.sd <- c(-0.16, 0.01)
WQ.mn.sd <- c(-0.68, 0.03)


# gone = read.csv("data.grid.csv", header=T) %>% 
#   dplyr::anti_join(PopData, by = c("lon", "lat")) %>% dplyr::pull(PIXEL_ID)
