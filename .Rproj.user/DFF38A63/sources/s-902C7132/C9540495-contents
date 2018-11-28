##########################
#  Script for modeling COTS outbreaks in the Great Barrier Reef
#    This script sets up the data structures necessary for modeling GOTS outbreaks in the GBR 
#      
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
# 
#  05 April 2015 -- started scripting

##########################


setwd(DATA_DIRECTORY)
PopData <- read.csv("Environmental/Environmental_Data.csv", header = TRUE) # we will just use a subset for testing
data.COTS <- read.csv("Disturbance/CoTS_data.csv", header = TRUE)
# cotsconn = read.csv("Connectivity/cotsconn_nodes.csv", header = T)
# colnames(cotsconn)[2] = "REEF_NAME"
# COTSMat = as.matrix(read.csv("Connectivity/cotsconn.csv", header = T))
reefs = read.csv("Connectivity/reefs.csv", header = F) %>% dplyr::pull(1)
# row.names(COTS.ConnMat) = reefs; colnames(COTS.ConnMat)= reefs
# COTS.ConnMat[1:10, 1:10]
load("Connectivity/COTS.ConnMat.Rdata")
# setwd(BASE_DIRECTORY)
# load(file = "R_Objects/ConnMatFull.Rdata")
#Remove rows and columns nolonger in dataset
# setwd(ENVDATA_DIRECTORY)
# gone = read.csv("data.grid.csv", header=T) %>% 
#   dplyr::anti_join(PopData, by = c("lon", "lat")) %>% dplyr::pull(PIXEL_ID)
# ConnMat.full = ConnMat.full[-gone, -gone]

# setwd(SPATIALDATA_DIRECTORY)
# shp = rgdal::readOGR("Indicative_Reef_boundary/Indicative_Reef_boundary.shp")
# zone = rgdal::readOGR("Great_Barrier_Reef_Marine_Park_Zoning/Great_Barrier_Reef_Marine_Park_Zoning.shp")
# shp2 = rgdal::readOGR("SDE_OWNER_crcgis_land250/SDE_OWNER_crcgis_land250.shp")

ourreefs = unique(as.character(PopData$REEF_NAME))
ourreefs = ourreefs[which(unique(as.character(PopData$REEF_NAME)) %in% reefs)]
PopData = dplyr::filter(PopData, REEF_NAME %in% ourreefs) %>% dplyr::mutate(REEF_NAME = factor(REEF_NAME))
COTS.data = dplyr::filter(COTS.data, REEF_NAME %in% ourreefs) %>% dplyr::mutate(REEF_NAME = factor(REEF_NAME))

Pixels = PopData %>% dplyr::group_by(REEF_NAME) %>% dplyr::summarise(Pixels = length(REEF_NAME))
Pixels = Pixels[match(ourreefs,Pixels$REEF_NAME),]
# npops = length(PopData$lon)
######!
# Convert Connectivity Matrix
######!

# The issue here is to distribute larvae among the reefs in Karlo's Connectivty Matrix
#  1. We need to combine all larvae produced at a Reef
#  2. These are then distributed to the other reefs via the connectivity matrix
#  3. Larvae distibuted evenly among cells within the reef

## need to reorder connectivity matrix to match our reefs
# Which reefs from Karlo's mtrix are present in ours?

        # ourreefs = unique(as.character(PopData$REEF_NAME))
        # ourreefs = ourreefs[which(unique(as.character(PopData$REEF_NAME)) %in% reefs)]
        # whichCONNreefs = which(reefs %in% ourreefs)
        # COTS.ConnMat = COTS.ConnMat[whichCONNreefs, whichCONNreefs]  
        # dim(COTS.ConnMat)
        # 
        # # check that we have the the same number of reefs and that there are no reefs missing from either dataset 
        # length(colnames(COTS.ConnMat)); length(whichourreefs)
        # which(!whichourreefs %in% colnames(COTS.ConnMat))
        # 
        # matchCONN = match(whichourreefs, colnames(COTS.ConnMat))
        # COTS.ConnMat = COTS.ConnMat[matchCONN, matchCONN]
        # head(whichourreefs); head(colnames(COTS.ConnMat))
        # setwd(BASE_DIRECTORY)
        # 
        # setwd(DATA_DIRECTORY)
        # save(COTS.ConnMat, file = "Connectivity/COTS.ConnMat.Rdata")
        # rm(COTS.ConnMat)
        # load("Connectivity/COTS.ConnMat.Rdata")


# Load Fertilisation vs Density Models Estimates
setwd(BASE_DIRECTORY)
load("R_Objects/FvDParams.Rdata")

# Load Environmental Data Frame
setwd(ENVDATA_DIRECTORY)
#PercentReef = read.csv("data.grid.csv", header=T)[,c(2,3,5)] %>% inner_join(PopData[,2:3], by = c("lon", "lat"))
data.grid = read.csv("data.grid.csv", header=T)[,c(2,3,5,44:50)] %>% 
  dplyr::inner_join(PopData, by = c("lon", "lat"))
data.grid = data.grid[c(11, 1:3, 12:51,4:10)]
WQ = PopData$Primary + PopData$Secondary + PopData$Tertiary
PopData = data.grid[1:7]

# PopData = PopData[1:npops,]
# COTS.data = COTS.data[1:npops,]
# data.grid = data.grid[1:npops,]

check = cbind(colnames(COTS.ConnMat), as.character(unique(COTS.data$REEF_NAME)), 
              as.character(unique(PopData$REEF_NAME)), as.character(unique(data.grid$REEF_NAME)))

# Gompertz parameters: disturbance & WQ effect sizes (as per MacNeil et al. in prep)
bleaching.mn.sd <- c(-0.19, 0.01)
COTS.mn.sd <- c(-0.54, 0.04)
disease.mn.sd <- c(-0.13, 0.01)
storms.mn.sd <- c(-0.64, 0.01)
unknown.mn.sd <- c(-0.16, 0.01)
WQ.mn.sd <- c(-0.68, 0.03)


# gone = read.csv("data.grid.csv", header=T) %>% 
#   dplyr::anti_join(PopData, by = c("lon", "lat")) %>% dplyr::pull(PIXEL_ID)

