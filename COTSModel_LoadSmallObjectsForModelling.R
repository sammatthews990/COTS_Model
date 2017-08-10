##########################
#  Script for modeling COTS outbreaks in the Great Barrier Reef
#    This script sets up the data structures necessary for modeling GOTS outbreaks in the GBR 
#      
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
# 
#  05 April 2015 -- started scripting

##########################


setwd(DATA_DIRECTORY)
 
PopData <- read.csv("reefs.csv", header = TRUE)[1:npops,] # we will just use a subset for testing
COTS.data <- read.csv("Disturbance/CoTS_data.csv", header = TRUE)[1:npops,]

setwd(BASE_DIRECTORY)
load("R_Objects/FvDParams.RData")
# # Load Connectivity Matrix
# setwd(BASE_DIRECTORY)
# #load("R_Objects/ProbDistance.Rdata")
# load("R_Objects/ConnMatFull.Rdata")
# ConnMat = ConnMat.full[1:npops, 1:npops]
# rm(ConnMat.full)

# Load Fertilisation vs Density Models Estimates
setwd(BASE_DIRECTORY)
load("R_Objects/FvDParams.Rdata")

# Load Environmental Data Frame
setwd(ENVDATA_DIRECTORY)
data.grid = read.csv("data.grid.csv", header=T)[1:npops,]
WQ = data.grid$Primary + data.grid$Secondary + data.grid$Tertiary

# Gompertz parameters: disturbance & WQ effect sizes (as per MacNeil et al. in prep)
bleaching.mn.sd <- c(-0.19, 0.01)
COTS.mn.sd <- c(-0.54, 0.04)
disease.mn.sd <- c(-0.13, 0.01)
storms.mn.sd <- c(-0.64, 0.01)
unknown.mn.sd <- c(-0.16, 0.01)
WQ.mn.sd <- c(-0.68, 0.03)


