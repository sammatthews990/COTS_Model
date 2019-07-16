##########################
#  Script for modeling COTS outbreaks in the Great Barrier Reef
#    This script sets up the data structures necessary for modeling GOTS outbreaks in the GBR 
#      
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
# 
#  05 April 2015 -- started scripting

##########################
#########################
# SET PROJECT DIRECTORIES (this should be the only place where local directories should be referenced)
#########################

USER = "SAM"

if(USER=="KEVIN") BASE_DIRECTORY <- "C:\\Users\\Kevin\\Dropbox\\CoTS_Model"             # NOTE: this should link to the Dropbox folder with shared project resources	                                                                        
if(USER=="KEVIN") CODE_DIRECTORY <- "C:\\Users\\Kevin\\GIT\\COTS_Model"              # NOTE: code directory should be your local copy of the GitHub repository

if(USER=="SAM") BASE_DIRECTORY <- "C:\\Users\\jc312264\\Dropbox\\CoTS_Model"
if(USER=="SAM") CODE_DIRECTORY <- "C:\\Users\\jc312264\\OneDrive - James Cook University\\COTS_Model"

if(USER=="SAM_UNI") BASE_DIRECTORY <- "C:\\Users\\jc312264\\Dropbox\\CoTS_Model"
if(USER=="SAM_UNI") CODE_DIRECTORY <- "C:\\Users\\jc312264\\OneDrive - James Cook University\\GitHub\\COTS_Model"

if(USER=="SAM_RENO") BASE_DIRECTORY <- "C:\\Users\\kshoemaker\\Dropbox\\CoTS_Model"
if(USER=="SAM_RENO") CODE_DIRECTORY <- "C:\\Users\\kshoemaker\\Documents\\GitHub\\COTS_Model"

SPATIALDATA_DIRECTORY <- paste(BASE_DIRECTORY,"\\Spatial Layers",sep="")                          # directory for storing relevant spatial data (ASC, SHP files)
if(is.na(file.info(SPATIALDATA_DIRECTORY)[1,"isdir"])) dir.create(SPATIALDATA_DIRECTORY)

DATA_DIRECTORY <- paste(BASE_DIRECTORY,"\\Data",sep="")                                 # directory for storing data (CSV files)
if(is.na(file.info(DATA_DIRECTORY)[1,"isdir"])) dir.create(DATA_DIRECTORY)

ENVDATA_DIRECTORY <- paste(DATA_DIRECTORY,"\\Environmental",sep="")                                 # directory for storing data (CSV files)
if(is.na(file.info(ENVDATA_DIRECTORY)[1,"isdir"])) dir.create(ENVDATA_DIRECTORY)

FIGURES_DIRECTORY <- paste(BASE_DIRECTORY,"\\Figures\\RawFigures",sep="")               # directory for storing raw figures generated from R
if(is.na(file.info(FIGURES_DIRECTORY)[1,"isdir"])) dir.create(FIGURES_DIRECTORY)

RESULTS_DIRECTORY <- paste(BASE_DIRECTORY,"\\results",sep="")                           # directory for storing relevant results
if(is.na(file.info(RESULTS_DIRECTORY)[1,"isdir"])) dir.create(RESULTS_DIRECTORY)

RDATA_DIRECTORY <- paste(BASE_DIRECTORY,"\\R_Workspaces",sep="")                        # directory for storing .RData files (R workspaces and data objects)
if(is.na(file.info(RDATA_DIRECTORY)[1,"isdir"])) dir.create(RDATA_DIRECTORY)

cat(sprintf("The current user is %s",USER))

setwd(DATA_DIRECTORY)
PopData = read.csv("Environmental/Environmental_Data.csv", header = TRUE) # we will just use a subset for testing
# COTS.data = read.csv("Disturbance/CoTS_data.csv", header = TRUE)
# COTS.ConnMat = as.matrix(read.csv("Connectivity/cotsconn.csv", header = T))
reefs = read.csv("Connectivity/reefs.csv", header = F) %>% dplyr::pull(1)
# row.names(COTS.ConnMat) = reefs; colnames(COTS.ConnMat)= reefs
# COTS.ConnMat[1:10, 1:10]ourreefs = unique(as.character(PopData$REEF_NAME))
ourreefs = unique(as.character(PopData$REEF_NAME))
ourreefs = ourreefs[which(unique(as.character(PopData$REEF_NAME)) %in% reefs)]
PopData = dplyr::filter(PopData, REEF_NAME %in% ourreefs) %>% dplyr::mutate(REEF_NAME = factor(REEF_NAME))
# COTS.data = dplyr::filter(COTS.data, REEF_NAME %in% ourreefs) %>% dplyr::mutate(REEF_NAME = factor(REEF_NAME))

Pixels = PopData %>% dplyr::group_by(REEF_NAME) %>% dplyr::summarise(Pixels = length(REEF_NAME))
Pixels = Pixels[match(ourreefs,Pixels$REEF_NAME),]
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

load("Connectivity/COTS.ConnMat.Rdata")


# LOAD DISTURBANCE DATA ----

setwd(DATA_DIRECTORY)

# Environmental Data + BRT/MRT Predictions of Coral Growth Parameters
data.grid = read.table("CoralModel/XYZ_BRTpred_MRTpred_GG.csv", header = TRUE, sep = ",")
data.grid = data.grid[order(data.grid$LONG, data.grid$LAT),]; colnames(data.grid)[2:3] = c("lon", "lat")
data.grid = PopData %>% dplyr::inner_join(dplyr::select(data.grid, 2:3, 51:62), by=c("lon", "lat"))

data.manta = read.table("CoralModel/Manta.csv", header = TRUE, sep = ",")
data.manta.env = read.table("CoralModel/Manta_ENV.txt", header = TRUE, sep = "\t")
data.manta$REEF_ID= gsub("s","",tolower(data.manta$FULLREEF_ID))
data.manta$REEF_ID = paste(substr(data.manta$REEF_ID,1,2), substr(data.manta$REEF_ID, 3 ,6), sep = "-")
data.manta$REEF_ID = gsub("22-104","22-104a",data.manta$REEF_ID)
data.manta$REEF_ID = gsub("23-052","23-052a",data.manta$REEF_ID)
data.manta$REEF_ID = gsub("17-012a","17-012",data.manta$REEF_ID)
data.manta$REEF_ID = gsub("23-082","23-082a",data.manta$REEF_ID)
data.manta$REEF_ID = gsub("23-055","23-055a",data.manta$REEF_ID)
data.manta = dplyr::left_join(data.manta, unique(data.grid[,4:5]), by="REEF_ID")
# Modelled Disturbance
data.bleaching = read.csv("Disturbance/Bleaching_data_98_02_16.csv", header = TRUE) %>%
  dplyr::filter(REEF_NAME %in% ourreefs) %>% dplyr::mutate(REEF_NAME = factor(REEF_NAME))
data.COTS = read.csv("Disturbance/CoTS_data.csv", header = TRUE) %>%
  dplyr::filter(REEF_NAME %in% ourreefs) %>% dplyr::mutate(REEF_NAME = factor(REEF_NAME))
data.storms =  read.csv("Disturbance/Cyclones_data.csv", header = TRUE) %>%
  dplyr::filter(REEF_NAME %in% ourreefs) %>% dplyr::mutate(REEF_NAME = factor(REEF_NAME))
data.disease = read.table("CoralModel/Disturb_disease.txt", header = TRUE, sep = "\t")




# Observed Disturbance
data.ltmp.bleaching <- dplyr::inner_join(data.COTS[1:5], 
                                         dplyr::select(read.table("CoralModel/LTMP_B_XYZ_updated.txt", 
                                                                  header = TRUE, sep = "\t"), -REEF_ID, -PIXEL_ID ), 
                                         by=c("lon", "lat"))
data.ltmp.COTS <- dplyr::inner_join(data.COTS[1:5], 
                                    dplyr::select(read.table("CoralModel/LTMP_C_XYZ_updated.txt", 
                                                             header = TRUE, sep = "\t"), -REEF_ID, -PIXEL_ID ), 
                                    by=c("lon", "lat"))
data.ltmp.disease <- dplyr::inner_join(data.COTS[1:5], 
                                       dplyr::select(read.table("CoralModel/LTMP_D_XYZ_updated.txt", 
                                                                header = TRUE, sep = "\t"), -REEF_ID, -PIXEL_ID ), 
                                       by=c("lon", "lat"))
data.ltmp.storms <-  dplyr::inner_join(data.COTS[1:5], 
                                       dplyr::select(read.table("CoralModel/LTMP_D_XYZ_updated.txt", 
                                                                header = TRUE, sep = "\t"), -REEF_ID, -PIXEL_ID ), 
                                       by=c("lon", "lat"))
data.ltmp.unknown <-  dplyr::inner_join(data.COTS[1:5], 
                                        dplyr::select(read.table("CoralModel/LTMP_U_XYZ_updated.txt", 
                                                                 header = TRUE, sep = "\t"), -REEF_ID, -PIXEL_ID ), 
                                        by=c("lon", "lat"))


# More Gompertz Parameters
data.WQ <- read.table("CoralModel/pst_grid.txt", header = TRUE, sep = "\t")

data.reef <- read.table("CoralModel/MRT_reef_clust_Gompertz.txt", header = TRUE, sep = "\t") 
names(data.reef)[5:7] <- c("clust.bent", "b0.mu", "b0.sd")

data.rap <- read.table("CoralModel/KERCOD_GP_Bent_Env3.txt", header = TRUE, sep = "\t")
rap.reefs <- data.rap$REEF[data.rap$P_CODE=="RAP"]

# Gompertz parameters: disturbance & WQ effect sizes (as per MacNeil et al. in prep)
bleaching.mn.sd <- c(-1.0182, 0.1028)
COTS.mn.sd <- c(-0.1713, 0.0097)
disease.mn.sd <- c(-0.1681, 0.0303)
storms.mn.sd <- c(-0.5232, 0.01173)
unknown.mn.sd <- c(-0.2681, 0.0331)

WQ_bleach <- c(0.851062345,0.096005095)
WQ_CoTS <- c(-2.034101163,0.134657287)
WQ_Cyclone <- c(0.649030968,0.056093201)
WQ_Disease <- c(-0.369621141,0.15480862)


# REMOVE NORTHERNMOST REEFS FROM SPATIAL TABLES
# data.ltmp.bleaching <- data.ltmp.bleaching[data.grid$lat < (-14),]
# data.ltmp.COTS <- data.ltmp.COTS[data.grid$lat < (-14),]
# data.ltmp.disease <- data.ltmp.disease[data.grid$lat < (-14),]
# data.ltmp.storms <- data.ltmp.storms[data.grid$lat < (-14),]
# data.ltmp.unknown <- data.ltmp.unknown[data.grid$lat < (-14),]
# 
# data.bleaching <- data.bleaching[data.grid$lat < (-14),]
# data.COTS <- data.COTS[data.grid$lat < (-14),]
# data.disease <- data.disease[data.grid$lat < (-14),]
# data.storms <- data.storms[data.grid$lat < (-14),]
# data.grid <- data.grid[data.grid$lat < (-14),]


WQ <- data.grid$Primary + data.grid$Secondary + data.grid$Tertiary

#### load Predicted chl from eReefs ----

load("eReefsPredictions.Rdata")
data.Chl = dat.predict.med.Chl %>% dplyr::inner_join(data.grid[,1:3], by = c("PIXEL_ID","lat", "lon"))
data.Salt = dat.predict.med.Salt %>% dplyr::inner_join(data.grid[,1:3], by = c("PIXEL_ID","lat", "lon"))
data.Temp = dat.predict.med.Temp %>% dplyr::inner_join(data.grid[,1:3], by = c("PIXEL_ID","lat", "lon"))

notin = dplyr::anti_join(dat.predict.med.Chl[,1:3], data.grid[,1:3])

data.chl.resid = dat.resid.chl[-notin$PIXEL_ID,,]
rm(dat.resid.chl, dat.predict.med.Chl, dat.predict.med.Salt, dat.predict.med.Temp, notin)

# Convert bleaching scores to mid points ------

#data.bleaching[,c("X1998","X2002","X2016")] <- round(data.bleaching[,c("X1998","X2002","X2016")])
data.bleaching.test = data.COTS; colnames(data.bleaching.test)[6:38] = paste0("bleach_", 1985:2017)
data.bleaching.test[6:38] = NA; data.bleaching.test[, c("bleach_1998", "bleach_2002", "bleach_2016")] = data.bleaching[6:8]

data.bleaching.test[6:38] = lapply(data.bleaching.test[6:38], as.numeric)
data.bleaching.bckp <- data.bleaching <- data.bleaching.test

data.bleaching[,-(1:5)][data.bleaching.bckp[,-(1:5)]==1] <- 0.05
data.bleaching[,-(1:5)][data.bleaching.bckp[,-(1:5)]==2] <- 0.15
data.bleaching[,-(1:5)][data.bleaching.bckp[,-(1:5)]==3] <- 0.45
data.bleaching[,-(1:5)][data.bleaching.bckp[,-(1:5)]==4] <- 0.80


# Create back up files of disturbance tables ----
data.bleaching[is.na(data.bleaching)] <- 0
data.COTS[is.na(data.COTS)] <- 0
data.storms[is.na(data.storms)] <- 0

data.bleaching.bckp <- data.bleaching
data.COTS.bckp <- data.COTS
data.storms.bckp <- data.storms

# ----


# Load Fertilisation vs Density Models Estimates
setwd(BASE_DIRECTORY)
load("R_Objects/FvDParams.Rdata")

# Add percent Reef
setwd(ENVDATA_DIRECTORY)
PercentReef = read.csv("data.grid.csv", header=T)[,c(2,3,5)] %>% dplyr::inner_join(PopData[,2:3], by = c("lon", "lat"))
data.grid = dplyr::inner_join(data.grid, PercentReef, by = c("lon", "lat"))
PopData = data.grid[1:7]






