########################
# Coral Model - Disturbance Layers
########################
#   Authors: Sam Matthews
#   
#   30 September 2016 - Started scripting
#
#   This script extracts distrurbance values (Interpolated: COTS, Bleaching; Modelled: Cyclones) for each 
#   LTMP reef from the AIMS reef monitoring program for each year
#
#
#   


##################
## 1. Load Data 
##################

source("C:/Users/jc312264/Documents/GitHub/COTS_popmodel/COTSModel_PrepareWorkspace.R")

setwd(DATA_DIRECTORY)
COTS <- read.table("COTS_Interpolated.txt", header=T, sep="\t")
Cyclones <- read.table("Cyclones.txt", header=T, sep="\t")
Bleaching <- read.table("Bleaching_Interpolated.txt", header=T, sep="\t")
LTMP <- read.csv("ltmp/benthos.csv", header=T)
##################
## 2. Write Out Standard grids for Disturbances
##################

write.table(Bleaching, "Bleaching.txt", row.names = F, sep = "\t")
write.table(COTS, "COTS.txt", row.names = F, sep = "\t")
write.table(Cyclones, "Cyclones.txt", row.names = F, sep = "\t")

##################
## 2. Assign Reef Monitoring Sites to closest grid points
##################


RM.Sites <- unique(LTMP[,c(1,7:8)])
colnames(RM.Sites)[2:3] <- c('y', 'x')
RM.Sites.round <- data.frame(REEF_NAME = RM.Sites[,1], x = round(RM.Sites[,3],2), y = round(RM.Sites[,2],2))
RM.Sites.round1 <- unique(RM.Sites.round[,1:3])

COTS.RM <- merge(COTS, RM.Sites.round1, by=c('x','y'), all.y = T)
Bleaching.RM <- merge(Bleaching, RM.Sites.round1, by=c('x','y'), all.y = T)
Cyclones.RM <- merge(Cyclones, RM.Sites.round1, by=c('x','y'), all.y = T)
COTS.RM <- COTS.RM[,c(36,1:35)]
Bleaching.RM <- Bleaching.RM[,c(5,1:4)]
Cyclones.RM <- Cyclones.RM[,c(34,1:33)]

setwd(DATA_DIRECTORY)
write.table(Bleaching.RM, "BleachingRM.txt", row.names = F, sep = "\t")
write.table(COTS.RM, "COTSRM.txt", row.names = F, sep = "\t")
write.table(Cyclones.RM, "CyclonesRM.txt", row.names = F, sep = "\t")

##### Need to add sites not included in original interpolations


# Now match to XYZ Grid





#########
## END SCRIPT
###########