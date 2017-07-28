########################
# CoTS Model -- AIMS LTMP Interpolation
########################
#   Authors: Sam Matthews
#   
#   28 September 2016 - Started scripting
#
#   This script interpolates The purpose of this script is to interpolate CoTS presence/abundance based 
#   on LTMP data for each year of collected data. Here I use the inverse distance weighted (IDW) 
#   interpolation method from the *gstat* package. The interpolation uses a fixed radius of 0.1 degree. 
#   Any cell in the grid which contained less than three points was assigned a ‘no data’ value.
#
#   The interpolation will only be applied to pixels that contain reef, and abundace values will be 
#   multiplied by the percentage of reef conatined within each cell.  
#
#
#   THINGS TO DO:
#         1. 


##################
## 1. Load Data 
##################

setwd(DATA_DIRECTORY)
manta <- read.csv("ltmp/manta.csv") #LTMP Manta Tow Data
#PopData <- read.table("PopData.txt", sep = "\t") #our environmental data grid 
PopData <- read.table("XYZGrid.txt", header = T, sep = "\t")

#we only need the ID variables for this process
#PopData <- PopData[,-c(6:38)]
colnames(PopData)[1:2] <- c('x','y')
PopData <- merge(PopData, RM.Sites.round1, by=c('x','y'),all=T)
##################
## 2. Summarise LTMP Data
##################

colnames(manta)[6] <- "LONG"
manta_year <- ddply(manta, c("REEF_NAME", "LONG", "LAT", "REPORT_YEAR"), 
                    summarise, mean_COTS = mean(COT_COUNT)) %>% 
                    dcast(REEF_NAME + LAT + LONG ~ REPORT_YEAR, value.var="mean_COTS")

##################
## 3. Convert to Spatial Objects for use in the gstat package
##################

coordinates(PopData) = ~ x + y
gridded(PopData) = TRUE

#for interpolation, each year has to have complete cases so I will have to write a function that 
#cuts up the data frame, removes incomplete cases and then passes it to idw function and eventually 
#spits out an interpolated data frame or Spatial Pixels Data Frame

##################
## 4. Interpolation Functions
##################

idw.interpolation <- function(Year, XYZ, maxdist, nmin) {
  #first we need to subset the manta_year data set
  Obs = na.omit(data.frame(manta_year[,2:3],manta_year[,Year-1983 +4])) #create observations of CoTS for year i
  colnames(Obs)[3] = "mean.COTS"
  coordinates(Obs) = ~LONG+LAT 
  COTS.idw = idw(mean.COTS~1, Obs, XYZ, maxdist = maxdist, nmin= nmin)
  #results = spplot(COTS.idw["var1.pred"], main = "CoTS interpolation")
  return(COTS.idw)
}

idw.interpolation.df <- function(Year, XYZ, maxdist, nmin) { ## returns an interpolated Data Frame
  #first we need to subset the manta_year data set
  Obs = na.omit(data.frame(manta_year[,2:3],manta_year[,Year-1983 +4])) #create observations of CoTS for year i
  colnames(Obs)[3] = "mean.COTS"
  coordinates(Obs) = ~LONG+LAT 
  COTS.idw = as.data.frame(idw(mean.COTS~1, Obs, XYZ, maxdist = maxdist, nmin= nmin))
  #results = spplot(COTS.idw["var1.pred"], main = "CoTS interpolation")
  return(COTS.idw)
}

##################
## 5. Apply Interpolations
##################

Years <- 1983:2015
COTS.interpolation.df <- lapply(Years, FUN = idw.interpolation.df, XYZ=PopData, maxdist=1, nmin=3)
names(COTS.interpolation.df) <- Years #name the list so they can be retrieved by year

##################
## 5. COmbine Interpolations into one data fram
##################

for(i in 1:length(COTS.interpolation.df)) names(COTS.interpolation.df[[i]])[3] <- paste0("X",1982+i) 
for(i in 1:length(COTS.interpolation.df)) COTS.interpolation.df[[i]][4] <- NULL
COTSInterp <- Reduce(function(x,y)merge(x,y, by=c("x", "y")), COTS.interpolation.df)  

#Add back in the PIXEL_ID, REEF_ID and Percetn reef
#COTSInterp <- left_join(as.data.frame(PopData[,1:2]), COTSInterp, by=c("x", "y"))

#rearragne order of columns
#COTSInterp <- COTSInterp[,c(1,4:5,2:3,6:38)]

#write file
setwd(DATA_DIRECTORY)
write.table(COTSInterp, file = "COTS_Interpolated.txt", row.names = F, sep="\t")

#############
# END SCRIPT
#############


