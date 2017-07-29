###############################
# RUN COTS MODEL
###############################

####################!
# 1 SET UP WORKSPACE ----
####################!


set.seed(1)
setwd(DATA_DIRECTORY)
PopData <- read.csv("reefs.csv", header = TRUE)[1:10,] # we will just use a subset for testing
COTS.data <- read.csv("Disturbance/CoTS_data.csv", header = TRUE)[1:10,]
setwd(ENVDATA_DIRECTORY)
data.grid = read.csv("data.grid.csv", header=T)[1:10,]
npops = 10

####################!
# 2 RUN MODEL ----
####################!

COTSabund = initializeCOTSabund(PopData = PopData, COTS.init = COTS.data, 
                                Year=1996, stagenames, COTS_StableStage = COTS_StableStage)   # initialize the COTS abundance object (for year 0) 
CoralCoverParams = intializeCoralCoverParams(data.grid = data.grid, nsims=10)

CoralCover = CoralCoverParams$HCINI[,1]

# Probably Change Storage to an array
Results = data.frame(sapply(PopData, rep.int, times=NYEARS*NSEASONS),
                     Year=rep(1996:2015,each=2*npops), Season=rep(c("summer", "winter"),each=npops), 
                     COTSJ1=NA, COTSJ2=NA, COTSA=NA, CoralCover=NA, DistCOTS=NA, DistCYCL=NA, DistBLCH=NA)
ConnMat=ConnMat[1:10,1:10]

for(year in 1996:2015){                  # loop through years
  for(season in SEASONS){               # loop through seasons
    
    COTSabund = doCOTSDispersal(season,COTSabund,SR,ConnMat)
    COTSabund = doCOTSDemography(season,COTSabund = COTSabund)
    CoralCover = doCoralConsumption(season,  COTSabund, CoralCover)
    #CoralCover = doCoralGrowth(season, CoralCover...)
    
    #CoralCover = doCoralDisturbance(season,CoralCover,...)           # coral disturbance processes, including from COTS
    Results[(Results$Year==year) & (Results$Season==season),
        c("COTSJ1", "COTSJ2", "COTSA", "CoralCover")] = cbind(COTSabund, CoralCover)
    
  }
}

test = Results %>% filter(PIXEL_ID==1) #something is wrong with the order of my operations start here tomorrow

