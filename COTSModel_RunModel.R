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
WQ <- data.grid$Primary + data.grid$Secondary + data.grid$Tertiary
# Gompertz parameters: disturbance & WQ effect sizes (as per MacNeil et al. in prep)
bleaching.mn.sd <- c(-0.19, 0.01)
COTS.mn.sd <- c(-0.54, 0.04)
disease.mn.sd <- c(-0.13, 0.01)
storms.mn.sd <- c(-0.64, 0.01)
unknown.mn.sd <- c(-0.16, 0.01)
WQ.mn.sd <- c(-0.68, 0.03)
####################!
# 2 RUN MODEL ----
####################!

COTSabund = initializeCOTSabund(PopData = PopData, COTS.init = COTS.data, 
                                Year=1996, stagenames, COTS_StableStage = COTS_StableStage)   # initialize the COTS abundance object (for year 0) 
CoralCoverParams = intializeCoralCoverParams(data.grid = data.grid, nsims=10)

CoralCover = CoralCoverParams$HCINI[,1][1:10]
B0=CoralCoverParams$B0[,1]

# Probably Change Storage to an array
Results = data.frame(sapply(PopData, rep.int, times=NYEARS*NSEASONS),
                     Year=rep(1996:2015,each=2*npops), Season=rep(c("summer", "winter"),each=npops), 
                     COTSJ1=NA, COTSJ2=NA, COTSA=NA, CoralCover=NA, DistCOTS=NA, DistCYCL=NA, DistBLCH=NA)
ConnMat=ConnMat[1:10,1:10]

for(year in 1996:2015){                  # loop through years
  for(season in SEASONS){               # loop through seasons
    COTSabund = doCOTSDispersal(season,COTSabund,SR,ConnMat)
    COTSabund = doCOTSDemography(season,COTSabund = COTSabund)
    CoralCover = doCoralConsumption(year, season,  COTSabund, CoralCover)
    CoralCover = doCoralGrowth(CoralCover, B0, WQ)
    #CoralCover = doCoralDisturbance(season,CoralCover,...)           # coral disturbance processes, including from COTS
    Results[(Results$Year==year) & (Results$Season==season),
        c("COTSJ1", "COTSJ2", "COTSA", "CoralCover")] = cbind(COTSabund, CoralCover)
    
  }
}

test = Results %>% filter(PIXEL_ID==1) #something is wrong with the order of my operations start here tomorrow

