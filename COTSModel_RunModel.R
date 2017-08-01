###############################
# RUN COTS MODEL
###############################

####################!
# 1 SET UP WORKSPACE ----
####################!


set.seed(1)
setwd(DATA_DIRECTORY)
npops = 10
PopData <- read.csv("reefs.csv", header = TRUE)[1:npops,] # we will just use a subset for testing
COTS.data <- read.csv("Disturbance/CoTS_data.csv", header = TRUE)[1:npops,]
ConnMat = load ###
setwd(ENVDATA_DIRECTORY)
data.grid = read.csv("data.grid.csv", header=T)[1:10,]
SexRatio =5
WQ <- data.grid$Primary + data.grid$Secondary + data.grid$Tertiary
NYEARS = 20
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

initializeModel = function(){

# Probably Change Storage to an array
Results = data.frame(sapply(PopData, rep.int, times=NYEARS*NSEASONS),
                     Year=rep(1996:2015,each=2*npops), Season=rep(c("summer", "winter"),each=npops), 
                     COTSJ1=NA, COTSJ2=NA, COTSA=NA, CoralCover=NA, DistCOTS=NA, DistCYCL=NA, DistBLCH=NA)

for(year in 1996){                  # loop through years
  for(season in SEASONS){               # loop through seasons
    COTSabund = doCOTSDispersal(season,COTSabund,SexRatio,ConnMat, avgPCF, sdPCF)
    COTSabund = doCOTSDemography(season, COTSabund, COTSmort, COTSremain)
    #COTSabund = doPredPreyDynamics(season, year, COTSabund)
    CoralCover = doCoralConsumption(year, season,  COTSabund, CoralCover, ConsRateS, ConsRateW)
    CoralCover = doCoralGrowth(CoralCover, B0, WQ)
    #CoralCover = doCoralDisturbance(season,CoralCover,...)           # coral disturbance processes, including from COTS
    Results[(Results$Year==year) & (Results$Season==season),
            c("COTSJ1", "COTSJ2", "COTSA", "CoralCover")] = cbind(COTSabund, CoralCover)
  }
}
return(Results)
}


masterDF = data.frame(SexRatio=5, ConsRateW=150,ConsRateS=250, avgPCF=avgPCF, sdPCF=sdPCF,
                    mortJ1=NA, mortJ2=NA, mortA=NA, remJ1=NA, remJ2=NA, remA=NA, cssJ1=NA, cssJ2=NA, cssA=NA)
masterDF[,6:14] = c(COTSmort, COTSremain, COTS_StableStage)

runModel = function(masterDF, PopData) {
  res = array(NA, dim=c(npops*NYEARS*NSEASONS, 4, nrow(masterDF)))
  for (i in 1:nrow(masterDF)) {
    SexRatio = masterDF[i, "SexRatio"]
    ConsRateW = masterDF[i, "ConsRateW"]
    ConsRateS = masterDF[i, "ConsRateS"]
    avgPCF = masterDF[i, "avgPCF"]
    sdPCF = masterDF[i, "sdPCF"]
    COTSmort = as.numeric(masterDF[i, c("mortJ1", "mortJ2", "mortA")])
    COTSremain = as.numeric(masterDF[i, c("remJ1", "remJ2", "remA")])
    COTS_StableStage = as.numeric(masterDF[i, c("cssJ1", "cssJ2", "cssA")])
    
    # Initialize
    CoralCoverParams = intializeCoralCoverParams(data.grid = data.grid, nsims=10)
    CoralCover = CoralCoverParams$HCINI[,1][1:10]
    B0=CoralCoverParams$B0[,1]
    COTSabund = initializeCOTSabund(PopData, COTS.data, 1996, stagenames, COTS_StableStage)   # initialize the COTS abundance object (for year 0) 
    Results = initializeModel()
    
    # year Loop
    for(year in 1997:2015){                  # loop through years
      for(season in SEASONS){               # loop through seasons
        COTSabund = doCOTSDispersal(season,COTSabund,SexRatio,ConnMat, avgPCF, sdPCF)
        COTSabund = doCOTSDemography(season, COTSabund, COTSmort, COTSremain)
        COTSabund = doPredPreyDynamics(season, year, COTSabund)
        CoralCover = doCoralConsumption(year, season,  COTSabund, CoralCover, ConsRateS, ConsRateW)
        CoralCover = doCoralGrowth(CoralCover, B0, WQ)
        #CoralCover = doCoralDisturbance(season,CoralCover,...)           # coral disturbance processes, including from COTS
        Results[(Results$Year==year) & (Results$Season==season),
                c("COTSJ1", "COTSJ2", "COTSA", "CoralCover")] = cbind(COTSabund, CoralCover)
      }
    }
    res[,,i] = Results[,c("COTSJ1", "COTSJ2", "COTSA", "CoralCover")] ##### PIck up Here, not sure why it wont load
  }
  list(Results[,1:5], res)
}

test=runModel(masterDF, PopData)

test = Results %>% filter(PIXEL_ID==1) #something is wrong with the order of my operations start here tomorrow

