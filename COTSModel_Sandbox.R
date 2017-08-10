########################
# COTSModel_Sandbox - Rough Working Out etc
########################
#   Authors: Kevin Shoemaker, Sam Matthews
rm(list=ls())
source("C://Users//jc312264//Documents//GitHub//COTS_Model//COTSModel_PrepareWorkspace.R")
# dirs = list(BASE=BASE_DIRECTORY, CODE=CODE_DIRECTORY, RESULTS=RESULTS_DIRECTORY, 
# DATA=DATA_DIRECTORY, ENVDATA=ENVDATA_DIRECTORY,SPATIALDATA=SPATIALDATA_DIRECTORY)
setwd(BASE_DIRECTORY)
load(file = "R_Objects/ConnMat.Rdata")
setwd(CODE_DIRECTORY)
source("COTSModel_LoadSmallObjectsForModelling.R")
setwd(CODE_DIRECTORY)
source("COTSModel_CoralFunctions.R")
source("COTSModel_COTSFunctions.R")
ConnMat=ConnMat
SEASONS=SEASONS

NREPS=10

masterDF = MakeLHSSamples(NREPS = 10)

rm(list = ls())
source("C://Users//jc312264//Documents//GitHub//COTS_Model//COTSModel_PrepareWorkspace.R")
# dirs = list(BASE=BASE_DIRECTORY, CODE=CODE_DIRECTORY, RESULTS=RESULTS_DIRECTORY, 
# DATA=DATA_DIRECTORY, ENVDATA=ENVDATA_DIRECTORY,SPATIALDATA=SPATIALDATA_DIRECTORY)
setwd(BASE_DIRECTORY)
load(file = "R_Objects/ConnMat.Rdata")
setwd(CODE_DIRECTORY)
source("COTSModel_LoadSmallObjectsForModelling.R")
setwd(CODE_DIRECTORY)
source("COTSModel_CoralFunctions.R")
source("COTSModel_COTSFunctions.R")
ConnMat=ConnMat

runModel = function(masterDF, PopData, COTS.data, data.grid, ConnMat, npops, rep, seasons) {
  SexRatio = masterDF[1, "SexRatio"]
  ConsRateW = masterDF[1, "ConsRateW"]
  ConsRateS = masterDF[1, "ConsRateS"]
  PCFParams = c(masterDF[1, "avgPCF"], masterDF[1,"sdPCF"])
  # avgPCF = masterDF[1, "avgPCF"]
  # sdPCF = masterDF[1, "sdPCF"]
  COTSmort = as.numeric(masterDF[1, c("mortJ1", "mortJ2", "mortA")])
  COTSremain = as.numeric(masterDF[1, c("remJ1", "remJ2", "remA")])
  COTS_StableStage = as.numeric(masterDF[1, c("cssJ1", "cssJ2", "cssA")])
  # avgAdultSize =
  # sdAdultSize = # These will change the fecundity estimatesC
  # need an Allee Effect
  # need to make stable stage vary by a scaling factor
  # make mortality and remain resource driven
  
  # Initialize
  PopData = PopData[1:npops, ]
  COTS.data = COTS.data[1:npops, ]
  data.grid = data.grid[1:npops, ]
  # CoralCoverParams = intializeCoralCoverParams(data.grid = data.grid, nsims=10, npops=npops)
  # CoralCover = CoralCoverParams$HCINI[,1][1:npops]
  # B0=CoralCoverParams$B0[,1][1:npops]
  # HC.asym = CoralCoverParams$HCMAX[,1][1:npops]
  ConnMat=ConnMat
  CoralCover=data.grid$pred.HCini.mean[1:npops]
  B0=data.grid$pred.b0.mean[1:npops]
  HC.asym=data.grid$pred.HCmax.mean[1:npops]
  WQ <- data.grid$Primary + data.grid$Secondary + data.grid$Tertiary
  #PCFParams = COTSPCF(npops, SexRatio = 5)
  K = setCarryingCapacity(npops)
  print(length(K$MinK.10A))
  COTSabund = initializeCOTSabund(PopData, COTS.data, 1996, stagenames, COTS_StableStage, npops)  # initialize the COTS abundance object (for year 0) 
  print(length(COTSabund[,3]))
  Results = initializeModel(PopData, COTSabund, CoralCover, SexRatio, 
                            ConsRateS, ConsRateW, B0, WQ, HC.asym, PCFParams, npops, ConnMat)
  
  
  # year Loop
  for(year in 1997:2015){                  # loop through years
    for(season in SEASONS){               # loop through seasons
      COTSabund = doCOTSDispersal(season,COTSabund,SexRatio,ConnMat, PCFParams, npops)
      COTSabund = doCOTSDemography(season, COTSabund, COTSmort, COTSremain)
      COTSabund = doPredPreyDynamics(season, year, COTSabund, Results,K)
      CoralCover = doCoralConsumption(year, season,  COTSabund, CoralCover, ConsRateS, ConsRateW)
      CoralCover = doCoralGrowth(CoralCover, B0, WQ, HC.asym)
      #CoralCover = doCoralDisturbance(season,CoralCover,...)           # coral disturbance processes, including from COTS
      Results[(Results$Year==year) & (Results$Season==season),
              c("COTSJ1", "COTSJ2", "COTSA", "CoralCover")] = cbind(COTSabund, CoralCover)
    }
  }
  setwd(RESULTS_DIRECTORY)
  name <- sprintf("Sample_%s.Rdata",rep)
  save(Results, file=name)
  
}

for (i in 1:20){
  source("C://Users//jc312264//Documents//GitHub//COTS_Model//COTSModel_PrepareWorkspace.R")
  masterDF = MakeLHSSamples(NREPS = 20)
  temp = MakeWorker(masterDF,npops=npops,seasons=SEASONS)(i)
  rm(list=ls())
}
