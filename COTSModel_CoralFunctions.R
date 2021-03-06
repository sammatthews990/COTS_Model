##########################
#  Functions for modeling coral dynamics in the great barrier reef as part of a framework for modeling COTS outbreaks  
#
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
# 
#  05 April 2015 -- started scripting
#  28 September 2016 -- intialize coral cover function
#  

  ## for coral dispersal: see kinlan and Gaines paper


#################!
# initializeCoralCover ----
#################!
# OBJECTIVE:
#    generate an matrix for storing the coral cover in each pixel. 
# PARAMS:
#    - reefmap: raster template for the study region: NA outside of reefs, reef ID value within reefs 
#    - initCoralCover: for every pixel in the study area, a vector of initial coral cover
#
# RETURNS:
#    - CoralCover: spatially-structured Coral Cover
#           CoralCover[,'F']: vector representing spatially structured cover for fast growing corals
#           CoralCover[,'S']: vector representing spatially structured abundance of Juvenile stage 2 individuals
#           COTSabund$A: vector representing spatially structured abundance of reproductive adult individuals
#           COTSabund$S: vector representing spatially structured abundance of senile adult individuals
#           NOTE: larvae are not considered explicitly here. 
###################


intializeCoralCoverParams = function(data.grid, nsims, npops){
  #browser()
  data.grid=data.grid[1:npops,]
  WQ <- data.grid$Primary + data.grid$Secondary + data.grid$Tertiary
  N = nsims
  HCINI <- HC.asym <- B0 <- matrix(NA, ncol = N, nrow = dim(data.grid)[1])
  #browser()
  for (j in 1:dim(data.grid)[1]) {
    #Define sigma (i.e. variance-covariance matrix) for ith grid cell 
    sigma <- matrix(c(data.grid$pred.HCini.sd[j]^2, 0.54*data.grid$pred.HCini.sd[j]*data.grid$pred.HCmax.sd[j], 0.14*data.grid$pred.HCini.sd[j]*data.grid$pred.b0.sd[j],
                      0.54*data.grid$pred.HCini.sd[j]*data.grid$pred.HCmax.sd[j], data.grid$pred.HCmax.sd[j]^2, -0.13*data.grid$pred.HCmax.sd[j]*data.grid$pred.b0.sd[j], 
                      0.14*data.grid$pred.HCini.sd[j]*data.grid$pred.b0.sd[j], -0.13*data.grid$pred.HCmax.sd[j]*data.grid$pred.b0.sd[j], data.grid$pred.b0.sd[j]^2),
                    3,3)
    
    #Pick N random parameters for jth grid cell
    #something going wrong when picking values
    pick <- MASS::mvrnorm(n=N, mu=c(data.grid$pred.HCini.mean[j], data.grid$pred.HCmax.mean[j], data.grid$pred.b0.mean[j]), Sigma = sigma)
    HCINI[j,] <- pick[,1]
    HC.asym[j,] <- pick[,2]
    B0[j,] <- pick[,3]
  }
  return(list(WQ=cbind(data.grid[,1:5], WQ), 
              HCINI=HCINI, HC.asym=HC.asym, B0=B0))
}

# CoralCoverParams = intializeCoralCoverParams(data.grid = data.grid, nsims=10, npops = npops)

#################!
# doCoralDistrurbances ----
#################!  



#################!
# doCoralGrowth ----
#################!  
# OBJECTIVE:
#    Allow Coral to grow follwinf both disturbance and COTS consumption 
# PARAMS:
#    - CoralCover: Coral cover from the previous timestep
#    - b0 Intrinsic growth parameter from Gompertz Model
#    - b1 Asymptotic growth parameter from Gompertz Model
#    - WQ Water quality parameter (not sure whether to predetermine this parameter or allow it to
#       change for the sensitivity analyses)
#    - HCMAX Estimated carrying capacity of coral cover
# RETURNS:
#    - CoralCover: spatially-structured Coral Cover
###################!

doCoralGrowth = function(CoralCover, B0, WQ, HC.asym) {
  b0.wq <- B0 + WQ * rnorm(length(WQ), mean=WQ.mn.sd[1], sd=WQ.mn.sd[2])
  b1.wq <- b0.wq / log(HC.asym)
  CoralCover <- log(CoralCover)
  CoralCover.t1 <- b0.wq + (1 - b1.wq)*CoralCover
  return(cbind(CoralCover=exp(CoralCover.t1), CoralGrowth=(exp(CoralCover.t1)-exp(CoralCover))))
}

#########################



####################!
# 3 INITIALIZE MODEL ----
####################!

initializeModel = function(PopData, data.grid, COTSabund,CoralCover, SexRatio, ConsRateS, 
                           ConsRateW, B0, WQ, HC.asym, PCFParams, npops, ConnMat, FvDParams){
  
  # Probably Change Storage to an array
  
  Results = data.frame(sapply(PopData, rep.int, times=NYEARS*NSEASONS),
                       Year=rep(1996:2015,each=2*npops), Season=rep(c("summer", "winter"),each=npops), 
                       COTSJ1=NA, COTSJ2=NA, COTSA=NA, CoralCover=NA, DistCOTS=NA, DistCYCL=NA, DistBLCH=NA)
  Results$CoralCover.Consum = NA
  Results$CoralCover.Growth = NA
  CoralCover=CoralCover
  for(year in 1996){                  # loop through years
    for(season in seasons){               # loop through seasons
      # browser()
      COTSabund = doCOTSDispersal(season,COTSabund,SexRatio,ConnMat, PCFParams, npops, FvDParams)
      COTSabund = doCOTSDemography(season, COTSabund, COTSmort, COTSremain)
      #COTSabund = doPredPreyDynamics(season, year, COTSabund,Reults, K)
      Consumption = doCoralConsumption(year, data.grid, season, COTSabund, CoralCover, ConsRateS, ConsRateW)
      CoralCover = Consumption[,'CRemaining']
      CoralConsum = Consumption['CChange']
      Growth = doCoralGrowth(CoralCover, B0, WQ, HC.asym)
      CoralCover = Growth[,'CoralCover']
      CoralGrowth = Growth[,'CoralGrowth']
      #CoralCover = doCoralDisturbance(season,CoralCover,...)           # coral disturbance processes, including from COTS
      Results[(Results$Year==year) & (Results$Season==season),
              c("COTSJ1", "COTSJ2", "COTSA", "CoralCover", "CoralCover.Consum", 'CoralCover.Growth')] = 
        cbind(COTSabund, CoralCover,CoralConsum, CoralGrowth)
    }
  }
  return(Results)
}


####################!
# 4 RUN MODEL ----
####################!

runModel = function(masterDF, PopData, COTS.data, Years = Years, data.grid, rep, Pred) {
  # browser()
  SexRatio = masterDF[rep, "SexRatio"]
  ConsRateW = masterDF[rep, "ConsRateW"]
  ConsRateS = masterDF[rep, "ConsRateS"]
  PCFParams = c(masterDF[rep, "avgPCF"], masterDF[rep,"sdPCF"])
  # avgPCF = masterDF[1, "avgPCF"]
  # sdPCF = masterDF[1, "sdPCF"]
  COTSmort = as.numeric(masterDF[rep, c("mortJ1", "mortJ2", "mortA")])
  COTSremain = as.numeric(masterDF[rep, c("remJ1", "remJ2", "remA")])
  COTS_StableStage = as.numeric(masterDF[rep, c("cssJ1", "cssJ2", "cssA")])
  # avgAdultSize =
  # sdAdultSize = # These will change the fecundity estimatesC
  # need an Allee Effect
  # need to make stable stage vary by a scaling factor
  # make mortality and remain resource driven
  
  #### FOR SOME REASON NONE OF THESE PARAMETERS ARE AVAILIABLE INSIDE THE FUNCTION
  
  # Initialize
  npops=npops
  seasons=seasons
  PopData = PopData[1:npops, ]
  COTS.data = COTS.data[1:npops, ]
  data.grid = data.grid[1:npops, ]
  # CoralCoverParams = intializeCoralCoverParams(data.grid = data.grid, nsims=10, npops=npops)
  # CoralCover = CoralCoverParams$HCINI[,1][1:npops]
  # B0=CoralCoverParams$B0[,1][1:npops]
  # HC.asym = CoralCoverParams$HCMAX[,1][1:npops]
  
  # Work out which reefs from our connectivity matrix are to be included
  # which reefs from npops are being used in the analysis
  whichreefs = unique(data.grid$REEF_NAME[1:npops])
  ConnMat = COTS.ConnMat[1:length(whichreefs), 1:length(whichreefs)]
  Pixels = Pixels[1:length(colnames(ConnMat)),]
  FvDParams=FvDParams
  CoralCover=data.grid$pred.HCini.mean[1:npops]
  B0=data.grid$pred.b0.mean[1:npops]
  HC.asym=data.grid$pred.HCmax.mean[1:npops]
  WQ <- data.grid$Primary + data.grid$Secondary + data.grid$Tertiary
  #PCFParams = COTSPCF(npops, SexRatio = 5)
  K = setCarryingCapacity(npops)
  print(length(K$MinK.10A))
  COTSabund = initializeCOTSabund(PopData, COTS.data, 1996, stagenames, COTS_StableStage, npops)  # initialize the COTS abundance object (for year 0) 
  print(length(COTSabund[,3]))
  # Results = initializeModel(PopData, data.grid, COTSabund, CoralCover=CoralCover, SexRatio, 
  #                           ConsRateS, ConsRateW, B0, WQ, HC.asym, PCFParams, npops, ConnMat, FvDParams)
  Results = data.frame(sapply(PopData[1:4], rep, times=NYEARS*NSEASONS),
                       sapply(PopData[5:7], rep, times=NYEARS*NSEASONS),
                       Year=rep(1996:2015,each=2*npops), Season=rep(c("summer", "winter"),each=npops), 
                       COTSJ1=NA, COTSJ2=NA, COTSA=NA, CoralCover=NA, DistCOTS=NA, DistCYCL=NA, DistBLCH=NA)
  Results$CoralCover.Consum = NA
  Results$CoralCover.Growth = NA
  
  # browser()
  # year Loop
  for(year in Years){  
    print(year)# loop through years
    for(season in seasons){             # loop through seasons
      # browser()
      COTSabund = doCOTSDispersal(season,COTSabund,SexRatio,ConnMat, PCFParams, Pred, FvDParams)
      COTSabund = doCOTSDemography(season, COTSabund, COTSmort, COTSremain)
      COTSabund = doPredPreyDynamics(season, year, COTSabund, Results,K, CoralCover, Crash)
      Consumption = doCoralConsumption(year, data.grid, season, COTSabund, CoralCover, ConsRateS, ConsRateW)
      CoralCover = Consumption[,'CRemaining']
      CoralConsum = round(Consumption[,'CChange'],4)
      Growth = doCoralGrowth(CoralCover, B0, WQ, HC.asym)
      CoralCover = Growth[,'CoralCover']
      CoralGrowth = round(Growth[,'CoralGrowth'],4)
      #CoralCover = doCoralDisturbance(season,CoralCover,...)           # coral disturbance processes, including from COTS
      # browser()
      Results[(Results$Year==year) & (Results$Season==season),
              c("COTSJ1", "COTSJ2", "COTSA", "CoralCover", "CoralCover.Consum", 'CoralCover.Growth')] = 
        cbind(COTSabund, CoralCover,CoralConsum, CoralGrowth)
    }
  }
  setwd(RESULTS_DIRECTORY)
  name <- sprintf("Sample_%s.Rdata",rep)
  save(Results, file=name)
}




