# COTSMod: Spatially Explicit COTS-Coral Metapopulation Model ----


# Load Functions and Modelling objects ----

rm(list=ls())
getwd()
load("RData/ModelWorkspace_2019-03-26.RData")
library(dplyr)
DIRECTORY = getwd()

masterDF = data.frame("SexRatio" = 5, 
                      "ConsRateW" = 150, 
                      "ConsRateS" = 250,
                      "avgPCF" = 10000000,
                      "sdPCF" = 1000000,
                      "mortJ1" =  0.99,
                      "mortJ2" = 0.95,
                      "mortA" = 0.6,
                      "remJ1" = 0,
                      "remJ2" = 0,
                      "remA" = 1,
                      "cssJ1" = 0.9803,
                      "cssJ2" = 0.0171,
                      "cssA" = 0.0026,
                      "Pred" = 0.98,
                      "p" = 0.25,
                      "Crash" = 3,
                      "OutbreakCrash" = 4)

NREPS = length(masterDF$OutbreakCrash)

# Results storage arrays (Pixel, Year/Season, Simulation)
res.cc = array(NA, dim=c(dim(data.grid)[1], nyears*2, nsimul))
# res.cots = array(NA, dim=c(dim(data.grid)[1], nyears*2, nsimul))

# Initialize Model ----

COTSfromCoralModel=T 
COTSfromSimul=F 
browse = FALSE 
inityear = 1995

COTSabund = initializeCOTSabund(data.grid, data.COTS, inityear, stagenames, COTS_StableStage, npops)
if(COTSfromSimul==F){
  COTSabund <- matrix(0,nrow=npops, ncol=3, dimnames = list(NULL, c("J_1", "J_2", "A")))
}

# Run Model ----  

  DIRECTORY = getwd()
  setwd(DIRECTORY)
  reps = 1
  # setwd("Results")
  # save(data.grid, file = "ReadWriteTest.Rdata")
  # save(data.grid, file = "Sample_1.Rdata")
  # setwd(DIRECTORY)
  # load(paste0("Rdata/ModelWorkspace_", Sys.Date(), ".Rdata"))
  # Replicate Loop
  SexRatio = masterDF[reps, "SexRatio"]
  ConsRate = as.vector(masterDF[reps, 2:3])
  PCFParams = c(masterDF[reps, "avgPCF"], masterDF[reps,"sdPCF"])
  COTSmort = as.numeric(masterDF[reps, c("mortJ1", "mortJ2", "mortA")])
  COTSremain = as.numeric(masterDF[reps, c("remJ1", "remJ2", "remA")])
  COTS_StableStage = as.numeric(masterDF[reps, c("cssJ1", "cssJ2", "cssA")])
  Pred = masterDF[reps,"Pred"]
  p = masterDF[reps,"p"]
  Crash = masterDF[reps,"Crash"]
  OutbreakCrash = masterDF[reps,"OutbreakCrash"]
  
  # Simulation loop
  for (j in 1:nsimul) {
    print(j)
    HC.1996 <- HCINI[,j]
    b0 <- B0[,j]
    b1 <- B0[,j]/log(HCMAX[,j])
    res.cc[,1,j] <- as.numeric(HC.1996)
    CoralCover = HC.1996
    # Year Loop
    for(i in 1:length(Years)){  
      print(i + 1995) # loop through years
      # Season Loop
      for(season in seasons){ 
        if(browse == TRUE) {
          browser()
        }
        # COTSabund = doPredPreyDynamics(COTSabund, CoralCover, p, Crash)
        # COTSabund = doCOTSDispersal(season,COTSabund,SexRatio,ConnMat, PCFParams, Pred, FvDParams) #Pruducing NAS
        # COTSabund = doCOTSDemography(season, COTSabund, COTSmort, COTSremain)
        # Consumption = doCoralConsumption(season, COTSabund, CoralCover, ConsRate) 
        # CoralCover = Consumption[,'CRemaining']
        # CoralConsum = round(Consumption[,'CChange'],4)
        CoralCover.Dist = doCoralDisturbance(i, j, season, CoralCover, 
                                             COTSfromCoralModel = COTSfromCoralModel, 
                                             storms.rsmpl, B.STORMS, WQ_Cyclone, 
                                             COTS.rsmpl, B.COTS, WQ_CoTS,
                                             bleaching.rsmpl, B.BLEACHING, WQ_bleach,
                                             disease.rsmpl, B.DISEASE, unknown.rsmpl, B.UNKNOWN)                                                    
        Disturbance = CoralCover.Dist - CoralCover
        CoralCover = CoralCover.Dist
        Growth = doCoralGrowth(season, CoralCover, b0, b1)
        CoralCover = Growth[,'CoralCover']
        CoralGrowth = round(Growth[,'CoralGrowth'],4)
        # Results[(Results$Year==i+1995) & (Results$Season==season),
        #         c("COTSJ1", "COTSJ2", "COTSA", "CoralCover", "CoralCover.DistLoss", "CoralCover.Consum", 'CoralCover.Growth')] =
        #   cbind(COTSabund, CoralCover, Disturbance, CoralConsum, CoralGrowth)
        # if(i>OutbreakCrash & season =="winter") {
        #   OutbreakCrasher = Results %>%
        #     dplyr::filter(Year > (i+1995-OutbreakCrash) & Year <= i+1995) %>%
        #     dplyr::group_by(REEF_NAME, Year) %>%
        #     dplyr::mutate(COTSA = (COTSA/100)*0.15) %>% # need to allow for detection
        #     dplyr::summarise(Mean.COTS = mean(COTSA)) %>%
        #     dplyr::mutate(is.outbreak = ifelse(Mean.COTS > 1, 1,0)) %>%
        #     dplyr::group_by(REEF_NAME) %>%
        #     dplyr::summarise(Crash = ifelse(sum(is.outbreak)==OutbreakCrash,1,0)) %>%
        #     dplyr::filter(Crash==1)
        #   matchit = match(data.grid$REEF_NAME,OutbreakCrasher$REEF_NAME)
        #   COTSabund[which(!is.na(matchit)),] = c(0,0,0)
        # Close outbreak crasher
        if (season=="summer") {
          res.cc[,2*i-1,j] = CoralCover
          # res.cots[,2*i-1,j] = COTSabund[,3]
        }
        res.cc[,2*i,j] = CoralCover
        # res.cots[,2*i,j] = COTSabund[,3]
      } # close season loop
    } # close Year loop
    
 }# Close simulation loop
  
  
  # Reef Level Summaries
  resCC.reef <- array(0, dim=c(length(unique(data.grid$REEF_ID)), nyears*2, nsimul))  
  resCC.cluster <- array(0, dim=c(length(unique(data.grid$bent.clust)), nyears*2, nsimul))
  # resCOTS.reef <- array(0, dim=c(length(unique(data.grid$REEF_ID)), nyears*2, nsimul))  
  # resCOTS.cluster <- array(0, dim=c(length(unique(data.grid$bent.clust)), nyears*2, nsimul))  
  for (i in 1:dim(res.cc)[3]) {
    resCC.reef[,,i] <- as.matrix(aggregate(res.cc[,,i], by=list(data.grid$REEF_ID), 
                                           FUN=mean, na.rm=T)[-1])
    resCC.cluster[,,i] <- as.matrix(aggregate(res.cc[,,i], by=list(data.grid$bent.clust), 
                                              FUN=mean, na.rm=T)[-1])
    # # resCOTS.reef[,,i] <- as.matrix(aggregate(res.cots[,,i], by=list(data.grid$REEF_ID), 
    #                                          FUN=mean, na.rm=T)[-1])
    # # resCOTS.cluster[,,i] <- as.matrix(aggregate(res.cots[,,i], by=list(data.grid$bent.clust), 
    #                                             FUN=mean, na.rm=T)[-1])
  }
  
  resCC.reef.mn <- apply(resCC.reef, c(1,2), mean, na.rm=T)
  resCC.reef.med <- apply(resCC.reef, c(1,2), median, na.rm=T)
  # resCC.reef.min <- apply(resCC.reef, c(1,2), quantile, probs=0.05, na.rm=T)
  # resCC.reef.max <- apply(resCC.reef, c(1,2), quantile, probs=0.95, na.rm=T)
  resCC.reef.25 <- apply(resCC.reef, c(1,2), quantile, probs=0.25, na.rm=T)
  resCC.reef.75 <- apply(resCC.reef, c(1,2), quantile, probs=0.75, na.rm=T)
  
  # resCOTS.reef.mn <- apply(resCOTS.reef, c(1,2), mean, na.rm=T)
  # resCOTS.reef.med <- apply(resCOTS.reef, c(1,2), median, na.rm=T)
  # resCOTS.reef.min <- apply(resCOTS.reef, c(1,2), quantile, probs=0.05, na.rm=T)
  # resCOTS.reef.max <- apply(resCOTS.reef, c(1,2), quantile, probs=0.95, na.rm=T)
  # resCOTS.reef.25 <- apply(resCOTS.reef, c(1,2), quantile, probs=0.25, na.rm=T)
  # resCOTS.reef.75 <- apply(resCOTS.reef, c(1,2), quantile, probs=0.75, na.rm=T)
  
  nReefs = length(unique(data.grid$REEF_ID))
  # Results df for Dashboard -- Reef Level
  ResultsDash = data.frame(sapply(unique(data.grid[4:5]), rep, times=nyears*NSEASONS),
                           Year=rep(Years,each=2*nReefs), 
                           Season=rep(c("summer", "winter"),each=nReefs), 
                           # COTS.mn=(as.vector(resCOTS.reef.mn)/100)*15, 
                           # COTS.Q50=(as.vector(resCOTS.reef.med)/100)*15, 
                           # # COTS.Q05=(as.vector(resCOTS.reef.min)/100)*15, 
                           # # COTS.Q95=(as.vector(resCOTS.reef.max)/100)*15, 
                           # COTS.Q25=(as.vector(resCOTS.reef.25)/100)*15, 
                           # COTS.Q75=(as.vector(resCOTS.reef.75)/100)*15,
                           CC.mn=as.vector(resCC.reef.mn), 
                           CC.Q50=as.vector(resCC.reef.med), 
                           # CC.Q05=as.vector(resCC.reef.min), 
                           # CC.Q95=as.vector(resCC.reef.max), 
                           CC.Q25=as.vector(resCC.reef.25), 
                           CC.Q75=as.vector(resCC.reef.75))
  
  setwd(DIRECTORY)
  setwd("Results")
  save(res.cc, res.cots, ResultsDash, bleaching.mn, storms.mn, disease.mn, unknown.mn, COTS.mn, file = "NoCOTS.Rdata") 
  

write.csv(ResultsDash, "ForDashboard_NoCOTS.csv")
















