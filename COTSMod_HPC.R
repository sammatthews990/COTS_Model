# COTSMod: Spatially Explicit COTS-Coral Metapopulation Model ----


# Load Functions and Modelling objects ----

rm(list=ls())

library(dplyr)
library(foreach)
library(doParallel)
DIRECTORY = getwd()

PRELOAD = F

if (PRELOAD == T) {
  load("RData/ModelWorkspace_2019-03-29_1.RData") # Choose an older workspace with pre run disturbance samples
  # Write Parameter Data Frame And Results Containers ----
  # Write Parameter Data Frame And Results Containers ----
  # Write Parameter Data Frame And Results Containers ----
  masterDF = data.frame("SexRatio" = 5, 
                        "ConsRateW" = 150, 
                        "ConsRateS" = 250,
                        "avgPCF" = 20000000,
                        "sdPCF" = 10000000,
                        "mortJ1" =  0.99,
                        "mortJ2" = 0.95,
                        "mortA" = 0.9,
                        "remJ1" = 0,
                        "remJ2" = 0,
                        "remA" = 1,
                        "cssJ1" = 0.9803,
                        "cssJ2" = 0.0171,
                        "cssA" = 0.0026,
                        "Pred" = 0.99,
                        "p" = 0.2,
                        "Crash" = 0,
                        "OutbreakCrash" = Inf,
                        "Fbase" = seq(0.05, 0.25, 0.05),
                        "CCRatioThresh" = 25)
  NREPS = length(masterDF$OutbreakCrash)
  masterDF$RUNNOCOTS = c(T, rep(F, NREPS-1))

} else {
  
  setwd(DIRECTORY)
  source("COTSModel_Utilityfunctions.R")  
  source("COTSModel_COTSFunctions.R")
  
  # source("COTSModel_LoadObjectsForModelling.R")
  # save.image(file = "RData/COTSMod_bckp.Rdata")
  load("RData/COTSMod_bckp.RData")
  
  # Set Global Parameters ----
  
  set.seed(123)
  NREPS <- 10 # How many rows from master DF    
  Years <- 1996:2017
  nyears <- length(Years)
  NSEASONS <- 2
  seasons <- c("summer","winter")
  npops <- 15802 #number of reefs we want to test
  nsimul <- 10
  
  VERBOSE <- TRUE        # flag whether functions should return detailed information
  DEBUG <- TRUE          # flag whether to output debug files etc. 
  
  projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   
  
  stagenames <- c('J_1', 'J_2', 'A')
  COTSmort <- c(0.8,0.7,0.2)
  names(COTSmort) <- stagenames
  COTSremain <- c(0.02,0.2,1) 
  names(COTSremain) <- stagenames
  COTS_StableStage <- c(0.9803, 0.0171, 0.0026)   
  
  # Resample Disturbance ----
  
  LHS <- lhs::randomLHS(nsimul,9)
  B.BLEACHING <- qnorm(LHS[,1], mean=bleaching.mn.sd[1], sd=bleaching.mn.sd[2])
  B.COTS <- qnorm(LHS[,2], mean=COTS.mn.sd[1], sd=COTS.mn.sd[2])
  B.DISEASE <- qnorm(LHS[,3], mean=disease.mn.sd[1], sd=disease.mn.sd[2])
  B.STORMS <- qnorm(LHS[,4], mean=storms.mn.sd[1], sd=storms.mn.sd[2])
  B.UNKNOWN <- qnorm(LHS[,5], mean=unknown.mn.sd[1], sd=unknown.mn.sd[2])
  # B.WQ <- qnorm(LHS[,6], mean=WQ.mn.sd[1], sd=WQ.mn.sd[2])
  
  
  A <- B0 <- HCINI <- HCMAX <- matrix(NA, ncol = nsimul, nrow = dim(data.grid)[1])
  
  for (i in 1:nsimul) {
    A[,i] <- qnorm(LHS[i,7], mean=data.grid$pred.a.mean, sd=data.grid$pred.a.sd)
    B0[,i] <- qnorm(LHS[i,8], mean=data.grid$pred.b0.mean, sd=data.grid$pred.b0.sd)
    HCINI[,i] <- qnorm(LHS[i,9], mean=data.grid$pred.HCini.mean, sd=data.grid$pred.HCini.sd)  
    HCMAX[,i] <- qnorm(LHS[i,9], mean=data.grid$pred.HCmax.mean, sd=data.grid$pred.HCmax.sd) 
  }
  
  HCMAX[HCINI > HCMAX] <- HCINI[HCINI > HCMAX]
  
  
  # res <- array(NA, dim=c(dim(data.grid)[1], nyears, nsimul)) ## Stores the coral cover values for each grid cell (rows), year (columns) and simulation (third dimension)
  bleaching.rsmpl <- COTS.rsmpl <- disease.rsmpl <- storms.rsmpl <- unknown.rsmpl <- array(NA, dim=c(dim(data.grid)[1], nyears, nsimul)) ## Stores the actual (i.e. resampled) disturbance intensity values for each grid cell (rows), year (columns) and simulation (third dimension)
  
  # Create back up files of disturbance tables 
  data.bleaching[is.na(data.bleaching)] <- 0
  data.COTS[is.na(data.COTS)] <- 0
  data.storms[is.na(data.storms)] <- 0
  
  data.bleaching.bckp <- data.bleaching
  data.COTS.bckp <- data.COTS
  data.storms.bckp <- data.storms 
  
  # Resample Distrubances for each simulation    
  for (j in 1:nsimul) {
  
    # Re-initialize disturbance data
    data.bleaching <- data.bleaching.bckp
    data.COTS <- data.COTS.bckp
    data.storms <-  data.storms.bckp
    
    # Resample disturbance data in each year
    data.unknown <- data.disease <- data.COTS
    colnames(data.unknown)[6:38] = paste0("UNKN_", 1985:2017)
    colnames(data.disease)[6:38] = paste0("DIS_", 1985:2017)
    data.unknown[,6:38] <- data.disease[,6:38] <- 0
    
    for (i in 1:nyears) {
      ## Simulate disease and unknown disturbance based on observed frequencies
      data.unknown[,i+16] <- sampling::srswor(round(length(data.unknown[,i+16])*0.01),length(data.unknown[,i+16]))
      data.disease[,i+16] <- sampling::srswor(round(length(data.disease[,i+16])*0.01),length(data.disease[,i+16]))
      
      ## Resample other disturbance based on P(Impact|Disturbance)
      count.cots <- length(data.COTS[,i+16][data.COTS[,i+16]>0])  
      count.storms <- length(data.storms[,i+16][data.storms[,i+16]>0])
      count.bleaching <- length(data.bleaching[,i+16][data.bleaching[,i+16]>0])
      if (count.cots>0)  data.COTS[,i+16][data.COTS[,i+16]>0][sample(count.cots, count.cots*.5)] <- 0
      if (count.storms>0)  data.storms[,i+16][data.storms[,i+16]>0][sample(count.storms, round(count.storms*.5))] <- 0     
      if (count.bleaching>0)  data.bleaching[,i+16][data.bleaching[,i+16]>0][sample(count.bleaching, count.bleaching*.1)] <- 0
    }
    data.disease[,-(1:5)][data.disease[,-(1:5)]>1] <- 1
    data.unknown[,-(1:5)][data.unknown[,-(1:5)]>1] <- 1
    
    data.storms[,"Hs4MW_2009"] <- data.storms[,"Hs4MW_2009"]*.5
    data.storms[,"Hs4MW_2013"] <- data.storms[,"Hs4MW_2013"]*.5
    data.storms[,"Hs4MW_2014"] <- data.storms[,"Hs4MW_2014"]*.5
    data.storms[,"Hs4MW_2015"] <- data.storms[,"Hs4MW_2015"]*.5
    
    # Add known disturbance for LTMP reefs
    data.bleaching[,-(1:16)][!is.na(data.ltmp.bleaching[,-(1:5)])] <- data.ltmp.bleaching[,-(1:5)][!is.na(data.ltmp.bleaching[,-(1:5)])]
    data.COTS[,-(1:17)][!is.na(data.ltmp.COTS[,-(1:6)])] <- data.ltmp.COTS[,-(1:6)][!is.na(data.ltmp.COTS[,-(1:6)])]
    data.storms[,-(1:16)][!is.na(data.ltmp.storms[,-(1:5)])] <- data.ltmp.storms[,-(1:5)][!is.na(data.ltmp.storms[,-(1:5)])]
    data.disease[,-(1:16)][!is.na(data.ltmp.disease[,-(1:5)])] <- data.ltmp.disease[,-(1:5)][!is.na(data.ltmp.disease[,-(1:5)])]
    data.unknown[,-(1:16)][!is.na(data.ltmp.unknown[,-(1:5)])] <- data.ltmp.unknown[,-(1:5)][!is.na(data.ltmp.unknown[,-(1:5)])]
    
    # store the resampled distrubances for the replicate j
    bleaching.rsmpl[,,j] <- as.matrix(data.bleaching[,-(1:16)])
    COTS.rsmpl[,,j] <- as.matrix(data.COTS[,-(1:16)])
    disease.rsmpl[,,j] <- as.matrix(data.disease[,-(1:16)])
    storms.rsmpl[,,j] <- as.matrix(data.storms[,-(1:16)])
    unknown.rsmpl[,,j] <- as.matrix(data.unknown[,-(1:16)])
    
    for (i in 1:nsimul) bleaching.rsmpl[,,i][is.na(bleaching.rsmpl[,,i])] <- 0
    
    bleaching.mn <- apply(bleaching.rsmpl, c(1,2), mean, na.rm=T)
    COTS.mn <- apply(COTS.rsmpl, c(1,2), mean, na.rm=T)
    disease.mn <- apply(disease.rsmpl, c(1,2), mean, na.rm=T)
    storms.mn <- apply(storms.rsmpl, c(1,2), mean, na.rm=T)
    unknown.mn <- apply(unknown.rsmpl, c(1,2), mean, na.rm=T)
  }
  # Write Parameter Data Frame And Results Containers ----
  masterDF = data.frame("SexRatio" = 5, 
                        "ConsRateW" = 150, 
                        "ConsRateS" = 250,
                        "avgPCF" = 20000000,
                        "sdPCF" = 10000000,
                        "mortJ1" =  0.99,
                        "mortJ2" = 0.95,
                        "mortA" = 0.9,
                        "remJ1" = 0,
                        "remJ2" = 0,
                        "remA" = 1,
                        "cssJ1" = 0.9803,
                        "cssJ2" = 0.0171,
                        "cssA" = 0.0026,
                        "Pred" = 0.99,
                        "p" = 0.2,
                        "Crash" = 0,
                        "OutbreakCrash" = Inf,
                        "Fbase" = seq(0.05, 0.25, 0.05),
                        "CCRatioThresh" = 25)
  
  NREPS = length(masterDF$OutbreakCrash)
  masterDF$RUNNOCOTS = c(T, rep(F, NREPS-1))
  
  # Results storage arrays (Pixel, Year/Season, Simulation)
  res.cc = array(NA, dim=c(dim(data.grid)[1], nyears*2, nsimul))
  res.cots = array(NA, dim=c(dim(data.grid)[1], nyears*2, nsimul))
  
  Results = data.frame(sapply(PopData[1:4], rep, times=nyears*NSEASONS),
                       sapply(PopData[5:7], rep, times=nyears*NSEASONS),
                       Year=rep(Years,each=2*npops), Season=rep(c("summer", "winter"),each=npops),
                       COTSJ1=NA, COTSJ2=NA, COTSA=NA, CoralCover=NA, CoralCover.DistLoss=NA)
  Results$CoralCover.Consum = NA
  Results$CoralCover.Growth = NA
  
  # Save Model Workspace for Prosperity ----
  setwd(DIRECTORY)
  nruns = length(list.files(path = "RData")[grep(Sys.Date(),list.files(path = "RData"))]) + 1
  save.image(file = paste0("RData/ModelWorkspace_", Sys.Date(),"_", nruns,".RData"))
}

# Run Model ----  

doParallel::registerDoParallel(cores=10)

foreach::foreach (reps = 1:NREPS) %dopar% {
  
  `%>%` <- magrittr::`%>%`
  DIRECTORY = getwd()
  setwd(DIRECTORY)

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
  Fbase = masterDF[reps,"Fbase"]
  CCRatioThresh = masterDF[reps,"CCRatioThresh"]
  RUNNOCOTS = masterDF[reps, "RUNNOCOTS"]
  # Initialize Model ----
  
  
  browse = FALSE 
  inityear = 1995
  COTSfromCoralModel=FALSE 
  COTSfromSimul=TRUE
 
  # Simulation loop
  for (j in 1:nsimul) {
    COTSabund = initializeCOTSabund(data.grid, data.COTS, inityear, 
                                    stagenames, COTS_StableStage, npops, j)
    if (RUNNOCOTS == T) {
      COTSfromCoralModel=T 
      COTSfromSimul=F
      COTSabund <- matrix(0,nrow=npops, ncol=3, dimnames = list(NULL, c("J_1", "J_2", "A")))
    }
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
        COTSabund = doPredPreyDynamics(COTSabund, CoralCover, p, Crash, CCRatioThresh)
        COTSabund = doCOTSDispersal(season,COTSabund,CoralCover,SexRatio,COTS.ConnMat, 
                                    PCFParams, Pred, FvDParams, Fbase, CCRatioThresh) #Pruducing NAS
        COTSabund = doCOTSDemography(season, COTSabund, COTSmort, COTSremain)
        Consumption = doCoralConsumption(season, COTSabund, CoralCover, ConsRate) 
        CoralCover = Consumption[,'CRemaining']
        CoralConsum = round(Consumption[,'CChange'],4)
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
        Results[(Results$Year==i+1995) & (Results$Season==season),
                c("COTSJ1", "COTSJ2", "COTSA", "CoralCover", "CoralCover.DistLoss", "CoralCover.Consum", 'CoralCover.Growth')] =
          cbind(COTSabund, CoralCover, Disturbance, CoralConsum, CoralGrowth)
        if(i>OutbreakCrash & season =="winter") {
          OutbreakCrasher = Results %>%
            dplyr::filter(Year > (i+1995-OutbreakCrash) & Year <= i+1995) %>%
            dplyr::group_by(REEF_NAME, Year) %>%
            dplyr::mutate(COTSA = (COTSA/100)*0.15) %>% # need to allow for detection
            dplyr::summarise(Mean.COTS = mean(COTSA)) %>%
            dplyr::mutate(is.outbreak = ifelse(Mean.COTS > 1, 1,0)) %>%
            dplyr::group_by(REEF_NAME) %>%
            dplyr::summarise(Crash = ifelse(sum(is.outbreak)==OutbreakCrash,1,0)) %>%
            dplyr::filter(Crash==1)
          matchit = match(data.grid$REEF_NAME,OutbreakCrasher$REEF_NAME)
          COTSabund[which(!is.na(matchit)),] = c(0,0,0)
        } # Close outbreak crasher
        if (season=="summer") {
          res.cc[,2*i-1,j] = CoralCover
          res.cots[,2*i-1,j] = COTSabund[,3]
        }
        res.cc[,2*i,j] = CoralCover
        res.cots[,2*i,j] = COTSabund[,3]
      } # close season loop
    } # close Year loop
    
  } # Close simulation loop
  
  
  # Reef Level Summaries
  resCC.reef <- array(0, dim=c(length(unique(data.grid$REEF_ID)), nyears*2, nsimul))  
  resCC.cluster <- array(0, dim=c(length(unique(data.grid$bent.clust)), nyears*2, nsimul))
  resCOTS.reef <- array(0, dim=c(length(unique(data.grid$REEF_ID)), nyears*2, nsimul))  
  resCOTS.cluster <- array(0, dim=c(length(unique(data.grid$bent.clust)), nyears*2, nsimul))  
  for (i in 1:dim(res.cc)[3]) {
    resCC.reef[,,i] <- as.matrix(aggregate(res.cc[,,i], by=list(data.grid$REEF_ID),
                                           FUN=mean, na.rm=T)[-1])
    resCC.cluster[,,i] <- as.matrix(aggregate(res.cc[,,i], by=list(data.grid$bent.clust),
                                              FUN=mean, na.rm=T)[-1])
    resCOTS.reef[,,i] <- as.matrix(aggregate(res.cots[,,i], by=list(data.grid$REEF_ID),
                                             FUN=mean, na.rm=T)[-1])
    resCOTS.cluster[,,i] <- as.matrix(aggregate(res.cots[,,i], by=list(data.grid$bent.clust),
                                                FUN=mean, na.rm=T)[-1])
  }
  
  resCC.reef.mn <- apply(resCC.reef, c(1,2), mean, na.rm=T)
  resCC.reef.med <- apply(resCC.reef, c(1,2), median, na.rm=T)
  # resCC.reef.min <- apply(resCC.reef, c(1,2), quantile, probs=0.05, na.rm=T)
  # resCC.reef.max <- apply(resCC.reef, c(1,2), quantile, probs=0.95, na.rm=T)
  resCC.reef.25 <- apply(resCC.reef, c(1,2), quantile, probs=0.25, na.rm=T)
  resCC.reef.75 <- apply(resCC.reef, c(1,2), quantile, probs=0.75, na.rm=T)
  
  resCOTS.reef.mn <- apply(resCOTS.reef, c(1,2), mean, na.rm=T)
  resCOTS.reef.med <- apply(resCOTS.reef, c(1,2), median, na.rm=T)
  # resCOTS.reef.min <- apply(resCOTS.reef, c(1,2), quantile, probs=0.05, na.rm=T)
  # resCOTS.reef.max <- apply(resCOTS.reef, c(1,2), quantile, probs=0.95, na.rm=T)
  resCOTS.reef.25 <- apply(resCOTS.reef, c(1,2), quantile, probs=0.25, na.rm=T)
  resCOTS.reef.75 <- apply(resCOTS.reef, c(1,2), quantile, probs=0.75, na.rm=T)
  
  nReefs = length(unique(data.grid$REEF_ID))
  # Results df for Dashboard -- Reef Level
  ResultsDash = data.frame(sapply(unique(data.grid[4:5]), rep, times=nyears*NSEASONS),
                           Year=rep(Years,each=2*nReefs), 
                           Season=rep(c("summer", "winter"),each=nReefs), 
                           COTS.mn=(as.vector(resCOTS.reef.mn)/667), 
                           COTS.Q50=(as.vector(resCOTS.reef.med)/667), 
                           # COTS.Q05=(as.vector(resCOTS.reef.min)/667), 
                           # COTS.Q95=(as.vector(resCOTS.reef.max)/667), 
                           COTS.Q25=(as.vector(resCOTS.reef.25)/667), 
                           COTS.Q75=(as.vector(resCOTS.reef.75)/667),
                           CC.mn=as.vector(resCC.reef.mn), 
                           CC.Q50=as.vector(resCC.reef.med), 
                           # CC.Q05=as.vector(resCC.reef.min), 
                           # CC.Q95=as.vector(resCC.reef.max), 
                           CC.Q25=as.vector(resCC.reef.25), 
                           CC.Q75=as.vector(resCC.reef.75)) 
  ResultsDash = dplyr::left_join(ResultsDash, unique(data.grid[5:7]), by="REEF_NAME") %>%
                dplyr::select(REEF_ID:Season, SECTOR:CROSS_SHELF, 5:12)
  
  setwd(DIRECTORY)
  setwd("Results")
  name <- sprintf("Sample_%s.Rdata", reps)
  save(res.cc, res.cots, ResultsDash, file = name) 
  
}

# Combine Files into summary 
setwd(DIRECTORY)
setwd("Results")
load("Sample_1.Rdata")

ForDashboard = ResultsDash
files = list.files()
files = files[grep("Sample",files)]

for (i in 2:NREPS) {
  load(sprintf("Sample_%s.Rdata", i))
  ForDashboard = cbind(ForDashboard, ResultsDash[7:14])
}

colnames(ForDashboard)[7:length(colnames(ForDashboard))] = paste0(rep(1:NREPS, each=8),"_", colnames(ResultsDash[7:14]))

nruns = length(list.files(path = "Archive")[grep(paste0("ForDashboard_",Sys.Date()),list.files(path = "Archive"))]) + 1
write.csv(ForDashboard, paste0("Archive/ForDashboard_", Sys.Date(),"_", nruns, ".csv"))
write.csv(masterDF, paste0("Archive/masterDF_", Sys.Date(),"_", nruns, ".csv"))
















