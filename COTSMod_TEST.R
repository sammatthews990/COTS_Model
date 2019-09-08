# COTSMod: Spatially Explicit COTS-Coral Metapopulation Model ----


# Load Functions and Modelling objects ----

rm(list=ls())
source("COTSModel_Utilityfunctions.R")
# loadPackages()
library(dplyr)
library(foreach)
library(doParallel)
library(ggplot2)
library(pscl)
# library(grid)

DIRECTORY = getwd()

PRELOAD = T
LHSPARAMS = T
SUBSET=F
RESUBSET=F

# DIRECTORY = "C:/Users/jc312264/OneDrive - James Cook University/GitHub/COTS_Model"

if (PRELOAD == T) {
  setwd(DIRECTORY)
  if (SUBSET==T) {load("RData/ModelWorkspace_TEST.RData")}
  else {load("RData/ModelWorkspace_FULL.RData")}
  # load("RData/ModelWorkspace_FULL.RData")
  DIRECTORY = getwd()
  # Choose an older workspace with pre run disturbance samples
  # Write Parameter Data Frame And Results Containers ----
  source("COTSModel_COTSFunctions.R")
  source("COTSModel_Utilityfunctions.R")
  
  LHSPARAMS = T
  
  if (LHSPARAMS==T) {
    # nsimul=100
    masterDF = MakeLHSSamples(80)
    setwd(DIRECTORY)
  } else {
    # nsimul=100
    masterDF = data.frame("SexRatio" = 5, 
                          "ConsRateW" = 0, 
                          "ConsRateS" = 0,
                          "avgPCF" = 20000000,
                          "sdPCF" = 10000000,
                          "mortJ1" =  0.99,
                          "mortJ2" = 0.7,
                          "mortA" = 0.1,
                          "remJ1" = 0,
                          "remJ2" = 0,
                          "remA" = 1,
                          "cssJ1" = 0.9803,
                          "cssJ2" = 0.0171,
                          "cssA" = 0.0026,
                          "Pred" = 0.985, #play with this
                          "Crash" = 0,
                          "OutbreakCrash" = Inf,
                          "Fbase" = 0,
                          "CCRatioThresh" = 30,
                          "CCRatioThresh2" = 10,
                          "maxmort" = 1, #and this
                          "selfseed" = 0.5, # and this
                          "chl.int" = 10, #and this
                          "Cbase" = 0.1,
                          "CMax" = 200,
                          "RUNNOCOTS"=F,
                          "J2M" = -0.1e5,
                          "J1M" = -1e7,
                          "J2R" = 0.000013,
                          "J1R" = 0.00000018, 
                          "AM" = c(15,15,20,20),
                          "AR" = 0.8,
                          "Linf" = 0.7,
                          "K" = 4e-04,
                          "t0" = 0)
  }
  NREPS = length(masterDF$OutbreakCrash)
  masterDF$RUNNOCOTS = c(T, rep(F, NREPS-1))
  masterDF$OutbreakCrash = rep(c(3,Inf))
  masterDF$DensDepend = TRUE
  
  
} else {
  
  setwd(DIRECTORY)
  
  source("COTSModel_LoadObjectsForModelling.R")
  setwd(DIRECTORY)
  # save.image(file = "RData/COTSMod_bckp.Rdata")
  #load("RData/COTSMod_bckp.RData")
  source("COTSModel_Utilityfunctions.R")
  setwd(DIRECTORY)
  source("COTSModel_COTSFunctions.R")
  data.grid$PIXEL_ID = PopData$PIXEL_ID = 1:length(data.grid[,1])
  
  
  # Set Global Parameters ----
  
  set.seed(123)
  NREPS <- 10 # How many rows from master DF    
  Years <- 1996:2017
  nyears <- length(Years)
  NSEASONS <- 2
  seasons <- c("summer","winter")
  npops <- 15802 #number of reefs we want to test
  nsimul <- 100
  
 
  
  VERBOSE <- TRUE        # flag whether functions should return detailed information
  DEBUG <- TRUE          # flag whether to output debug files etc. 
  
  projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   
  
  stagenames <- c('J_1', 'J_2', 'A')
  # COTSmort <- c(0.8,0.7,0.2)
  # names(COTSmort) <- stagenames
  # COTSremain <- c(0.02,0.2,1) 
  # names(COTSremain) <- stagenames
  # COTS_StableStage <- c(0.9803, 0.0171, 0.0026)   
  
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
    # browser()
    # Add known disturbance for LTMP reefs
    # data.ltmp.COTS$X1996 = data.ltmp.COTS$X1997
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
  if (LHSPARAMS==T) {
    masterDF = MakeLHSSamples(20)
    masterDF$OutbreakCrash = Inf
    setwd(DIRECTORY)
  } else {
    masterDF = data.frame("SexRatio" = 5, 
                          "ConsRateW" = 0, 
                          "ConsRateS" = 0,
                          "avgPCF" = 20000000,
                          "sdPCF" = 10000000,
                          "mortJ1" =  0.99,
                          "mortJ2" = 0.7,
                          "mortA" = 0.1,
                          "remJ1" = 0,
                          "remJ2" = 0,
                          "remA" = 1,
                          "cssJ1" = 0.9803,
                          "cssJ2" = 0.0171,
                          "cssA" = 0.0026,
                          "Pred" = 0.985, #play with this
                          "Crash" = 0,
                          "OutbreakCrash" = Inf,
                          "Fbase" = 0,
                          "CCRatioThresh" = 30,
                          "CCRatioThresh2" = 10,
                          "maxmort" = 1, #and this
                          "selfseed" = 1, # and this
                          "chl.int" = seq(-0.04, 1, length.out = 4), #and this
                          "Cbase" = 0.1,
                          "CMax" = 350,
                          "RUNNOCOTS"=F,
                          "J2M" = 1,
                          "J1M" = 2.5,
                          "J2R" = 0.00002,
                          "J1R" = 0.0000002,
                          "Linf" = 0.7,
                          "K" = 7e-04,
                          "t0" = 0)
  }
  
  NREPS = length(masterDF$OutbreakCrash)
  masterDF$RUNNOCOTS = c(T, rep(F, NREPS-1))
  masterDF$OutbreakCrash = Inf
  masterDF$DensDepend = TRUE
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
  save.image(file = "RData/ModelWorkspace_FULL.RData", version = 2)
  # if (SUBSET==T) {save.image(file = "RData/ModelWorkspace_TEST.RData")}else{save.image(file = "RData/ModelWorkspace_FULL.RData")}
  save.image(file = paste0("RData/ModelWorkspace_", Sys.Date(),"_", nruns,".RData"),version = 2)
}

# browser()

WHERE = list(c("CA", "CL", "IN","TO"), c("M"))
RESUBSET = F

if (RESUBSET == T) {
  load("RData/ModelWorkspace_FULL.RData")
  REEFSUB = data.grid %>% dplyr::select(1:7) %>% 
    dplyr::filter(SECTOR %in% WHERE[[1]] & CROSS_SHELF %in% WHERE[[2]]) %>% dplyr::arrange(PIXEL_ID)
  # List of objects to reduce in size
  reduce = c("A", "B0", "bleaching.mn", "COTS.mn", "disease.mn", "HCINI", "HCMAX", "PercentReef", "PopData", "storms.mn", "unknown.mn",
             "bleaching.rsmpl", "COTS.rsmpl", "data.chl.resid", "disease.rsmpl", "storms.rsmpl", "unknown.rsmpl", "WQ", "res.cc", "res.cots", ls(pattern = "data."))
  reduce = reduce[-which(reduce %in% c("data.manta","data.manta.env","data.rap", "data.reef", "data.WQ", "data.COTSPred"))]
  dfs = mget(list(reduce)[[1]])
  for (i in 1:length(dfs)) {
    ndim = length(dim(dfs[[i]]))
    if (ndim==0)  {dfs[[i]] = dfs[[i]][REEFSUB$PIXEL_ID]}
    if (ndim==1)  {dfs[[i]] = dfs[[i]][REEFSUB$PIXEL_ID]}
    if (ndim==2)  {dfs[[i]] = dfs[[i]][REEFSUB$PIXEL_ID,]}
    if (ndim==3)  {dfs[[i]] = dfs[[i]][REEFSUB$PIXEL_ID,,]}
  }
  list2env(dfs, .GlobalEnv)
  REEFSUB.index = which(Pixels$REEF_NAME %in% REEFSUB$REEF_NAME)
  data.COTSPred = data.COTSPred %>% dplyr::filter(REEF_NAME %in% REEFSUB$REEF_NAME)
  Pixels = Pixels[REEFSUB.index,]
  COTS.ConnMat = COTS.ConnMat[REEFSUB.index, REEFSUB.index]
  npops = length(REEFSUB$PIXEL_ID)
  rm(dfs)
  save.image(file = "RData/ModelWorkspace_TEST.RData")
  }

# Parallel Loop ----


cl = parallel::makeCluster(2)
doParallel::registerDoParallel(cores = 20)

foreach::foreach (reps = 1:NREPS) %dopar% {
  
  `%>%` <- magrittr::`%>%`
  DIRECTORY = getwd()
  rm(FvDParams)
  # DIRECTORY = "C:/Users/jc312264/OneDrive - James Cook University/GitHub/COTS_Model"
  nsimul=10
  SexRatio = masterDF[reps, "SexRatio"]
  ConsRate = as.vector(masterDF[reps, 2:3])
  PCFParams = c(masterDF[reps, "avgPCF"], masterDF[reps,"sdPCF"])
  COTSmort = as.numeric(masterDF[reps, c("mortJ1", "mortJ2", "mortA")])
  FvDParams = as.numeric(masterDF[reps, c("Linf", "K", "t0")])
  COTSremain = as.numeric(masterDF[reps, c("remJ1", "remJ2", "remA")])
  COTS_StableStage = as.numeric(masterDF[reps, c("cssJ1", "cssJ2", "cssA")])
  Pred = masterDF[reps,"Pred"]
  Crash = masterDF[reps,"Crash"]
  OutbreakCrash = masterDF[reps,"OutbreakCrash"]
  Fbase = masterDF[reps,"Fbase"]
  CCRatioThresh = masterDF[reps,"CCRatioThresh"]
  CCRatioThresh2 = masterDF[reps,"CCRatioThresh2"]
  maxmort = masterDF[reps,"maxmort"]
  selfseed = masterDF[reps,"selfseed"]
  chl.int = masterDF[reps, "chl.int"]
  Cbase = masterDF[reps, "Cbase"]
  CMax = masterDF[reps, "CMax"]
  J2M = masterDF[reps, "J2M"]
  J1M = masterDF[reps, "J1M"]
  J2R = masterDF[reps, "J2R"]
  J1R = masterDF[reps, "J1R"]
  AM = masterDF[reps, "AM"]
  AR = masterDF[reps, "AR"]
  DensDepend = masterDF[reps, "DensDepend"]
  RUNNOCOTS = masterDF[reps, "RUNNOCOTS"]
  
  # Initialize Model ----
  chl.lm$coefficients[1] = chl.int
  seasons <- c("summer","winter")
  browse = FALSE 
  inityear = 1996
  COTSfromCoralModel=FALSE 
  COTSfromSimul=TRUE
  
  # Simulation loop
  for (j in 1:nsimul) {
    COTSabund = initializeCOTSabund(data.grid, data.COTS, inityear, 
                                    stagenames, COTS_StableStage, npops, j)
    if (RUNNOCOTS == T) {
      COTSfromCoralModel=T 
      COTSfromSimul=F
      #COTSabund <- matrix(0,nrow=npops, ncol=3, dimnames = list(NULL, c("J_1", "J_2", "A")))
    }
    print(j)
    HC.1996 <- HCINI[,j]
    b0 <- B0[,j]
    b1 <- B0[,j]/log(HCMAX[,j])
    res.cc[,1,j] <- as.numeric(HC.1996)
    CoralCover = HC.1996
    # Year Loop
    for(i in 1:length(Years)){  
      (Year = i + 1995) # loop through years
      # Season Loop
      for(season in seasons){ 
        if(browse == TRUE) {
          browser()
        }
        # browser()
        COTSabund = doPredPreyDynamics(COTSabund, CoralCover, Crash, CCRatioThresh, CCRatioThresh2, maxmort, 
                                       J2M, J1M,J2R, J1R,AM, AR, season, COTSmort, DensDepend)
        COTSabund = doCOTSDemography(season, COTSabund, COTSmort, COTSremain)
        Consumption = doCoralConsumption(season, COTSabund, CoralCover, ConsRate, COTSfromCoralModel, Cbase,CMax, CCRatioThresh)
        CoralCover = Consumption[,'CRemaining']
        CoralConsum = round(Consumption[,'CChange'],4)
        COTSabund = doCOTSDispersal(season,COTSabund,CoralCover,SexRatio,COTS.ConnMat, 
                                    PCFParams, Pred, FvDParams, Fbase, CCRatioThresh,Year, data.chl, data.chl.resid, j, selfseed) 
        
        #Pruducing NAS
         
        
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
        # Convert to manta tow abundance
        Results$COTSA = predict(MTCalib.gaminv,newdata=data.frame(DENS=Results$COTSA))
        if(i>=OutbreakCrash & season =="winter") {
          # browser()
          OutbreakCrasher = Results %>%
            dplyr::filter(Year >= (i+1995-OutbreakCrash) & Year <= i+1995) %>%
            dplyr::group_by(REEF_NAME, Year) %>% # need to allow for detection
            dplyr::summarise(Mean.COTS = mean(COTSA)) %>%
            dplyr::mutate(is.outbreak = ifelse(Mean.COTS > 0.2, 1,0)) #%>%
            dplyr::group_by(REEF_NAME) %>%
            dplyr::summarise(sumcrash = sum(is.outbreak),
              Crash = ifelse(sum(is.outbreak)>=OutbreakCrash,1,0)) #%>%
            dplyr::filter(Crash>0)
            # browser()
          # if (nrow(OutbreakCrasher) >0) {browser()}  
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
  
  
  data.grid$REEF_ID = factor(data.grid$REEF_ID)
  # Reef Level Summaries
  resCC.reef <- array(0, dim=c(length(unique(data.grid$REEF_ID)), nyears*2, nsimul))  
  resCC.cluster <- array(0, dim=c(length(unique(data.grid$bent.clust)), nyears*2, nsimul))
  resCOTS.reef <- array(0, dim=c(length(unique(data.grid$REEF_ID)), nyears*2, nsimul))  
  resCOTS.cluster <- array(0, dim=c(length(unique(data.grid$bent.clust)), nyears*2, nsimul))  
  for (i in 1:nsimul) {
    resCC.reef[,,i] <- as.matrix(aggregate(res.cc[,,i], by=list(data.grid$REEF_ID),
                                           FUN=mean, na.rm=T)[-1])
    resCC.cluster[,,i] <- as.matrix(aggregate(res.cc[,,i], by=list(data.grid$bent.clust),
                                              FUN=mean, na.rm=T)[-1])
    resCOTS.reef[,,i] <- as.matrix(aggregate(res.cots[,,i], by=list(data.grid$REEF_ID),
                                             FUN=mean, na.rm=T)[-1])
    resCOTS.cluster[,,i] <- as.matrix(aggregate(res.cots[,,i], by=list(data.grid$bent.clust),
                                                FUN=mean, na.rm=T)[-1])
  }
  
  # resCC.mn <- apply(res.cc, c(1,2), mean, na.rm=T)
  # resCC.med <- apply(res, c(1,2), median, na.rm=T)
  # resCC.min <- apply(res, c(1,2), quantile, probs=0.05, na.rm=T)
  # resCC.max <- apply(res, c(1,2), quantile, probs=0.95, na.rm=T)
  # resCC.25 <- apply(res, c(1,2), quantile, probs=0.25, na.rm=T)
  # resCC.75 <- apply(res, c(1,2), quantile, probs=0.75, na.rm=T)
  
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
  
  ResultsDash = data.frame(sapply(data.frame(REEF_ID =levels(data.grid$REEF_ID)), rep, times=nyears*NSEASONS),
                           Year=rep(Years,each=2*nReefs), 
                           Season=rep(c("summer", "winter"),each=nReefs), 
                           COTS.mn=(as.vector(resCOTS.reef.mn)), 
                           COTS.Q50=(as.vector(resCOTS.reef.med)), 
                           # COTS.Q05=(as.vector(resCOTS.reef.min)/667), 
                           # COTS.Q95=(as.vector(resCOTS.reef.max)/667), 
                           COTS.Q25=(as.vector(resCOTS.reef.25)), 
                           COTS.Q75=(as.vector(resCOTS.reef.75)),
                           CC.mn=as.vector(resCC.reef.mn), 
                           CC.Q50=as.vector(resCC.reef.med), 
                           # CC.Q05=as.vector(resCC.reef.min), 
                           # CC.Q95=as.vector(resCC.reef.max), 
                           CC.Q25=as.vector(resCC.reef.25), 
                           CC.Q75=as.vector(resCC.reef.75)) 
  ResultsDash = dplyr::left_join(ResultsDash, unique(data.grid[4:7]), by="REEF_ID") %>%
    dplyr::select(REEF_ID, REEF_NAME,SECTOR:CROSS_SHELF, Year, Season, 4:12)
  
 
  #Save Results ----
  # setwd("Results")
  name <- sprintf("Results/Sample_%s.Rdata", reps)
  save(
    # res.cc, res.cots, 
    ResultsDash, file = name) 
  
}

sort(sapply(ls(), function(x){object.size(get(x))}))
rm(dfs)

parallel::stopCluster(cl)

# Combine Files into summary ----
setwd(DIRECTORY)
# rm(list=ls())
# source("COTSModel_Utilityfunctions.R")
# DIRECTORY=getwd()
setwd("Results")
load("Sample_1.Rdata")

ForGraph = ResultsDash
ForGraph$REP=1
files = list.files()
files = files[grep("Sample",files)]
NREPS = nrow(masterDF)

for (i in 2:NREPS) {
  load(sprintf("Sample_%s.Rdata", i))
  ResultsDash$REP=i
  # ForDashboard = cbind(ForDashboard, ResultsDash[7:14])
  ForGraph = rbind(ForGraph, ResultsDash)
}

setwd(DIRECTORY)

ForGraph = ForGraph %>% 
  mutate(COTS.mn = round(predict(MTCalib.gaminv, newdata=data.frame(DENS=COTS.mn))^2,3),
         COTS.Q50 = round(predict(MTCalib.gaminv, newdata=data.frame(DENS=COTS.Q50))^2,3),
         COTS.Q25 = round(predict(MTCalib.gaminv, newdata=data.frame(DENS=COTS.Q25))^2,3),
         COTS.Q75 = round(predict(MTCalib.gaminv, newdata=data.frame(DENS=COTS.Q75))^2,3))


ForGraph$YearDec = as.numeric(ifelse(ForGraph$Season=="summer", ForGraph$Year, c(ForGraph$Year, ".5")))
# colnames(ForDashboard)[7:length(colnames(ForDashboard))] = paste0(rep(1:NREPS, each=8),"_", colnames(ResultsDash[7:14]))



#### Validation ----
manta.reefs = data.manta %>% 
  dplyr::filter(REEF_NAME.y %in% unique(data.grid$REEF_NAME) & REPORT_YEAR > 1995) %>% 
  mutate(REEF_NAME = as.factor(REEF_NAME.y)) %>%
  group_by(A_SECTOR, REEF_NAME) %>% 
  dplyr::summarise(n=length(REPORT_YEAR)) %>% filter(n>10) %>% 
  pull(REEF_NAME) %>% as.character()

manta.reefs.all = data.manta %>% 
  dplyr::filter(REEF_NAME.y %in% unique(data.grid$REEF_NAME) & REPORT_YEAR > 1995) %>% 
  mutate(REEF_NAME = as.factor(REEF_NAME.y)) %>%
  group_by(A_SECTOR, REEF_NAME) %>% 
  dplyr::summarise(n=length(REPORT_YEAR)) %>% filter(n>4) %>%
  pull(REEF_NAME) %>% as.character()

manta.reefs.valid = data.manta %>% 
  dplyr::filter(REEF_NAME.y %in% unique(data.grid$REEF_NAME) & REPORT_YEAR > 1995) %>% 
  mutate(REEF_NAME = as.factor(REEF_NAME.y)) %>%
  group_by(A_SECTOR, REEF_NAME) %>% 
  dplyr::summarise(n=length(REPORT_YEAR)) %>% filter(n>4) %>%
  dplyr::filter(!REEF_NAME %in% manta.reefs) %>%
  pull(REEF_NAME) %>% as.character()

res.plot = ForGraph %>% filter(REEF_NAME %in% data.grid$REEF_NAME) %>% 
  dplyr::group_by(REEF_NAME,REP,Year,SECTOR, CROSS_SHELF) %>%
  dplyr::summarise(COTS.mn = mean(COTS.mn), # add in detection
                   COTS.25 = mean(COTS.Q25),
                   COTS.75 = mean(COTS.Q75),
                   CC.mn = mean(CC.mn),
                   CC.25 = mean(CC.Q25),
                   CC.75 = mean(CC.Q75))
res.plot.manta = res.plot[which(res.plot$REEF_NAME %in% manta.reefs),]
res.plot.manta.all = res.plot[which(res.plot$REEF_NAME %in% manta.reefs.all),]
res.plot.manta.valid = res.plot[which(res.plot$REEF_NAME %in% manta.reefs.valid),]

data.manta.valid = data.manta %>% filter(REEF_NAME.y %in% manta.reefs & REPORT_YEAR > 1995) %>% 
  mutate(REEF_NAME = REEF_NAME.y,
         Year = REPORT_YEAR - 1) %>% 
  right_join(res.plot.manta, by=c("REEF_NAME", "Year"))

data.manta.valid.all = data.manta %>% filter(REEF_NAME.y %in% manta.reefs.all & REPORT_YEAR > 1995) %>% 
  mutate(REEF_NAME = REEF_NAME.y,
         Year = REPORT_YEAR - 1) %>% 
  right_join(res.plot.manta.all, by=c("REEF_NAME", "Year"))

data.manta.valid.valid = data.manta %>% filter(REEF_NAME.y %in% manta.reefs.all & REPORT_YEAR > 1995) %>% 
  mutate(REEF_NAME = REEF_NAME.y,
         Year = REPORT_YEAR - 1) %>% 
  right_join(res.plot.manta.valid, by=c("REEF_NAME", "Year"))


SECTORS = unique(as.character(data.manta.valid$SECTOR))
MPE = MPE.valid = data.frame(SECTOR= rep(SECTORS, each=NREPS), 
                 REP = 1:NREPS,
                 MPE.CC = vector(mode = "numeric", NREPS*length(SECTORS)),
                 MPE.CC.c = vector(mode = "numeric", NREPS*length(SECTORS)),
                 MPE.CC.R2 = vector(mode = "numeric", NREPS*length(SECTORS)),
                 MPE.CC.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 MPE.CC.ALL.c = vector(mode = "numeric", NREPS*length(SECTORS)),
                 MPE.CC.ALL.R2 = vector(mode = "numeric", NREPS*length(SECTORS)),
                 MPE.COTS = vector(mode = "numeric", NREPS*length(SECTORS)),
                 MPE.COTS.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 ACC.Pres = vector(mode = "numeric", NREPS*length(SECTORS)),
                 ACC.Pres.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 KAP.Pres = vector(mode = "numeric", NREPS*length(SECTORS)),
                 KAP.Pres.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 P.Pres = vector(mode = "numeric", NREPS*length(SECTORS)),
                 P.Pres.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 KAP.Out = vector(mode = "numeric", NREPS*length(SECTORS)),
                 KAP.Out.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 P.Out = vector(mode = "numeric", NREPS*length(SECTORS)),
                 P.Out.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 ACC.Out = vector(mode = "numeric", NREPS*length(SECTORS)),
                 ACC.Out.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 pR2.COTS = vector(mode = "numeric", NREPS*length(SECTORS)),
                 pR2.COTS.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 pR2sqrt.COTS = vector(mode = "numeric", NREPS*length(SECTORS)),
                 pR2sqrt.COTS.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 OUT.R2 = vector(mode = "numeric", NREPS*length(SECTORS)),
                 OUT.C = vector(mode = "numeric", NREPS*length(SECTORS)),
                 ACC.CLASS = vector(mode = "numeric", NREPS*length(SECTORS)),
                 KAP.CLASS = vector(mode = "numeric", NREPS*length(SECTORS)),
                 P.CLASS = vector(mode = "numeric", NREPS*length(SECTORS)))
                 

# Calculate Calibration Metrics ----
#Compute validation table at sector level ----
for (i in 1:NREPS) {
  # browser()
  print(i)
  manta.i = dplyr::filter(data.manta.valid, REP==i)
  lm = lm(CC.mn~MEAN_LIVE.corr, data=manta.i)
  MPE[MPE$REP==i,"MPE.CC.ALL"] = mean(abs(summary(lm)$residuals)) ## Mean prediction error = 7.7 %
  MPE[MPE$REP==i,"MPE.CC.ALL.c"] = lm$coefficients[2]
  MPE[MPE$REP==i,"MPE.CC.ALL.R2"] = summary(lm)$r.squared
  lm = lm(COTS.mn~MEAN_COTS, data=manta.i)
  summary(lm)$r.squared  ## R2 = 0.637
  MPE[MPE$REP==i,"MPE.COTS.ALL"] = mean(abs(summary(lm)$residuals))
  # make confusion matrix --> get KAPPA and Accuracy
  obs.bin = factor(ifelse(manta.i$MEAN_COTS>0.01,1,0),levels = c(0,1)) ; pred.bin = factor(ifelse(manta.i$COTS.mn > 0.01, 1,0),levels = c(0,1))
  obs.binOUT = factor(ifelse(manta.i$MEAN_COTS >0.11,1,0), levels = c(0,1)) ; pred.binOUT = factor(ifelse(manta.i$COTS.mn > 0.22, 1,0),levels = c(0,1))
  pres.conf = caret::confusionMatrix(pred.bin, obs.bin)
  out.conf = caret::confusionMatrix(pred.binOUT, obs.binOUT) #### Confusion Matrix
  # Multi-Class confusion matrix
  obs.clas = factor(ifelse(manta.i$MEAN_COTS<0.01,"NC",
                        ifelse(manta.i$MEAN_COTS<0.11,"NO",
                                  ifelse(manta.i$MEAN_COTS<0.22,"PO",
                                         ifelse(manta.i$MEAN_COTS<1,"EO", "SO")))),
                                         levels = c("NC", "PO", "EO", "SO")) 
  pred.clas = factor(ifelse(manta.i$COTS.mn<0.01,"NC",
                        ifelse(manta.i$COTS.mn<0.11,"NO",
                                  ifelse(manta.i$COTS.mn<0.22,"PO",
                                         ifelse(manta.i$COTS.mn<1,"EO", "SO")))),
                    levels = c("NC", "PO", "EO", "SO")) 
  clas.conf = caret::confusionMatrix(pred.clas, obs.clas)
  MPE[MPE$REP==i,"ACC.CLASS"] = clas.conf$overall[1]
  MPE[MPE$REP==i,"KAP.CLASS"] = clas.conf$overall[2]
  MPE[MPE$REP==i,"P.CLASS"] = clas.conf$overall[6]
  MPE[MPE$REP==i,c("ACC.Pres.ALL", "KAP.Pres.ALL", "P.Pres.ALL")] = round(rep(pres.conf$overall[c(1:2, 6)], each=length(SECTORS)),3)
  MPE[MPE$REP==i,c("ACC.Out.ALL", "KAP.Out.ALL", "P.Out.ALL")] = round(rep(out.conf$overall[c(1:2, 6)], each=length(SECTORS)),3)
  df = manta.i %>%  
    dplyr::select(REEF_NAME, Year, REP, COTS.mn, MEAN_COTS) %>% 
    mutate(COTS.mn = predict(MTCalib.gam, newdata=data.frame(MT=COTS.mn)),
           MEAN_COTS = predict(MTCalib.gam, newdata=data.frame(MT=MEAN_COTS)),
           pred.bin = pred.bin) %>%
    mutate(COTS.mn = ifelse(COTS.mn <=0.01, 0, COTS.mn),
           MEAN_COTS = round(ifelse(MEAN_COTS <0, 0, MEAN_COTS),0))
  df.out = manta.i %>%  
    dplyr::select(REEF_NAME, SECTOR, Year, REP, COTS.mn, MEAN_COTS) %>% 
    dplyr::group_by(SECTOR,REEF_NAME, REP) %>% dplyr::filter(MEAN_COTS > 0.11) %>%
      dplyr::summarise(MaxMEAN = mean(MEAN_COTS))
  df.out.mod = manta.i %>%  
    dplyr::select(REEF_NAME, Year, REP, COTS.mn, MEAN_COTS) %>% 
    dplyr::group_by(REEF_NAME, REP) %>% dplyr::filter(REEF_NAME %in% df.out$REEF_NAME) %>%
    dplyr::filter(COTS.mn > 0.11) %>%
    dplyr::summarise(MaxMEAN.mod = mean(COTS.mn)) %>% dplyr::right_join(df.out)
  lm = lm(MaxMEAN~MaxMEAN.mod, data=df.out.mod)
  summary(lm)$r.squared  ## R2 = 0.637
  MPE[MPE$REP==i,"OUT.R2"] = summary(lm)$r.squared
  MPE[MPE$REP==i,"OUT.C"] = lm$coefficients[2]
  # run Zeron-infl neg bin on results using pscl to get pseudo R2 --> convert to integer value
  try({
    m1 = pscl::zeroinfl(MEAN_COTS ~ COTS.mn | pred.bin,
                        data = df, dist = "negbin", EM = TRUE)
    m2 = pscl::zeroinfl(round(sqrt(MEAN_COTS),0) ~ round(sqrt(COTS.mn),0) | pred.bin,
                      data = df, dist = "negbin", EM = TRUE)
    MPE[MPE$REP==i,"pR2.COTS.ALL"] = round(r2_zeroinflated(m1)$R2,3)
    MPE[MPE$REP==i,"pR2sqrt.COTS.ALL"] = round(r2_zeroinflated(m2)$R2,3)
  }, silent=T)
  
  
    for (j in 1:length(SECTORS)) {
    # browser()
    print(j)
    manta.i = dplyr::filter(data.manta.valid, REP==i & SECTOR == SECTORS[j])
    lm = lm(CC.mn~MEAN_LIVE.corr, data=manta.i)
    summary(lm)$r.squared  ## R2 = 0.637
    MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,"MPE.CC"] = mean(abs(summary(lm)$residuals)) ## Mean prediction error = 7.7 %
    MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,"MPE.CC.c"] = lm$coefficients[2]
    MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,"MPE.CC.R2"] = summary(lm)$r.squared
    lm = lm(COTS.mn~MEAN_COTS, data=manta.i)
    summary(lm)$r.squared  ## R2 = 0.637
    MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,"MPE.COTS"] = mean(abs(summary(lm)$residuals)) ##
    # make confusion matrix --> get KAPPA and Accuracy
    obs.bin = factor(ifelse(manta.i$MEAN_COTS>0.01,1,0),levels = c(0,1)) ; pred.bin = factor(ifelse(manta.i$COTS.mn > 0.01, 1,0),levels = c(0,1))
    obs.binOUT = factor(ifelse(manta.i$MEAN_COTS >0.11,1,0), levels = c(0,1)) ; pred.binOUT = factor(ifelse(manta.i$COTS.mn > 0.22, 1,0),levels = c(0,1))
    pres.conf = caret::confusionMatrix(pred.bin, obs.bin)
    out.conf = caret::confusionMatrix(pred.binOUT, obs.binOUT)#### Confusion Matrix
    MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,c("ACC.Pres", "KAP.Pres", "P.Pres")] = pres.conf$overall[c(1:2, 6)]
    MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,c("ACC.Out", "KAP.Out", "P.Out")] = out.conf$overall[c(1:2, 6)]
    df = manta.i %>%  
      dplyr::select(REEF_NAME, Year, REP, COTS.mn, MEAN_COTS) %>% 
      mutate(COTS.mn = predict(MTCalib.gam, newdata=data.frame(MT=COTS.mn)),
             MEAN_COTS = predict(MTCalib.gam, newdata=data.frame(MT=MEAN_COTS)),
             pred.bin = pred.bin) %>%
      mutate(COTS.mn = ifelse(COTS.mn <=0.01, 0, COTS.mn),
             MEAN_COTS = round(ifelse(MEAN_COTS <0, 0, MEAN_COTS),0))
    
    # run Zeron-infl neg bin on results using pscl to get pseudo R2 --> convert to integer value
    try({
      m1 = pscl::zeroinfl(MEAN_COTS ~ COTS.mn | pred.bin,
                        data = df, dist = "negbin", EM = TRUE)
      m2 = pscl::zeroinfl(round(sqrt(MEAN_COTS),0) ~ round(sqrt(COTS.mn),0) | pred.bin,
                        data = df, dist = "negbin", EM = TRUE)
      MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,"pR2.COTS"] = round(r2_zeroinflated(m1)$R2,3)
      MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,"pR2sqrt.COTS"] = round(r2_zeroinflated(m2)$R2,3)
    }, silent = T)
    
  }
}  

# Calculate Validation Metrics ----
#Compute validation table at sector level ----
for (i in 1:NREPS) {
  # browser()
  print(i)
  manta.i = dplyr::filter(data.manta.valid.valid, REP==i)
  lm = lm(CC.mn~MEAN_LIVE.corr, data=manta.i)
  MPE.valid[MPE.valid$REP==i,"MPE.CC.ALL"] = mean(abs(summary(lm)$residuals)) ## Mean prediction error = 7.7 %
  MPE.valid[MPE.valid$REP==i,"MPE.CC.ALL.c"] = lm$coefficients[2]
  MPE.valid[MPE.valid$REP==i,"MPE.CC.ALL.R2"] = summary(lm)$r.squared
  lm = lm(COTS.mn~MEAN_COTS, data=manta.i)
  summary(lm)$r.squared  ## R2 = 0.637
  MPE.valid[MPE.valid$REP==i,"MPE.COTS.ALL"] = mean(abs(summary(lm)$residuals))
  # make confusion matrix --> get KAPPA and Accuracy
  obs.bin = factor(ifelse(manta.i$MEAN_COTS>0.01,1,0),levels = c(0,1)) ; pred.bin = factor(ifelse(manta.i$COTS.mn > 0.01, 1,0),levels = c(0,1))
  obs.binOUT = factor(ifelse(manta.i$MEAN_COTS >0.11,1,0), levels = c(0,1)) ; pred.binOUT = factor(ifelse(manta.i$COTS.mn > 0.22, 1,0),levels = c(0,1))
  pres.conf = caret::confusionMatrix(pred.bin, obs.bin)
  out.conf = caret::confusionMatrix(pred.binOUT, obs.binOUT) #### Confusion Matrix
  # Multi-Class confusion matrix
  obs.clas = factor(ifelse(manta.i$MEAN_COTS<0.01,"NC",
                         ifelse(manta.i$MEAN_COTS<0.11,"NO",
                           ifelse(manta.i$MEAN_COTS<0.22,"PO",
                                  ifelse(manta.i$MEAN_COTS<1,"EO", "SO")))),
                    levels = c("NC", "PO", "EO", "SO")) 
  pred.clas = factor(ifelse(manta.i$COTS.mn<0.01,"NC",
                        ifelse(manta.i$COTS.mn<0.11,"NO",
                            ifelse(manta.i$COTS.mn<0.22,"PO",
                                   ifelse(manta.i$COTS.mn<1,"EO", "SO")))),
                     levels = c("NC", "PO", "EO", "SO")) 
  clas.conf = caret::confusionMatrix(pred.clas, obs.clas)
  MPE.valid[MPE.valid$REP==i,"ACC.CLASS"] = clas.conf$overall[1]
  MPE.valid[MPE.valid$REP==i,"KAP.CLASS"] = clas.conf$overall[2]
  MPE.valid[MPE.valid$REP==i,"P.CLASS"] = clas.conf$overall[6]
  MPE.valid[MPE.valid$REP==i,c("ACC.Pres.ALL", "KAP.Pres.ALL", "P.Pres.ALL")] = round(rep(pres.conf$overall[c(1:2, 6)], each=length(SECTORS)),3)
  MPE.valid[MPE.valid$REP==i,c("ACC.Out.ALL", "KAP.Out.ALL", "P.Out.ALL")] = round(rep(out.conf$overall[c(1:2, 6)], each=length(SECTORS)),3)
  df = manta.i %>%  
    dplyr::select(REEF_NAME, Year, REP, COTS.mn, MEAN_COTS) %>% 
    mutate(COTS.mn = predict(MTCalib.gam, newdata=data.frame(MT=COTS.mn)),
           MEAN_COTS = predict(MTCalib.gam, newdata=data.frame(MT=MEAN_COTS)),
           pred.bin = pred.bin) %>%
    mutate(COTS.mn = ifelse(COTS.mn <=0.01, 0, COTS.mn),
           MEAN_COTS = round(ifelse(MEAN_COTS <0, 0, MEAN_COTS),0))
  df.out = manta.i %>%  
    dplyr::select(REEF_NAME, SECTOR, Year, REP, COTS.mn, MEAN_COTS) %>% 
    dplyr::group_by(SECTOR,REEF_NAME, REP) %>% dplyr::filter(MEAN_COTS > 0.11) %>%
    dplyr::summarise(MaxMEAN = mean(MEAN_COTS))
  df.out.mod = manta.i %>%  
    dplyr::select(REEF_NAME, Year, REP, COTS.mn, MEAN_COTS) %>% 
    dplyr::group_by(REEF_NAME, REP) %>% dplyr::filter(REEF_NAME %in% df.out$REEF_NAME) %>%
    dplyr::filter(COTS.mn > 0.11) %>%
    dplyr::summarise(MaxMEAN.mod = mean(COTS.mn)) %>% dplyr::right_join(df.out)
  lm = lm(MaxMEAN~MaxMEAN.mod, data=df.out.mod)
  summary(lm)$r.squared  ## R2 = 0.637
  MPE.valid[MPE.valid$REP==i,"OUT.R2"] = summary(lm)$r.squared
  MPE.valid[MPE.valid$REP==i,"OUT.C"] = lm$coefficients[2]
  # run Zeron-infl neg bin on results using pscl to get pseudo R2 --> convert to integer value
  try({
    m1 = pscl::zeroinfl(MEAN_COTS ~ COTS.mn | pred.bin,
                        data = df, dist = "negbin", EM = TRUE)
    m2 = pscl::zeroinfl(round(sqrt(MEAN_COTS),0) ~ round(sqrt(COTS.mn),0) | pred.bin,
                        data = df, dist = "negbin", EM = TRUE)
    MPE.valid[MPE.valid$REP==i,"pR2.COTS.ALL"] = round(r2_zeroinflated(m1)$R2,3)
    MPE.valid[MPE.valid$REP==i,"pR2sqrt.COTS.ALL"] = round(r2_zeroinflated(m2)$R2,3)
  }, silent=T)
  
  
  for (j in 1:length(SECTORS)) {
    # browser()
    print(j)
    manta.i = dplyr::filter(data.manta.valid.valid, REP==i & SECTOR == SECTORS[j])
    lm = lm(CC.mn~MEAN_LIVE.corr, data=manta.i)
    summary(lm)$r.squared  ## R2 = 0.637
    MPE.valid[MPE.valid$SECTOR==SECTORS[j] & MPE.valid$REP==i,"MPE.CC"] = mean(abs(summary(lm)$residuals)) ## Mean prediction error = 7.7 %
    MPE.valid[MPE.valid$SECTOR==SECTORS[j] & MPE.valid$REP==i,"MPE.CC.c"] = lm$coefficients[2]
    MPE.valid[MPE.valid$SECTOR==SECTORS[j] & MPE.valid$REP==i,"MPE.CC.R2"] = summary(lm)$r.squared
    lm = lm(COTS.mn~MEAN_COTS, data=manta.i)
    summary(lm)$r.squared  ## R2 = 0.637
    MPE.valid[MPE.valid$SECTOR==SECTORS[j] & MPE.valid$REP==i,"MPE.COTS"] = mean(abs(summary(lm)$residuals)) ##
    # make confusion matrix --> get KAPPA and Accuracy
    obs.bin = factor(ifelse(manta.i$MEAN_COTS>0.01,1,0),levels = c(0,1)) ; pred.bin = factor(ifelse(manta.i$COTS.mn > 0.01, 1,0),levels = c(0,1))
    obs.binOUT = factor(ifelse(manta.i$MEAN_COTS >0.11,1,0), levels = c(0,1)) ; pred.binOUT = factor(ifelse(manta.i$COTS.mn > 0.22, 1,0),levels = c(0,1))
    pres.conf = caret::confusionMatrix(pred.bin, obs.bin)
    out.conf = caret::confusionMatrix(pred.binOUT, obs.binOUT)#### Confusion Matrix
    MPE.valid[MPE.valid$SECTOR==SECTORS[j] & MPE.valid$REP==i,c("ACC.Pres", "KAP.Pres", "P.Pres")] = pres.conf$overall[c(1:2, 6)]
    MPE.valid[MPE.valid$SECTOR==SECTORS[j] & MPE.valid$REP==i,c("ACC.Out", "KAP.Out", "P.Out")] = out.conf$overall[c(1:2, 6)]
    df = manta.i %>%  
      dplyr::select(REEF_NAME, Year, REP, COTS.mn, MEAN_COTS) %>% 
      mutate(COTS.mn = predict(MTCalib.gam, newdata=data.frame(MT=COTS.mn)),
             MEAN_COTS = predict(MTCalib.gam, newdata=data.frame(MT=MEAN_COTS)),
             pred.bin = pred.bin) %>%
      mutate(COTS.mn = ifelse(COTS.mn <=0.01, 0, COTS.mn),
             MEAN_COTS = round(ifelse(MEAN_COTS <0, 0, MEAN_COTS),0))
    
    # run Zeron-infl neg bin on results using pscl to get pseudo R2 --> convert to integer value
    try({
      m1 = pscl::zeroinfl(MEAN_COTS ~ COTS.mn | pred.bin,
                          data = df, dist = "negbin", EM = TRUE)
      m2 = pscl::zeroinfl(round(sqrt(MEAN_COTS),0) ~ round(sqrt(COTS.mn),0) | pred.bin,
                          data = df, dist = "negbin", EM = TRUE)
      MPE.valid[MPE.valid$SECTOR==SECTORS[j] & MPE.valid$REP==i,"pR2.COTS"] = round(r2_zeroinflated(m1)$R2,3)
      MPE.valid[MPE.valid$SECTOR==SECTORS[j] & MPE.valid$REP==i,"pR2sqrt.COTS"] = round(r2_zeroinflated(m2)$R2,3)
    }, silent = T)
    
  }
}  



range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
MPE.overall = MPE %>% group_by(REP) %>% 
  dplyr::select(contains("ALL"), 27:31 ) %>%
  # dplyr::select(MPE.CC.ALL, MPE.COTS.ALL, ACC.Pres.ALL, KAP.Pres.ALL, KAP.Out.ALL, 
  #               ACC.Out.ALL, pR2sqrt.COTS.ALL,P.Pres.ALL, P.Out.ALL) %>%
  summarise_all(mean) %>%
  mutate(Overall = (range01(ACC.Pres.ALL) + range01(KAP.Pres.ALL) + range01(ACC.CLASS) +
                      range01(OUT.R2) +
                      range01(ACC.Out.ALL) +range01(KAP.Out.ALL) + (1-range01(MPE.CC.ALL)))/7)

MPE.valid.overall = MPE.valid %>% group_by(REP) %>% 
  dplyr::select(contains("ALL"), 27:31) %>%
  # dplyr::select(MPE.CC.ALL, MPE.COTS.ALL, ACC.Pres.ALL, KAP.Pres.ALL, KAP.Out.ALL, 
  #               ACC.Out.ALL, pR2sqrt.COTS.ALL,P.Pres.ALL, P.Out.ALL) %>%
  summarise_all(mean) %>%
  mutate(Overall = (range01(ACC.Pres.ALL) + range01(KAP.Pres.ALL) + range01(ACC.CLASS) +
                      range01(OUT.R2) +
                      range01(ACC.Out.ALL) +range01(KAP.Out.ALL) + (1-range01(MPE.CC.ALL)))/7)
write.csv(MPE.overall, file = "Results/Archive/MPECalib_2019-08-29_5.csv", row.names = F)
write.csv(MPE.valid.overall, file = "Results/Archive/MPEValid_2019-08-29_5.csv", row.names = F)

BESTREPS = unique(MPE.overall[order(MPE.overall$ACC.Out.ALL, decreasing = T),"REP"])[1:10,1] %>% pull(REP)
# MPE.BEST = MPE %>% filter(REP %in% BESTREPS)
# BESTREPS.CC = unique(MPE[order(MPE$MPE.CC.ALL,decreasing = T),"REP"])[1:5]
# BESTPARAMS = masterDF[BESTREPS,]
# BESTPARAMS.CC = masterDF[BESTREPS.CC,]
# MPE.T = MPE[MPE$SECTOR=="TO",]
# BESTREPS.T = MPE.T[order(MPE.T$ACC.Out)[1:5],"REP"]
# BESTPARAMS.T = masterDF[BESTREPS.T,]


# Plot Validation Reefs----
GG.COTS = ggplot(data.manta.valid %>% filter(REP %in% 2:4 & REEF_NAME %in% manta.reefs), aes(x=Year, y=COTS.mn)) +
  geom_point(aes(y=MEAN_COTS)) + geom_line(aes(x=Year, y=MEAN_COTS)) +
  geom_line(aes(colour=as.factor(REP))) +  geom_ribbon(aes(ymin = COTS.25, ymax=COTS.75, fill=factor(REP)),alpha=0.2) +
  facet_wrap(~REEF_NAME, scales = "free_y", ncol=3 )
GG.COTS
GG.CORAL = ggplot(data.manta.valid  %>% filter(REP %in% 2:4 & REEF_NAME %in% manta.reefs), aes(x=Year, y=CC.mn)) +
    geom_point(aes(y=MEAN_LIVE.corr)) + geom_line(aes(x=REPORT_YEAR-1, y=MEAN_LIVE.corr)) +
    geom_line(aes(colour=as.factor(REP))) + geom_ribbon(aes(ymin = CC.25, ymax=CC.75, fill=factor(REP)),alpha=0.2) +
    facet_wrap(~SECTOR+REEF_NAME, scales = "free_y") + ylim(c(0,50))
GG.CORAL


cl = parallel::makeCluster(2)
doParallel::registerDoParallel(cores = 20)

foreach::foreach (i = 1:NREPS) %dopar% {
# # dir.create("Results/Figures/Replicates")
# for (i in 1:NREPS) {
  # browser()
  `%>%` <- magrittr::`%>%`
  library(ggplot2)
  library(dplyr)
  gg=ggplot(data.manta.valid %>% filter(REP == i & REEF_NAME %in% manta.reefs), aes(x=Year, y=COTS.mn)) +
    geom_point(aes(y=MEAN_COTS)) + geom_line(aes(x=Year, y=MEAN_COTS)) +
    geom_line(aes(colour=as.factor(REP))) +  geom_ribbon(aes(ymin = COTS.25, ymax=COTS.75, fill=factor(REP)),alpha=0.2) +
    facet_wrap(~SECTOR+REEF_NAME, scales = "free_y")
  ggsave(paste0("Results/Figures/Replicates/COTS_",i,"_",Sys.Date(),".png"),gg, device="png", width=14, height=8, dpi = 300)

# }
#
# for (i in 1:NREPS) {
  gg = ggplot(data.manta.valid  %>% filter(REP ==i & REEF_NAME %in% manta.reefs), aes(x=Year, y=CC.mn)) +
    geom_point(aes(y=MEAN_LIVE.corr)) + geom_line(aes(x=REPORT_YEAR-1, y=MEAN_LIVE.corr)) +
    geom_line(aes(colour=as.factor(REP))) + geom_ribbon(aes(ymin = CC.25, ymax=CC.75, fill=factor(REP)),alpha=0.2) +
    facet_wrap(~SECTOR+REEF_NAME, scales = "free_y") + ylim(c(0,50))
  ggsave(paste0("Results/Figures/Replicates/Coral_",i,"_",Sys.Date(),".png"),gg, device="png", width=14, height=8, dpi = 300)
}
parallel::stopCluster(cl)
# GG.COTS.T = ggplot(data.manta.valid %>% filter(REP %in% BESTREPS.T & REEF_NAME %in% manta.reefs), aes(x=Year, y=COTS.mn)) + 
#   geom_point(aes(y=MEAN_COTS)) + geom_line(aes(x=Year, y=MEAN_COTS)) +
#   geom_line(aes(colour=as.factor(REP))) +  geom_ribbon(aes(ymin = COTS.25, ymax=COTS.75, fill=factor(REP)),alpha=0.2) +
#   facet_wrap(~SECTOR+REEF_NAME,scales = "free_y")
# GG.COTS.T
# 
# 
# # Check SECTORS
# ggplot(data.manta.valid %>% 
#          filter(REP %in% 6:10 & REEF_NAME %in% manta.reefs & A_SECTOR %in% c("TO")), 
#        aes(x=Year, y=COTS.mn)) + 
#   geom_point(aes(y=MEAN_COTS)) + geom_line(aes(x=Year, y=MEAN_COTS)) +
#   geom_line(aes(colour=as.factor(REP))) +  geom_ribbon(aes(ymin = COTS.25, ymax=COTS.75, fill=factor(REP)),alpha=0.2) +
#   facet_grid(cols=vars(REEF_NAME), rows = vars(REP),scales = "free_y") + ylim(c(0,5))
# ggplot(data.manta.valid %>% 
#          filter(REP %in% 1:5 & REEF_NAME %in% manta.reefs & A_SECTOR %in% c("TO")), 
#        aes(x=Year, y=CC.mn)) + 
#   geom_point(aes(y=MEAN_LIVE.corr)) + geom_line(aes(x=REPORT_YEAR-1, y=MEAN_LIVE.corr)) +
#   geom_line(aes(colour=as.factor(REP))) + geom_ribbon(aes(ymin = CC.25, ymax=CC.75, fill=factor(REP)),alpha=0.2) +
#   facet_wrap(~SECTOR+REEF_NAME, scales = "free_y") + ylim(c(0,50))
# 
# ggplot(data.manta.valid.all %>% 
#          filter(REP %in% 16:20  & SECTOR %in% c("SW")), 
#        aes(x=Year, y=COTS.mn)) + 
#   geom_point(aes(y=MEAN_COTS)) + geom_line(aes(x=Year, y=MEAN_COTS)) +
#   geom_line(aes(colour=as.factor(REP))) +  geom_ribbon(aes(ymin = COTS.25, ymax=COTS.75, fill=factor(REP)),alpha=0.2) +
#   facet_grid(cols=vars(REEF_NAME), rows = vars(REP),scales = "free_y") + ylim(c(0,5))
# 
# ggsave("Results/Figures/COTS.tiff",GG.COTS, device="tiff", width=14, height=8, dpi = 300)
# ggsave("Results/Figures/COTS.png",GG.COTS, device="png", width=14, height=8, dpi = 300)
# 
# GG.CORAL = ggplot(data.manta.valid  %>% filter(REP %in% 7:10 & REEF_NAME %in% manta.reefs), aes(x=Year, y=CC.mn)) +
#   geom_point(aes(y=MEAN_LIVE.corr)) + geom_line(aes(x=REPORT_YEAR-1, y=MEAN_LIVE.corr)) +
#   geom_line(aes(colour=as.factor(REP))) + geom_ribbon(aes(ymin = CC.25, ymax=CC.75, fill=factor(REP)),alpha=0.2) +
#   facet_wrap(~SECTOR+REEF_NAME, scales = "free_y") + ylim(c(0,50))
# GG.CORAL



# ggsave("Results/Figures/Coral.tiff", GG.CORAL, device="tiff", width=14, height=8, dpi = 300)
# ggsave("Results/Figures/Coral.png", GG.CORAL, device="png", width=14, height=8, dpi = 300)

# Best Model Average
# res.plot.manta.av = res.plot.manta  %>% dplyr::filter(REP %in% BESTREPS & REEF_NAME %in% manta.reefs) %>%
#   dplyr::group_by(REEF_NAME, Year, SECTOR, CROSS_SHELF) %>%
#   dplyr::summarise(COTS.md = median(COTS.mn),
#             COTS.25 = quantile(COTS.mn, probs = 0.25),
#             COTS.75 = quantile(COTS.mn, probs = 0.75),
#             CC.md =  median(CC.mn),
#             CC.25 = quantile(CC.mn, probs = 0.25),
#             CC.75 = quantile(CC.mn, probs = 0.75))

# data.manta.valid.av = data.manta %>% filter(REEF_NAME.y %in% manta.reefs & REPORT_YEAR > 1995) %>% 
#   mutate(REEF_NAME = REEF_NAME.y,
#          Year = REPORT_YEAR - 1) %>% 
#   inner_join(res.plot.manta.av)
# 
# GG.COTS.AV = ggplot(data.manta.valid.av, aes(x=Year, y=COTS.md)) + 
#   geom_point(aes(y=MEAN_COTS)) + geom_line(aes(x=Year, y=MEAN_COTS)) +
#   geom_line() +  geom_ribbon(aes(ymin = COTS.25, ymax=COTS.75),alpha=0.2) +
#   facet_wrap(~SECTOR+REEF_NAME,scales = "free_y")
# GG.COTS.AV

# Save PLots and Results ----
# ggsave("Results/Figures/COTS.tiff",GG.COTS.AV, device="tiff", width=14, height=8, dpi = 300)
# ggsave("Results/Figures/COTS.png",GG.COTS.AV, device="png", width=14, height=8, dpi = 300)

# GG.CORAL.AV = ggplot(data.manta.valid  %>% filter(REP %in% BESTREPS & REEF_NAME %in% manta.reefs), aes(x=Year, y=CC.mn)) + 
#   geom_point(aes(y=MEAN_LIVE.corr)) + geom_line(aes(x=REPORT_YEAR-1, y=MEAN_LIVE.corr)) +
#   geom_line(aes(colour=as.factor(REP))) +
#   facet_wrap(~SECTOR+REEF_NAME, scales = "free_y") + ylim(c(0,50))
# GG.CORAL.AV
# 
# ggsave("Results/Figures/Coral.tiff", GG.CORAL.AV, device="tiff", width=14, height=8, dpi = 300)
# ggsave("Results/Figures/Coral.png", GG.CORAL.AV, device="png", width=14, height=8, dpi = 300)

nruns = length(list.files(path = "Results/Archive")[grep(paste0("Validation_",Sys.Date()),list.files(path = "Archive"))]) + 1
# write.csv(ForDashboard, paste0("Results/Archive/ForDashboard_", Sys.Date(),"_", nruns, ".csv"))
# write.csv(ForDashboard, "Results/ResultsDashboard/ForDashboard.csv")
write.csv(masterDF, paste0("Results/Archive/masterDF_", Sys.Date(),"_", nruns, ".csv"))
write.csv(masterDF, "Results/Archive/masterDF.csv")

save(res.plot, MPE, MPE.overall,  data.manta.valid, data.manta.valid.all, masterDF,
     file = paste0("Results/Archive/Validation_", Sys.Date(),"_", nruns, ".RData"))


# Make Yearly COTS-Coral Graphs ----
source("COTSMod_PlotSpatial.R")
dir.create("Results/Figures/YearlyPlots")

# Set Up plotting Window ---
BESTREP=53
res.fig = res.plot %>% dplyr::filter(REP %in% BESTREP)
res.fig.GBR = res.fig %>% group_by(Year) %>%
  dplyr::summarise(CC.md = median(CC.mn),
            CC.Q05 = quantile(CC.mn, probs=0.05),
            CC.Q25 = quantile(CC.mn, probs=0.25),
            CC.Q75 = quantile(CC.mn, probs=0.75),
            CC.Q95 = quantile(CC.mn, probs=0.95),
            COTS.md = median(COTS.mn),
            COTS.Q05 = quantile(COTS.mn, probs=0.05),
            COTS.Q25 = quantile(COTS.mn, probs=0.25),
            COTS.Q75 = quantile(COTS.mn, probs=0.75),
            COTS.Q95 = quantile(COTS.mn, probs=0.95))

# cl = parallel::makeCluster(1)
# doParallel::registerDoParallel(cores = 3)

# foreach::foreach (i = 1996:2017) %dopar% {# Plot map (left panel)
for(i in 1996:2017){
  library(dplyr)
  library(ggplot2)
  library(PBSmapping)
  library(grid)
  i=2017
  res.fig.i = dplyr::select(data.grid, REEF_NAME, SECTOR, CROSS_SHELF, lat, lon) %>% inner_join(res.fig %>% filter(Year %in% i))

  # dev.off()
  f = colorRamp(c("yellow", "red"))
  layout(rbind(c(1,1,2,2),
                      c(1,1,2,2),
                      c(1,1,2,2),
                      c(3,3,3,4)))

  # PLOT CORAL COVER
  color = rgb(f(res.fig.i$CC.mn/60)/255)
  par(mai=c(0,0,0,0))
  plotPolys(shape45[1:1000,], col="gray95", border="gray70", xlim=c(148,157), ylim=c(-36, -21),
            bg="white", cex=1.2, axes=F, ylab=NA, plt=c(0,1,0,1), colHoles=NA, projection = 1)
  lines(gridln.elide, col="lightgrey")
  addPolys(shape45, col="gray95", border="gray70", bg="white", xlim=c(149,156), ylim=c(-36, -20),
           colHoles=NA)
  points(grid45$x, grid45$y, pch=19, cex=0.5, col=addalpha(color, .5))


  text(149.9,-22.5,i, cex=4, font=2)

  text(155.3, -22.5, "150E", col="lightgrey", srt=40, cex=1)
  text(155.3, -29.4, "155E", col="lightgrey", srt=40, cex=1)
  text(155.3, -26, "15S", col="lightgrey", srt=-40, cex=1)
  text(155.3, -33.2, "20S", col="lightgrey", srt=-40, cex=1)

  # pnts = cbind(x =c(149.3,150.2, 150.2,149.3), y =c(-35.6, -35.6,-34,-34))
  pnts = cbind(x =c(153.2,154.2, 154.2,153.2), y =c(-28.8, -28.8,-27,-27))
  # rect(149.2, -35.8, 151.1, -33.2, col="white", border = "black")
  # rect(153, -29.1, 155.3, -26.2, col="white", border = "black")
  color.grad = rgb(f(seq(0,1,0.01))/255)
  SDMTools::legend.gradient(pnts,color.grad,
                  c(0,60), title = "Coral Cover", cex=1.4)

  # Make GGPLOT FIGS
  res.plot.i = data.manta.valid %>%
    filter(REP %in% BESTREP &
           Year %in% 1996:i &
           REEF_NAME %in% c("MacGillivray Reef (14-114)", "Lady Musgrave Reef (23-082a)",
                            "Wardle Reef (17-032)"))

  GG.CORAL = ggplot(res.plot.i, aes(x=Year, y=CC.mn, fill=as.factor(REEF_NAME))) +
    geom_point(aes(y=MEAN_LIVE.corr)) + geom_line(aes(x=REPORT_YEAR-1, y=MEAN_LIVE.corr)) +
    geom_line(aes(colour=as.factor(REEF_NAME))) + geom_ribbon(aes(ymin = CC.25, ymax=CC.75,colour=as.factor(REEF_NAME)),alpha=0.4) +
    facet_wrap(~REEF_NAME, scales = "free_y", ncol = 1) + ylim(c(0,60)) + ylab("Mean % Coral Cover") +
    theme_bw(base_size=11) + theme(legend.position = "none") + xlim(c(1996,2017))

  GG.GBR = ggplot(res.fig.GBR %>% filter(Year %in% 1996:i), aes(x=Year, y=CC.md)) +
    geom_line(aes(x=Year, y=CC.md)) +
    geom_ribbon(aes(ymin = CC.Q05, ymax=CC.Q95),alpha=0.2) +
    geom_ribbon(aes(ymin = CC.Q25, ymax=CC.Q75),alpha=0.5) +
    ylab("Mean % Coral Cover") +
    theme_bw(base_size=11) + theme(legend.position = "none") + xlim(c(1996,2017))


  # GG.CORAL
  ## PLot GGPLOTs using viewport
  vp<-grid::viewport(x=0.7,y=0.6,width=0.6,height = 0.8)
  print(GG.CORAL, vp = vp)
  vp<-grid::viewport(x=0.4,y=0.1,width=0.72,height = 0.2)
  print(GG.GBR, vp = vp)
  vp<-grid::viewport(x=0.75,y=0.1,width=0.3,height = 0.2)
  print(grid::grid.text(as.character(i), x=0.85, y=0.15, gp=gpar(size=70, fontface="bold")), vp = vp)
  vp<-grid::viewport(x=0.75,y=0.1,width=0.3,height = 0.2)
  print(grid::grid.text("", x=0.85, y=0.1, gp=gpar(size=36, fontface="bold")), vp = vp)

  export::graph2tif(file = paste0("Results/Figures/YearlyPlots_",i), width =11, height=9, dpi=300)
}

# cl = parallel::makeCluster(1)
# doParallel::registerDoParallel(cores = 3)

# foreach::foreach (i = 1996:2017) %dopar% {
  library(dplyr)
  library(ggplot2)
  library(PBSmapping)
  library(grid)
# COTS Figs
for(i in 1996:2017){
  # i=2017
  res.fig.i = dplyr::select(data.grid, REEF_NAME, SECTOR, CROSS_SHELF, lat, lon) %>% inner_join(res.fig %>% filter(Year %in% i))

  # dev.off()
  f = colorRamp(c("green","yellow", "red"))
  layout(rbind(c(1,1,2,2),
               c(1,1,2,2),
               c(1,1,2,2),
               c(3,3,3,4)))

  # PLOT CORAL COVER
  # res.fig.i$COTS.mn2 = ifelse(res.fig.i$COTS.mn >1,1,res.fig.i$COTS.mn)
  # color = rgb(f(res.fig.i$COTS.mn/max(round(res.fig$COTS.mn,1)))/255)
  color = ifelse(res.fig.i$COTS.mn<0.01,"grey",
                 ifelse(res.fig.i$COTS.mn <0.11,"green",
                        ifelse(res.fig.i$COTS.mn <0.22, "yellow",
                               ifelse(res.fig.i$COTS.mn <0.5, "orange", "red"))))
  par(mai=c(0,0,0,0))
  plotPolys(shape45[1:1000,], col="gray95", border="gray70", xlim=c(148,157), ylim=c(-36, -21),
            bg="white", cex=1.2, axes=F, ylab=NA, plt=c(0,1,0,1), colHoles=NA, projection = 1)
  lines(gridln.elide, col="lightgrey")
  addPolys(shape45, col="gray95", border="gray70", bg="white", xlim=c(149,156), ylim=c(-36, -20),
           colHoles=NA)
  points(grid45$x, grid45$y, pch=19, cex=0.5, col=addalpha(color, .5))


  text(149.9,-22.5,i, cex=4, font=2)

  text(155.3, -22.5, "150E", col="lightgrey", srt=40, cex=1)
  text(155.3, -29.4, "155E", col="lightgrey", srt=40, cex=1)
  text(155.3, -26, "15S", col="lightgrey", srt=-40, cex=1)
  text(155.3, -33.2, "20S", col="lightgrey", srt=-40, cex=1)

  # pnts = cbind(x =c(149.3,150.2, 150.2,149.3), y =c(-35.6, -35.6,-34,-34))
  pnts = cbind(x =c(153.2,154.2, 154.2,153.2), y =c(-28.8, -28.8,-27,-27))
  # rect(149.2, -35.8, 151.1, -33.2, col="white", border = "black")
  # rect(153, -29.1, 155.3, -26.2, col="white", border = "black")
  # color.grad = rgb(f(seq(0,1,0.01))/255)
  # SDMTools::legend.gradient(pnts,color.grad,
  #                 c(0,max(round(res.fig$COTS.mn,1))), title = "COTS", cex=1.4)
  legend(x='bottomleft', legend = c("NC", "NO", "PO", "EO", "SO"), 
         fill = c("grey", "green", "yellow", "orange", "red"))

  # Make GGPLOT FIGS
  res.plot.i = data.manta.valid %>%
    filter(REP %in% BESTREP &
             Year %in% 1996:i &
             REEF_NAME %in% c("MacGillivray Reef (14-114)", "Lady Musgrave Reef (23-082a)",
                              "Wardle Reef (17-032)"))

  GG.COTS = ggplot(res.plot.i, aes(x=Year, y=COTS.mn, fill=as.factor(REEF_NAME))) +
    geom_point(aes(y=MEAN_COTS)) + geom_line(aes(x=REPORT_YEAR-1, y=MEAN_COTS)) +
    geom_line(aes(colour=as.factor(REEF_NAME))) + geom_ribbon(aes(ymin = COTS.25, ymax=COTS.75,colour=as.factor(REEF_NAME)),alpha=0.4) +
    facet_wrap(~REEF_NAME, scales = "free_y", ncol = 1)  + ylab("Mean COTS/Manta Tow") +
    theme_bw(base_size=11) + theme(legend.position = "none") + xlim(c(1996,2017))

  GG.GBR.COTS = ggplot(res.fig.GBR %>% filter(Year %in% 1996:i), aes(x=Year, y=COTS.md)) +
    geom_line(aes(x=Year, y=COTS.md)) +
    geom_ribbon(aes(ymin = COTS.Q05, ymax=COTS.Q95),alpha=0.2) +
    geom_ribbon(aes(ymin = COTS.Q25, ymax=COTS.Q75),alpha=0.5) +
    ylab("Mean COTS/Manta Tow") +
    theme_bw(base_size=11) + theme(legend.position = "none") + xlim(c(1996,2017))


  # GG.CORAL
  ## PLot GGPLOTs using viewport
  vp<-grid::viewport(x=0.7,y=0.6,width=0.6,height = 0.8)
  print(GG.COTS, vp = vp)
  vp<-grid::viewport(x=0.4,y=0.1,width=0.72,height = 0.2)
  print(GG.GBR.COTS, vp = vp)
  vp<-grid::viewport(x=0.75,y=0.1,width=0.3,height = 0.2)
  print(grid::grid.text(as.character(i), x=0.85, y=0.15, gp=gpar(size=70, fontface="bold")), vp = vp)
  vp<-grid::viewport(x=0.75,y=0.1,width=0.3,height = 0.2)
  print(grid::grid.text("", x=0.85, y=0.1, gp=gpar(size=36, fontface="bold")), vp = vp) # Put pred error here

  export::graph2tif(file = paste0("Results/Figures/YearlyPlotsCOTS_",i), width =11, height=9, dpi=300)
}


