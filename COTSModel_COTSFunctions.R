########################@
# Functions for modeling COTS demography and dispersal -----

#   Authors: Kevin Shoemaker, Sam Matthews
#
#   28 September 2016 - Updated initializeCoTSabund to matrix form
#                     - started scripting mortality, fecundity, dispersal functions
#
#   THINGS TO DO:
#         1. Add coral cover to dispersal probabilities
#         2. Incoroporate the [chl] conditions on the dispersal path
#         3. Write function to determine Carrying Capacity for Pred Prey Dynamics
#         4. DENSITY DEPENDENCE ON FECUNDITY AND TRANSITION



###################!
# BUILD DISPERSAL MATRIX ----
###################!

# setwd(DATA_DIRECTORY)
# Coords <- PopData[,c('lon','lat')]
# coordinates(Coords) <- ~lon+lat
# Gdist <- gDistance(Coords, Coords, byid = T)
# Gdist[1:10,1:10]
# # Convert to km
# Gdist <- Gdist*100
# # Limit the distance matrix to ~500km. i.e set >500km to NA
# Gdist[Gdist>500] <- 0
# Gdist.Sp <- Matrix::Matrix(Gdist, sparse=T)
# # Convert to probability matrix where probability of dispersal is inverse of distance
# Pdist = signif(1/Gdist.Sp,3)
# Pdist.Sp <- Matrix::Matrix(Pdist, sparse=T)

###################!
# BASIC VITAL RATE AND GROWTH FUNCTIONS ----
###################!

# converts diameter (in mm) to mass (in grams)

COTS_MassFromDiam <- function(Diam){
  Mass <- signif(6.29*10^-5*Diam^2.929,3)
  return(Mass)
}


# converts female mass to total larval fecundity 

COTS_FecFromMass <- function(Mass){
  Fec <- signif(558*Mass^1.439,4)
  return(Fec)
}


###################!

###################!
# initializeCOTSabund ----
###################!
# OBJECTIVE:
#    generate an object for storing the COTS abundance in each pixel. 
# PARAMS:
#    - reefmap: raster template for the study region: NA outside of reefs, reef ID value within reefs 
#    - PopData: data frame containiing PIXEL ID's, percent reef cover and environ vriable
#    - COTSInterp: txt file containing interpolated values of COTS Manta tow, giving a value for CoTS density
#    - Year: Which year we are using as our starting values
#    - Detectability: detectability of adult CoTS from MAnata tow surveys
#    - stagenames: vector of stagenames eg J1, J2, A1
#    - nstages: number of stages
#    - nreefs: number of reefs in simulation
#    - npops: number of separate populations - initially using every reef pixel as a reef
# RETURNS:
#    - COTSabund: spatially-structured and stage-structured COTS abundance
#           COTSabund$J_1: vector representing spatially structured abundance of Juvenile stage 1 individuals
#           COTSabund$J_2: vector representing spatially structured abundance of Juvenile stage 2 individuals
#           COTSabund$A: vector representing spatially structured abundance of reproductive adult individuals
#           COTSabund$S: vector representing spatially structured abundance of senile adult individuals
#           NOTE: larvae are not considered explicitly here. 
###################!

# De'ath Calibration
df1 = data.frame(MT = c(0,0.1,0.22,1), DENS = c(0,3500,4900,11000))

MTCalib.gam = lm(DENS~sqrt(MT), data=df1)
MTCalib.gaminv = lm(sqrt(MT)~DENS, data=df1)
# plot( df1$DENS, sqrt(df1$MT))
# predict(MTCalib.gam, newdata = data.frame(MT=0.06))
# predict(MTCalib.gaminv, newdata = data.frame(DENS=2638))^2
# ggplot(df1, aes(x=MT, y=DENS)) +geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=1))
initializeCOTSabund <- function(data.grid, data.COTS, Year, stagenames, 
                                COTS_StableStage, npops, j){
  # browser()
  ### set NA Values in interpolation CoTS.init to 0
  data.COTS[is.na(data.COTS)] <- 0
  
  nstages <- length(stagenames)
  
  ### Set up the COTS abundance object
  COTSabund <- matrix(0,nrow=npops, ncol=nstages)
  colnames(COTSabund) <- stagenames
  
  ### Set up reference for year
  colname <- paste('COTS_', Year, sep="")
  # browser()
  preds = predict(MTCalib.gam, newdata=data.frame(MT=COTS.rsmpl[,Year-1995, j]))
  COTSabund[,'A'] <- round(ifelse(preds<0.005,0,preds),0)
  ### Update abundances based from interpolated manta tow data
  # COTSabund[,'A'] <- round(data.COTS[,colname] * (666/0.7) * (data.grid$PercentReef/100),0)   # 666 converts manta tow to 1kmx1km
  # COTSabund[,'A'] <- round(COTS.rsmpl[,1997-Year, j] * (666/0.7) * (data.grid$PercentReef/100),0) 
  COTSabund[,'J_2'] <- round(COTSabund[,'A'] * as.numeric(COTS_StableStage[2]/COTS_StableStage[3]),0)
  COTSabund[,'J_1'] <- round(COTSabund[,'A'] * as.numeric(COTS_StableStage[1]/COTS_StableStage[3]),0)
  return(COTSabund)
}

# COTSabund <- initializeCOTSabund(data.grid, data.COTS, 1995, stagenames, 
#                                  COTS_StableStage, npops, 1)

    
###################!
# CoTS_StageTransition ----
###################!
# OBJECTIVE: Transition all individuals through life stages
#    
# PARAMS: 
#     - COTSabund: Matrix of CoTS abundances for which to transition
#     - COTSmort: named Vector of natural mortality rates for each stage
#     - COTSremain: named Vector of proportions of individuals that remain in current life stage
#    
# RETURNS: 
#     - newCOTSabund: COTS abund updated after 6 month time step
#     
###################!

COTS_StageTransition <- function(COTSabund, COTSmort, COTSremain) {
  
  #Set up matrices
  newCOTSabund <- matrix(0,nrow=nrow(COTSabund), ncol=ncol(COTSabund))
  colnames(newCOTSabund) <- colnames(COTSabund)
  COTS_Mort <- matrix(0,nrow=nrow(COTSabund), ncol=ncol(COTSabund))
  colnames(COTS_Mort) <- colnames(COTSabund)
  COTS_Remain <- matrix(0,nrow=nrow(COTSabund), ncol=ncol(COTSabund))
  colnames(COTS_Remain) <- colnames(COTSabund)
  COTS_Trans <- matrix(0,nrow=nrow(COTSabund), ncol=ncol(COTSabund))
  colnames(COTS_Trans) <- colnames(COTSabund)
  
  # apply mortality
  COTS_Mort <- sweep(COTSabund,MARGIN=2,COTSmort,`*`)
  # update abundance
  # newCOTSabund <- COTSabund - COTS_Mort
  newCOTSabund <- COTSabund # this removes the mortality from this function so it is all handled in the pred prey dynamics
  
  # number of COTS remaining and transitioning for each stage based on post-mortality abundaces
  COTS_Remain <- sweep(newCOTSabund,MARGIN=2,COTSremain,`*`)
  COTS_Trans <- sweep(newCOTSabund,MARGIN=2,1-COTSremain,`*`)
  
  # update newCOTSabund
  newCOTSabund[, 'J_1'] <- round(COTS_Remain[,'J_1'],0)
  newCOTSabund[, 'J_2'] <- round(COTS_Remain[,'J_2'] + COTS_Trans[,'J_1'],0)
  newCOTSabund[, 'A'] <- round(COTS_Remain[,'A'] + COTS_Trans[,'J_2'],0)
  return(newCOTSabund)
}

# COTSabund <- COTS_StageTransition(COTSabund = initCOTS, COTSmort = COTSmort, COTSremain = COTSremain)




####################!
# CoTS_Fecundity ----
####################!
# OBJECTIVE: 1. Assume a size distribution amongst adults and then sample diameters
#            2. Convert diameter to mass 
#            3. Convert mass to eggs
#            4. Include maternal condition 
#            
#    
# PARAMS: 
#     - COTSabund: standard abundance matrix for the time step
#     - mean: mean of the size distribution of COTS adults
#     - sd: standard deviation of size distribution
#     - npops: number of pixels
#     - mat.cond: maternal condition is a function of resources on source reef
# RETURNS:
#     - TotalLarvae: Total larvae produced per pixel
#     - mat.cond: maternal condition is a function of resources on source reef
#         to be used in the dispersal/settlment fucntion
#     
###################!

COTS_Fecundity <- function(COTSabund, mean, sd, npops, SR) {
    
  ### Intitialize matrix to store total eggs
  COTS_Eggs <- vector(mode = "numeric", length = npops)
    for (r in 1:npops) {
        Sizes <- rnorm(COTSabund[r,'A']*(10-SR)/10, mean, sd)
        Sizes[Sizes<0] = 0
        COTS_Eggs[r] <- sum(COTS_FecFromMass(COTS_MassFromDiam(Sizes)))
    }
    return(COTS_Eggs)
}

COTSPCF = function (npops, SexRatio) {
  COTSabund = matrix(NA, nrow=npops, ncol=3)
  COTSabund[,3] = seq(1,100000, length.out = npops)
  colnames(COTSabund) = stagenames
  nEggs <- COTS_Fecundity(COTSabund, 250, 70, SR = SexRatio, npops=npops)
  avgPCF = mean(nEggs[-1]/COTSabund[,'A'][-1])
  sdPCF = sd(nEggs[-1]/COTSabund[,'A'][-1])
  return(c(avgPCF,sdPCF))
}



COTS_Fecundity.PCF = function(COTSabund, PCFParams, SexRatio) {
  nEggs = COTSabund[,'A']*rnorm(1, PCFParams[1], PCFParams[2])*((10-SexRatio)/10)
  #nEggs.Pois = vector("numeric", length = npops)
  #for (i in 1:length(nEggs)){
   # nEggs.Pois[2] = rpois(1,1e10)
  #}
  return(nEggs)
}

# nEggs = COTS_Fecundity.PCF(COTSabund, PCFParams, SexRatio=5)

###################!
# CoTS_Fertilisation ----
###################!
# OBJECTIVE: Determine the fertilisation success(%) at different densities of COTS (between 0-150,000)
#    
# PARAMS:
#     - nLarvae: vector of number of larvae produced for each pixel
#     - SexRatio: integer from 1:9 indicating the male:female sex ratio ie
#     -   1 = 0.1M:0.9F
#    
# RETURNS:
#     - 
#     
###################!

COTS_Fertilisation = function(nEggs, nCOTS, SexRatio, FvDParams) {
  # Von Bertanlanffy Growth = Proportion of Eggs fertilised
  # browser()
  fEggs = FvDParams[SexRatio,"Linf"] * (1 - exp(-FvDParams[SexRatio,"K"] * (nCOTS - FvDParams[SexRatio,"t0"])))
  nLarvae = nEggs * fEggs
  return(nLarvae)
}

# nLarvae = COTS_Fertilisation(nEggs = nEggs, nCOTS = COTSabund[,'A'], SexRatio =5, FvDParams = FvDParams)
# range(nLarvae, na.rm = T)
##### PRODUCING NEGATIVE LARVAE



###################!


##############!
# CoTS_Dispersal ----
#############!  
# OBJECTIVE: Disperse larvae throughout the system based on the coralcover (maternal
#             condition), a connectivity matrix, and a vector containing self recruitment
# PARAMS: 
#   - ConnMat: 16035x16035 connectivity matrix
#   - SelfRecruit: 16035 vector containing estimates of self recruitment
#   - CoralCover: 16035 vector containing coral cover for the preceding year to determine
#               the maternal condition (needs to be looked into further)
# RETURNS: 
#   - LarveMax: 16035 vector containing the maximum larve received by each reef
  
# SelfRecruit = rep(0.1, npops)
# Pdist[Pdist==Inf] = 1
# Pdist.test = Pdist[1:1000, 1:1000]  

CoTS_Dispersal <- function(ConnMat, COTSabund, nLarvae, CoralCover, Pred, npops){
  # PsuedoCode:
  # 1. Loop through every population of fertilised eggs
  # 2. Discard 80% of eggs, then distribute the remainder as a function of distance
  # 3. Multiply nLarvae * connectivity probabilities to get larvae max
  # 5. Determine path of larvae to estimate survival (leave this for now)
  # 4. add vector of larvae produced to COTSabund
  # browser()
  nLarvae = ifelse(nLarvae < 0 , 0, nLarvae) ### FIx this
  nLarvae = as.matrix((1-Pred)*nLarvae) # assume 80% lost to the sea
  row.names(nLarvae) = data.grid$REEF_NAME
  nLarvae_Reef = rowsum(nLarvae, row.names(nLarvae), reorder = F) 
  nArriving_Reef = ConnMat
  nArriving_Reef[nArriving_Reef > 0] = 0
  
  for (i in 1:length(nLarvae_Reef)) {
    nArriving_Reef[i,] = signif(nLarvae_Reef[i]*(ConnMat[i,]),3)
  }
  nArriving_Reef = base::colSums(nArriving_Reef, na.rm = T) # create total juveniles arriving at a reef
  nArriving_Reef_PerPix = nArriving_Reef/Pixels$Pixels
  names(nArriving_Reef_PerPix) = colnames(ConnMat)
  match(data.grid$REEF_NAME, colnames(ConnMat))
  nArriving = nArriving_Reef_PerPix[match(data.grid$REEF_NAME, colnames(ConnMat))]
  COTSabund[,'J_1'] = COTSabund[,'J_1'] + nArriving # add these to our abundance
  # need to add a line that increases survival based on maternal resources
  return(COTSabund)
}

# COTSabund.t1 = CoTS_Dispersal(ConnMat = COTS.ConnMat, COTSabund = COTSabund, npops=npops,
#                               nLarvae=nLarvae, Pred=0.80)


##############!
# CoTS_Dispersal_Conn ----
#############!  
# OBJECTIVE: Disperse larvae throughout the system based on the connectivity matrix, 
# PARAMS: 
#   - ConnMat: 3806x3806 connectivity matrix
#   - COTSabund: 15928 x 3 matrix containing COTS abundances from current timestep
#   - CoralCover: 15928 vector containing coral cover for the current year to determine
#               the maternal condition (needs to be looked into further)
# RETURNS: 
#   - COTSabund: updated 15928 x 3 matrix containing COTS abundances from current timestep

# START HERE ----

CoTS_Dispersal_Conn <- function(COTSConnMat, COTSabund, nLarvae, CoralCover, Loss, npops){
  # PsuedoCode:
  # 1. Loop through every population of fertilised eggs
  # 2. Sum eggs up to reef level
  # 3. Multiply nLarvae * connectivity probabilities to get larvae ditributed to each reef
  # 5. Split up larve to each of the cells within the reef based on coral percentage
  # 4. add vector of larvae produced to COTSabund
  # browser()
  nLarvae = (1-Loss)*nLarvae # assume 80% lost to the sea
  nArriving = matrix(nrow=npops, ncol=npops)
  for (i in 1:npops) {
    nArriving[i,] = signif(nLarvae[i]*(ConnMat[i,]/sum(ConnMat[i,])),3)
  }
  nArriving = base::colSums(nArriving, na.rm = T) # create total juveniles arriving
  COTSabund[,'J_1'] = COTSabund[,'J_1'] + nArriving # add these to our abundance
  # need to add a line that increases survival based on maternal resources
  return(COTSabund)
}



##############!
# doCOTSDispersal ----
#############! 

# This is the master fucntion that incorporates the above functions

doCOTSDispersal = function(season, COTSabund, CoralCover, SexRatio, ConnMat, PCFParams, Pred, 
                           FvDParams, Fbase, CCRatioThresh, Year, data.chl, data.chl.resid, j, selfseed){
  # browser()
  COTSabund = COTSabund
  b = (1-Fbase)/CCRatioThresh # Make Ratio a parameter for tuning
  
  COTSMT = predict(MTCalib.gaminv, newdata=data.frame(DENS=COTSabund[,3]))^2
  COTSMT = ifelse(COTSMT <0.00001, 0,COTSMT)
  Ratio = (CoralCover*data.grid$PercentReef/100)/COTSMT # Do i need % Reef here?
  if (season=="summer"){
    # Per Capita Fecundity
    nEggs = COTSabund[,'A']*rnorm(1, PCFParams[1], PCFParams[2])*((10-SexRatio)/10)
    # Add in density dependent fecundity
    nEggs[which(Ratio<CCRatioThresh)] = (nEggs*(Fbase + (b*Ratio)))[which(Ratio<CCRatioThresh)]
    # Fertilisation by Density
    fEggs = FvDParams["Linf"] * (1 - exp(-FvDParams["K"] * (COTSabund[,'A'] - FvDParams["t0"])))
    nLarvae = round(nEggs * fEggs,0)
    # Do Dispersal ----
    nLarvae = ifelse(nLarvae < 0 , 0, nLarvae) ### FIx this
    nLarvae = round(as.matrix((1-Pred)*nLarvae),0) # assume 98% Larval Predation
    
    # Add in Chlorophyll relationship for survival
    chl = data.Chl[,as.character(Year)] + data.chl.resid[,(Year-1990+1),j]
    chl = ifelse(chl<0,0,chl)
    surv.pred = inv.logit(predict(chl.lm, newdata = data.frame(chl=chl)))
    
    nLarvae = nLarvae*surv.pred
    # Add in Larval Nutrition at home reef
    
    row.names(nLarvae) = data.grid$REEF_NAME
    nLarvae_Reef = rowsum(nLarvae, row.names(nLarvae), reorder = F) 
    nArriving_Reef = ConnMat
    nArriving_Reef[nArriving_Reef > 0] = 0
    
    for (i in 1:length(nLarvae_Reef)) {
      nArriving_Reef[i,] = round(nLarvae_Reef[i]*(ConnMat[i,]),0)
    }
    # browser()
    diag(nArriving_Reef) = diag(nArriving_Reef)*selfseed
    nArriving_Reef = base::colSums(nArriving_Reef, na.rm = T) # create total juveniles arriving at a reef
    nArriving_Reef_PerPix = round(nArriving_Reef/Pixels$Pixels,0)
    names(nArriving_Reef_PerPix) = colnames(ConnMat)
    # match(data.grid$REEF_NAME, colnames(ConnMat))
    nArriving = nArriving_Reef_PerPix[match(data.grid$REEF_NAME, colnames(ConnMat))]
    # Add in Larval Nurtrition at destination reef
    
    COTSabund[,'J_1'] = COTSabund[,'J_1'] + nArriving # add these to our abundance
   
  }
  return(COTSabund)
}


# test.winter = doCOTSDispersal("winter", initCOTS, SR=5, ConnMat = Pdist.test)

##############!
# doCOTSdemography ----
#############! 

# This is the master fucntion that incorporates the above functions
doCOTSDemography = function(season, COTSabund, COTSmort, COTSremain){
  newCOTSabund = COTSabund
  if (season=="winter"){ 
    #Set up matrices
    newCOTSabund <- matrix(0,nrow=nrow(COTSabund), ncol=ncol(COTSabund))
    colnames(newCOTSabund) <- colnames(COTSabund)
    COTS_Mort <- matrix(0,nrow=nrow(COTSabund), ncol=ncol(COTSabund))
    colnames(COTS_Mort) <- colnames(COTSabund)
    COTS_Remain <- matrix(0,nrow=nrow(COTSabund), ncol=ncol(COTSabund))
    colnames(COTS_Remain) <- colnames(COTSabund)
    COTS_Trans <- matrix(0,nrow=nrow(COTSabund), ncol=ncol(COTSabund))
    colnames(COTS_Trans) <- colnames(COTSabund)
    
    # apply mortality
    # COTS_Mort <- sweep(COTSabund,MARGIN=2,COTSmort,`*`)
    # update abundance
    newCOTSabund <- COTSabund - COTS_Mort
    
    # number of COTS remaining and transitioning for each stage based on post-mortality abundaces
    COTS_Remain <- sweep(newCOTSabund,MARGIN=2,COTSremain,`*`)
    COTS_Trans <- sweep(newCOTSabund,MARGIN=2,1-COTSremain,`*`)
    
    # update newCOTSabund
    newCOTSabund[, 'J_1'] <- round(COTS_Remain[,'J_1'],0)
    newCOTSabund[, 'J_2'] <- round(COTS_Remain[,'J_2'] + COTS_Trans[,'J_1'],0)
    newCOTSabund[, 'A'] <- round(COTS_Remain[,'A'] + COTS_Trans[,'J_2'],0)
  }
  return(newCOTSabund)
}

# test.summer = doCOTSDemography("summer", initCOTS)
# test.winter = doCOTSDemography("winter", initCOTS)





####################!
# doCoralConsumption ----
####################!
# OBJECTIVE:
#    - Determine the percentage Coral Cover consumed 
# PARAMS:
#    - season: "su,,er, or "winter" to determine the leve of coral consumption
#    - COTSabund: spatially-structured and stage-structured COTS abundance
#    - CoralCover: Spatially structured coral cover
# RETURNS:
#    - CoralCover: Spatially structured coral cover
###################!


doCoralConsumption = function(season, COTSabund, CoralCover, ConsRate, COTSfromCoralModel, Cbase, CMax) {
  if (COTSfromCoralModel==T) {
    CRemaining = CoralCover
    CChange = rep(0, length(CoralCover))
  } else {
  #browser() 
  COTSMT = predict(MTCalib.gaminv, newdata=data.frame(DENS=COTSabund[,3]))^2
  COTSMT = ifelse(COTSMT <0.00001, 0,COTSMT)
  Ratio = (CoralCover*data.grid$PercentReef/100)/COTSMT
  b = (1-Cbase)/CCRatioThresh
  #######STARRT HERE!!!!!!!!!!!!!!!!
  ConsRate.D = rep(CMax, length(CoralCover))
  ConsRate.D[which(Ratio<CCRatioThresh)] = (CMax*(Cbase + (b*Ratio)))[which(Ratio<CCRatioThresh)]
  if (season =="summer") {
    # browser()
    #CoralCover= Results[(Results$Year==year-1) & (Results$Season=="winter"),"CoralCover"]
    CAvailable = (CoralCover*data.grid$PercentReef/10000)*1e6*1e4 # in cm2
    # Do I want consumption to be density related?
    #CConsumed = ConsRate[,1]*COTSabund[,"A"]*182
    CConsumed = ConsRate.D*COTSabund[,"A"]*182
    CRemaining=((CAvailable-CConsumed)/1e10)*(10000/data.grid$PercentReef)
    CChange = CRemaining-CoralCover
    CRemaining[CRemaining < 0.5] <- 0.5
  } 
  if (season =="winter") {
    #CoralCover= Results[(Results$Year==year) & (Results$Season=="summer"),"CoralCover"]
    CAvailable = (CoralCover*data.grid$PercentReef/10000)*1e6*1e4 # in cm2
    # CConsumed = ConsRate[,2]*COTSabund[,"A"]*182
    CConsumed = ConsRate.D*COTSabund[,"A"]*182
    CRemaining=((CAvailable-CConsumed)/1e10)*(10000/data.grid$PercentReef)
    CChange = CRemaining-CoralCover
    CRemaining[CRemaining < 0.5] <- 0.5
  }
  }
  return(cbind(CRemaining, CChange))
}

# CoralCover = doCoralConsumption("winter", COTSabund = COTSabund, CoralCover = CoralCover[1:1000])


####################!
# doPredPreyDynamics ----
####################!
# OBJECTIVE:
#    - Restrict adult survival using Lotka Volterra Dynamics 
# PARAMS:
#    - season: "summer, or "winter" to determine the leve of coral consumption
#    - COTSabund: spatially-structured and stage-structured COTS abundance
#    - CoralCover: Spatially structured coral cover
# RETURNS:
#    - CoTS abund: spatially-structured and stage-structured COTS abundance
###################!

setCarryingCapacity = function(npops) {
  # browser()
  data.grid=data.grid[1:npops,]
  CC.10=rep(10,npops)
  # Calculate percent growth
  MinG.10 = doCoralGrowth(CC.10, B0=data.grid$pred.b0.mean[1:npops], WQ[1:npops], HC.asym = data.grid$pred.HCmax.mean[1:npops])[,2]
  # Number of COTS sustaianble at this level
  ConsRate = 250 # coral consumption rate
  CAvailable = (MinG.10*data.grid$PercentReef/10000)*1e6*1e4 # in cm2
  MinK.10A = CAvailable/(ConsRate*365)
  MinK.10J1 = MinK.10A*COTS_StableStage[1]/COTS_StableStage[3]
  MinK.10J2 = MinK.10A*COTS_StableStage[2]/COTS_StableStage[3]
  MinK.10 = as.matrix(cbind(MinK.10J1, MinK.10J2, MinK.10A)) ##THERE IS AN ISSUE WITH THE MAGNITUDE OF COTS CARRYING CAPACITY
  
  CC.0.5=rep(0.5,npops)
  # Calculate percent growth
  MinG.0.5 = doCoralGrowth(CC.0.5, B0=data.grid$pred.b0.mean[1:npops], WQ[1:npops], HC.asym = data.grid$pred.HCmax.mean[1:npops])[,2]
  # Number of COTS sustaianble at this level
  ConsRate = 250 # coral consumption rate
  CAvailable = (MinG.0.5*data.grid$PercentReef/10000)*1e6*1e4 # in cm2
  MinK.0.5A = CAvailable/(ConsRate*365)
  MinK.0.5J1 = MinK.0.5A*COTS_StableStage[1]/COTS_StableStage[3]
  MinK.0.5J2 = MinK.0.5A*COTS_StableStage[2]/COTS_StableStage[3]
  MinK.0.5 = as.matrix(cbind(MinK.0.5J1, MinK.0.5J2, MinK.0.5A))
return(list=c(as.data.frame(MinK.0.5), data.frame(MinK.10)))
}


CoralCOTSMort = function(p,CoralCover) {
  (1 - (p*CoralCover/(10+CoralCover)))
}
logistic.mort = function(phi1, phi2,phi3, x) {
y <-phi1/(1+exp(-(phi2+phi3*x)))
return(y)
}
logistic.mort = function(phi1, phi2,phi3, x) {
  y <-phi1/(1+exp(-(phi2+phi3*x)))
  return(y)
}
# threshdf = data.frame(MT = c(0.22,1,3,6,10,20),
#                       A=predict(MTCalib.gam, newdata=data.frame(MT=c(0.22,1,3,6,10,20))),
#                       J2 = predict(MTCalib.gam, newdata=data.frame(MT=c(0.22,1,3,6,10,20)))*COTS_StableStage[2]/COTS_StableStage[3],
#                       J1 = predict(MTCalib.gam, newdata=data.frame(MT=c(0.22,1,3,6,10,20)))*COTS_StableStage[1]/COTS_StableStage[3])
# COTS_StableStage
# plot(seq(0,20000000, by=10000), logistic.mort(1, 2, 0.00000015, seq(0,20000000, by=10000)))
# points(seq(0,20000000, by=10000), logistic.mort(1, 2, 0.00000005, seq(0,20000000, by=10000)))
# points(seq(0,20000000, by=10000), logistic.mort(1, 1.8, 0.0000001, seq(0,20000000, by=10000)))
# plot(seq(0,200000, by=10000), logistic.mort(1, 0.5, 0.00003, seq(0,200000, by=10000)))
# points(seq(0,200000, by=10000), logistic.mort(1, 0.4, 0.00001, seq(0,200000, by=10000)))
# points(seq(0,200000, by=10000), logistic.mort(1, 0.4, 0.000005, seq(0,200000, by=10000)))
# plot(CoralCOTSMort(0.2,seq(0,100, 1)))

doPredPreyDynamics = function(COTSabund, CoralCover, Crash, CCRatioThresh, CCRatioThresh2, maxmort, J2M, J1M, J2R, J1R) {
  # Implement ratio dependent mortality on Adults and J_2
  if (season=="winter") {
  # for (i in 1:(length(WhichPopCrash)-1)){
  #   WhichPopCrash[[i+1]] = WhichPopCrash[[i]]
  # }
  # browser()
  b = (COTSmort[3]-1)/(CCRatioThresh-CCRatioThresh2) # Make Ratio a parameter for tuning
  b2 = (COTSmort[2]-1)/(CCRatioThresh-CCRatioThresh2)
  b3 = (COTSmort[1]-1)/(CCRatioThresh-CCRatioThresh2)
  COTSMT = predict(MTCalib.gaminv, newdata=data.frame(DENS=COTSabund[,3]))^2
  COTSMT = ifelse(COTSMT <0.001, 0,COTSMT)
  Ratio = (CoralCover*data.grid$PercentReef/100)/COTSMT
  J2Mort = logistic.mort(1, J2M, J2R, COTSabund[,2])
  J1Mort = logistic.mort(1, J1M, J1R, COTSabund[,1])
  # WhichPopCrash[[1]] = which(Ratio<CCRatioThresh2) # find pops to crash
  # COTSabund[WhichPopCrash[[1]],2:3] = COTSabund[WhichPopCrash[[1]],2:3]*(1-maxmort)
  COTSabund[,"A"][which(Ratio<CCRatioThresh2)] = round((COTSabund[,"A"]*(1-maxmort))[which(Ratio<CCRatioThresh2)],0)
  COTSabund[,"A"][which(Ratio<CCRatioThresh & Ratio >= CCRatioThresh2)] = 
    round((COTSabund[,"A"]*(1 + (b*Ratio)))[which(Ratio<CCRatioThresh & Ratio >= CCRatioThresh2)],0)
  COTSabund[,"A"][which(Ratio>=CCRatioThresh)] = round((COTSabund[,"A"]*(1-COTSmort[3]))[which(Ratio>=CCRatioThresh)],0)
  COTSabund[,"J_2"] = COTSabund[,"J_2"]*(1-J2Mort)
  COTSabund[,"J_1"] = COTSabund[,"J_1"]*(1-J1Mort)
  COTSabund[,"J_2"][which(Ratio<CCRatioThresh2)] = round((COTSabund[,"J_2"]*(1-maxmort))[which(Ratio<CCRatioThresh2)],0)
  # COTSabund[,"J_2"][which(Ratio<CCRatioThresh & Ratio >= CCRatioThresh2)] = 
  #   round((COTSabund[,"J_2"]*(1 + (b2*Ratio)))[which(Ratio<CCRatioThresh & Ratio >= CCRatioThresh2)],0)
  # COTSabund[,"J_2"][which(Ratio>=CCRatioThresh)] = round((COTSabund[,"J_2"]*(1-COTSmort[2]))[which(Ratio>=CCRatioThresh)],0)
  COTSabund[,"J_1"][which(Ratio<CCRatioThresh2)] = round((COTSabund[,"J_1"]*(1-maxmort))[which(Ratio<CCRatioThresh2)],0)
  # COTSabund[,"J_1"][which(Ratio<CCRatioThresh & Ratio >= CCRatioThresh2)] =
  #   round((COTSabund[,"J_1"]*(1 + (b3*Ratio)))[which(Ratio<CCRatioThresh & Ratio >= CCRatioThresh2)],0)
  # COTSabund[,"J_1"][which(Ratio>=CCRatioThresh)] = round((COTSabund[,"J_1"]*(1-COTSmort[1]))[which(Ratio>=CCRatioThresh)],0)
  }
  # COTS.m.CC = (1 - (p*CoralCover/(10+CoralCover)))
  # COTSabund[,"A"] = COTSabund[,"A"]*exp(-COTS.m.CC*COTSmort[3])
  # COTSabund[,"J_2"] = COTSabund[,"J_2"]*exp(-COTS.m.CC*COTSmort[2])
  # COTSabund[,"J_1"] = COTSabund[,"J_1"]*exp(-COTS.m.CC*COTSmort[1])
  return(COTSabund)
} 

# Exponential Model of NAtural Mortality for Juveniles Keesing and Halford ----
# df = data.frame(size = c(1,3,5), M = c(6.5,1.24,0.45))
# mortlm = lm(log(M)~size, df)
# juv_growth = function(days){
#   7.08*exp(0.0045*days) - 6.60
# }
# size = juv_growth(1:365)
# mort = exp(predict(mortlm,newdata = data.frame(size = size))) + 0.13 # 0.13 is fish pred
# 
# mortdf = data.frame(Day = 1:365, size = size, mort = 1 -(mort/100))
# prod(mortdf$mort)

# Results[(Results$Year==1999),"CoralCover"] = rep(2,6000)
# 
# CoralCovertest = rep(2,6000)
# doPredPreyDynamics("summer", 2000, Results = Results, K=K, COTSabund = COTSabund, CoralCover = CoralCovertest )


#### Do Coral Distrubance ----
doCoralDisturbance = function (i, j, season, CoralCover, COTSfromCoralModel = F, 
                               storms.rsmpl, B.STORMS, WQ_Cyclone, 
                               COTS.rsmpl, B.COTS, WQ_CoTS,
                               bleaching.rsmpl, B.BLEACHING, WQ_bleach,
                               disease.rsmpl, B.DISEASE, unknown.rsmpl, B.UNKNOWN) {
  CoralCover = log(CoralCover)
  if (season =="summer"){
    # browser()
    ### Apply disturbances (bleaching, CoTS, disease, storms) in year 1994+i (starting from 1996)
    storms.rsmpl[,i,j][WQ > (-1)*B.STORMS[j]/WQ_Cyclone[1]] <- 0
    CoralCover[storms.rsmpl[,i,j]>0] <- CoralCover[storms.rsmpl[,i,j]>0] + 
      storms.rsmpl[,i,j][storms.rsmpl[,i,j]>0] * (B.STORMS[j] + WQ[storms.rsmpl[,i,j]>0] * WQ_Cyclone[1])
    if (COTSfromCoralModel == TRUE) {
      CoralCover[COTS.rsmpl[,i,j]>0] <- CoralCover[COTS.rsmpl[,i,j]>0] +
        COTS.rsmpl[,i,j][COTS.rsmpl[,i,j]>0] * (B.COTS[j] + WQ[COTS.rsmpl[,i,j]>0] * WQ_CoTS[1])
    }
    CoralCover[bleaching.rsmpl[,i,j]>0] <- CoralCover[bleaching.rsmpl[,i,j]>0] + 
      bleaching.rsmpl[,i,j][bleaching.rsmpl[,i,j]>0] * (B.BLEACHING[j] + WQ[bleaching.rsmpl[,i,j]>0] * WQ_bleach[1])
    CoralCover[disease.rsmpl[,i,j]>0] <- CoralCover[disease.rsmpl[,i,j]>0] + B.DISEASE[j] 
    CoralCover[unknown.rsmpl[,i,j]>0] <- CoralCover[unknown.rsmpl[,i,j]>0] + B.UNKNOWN[j] 
    CoralCover[CoralCover < log(0.1)] <- log(0.1) # sets minimal value to 10% (as 0% does not allow for recovery. 10% is the minimum HC cover observed in the LTMP data)
    # CoralCover <- b0 + (1 - b1)* CoralCover
  }
  return(exp(CoralCover))
  ### Store results
  #res[,year,j] <- exp(CoralCover)
  
}

doCoralGrowth = function(season, CoralCover, b0, b1) {
  if(season == "summer") {
    # b0.wq <- B0 + WQ * rnorm(length(WQ), mean=WQ.mn.sd[1], sd=WQ.mn.sd[2])
    # b1.wq <- b0.wq / log(HC.asym)
    CoralCover <- log(CoralCover)
    # CoralCover.t1 <- b0.wq + (1 - b1.wq)*CoralCover
    CoralCover.t1 = b0 + (1 - b1)*CoralCover
    return(cbind(CoralCover=exp(CoralCover.t1), CoralGrowth=(exp(CoralCover.t1)-exp(CoralCover))))
  }
  return(cbind(CoralCover=CoralCover, CoralGrowth=NA))
}

# doCoralGrowth("winter", 35,0.92,0.25)

####################
### COTS SANDBOX: for testing, etc.
####################
# 
#   
# 
# (CoralCover = CoralCoverParams$HCINI[,1][1:10])
# doCoralConsumption(1996,"summer",  COTSabund, CoralCover)
# doCoralConsumption(1996,"winter",  COTSabund, CoralCover)
# 
# 
#   
# # Step 1 - Create Geographic Distance Matrix
# # Step 2 - Convert to Sparse
#   
#   #############
#   # Build Dispersal Matrix (Pdist)----
#   
#   # Import coords of all of our sites ..NB for now this is the Env_Data
#   # setwd(DATA_DIRECTORY)
#   # Gdist <- load("Reefs.csv", header = TRUE) 
#   # Coords <- PopData[,c('lon','lat')]
#   # coordinates(Coords) <- ~lon+lat
#   # Gdist <- gDistance(Coords, Coords, byid = T)
#   
#   
#   # Convert to km
#   Gdist <- Gdist*100
#   
#   # Limit the distance matrix to ~500km. i.e set >500km to NA
#   Gdist[Gdist>500] <- 0
#   Gdist.Sp <- Matrix::Matrix(Gdist, sparse=T)
#   Gdist.Sp[1:10,1:10]
#   
#   #function to ignore NA's when summing
#   plus <- function(x) {
#     if(all(is.na(x))){
#       c(x[0],NA)} else {
#         sum(x,na.rm = TRUE)}
#   }
#   # Assume probability to be the inverse of distance --> need cumulative distrubtion function
#   Pdist <- apply(Gdist.Sp[1:10,1:10], 2, function(x) ((1/x)/plus(1/x))) 
#   signif(1/Gdist.Sp[1:10,1:10],3)
#   
#   
# 
#   
#     
# # ----------------------------------------------------   
# 
# x = matrix(rnorm(20), ncol=4)
# rownames(x) = paste("X", 1:nrow(x), sep=".")
# y = matrix(rnorm(12), ncol=4)
# rownames(y) = paste("Y", 1:nrow(y), sep=".")
# 
# 
# #find geographic distances between all sites
# 
# Pop1 <- PopData
# coordinates(Pop1) <- ~x+y
# Gdist <- gDistance(Pop1, Pop1, byid = T)
# 
# Gdist[1:10,1:10]
# # scale geographic distances between 0-1
# 
# # 1/((1-GDistNorm) - 0.5)
# 
# 
# m <- matrix(1, 3,3)
# v <- 1:3
# m*v
   ## build a very rough population transition matrix...

#  typical reproductive female is 300 mm in diameter

# Mass <- COTS_MassFromDiam(300)    # 1132 grams
# Fec <- COTS_FecFromMass(Mass)     # each female produces approx 14 million larvae!
# 
#   # assume that maybe 0.0001 of these larvae establish on a reef
# Fec <- Fec*0.0001
# 
# TransMat <- matrix(c(0,0.03,0,0,0,0,0.2,0,Fec/2,0,0.3,0.1,Fec/(2*6),0,0,0.6),nrow=4)
# 
# loadPackage("popbio")
# library(popbio)
# lambda(TransMat)         # strong positive growth rate: 3.57
# stable.stage(TransMat)   #stable age distribution

#### take away: stable stage distribution is heavily biased towards juveniles

# ##### Testing VBG ----
# 
# VBG.Models[[1]]
# 
# fert <- FvD[[7]][[2]][,2]
# dens <- FvD[[7]][[2]][,1]
# rm(dens, fert)
# head(FvD[[1]][[2]])
# 
# theta <- c(1, 0.1, 0.1)
# 
# out7 <- optim(theta, fn = SSQ, method = "BFGS", dens = na.omit(dens), fert=fert, hessian = TRUE)
# out7$par

