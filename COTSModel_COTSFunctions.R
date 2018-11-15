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


initializeCOTSabund <- function(PopData, COTS.data, Year, stagenames, COTS_StableStage, npops){
  # browser()
  ### set NA Values in interpolation CoTS.init to 0
  COTS.data[is.na(COTS.data)] <- 0
  
  nstages <- length(stagenames)
  
  ### Set up the COTS abundance object
  COTSabund <- matrix(0,nrow=npops, ncol=nstages)
  colnames(COTSabund) <- stagenames
  
  ### Set up reference for year
  colname <- paste('COTS_', Year, sep="")
  
  ### Update abundances based from interpolated manta tow data
  COTSabund[,'A'] <- COTS.data[,colname] * 1500 * (PopData$PercentReef/100)   #need to multiply by function from observations to density
  COTSabund[,'J_2'] <- COTSabund[,'A'] * as.numeric(COTS_StableStage[2]/COTS_StableStage[3])
  COTSabund[,'J_1'] <- COTSabund[,'A'] * as.numeric(COTS_StableStage[1]/COTS_StableStage[3])
  return(COTSabund)
}

# initCOTS <- initializeCOTSabund(PopData, COTS.data, 1996, stagenames, COTS_StableStage, npops)

    
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
  newCOTSabund <- COTSabund - COTS_Mort
  
  # number of COTS remaining and transitioning for each stage based on post-mortality abundaces
  COTS_Remain <- sweep(newCOTSabund,MARGIN=2,COTSremain,`*`)
  COTS_Trans <- sweep(newCOTSabund,MARGIN=2,1-COTSremain,`*`)
  
  # update newCOTSabund
  newCOTSabund[, 'J_1'] <- COTS_Remain[,'J_1']
  newCOTSabund[, 'J_2'] <- COTS_Remain[,'J_2'] + COTS_Trans[,'J_1']
  newCOTSabund[, 'A'] <- COTS_Remain[,'A'] + COTS_Trans[,'J_2']
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
nEggs <- COTS_Fecundity(COTSabund, 350, 100, SR = SexRatio, npops=npops)
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

# COTS_Fecundity.PCF(COTSabund, PCFParams, SexRatio=5)

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
  fEggs <- FvDParams[SexRatio,"Linf"] * (1 - exp(-FvDParams[SexRatio,"K"] * (nCOTS - FvDParams[SexRatio,"t0"])))
  nLarvae = nEggs * fEggs
  return(nLarvae)
}

# nLarvae = COTS_Fertilisation(nEggs = nEggs, nCOTS = COTSabund[,'A'], SR=5, Params = Params)





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

CoTS_Dispersal <- function(ConnMat, COTSabund, nLarvae, CoralCover, Loss, npops){
  # PsuedoCode:
  # 1. Loop through every population of fertilised eggs
  # 2. Discard 80% of eggs, then distribute the remainder as a function of distance
  # 3. Multiply nLarvae * connectivity probabilities to get larvae max
  # 5. Determine path of larvae to estimate survival (leave this for now)
  # 4. add vector of larvae produced to COTSabund
  
  nLarvae = (1-Loss)*nLarvae # assume 80% lost to the sea
  nArriving = matrix(nrow=npops, ncol=npops)
  for (i in 1:npops) {
    nArriving[i,] = signif(nLarvae[i]*(ConnMat[i,]/sum(ConnMat[i,])),3)
  }
  nArriving = base::colSums(nArriving) # create total juveniles arriving
  COTSabund[,'J_1'] = COTSabund[,'J_1'] + nArriving # add these to our abundance
  # need to add a line that increases survival based on maternal resources
  return(COTSabund)
}

# COTSabund.t1 = CoTS_Dispersal(ConnMat = Pdist.test, COTSabund = COTSabund,
#                               nLarvae=nLarvae, Loss=0.8)


##############!
# doCOTSDispersal ----
#############! 

# This is the master fucntion that incorporates the above functions

doCOTSDispersal = function(season, COTSabund, SexRatio, ConnMat, PCFParams, npops, FvDParams){
  #browser()
  COTSabund.t1 = COTSabund
  if (season=="summer"){
    nEggs = COTS_Fecundity.PCF(COTSabund, PCFParams = PCFParams, SexRatio = SexRatio)
    nLarvae = COTS_Fertilisation(nEggs = nEggs, SexRatio = SexRatio, FvDParams = FvDParams, nCOTS = COTSabund[,'A'])
    COTSabund.t1 = CoTS_Dispersal(ConnMat = ConnMat, COTSabund = COTSabund,
                                  nLarvae=nLarvae, Loss=0.99, npops=npops)
  }
  return(COTSabund.t1)
}

# test.summer = doCOTSDispersal("summer", initCOTS, SR=5, ConnMat = Pdist.test)
# test.winter = doCOTSDispersal("winter", initCOTS, SR=5, ConnMat = Pdist.test)

##############!
# doCOTSdemography ----
#############! 

# This is the master fucntion that incorporates the above functions
doCOTSDemography = function(season, COTSabund, COTSmort, COTSremain){
  COTSabund.t1 = COTSabund
  if (season=="winter"){ 
    COTSabund.t1 = COTS_StageTransition(COTSabund, COTSmort, COTSremain)
  }
  return(COTSabund.t1)
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

doCoralConsumption = function(year, season, COTSabund, CoralCover, ConsRateS, ConsRateW) {
  if (season =="summer") {
    #CoralCover= Results[(Results$Year==year-1) & (Results$Season=="winter"),"CoralCover"]
    CAvailable = (CoralCover*data.grid$PercentReef/10000)*1e6*1e4 # in cm2
    CConsumed = ConsRateS*COTSabund[,"A"]*182
    CRemaining=((CAvailable-CConsumed)/1e10)*(10000/data.grid$PercentReef)
    CRemaining[CRemaining < 0.5] <- 0.5
  } 
  if (season =="winter") {
    #CoralCover= Results[(Results$Year==year) & (Results$Season=="summer"),"CoralCover"]
    CAvailable = (CoralCover*data.grid$PercentReef/10000)*1e6*1e4 # in cm2
    CConsumed = ConsRateW*COTSabund[,"A"]*182
    CRemaining=((CAvailable-CConsumed)/1e10)*(10000/data.grid$PercentReef)
    CRemaining[CRemaining < 0.5] <- 0.5
  }
  return(CRemaining)
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
data.grid=data.grid[1:npops,]
CC.10=rep(10,npops)
# Calculate percent growth
MinG.10 = doCoralGrowth(CC.10, B0=data.grid$pred.b0.mean[1:npops], WQ[1:npops], HCMAX = data.grid$pred.HCmax.mean[1:npops])-CC.10
# Number of COTS sustaianble at this level
ConsRate = 250 # coral consumption rate
CAvailable = (MinG.10*data.grid$PercentReef/100)*1e6*1e4 # in cm2
MinK.10A = CAvailable/(ConsRate*182)
MinK.10J1 = MinK.10A*COTS_StableStage[1]/COTS_StableStage[3]
MinK.10J2 = MinK.10A*COTS_StableStage[2]/COTS_StableStage[3]
MinK.10 = as.matrix(cbind(MinK.10J1, MinK.10J2, MinK.10A)) ##THERE IS AN ISSUE WITH THE MAGNITUDE OF COTS CARRYING CAPACITY

CC.0.5=rep(0.5,npops)
# Calculate percent growth
MinG.0.5 = doCoralGrowth(CC.0.5, B0=data.grid$pred.b0.mean[1:npops], WQ[1:npops], HCMAX = data.grid$pred.HCmax.mean[1:npops])-CC.0.5
# Number of COTS sustaianble at this level
ConsRate = 250 # coral consumption rate
CAvailable = (MinG.0.5*data.grid$PercentReef/100)*1e6*1e4 # in cm2
MinK.0.5A = CAvailable/(ConsRate*182)
MinK.0.5J1 = MinK.0.5A*COTS_StableStage[1]/COTS_StableStage[3]
MinK.0.5J2 = MinK.0.5A*COTS_StableStage[2]/COTS_StableStage[3]
MinK.0.5 = as.matrix(cbind(MinK.0.5J1, MinK.0.5J2, MinK.0.5A))
return(list=c(as.data.frame(MinK.0.5), data.frame(MinK.10)))
}


doPredPreyDynamics = function(season, year, COTSabund,Results,K) {
  # after 1 year at low levels COTS densities get brought down to levels supported by growth
  # work out the %CC growth from 0.5% and set that as the number of COTS
  # MinK is the Maximum CoTS that can be supported at depleted coral cover
  if(season=="summer"){
    prevCC = dplyr::filter(Results, Year==year-1 & Season=='winter') %>% dplyr::select(CoralCover)
    prevCC = as.matrix(cbind(prevCC, prevCC,prevCC))
    ####THis is the problem
    newCOTS = ifelse(prevCC < 10 & prevCC > 1 & COTSabund[,'A']> K$MinK.10A, K$MinK.10A/10000,
           ifelse(prevCC < 1 & COTSabund[,'A']> K$MinK.0.5A, K$MinK.0.5A/10000, COTSabund))
    colnames(newCOTS) = stagenames
  }
  if(season=="winter"){
    prevCC = dplyr::filter(Results, Year==year & Season=='summer') %>% dplyr::select(CoralCover)
    prevCC = as.matrix(cbind(prevCC, prevCC,prevCC))
    newCOTS = ifelse(prevCC < 10 & prevCC > 1 & COTSabund[,'A']> K$MinK.10A, K$MinK.10A/10000,
              ifelse(prevCC < 1 & COTSabund[,'A']> K$MinK.0.5A, K$MinK.0.5A/10000, COTSabund))
    colnames(newCOTS) = stagenames
  }
  return(newCOTS)
} 




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

Mass <- COTS_MassFromDiam(300)    # 1132 grams
Fec <- COTS_FecFromMass(Mass)     # each female produces approx 14 million larvae!

  # assume that maybe 0.0001 of these larvae establish on a reef
Fec <- Fec*0.0001

TransMat <- matrix(c(0,0.03,0,0,0,0,0.2,0,Fec/2,0,0.3,0.1,Fec/(2*6),0,0,0.6),nrow=4)

loadPackage("popbio")
library(popbio)
lambda(TransMat)         # strong positive growth rate: 3.57
stable.stage(TransMat)   #stable age distribution

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

