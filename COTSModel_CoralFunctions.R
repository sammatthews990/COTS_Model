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

setwd(ENVDATA_DIRECTORY)
data.grid = read.csv("data.grid.csv", header=T) 

intializeCoralCoverParams = function(data.grid, nsims){
  WQ <- data.grid$Primary + data.grid$Secondary + data.grid$Tertiary
  N <- nsims
  HCINI <- HCMAX <- B0 <- matrix(NA, ncol = N, nrow = dim(data.grid)[1])
  
  for (i in 1:dim(data.grid)[1]) {
    #Define sigma (i.e. variance-covariance matrix) for ith grid cell 
    sigma <- matrix(c(I(data.grid$pred.HCini.sd[i])^2, 0.54*data.grid$pred.HCini.sd[i]*data.grid$pred.HCmax.sd[i], 0.14*data.grid$pred.HCini.sd[i]*data.grid$pred.b0.sd[i],
                      0.54*data.grid$pred.HCini.sd[i]*data.grid$pred.HCmax.sd[i], I(data.grid$pred.HCmax.sd[i])^2, -0.13*data.grid$pred.HCmax.sd[i]*data.grid$pred.b0.sd[i], 
                      0.14*data.grid$pred.HCini.sd[i]*data.grid$pred.b0.sd[i], -0.13*data.grid$pred.HCmax.sd[i]*data.grid$pred.b0.sd[i], I(data.grid$pred.b0.sd[i])^2),
                    3,3)
    
    #Pick N random parameters for ith grid cell
    pick <- mvrnorm(n=N, mu=c(data.grid$pred.HCini.mean[i], data.grid$pred.HCmax.mean[i], data.grid$pred.b0.mean[i]), Sigma = sigma)
    HCINI[i,] <- pick[,1]
    HCMAX[i,] <- pick[,2]
    B0[i,] <- pick[,3]
  }
  return(list(WQ=cbind(data.grid[,1:5], WQ), 
              HCINI=HCINI, HCMAX=HCMAX, B0=B0))
}

CoralCoverParams = intializeCoralCoverParams(data.grid = data.grid, nsims=10)

#################!
# doCoralDistrurbances ----
#################!  

#########################