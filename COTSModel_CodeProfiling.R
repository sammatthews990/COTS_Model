###### PROFILING ----

library(profvis)
profvis({
  COTSabund = initializeCOTSabund(PopData = PopData, COTS.init = COTS.data, 
                                  Year=1996, stagenames, COTS_StableStage = COTS_StableStage)   # initialize the COTS abundance object (for year 0) 
  CoralCover = CoralCoverParams$HCINI[,1][1:10]
  
  # Probably Change Storage to an array
  Results = data.frame(sapply(PopData, rep.int, times=NYEARS*NSEASONS),
                       Year=rep(1996:2015,each=2*npops), Season=rep(c("summer", "winter"),each=npops), 
                       COTSJ1=NA, COTSJ2=NA, COTSA=NA, CoralCover=NA, DistCOTS=NA, DistCYCL=NA, DistBLCH=NA)
  
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
})

profvis({
COTS_Eggs <- vector(mode = "numeric", length = npops)
for (r in 1:npops) {
  Sizes <- rnorm(COTSabund[r,'A']*5/10, 35, 10)
  Sizes[Sizes<0] = 0
  COTS_Eggs[r] <- sum(COTS_FecFromMass(COTS_MassFromDiam(Sizes)))
}
})
