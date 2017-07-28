###############################
# RUN COTS MODEL
###############################


COTSabund = initializeCOTSabund(PopData = PopData, COTS.init = COTS.data, 
                                Year=1996, stagenames, COTS_StableStage = COTS_StableStage)   # initialize the COTS abundance object (for year 0) 
initializeCoralCover(,...)    # initialize the coral cover object (for year 0)

Results = data.frame(sapply(PopData, rep.int, times=NYEARS*NSEASONS),
                     Year=rep(1996:2005,each=2*npops), Season=rep(c("summer", "winter"),each=npops), 
                     COTSJ1=NA, COTSJ2=NA, COTSA=NA, CoralCover=NA, DistCOTS=NA, DistCYCL=NA, DistBLCH=NA)

for(year in 1996:1998){                  # loop through years
  for(season in SEASONS){               # loop through seasons
    COTSabund = doCOTSDispersal(season,COTSabund,SR,ConnMat)
    COTSabund = doCOTSDemography(season,COTSabund = COTSabund.t1)
    # CoralCover.t1 = doCoralDispersal(season, CoralCover...)
    # CoralCover.t1 = doCoralDisturbance(season,CoralCover,...)           # coral disturbance processes, including from COTS
    Results[(Results$Year==year) & (Results$Season==season),
        c("COTSJ1", "COTSJ2", "COTSA")] = COTSabund
    
  }
}

