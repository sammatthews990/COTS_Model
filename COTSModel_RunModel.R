###############################
# RUN COTS MODEL
###############################


# THINGS TO DO
# 1. Fix estimates of percapita fecundity.. at the moment the mean and sd vary independently
# 2. Analyse Results with BRT's
# 3. Graph Coral Cover Alongside COTS

####################!
# 1 SET UP WORKSPACE ----
####################!

source("C://Users//jc312264//Documents//GitHub//COTS_Model//COTSModel_PrepareWorkspace.R")
setwd(CODE_DIRECTORY)
# optional to change the number of populations we will test for
npops = 20
source("COTSModel_LoadObjectsForModelling.R")
setwd(CODE_DIRECTORY)
source("COTSModel_CoralFunctions.R")
source("COTSModel_COTSFunctions.R")



####################!
# 2 INITIALIZE MODEL ----
####################!

initializeModel = function(PopData,COTSabund,CoralCover, SexRatio, ConsRateS, ConsRateW, B0, WQ, HC.asym){

# Probably Change Storage to an array
Results = data.frame(sapply(PopData, rep.int, times=NYEARS*NSEASONS),
                     Year=rep(1996:2015,each=2*npops), Season=rep(c("summer", "winter"),each=npops), 
                     COTSJ1=NA, COTSJ2=NA, COTSA=NA, CoralCover=NA, DistCOTS=NA, DistCYCL=NA, DistBLCH=NA)

for(year in 1996){                  # loop through years
  for(season in SEASONS){               # loop through seasons
    COTSabund = doCOTSDispersal(season,COTSabund,SexRatio,ConnMat, PCFParams)
    COTSabund = doCOTSDemography(season, COTSabund, COTSmort, COTSremain)
    #COTSabund = doPredPreyDynamics(season, year, COTSabund,Reults)
    CoralCover = doCoralConsumption(year, season,  COTSabund, CoralCover, ConsRateS, ConsRateW)
    CoralCover = doCoralGrowth(CoralCover, B0, WQ, HC.asym)
    #CoralCover = doCoralDisturbance(season,CoralCover,...)           # coral disturbance processes, including from COTS
    Results[(Results$Year==year) & (Results$Season==season),
            c("COTSJ1", "COTSJ2", "COTSA", "CoralCover")] = cbind(COTSabund, CoralCover)
  }
}
return(Results)
}


####################!
# 3 GENERATE LATIN HYPERCUBESAMPLE ----
####################!


masterDF = MakeLHSSamples(NREPS = 10)

####################!
# 3 RUN MODEL ----
####################!

runModel = function(masterDF, PopData, rep) {
    SexRatio = masterDF[rep, "SexRatio"]
    ConsRateW = masterDF[rep, "ConsRateW"]
    ConsRateS = masterDF[rep, "ConsRateS"]
    PCFParams = c(masterDF[rep, "avgPCF"], masterDF[rep,"sdPCF"])
    # avgPCF = masterDF[rep, "avgPCF"]
    # sdPCF = masterDF[rep, "sdPCF"]
    COTSmort = as.numeric(masterDF[rep, c("mortJ1", "mortJ2", "mortA")])
    COTSremain = as.numeric(masterDF[rep, c("remJ1", "remJ2", "remA")])
    COTS_StableStage = as.numeric(masterDF[rep, c("cssJ1", "cssJ2", "cssA")])
    # avgAdultSize =
    # sdAdultSize = # These will change the fecundity estimatesC
    # need an Allee Effect
    # need to make stable stage vary by a scaling factor
    # make mortality and remain resource driven
    
    # Initialize
    CoralCoverParams = intializeCoralCoverParams(data.grid = data.grid, nsims=10)
    CoralCover = CoralCoverParams$HCINI[,1][1:10]
    B0=CoralCoverParams$B0[,1]
    HC.asym = CoralCoverParams$HCMAX[,1]
    COTSabund = initializeCOTSabund(PopData, COTS.data, 1996, stagenames, COTS_StableStage)   # initialize the COTS abundance object (for year 0) 
    Results = initializeModel(PopData, COTSabund, CoralCover, SexRatio, ConsRateS, ConsRateW, B0, WQ, HC.asym)
       
    # year Loop
    for(year in 1997:2015){                  # loop through years
      for(season in SEASONS){               # loop through seasons
        COTSabund = doCOTSDispersal(season,COTSabund,SexRatio,ConnMat, PCFParams)
        COTSabund = doCOTSDemography(season, COTSabund, COTSmort, COTSremain)
        COTSabund = doPredPreyDynamics(season, year, COTSabund, Results)
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

for (i in 1:npops){
  runModel(masterDF, PopData,i)
}
# setwd(RESULTS_DIRECTORY)
# load("Sample_3.Rdata")


####################!
# 3 PLOT RESULTS ----
####################!

# 1 Make Coral Cover and Adult CoTS abundance matrices

# 2 Load in ResultsDataFrames sequentially

CoralMat = matrix(NA, nrow = (npops*NYEARS*NSEASONS), ncol = NREPS+3)
COTSMat = matrix(NA, nrow = (npops*NYEARS*NSEASONS), ncol = NREPS+3)
load("Sample_1.Rdata")
CoralMat[,1]=COTSMat[,1] = Results[, "PIXEL_ID"]
CoralMat[,2]=COTSMat[,2] = Results[, "Year"]
CoralMat[,3]=COTSMat[,3] = Results[, "Season"]
colnames(CoralMat) = c("PIXEL_ID", "Year", "Season", sprintf("CoralLHS%s", 1:NREPS))
colnames(COTSMat) = c("PIXEL_ID", "Year", "Season", sprintf("COTSLHS%s", 1:NREPS))

for (i in 1:NREPS){
  load(sprintf("Sample_%s.Rdata", i))
  CoralMat[,i+3] = Results[,"CoralCover"]
  COTSMat[,i+3] = Results[,"COTSA"]
}

# Select Reefs to Plot

df = as.data.frame(CoralMat) %>% dplyr::filter(PIXEL_ID==10) %>% tidyr::gather(LHS, CoralCover, 4:(NREPS+3))
df$rownum = rep((1:(NYEARS*NSEASONS)), NREPS)

Years = as.character(1996:2015)
# plot
ggplot(data = df, aes(x=rownum, y=CoralCover)) + geom_line(aes(colour=LHS))+
  scale_x_continuous("Year", breaks = seq(1,40,4), labels = Years[c(TRUE,FALSE)]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####################!
# 3 ANALYSE BRT ----
####################!