###############################!
# RUN COTS MODEL ----
###############################!


# THINGS TO DO
# 1. Fix estimates of percapita fecundity.. at the moment the mean and sd vary independently
# 2. Analyse Results with BRT's
# 3. Graph Coral Cover Alongside COTS
# 4. Find way to link Mortality and Remain rates to Coral Cover
# 5. Make stable stage vary using a sinlge parameter

####################!
# 1 SET UP WORKSPACE ----
####################!

source("C://Users//jc312264//Documents//GitHub//COTS_Model//COTSModel_PrepareWorkspace.R")
dirs = list(BASE=BASE_DIRECTORY, CODE=CODE_DIRECTORY, RESULTS=RESULTS_DIRECTORY, 
            DATA=DATA_DIRECTORY, ENVDATA=ENVDATA_DIRECTORY,SPATIALDATA=SPATIALDATA_DIRECTORY)
setwd(CODE_DIRECTORY)
# optional to change the number of populations we will test for
npops = 10
source("COTSModel_LoadObjectsForModelling.R")
setwd(BASE_DIRECTORY)
save(ConnMat, file = "R_Objects/ConnMat.Rdata")
setwd(CODE_DIRECTORY)
source("COTSModel_CoralFunctions.R")
source("COTSModel_COTSFunctions.R")

####################!
# 2 GENERATE LATIN HYPERCUBESAMPLE ----
####################!

NREPS=10

masterDF = MakeLHSSamples(NREPS = 10)

####################!
# 2 SET UP CLUSTERS ----
####################!

num_cores <- parallel::detectCores() - 2

# Parallel for loop
cl <- parallel::makeCluster(num_cores,outfile="LOG.TXT")
doParallel::registerDoParallel(cl=cl)    # make the cluster

####################!
# 2 RUN PARALLEL ----
####################!

library(foreach)

allsamples <- foreach::foreach(i = 1:nrow(masterDF) 

)%dopar% { 
                      
RunMaster = function (masterDF, PopData, dirs){
  
  force(masterDF)
  force(PopData)
  force(dirs)
  # source("C://Users//jc312264//Documents//GitHub//COTS_Model//COTSModel_PrepareWorkspace.R")
  # dirs = list(BASE=BASE_DIRECTORY, CODE=CODE_DIRECTORY, RESULTS=RESULTS_DIRECTORY, 
  #             DATA=DATA_DIRECTORY, ENVDATA=ENVDATA_DIRECTORY,SPATIALDATA=SPATIALDATA_DIRECTORY)
  setwd(dirs$BASE)
  load(file = "R_Objects/ConnMat.Rdata")
  setwd(dirs$CODE)
  source("COTSModel_LoadSmallObjectsForModelling.R")
  setwd(dirs$CODE)
  source("COTSModel_CoralFunctions.R")
  source("COTSModel_COTSFunctions.R")
  

  ####################!
  # 3 INITIALIZE MODEL ----
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
  # 4 RUN MODEL ----
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
      setwd(dirs$RESULTS)
      name <- sprintf("Sample_%s.Rdata",rep)
      save(Results, file=name)
      
  }
  
  # for (i in 1:NREPS){
  #   runModel(masterDF, PopData,i)
  # }

runModel(masterDF, PopData, rep=i)

}

RunMaster(masterDF, PopData, dirs)(i)

}

###################!
# CLOSE CLUSTER ----
###################!

if(!is.null(cl)) {
  parallel::stopCluster(cl)
  cl <- c()
}

####################!
# 3 HARVEST RESULTS ----
####################!

#### MAKE RESULTS DATA FRAMES

NREPS=10

HarvestData = function(RESULTS_DIRECTORY) {
  
setwd(RESULTS_DIRECTORY)

CoralMat = matrix(NA, nrow = (npops*NYEARS*NSEASONS), ncol = NREPS+3)
COTSMat = matrix(NA, nrow = (npops*NYEARS*NSEASONS), ncol = NREPS+3)
load("Sample_1.Rdata")
CoralMat[,1]=COTSMat[,1] = Results[, "PIXEL_ID"]
CoralMat[,2]=COTSMat[,2] = Results[, "Year"]
CoralMat[,3]=COTSMat[,3] = Results[, "Season"]
colnames(CoralMat) = c("PIXEL_ID", "Year", "Season",  1:NREPS)
colnames(COTSMat) = c("PIXEL_ID", "Year", "Season", 1:NREPS)

for (i in 1:NREPS){
  load(sprintf("Sample_%s.Rdata", i))
  CoralMat[,i+3] = Results[,"CoralCover"]
  COTSMat[,i+3] = Results[,"COTSA"]
}

# # Select Reefs to Plot
# 
# df = as.data.frame(CoralMat) %>% dplyr::filter(PIXEL_ID==10) %>% tidyr::gather(LHS, CoralCover, 4:(NREPS+3))
# df$rownum = rep((1:(NYEARS*NSEASONS)), NREPS)
# df1 = as.data.frame(COTSMat) %>% dplyr::filter(PIXEL_ID==10) %>% tidyr::gather(LHS, COTS, 4:(NREPS+3))
# df1$rownum = rep((1:(NYEARS*NSEASONS)), NREPS)
# 
# #convert COTS abundance back to obs/manta tow
# 
# Years = as.character(1996:2015)
# # plot Coral Cover
# ggplot(data = df, aes(x=rownum, y=CoralCover)) + geom_line(aes(colour=LHS))+
#   scale_x_continuous("Year", breaks = seq(1,40,4), labels = Years[c(TRUE,FALSE)]) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   # geom_line(data = df1, aes(x=rownum, y=COTS/1e9, colour=LHS))+
#   # scale_y_continuous(sec.axis = sec_axis(~.^1e9, name = "CoTS Abundance"))+
#   theme_classic()
# # Plot CoTS abundance
# ggplot(data = df1, aes(x=rownum, y=COTS)) + geom_line(aes(colour=LHS))+
#   scale_x_continuous("Year", breaks = seq(1,40,4), labels = Years[c(TRUE,FALSE)]) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   # geom_line(data = df1, aes(x=rownum, y=COTS/1e9, colour=LHS))+
#   # scale_y_continuous(sec.axis = sec_axis(~.^1e9, name = "CoTS Abundance"))+
#   theme_classic()

# Add Columns to master df to show max COTS and time between outbreaks (>2500)

CreateResponseVars = function(COTSMat) {

df1 = as.data.frame(COTSMat) %>% tidyr::gather(LHS, COTS, 4:(NREPS+3))
ResponseVars = df1 %>% dplyr::group_by(as.numeric(LHS)) %>% 
        summarize(COTSmax = max(COTS), COTSsd = sd(COTS))
ResponseVars$meanrecovery = ResponseVars$meanoutbreak = ResponseVars$totoutbreaks = ResponseVars$time2outbreak = NA

    Outbreak.interval = function(LHSSubset) {
      time2outbreak = vector("numeric", npops)
      totoutbreaks = vector("numeric", npops)
      meanoutbreak = vector("numeric", npops)
      meanrecovery  = vector("numeric", npops)
      for (i in 1:npops){
        sub = LHSSubset %>% dplyr::filter(PIXEL_ID==i)
        time2outbreak[i]=sum(cumprod(!sub$COTS>2500))
        totoutbreaks[i]=sum(sub$COTS>2500)
        OutbreakRLE=rle(sub$COTS>2500)  
        meanoutbreak[i] = mean(OutbreakRLE$lengths[OutbreakRLE$values][-1])
        meanrecovery[i] = mean(OutbreakRLE$lengths[!OutbreakRLE$values][-1])
      }
      
      return(c(mean(time2outbreak), mean(totoutbreaks), mean(meanoutbreak), mean(meanrecovery)))
    }
for (i in 1:NREPS){
  LHSSubset= df1 %>% dplyr::filter(as.numeric(LHS)==i)
  ResponseVars[i, 4:7] = Outbreak.interval(LHSSubset)
  }
return(ResponseVars)
}

ResponseVars = CreateResponseVars(COTSMat) 

return(list(ResponseVars=ResponseVars,CoralMat=CoralMat, COTSMat=COTSMat))

}

myresults = HarvestData(RESULTS_DIRECTORY)$ResponseVars

####################!
# 3 ANALYSE BRT ----
####################!

masterDF$SexRatio <- as.factor(masterDF$SexRatio)  

masterDF = na.omit(cbind(masterDF, myresults))

CoTS.gbmRec <- dismo::gbm.step(data=as.data.frame(masterDF),
                       gbm.x = c(1:10), gbm.y = "meanrecovery", # need to add sector and shelf
                       family = "gaussian", #what is the family of distribution of our response 
                       learning.rate = 0.001,
                       tree.complexity = 3, #how "deep" are each tree, i.e how complex are the interactions 3 is for two
                       bag.fraction=0.5, # how much of the train data do we use for each tree
                       n.trees = 100,
                       max.trees = 100000)
summary(CoTS.gbmRec)

