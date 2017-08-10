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

rm(list=ls())
source("C://Users//jc312264//Documents//GitHub//COTS_Model//COTSModel_PrepareWorkspace.R")
# dirs = list(BASE=BASE_DIRECTORY, CODE=CODE_DIRECTORY, RESULTS=RESULTS_DIRECTORY, 
#             DATA=DATA_DIRECTORY, ENVDATA=ENVDATA_DIRECTORY,SPATIALDATA=SPATIALDATA_DIRECTORY)
# setwd(CODE_DIRECTORY)
# # optional to change the number of populations we will test for
# # source("COTSModel_LoadObjectsForModelling.R")
# setwd(BASE_DIRECTORY)
# # save(ConnMat, file = "R_Objects/ConnMat.Rdata")
# load("R_Objects/ConnMat.Rdata")
# setwd(CODE_DIRECTORY)
# source("COTSModel_CoralFunctions.R")
# source("COTSModel_COTSFunctions.R")


####################!
# 4 Make Worker Function For Parallelization ----
####################!

MakeWorker = function (masterDF, npops){
  
  force(masterDF)
  force(npops)
  # force(seasons)
  
  source("C://Users//jc312264//Documents//GitHub//COTS_Model//COTSModel_PrepareWorkspace.R")
  # dirs = list(BASE=BASE_DIRECTORY, CODE=CODE_DIRECTORY, RESULTS=RESULTS_DIRECTORY, 
  # DATA=DATA_DIRECTORY, ENVDATA=ENVDATA_DIRECTORY,SPATIALDATA=SPATIALDATA_DIRECTORY)
  setwd(BASE_DIRECTORY)
  load(file = "R_Objects/ConnMat.Rdata")
  setwd(CODE_DIRECTORY)
  source("COTSModel_LoadSmallObjectsForModelling.R")
  setwd(CODE_DIRECTORY)
  source("COTSModel_CoralFunctions.R")
  source("COTSModel_COTSFunctions.R")
  ConnMat=ConnMat
  # SEASONS=seasons
  # PopData = PopData[1:npops,]
  # COTS.data = COTS.data[1:npops,]
  # data.grid = data.grid[1:npops,]
  #force(PopData);force(COTS.data); force(data.grid)
  
  Worker = function(i){
    runModel(masterDF=masterDF, PopData=PopData,COTS.data = COTS.data, 
             data.grid = data.grid, ConnMat = ConnMat, npops=npops, rep=i, seasons = SEASONS)
  }
  return(Worker)
}

#MakeWorker(masterDF, npops=50)(5)

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

) %dopar% { 
          
  #dirs = list(BASE=BASE_DIRECTORY, CODE=CODE_DIRECTORY, RESULTS=RESULTS_DIRECTORY, 
  #           DATA=DATA_DIRECTORY, ENVDATA=ENVDATA_DIRECTORY,SPATIALDATA=SPATIALDATA_DIRECTORY)
  temp = MakeWorker(masterDF=masterDF,npops=npops, seasons=SEASONS)(i)

}

# for (i in 1:NREPS){
#   runModel(masterDF, PopData,i)
# }

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
        dplyr::summarize(COTSmax = max(COTS), COTSsd = sd(COTS))
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

