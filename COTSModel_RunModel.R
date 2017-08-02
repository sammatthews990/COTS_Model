###############################
# RUN COTS MODEL
###############################

####################!
# 1 SET UP WORKSPACE ----
####################!

source("C://Users//jc312264//Documents//GitHub//COTS_Model//COTSModel_PrepareWorkspace.R")
load # Rdata containing connectivity matrices and precalculations

set.seed(1)
setwd(DATA_DIRECTORY)
npops = 10 #number of reefs we want to test
PopData <- read.csv("reefs.csv", header = TRUE)[1:npops,] # we will just use a subset for testing
COTS.data <- read.csv("Disturbance/CoTS_data.csv", header = TRUE)[1:npops,]

# Load Connectivity Matrix
setwd(BASE_DIRECTORY)
load("R_Objects/ProbDistance.Rdata")
ConnMat.full = as.matrix(Pdist.Sp)
ConnMat = ConnMat.full[1:npops, 1:npops]
# ConnMat[ConnMat>1]=

sum(ConnMat.full[1,][-1])
  
setwd(ENVDATA_DIRECTORY)
data.grid = read.csv("data.grid.csv", header=T)[1:npops,]
WQ = data.grid$Primary + data.grid$Secondary + data.grid$Tertiary

# Gompertz parameters: disturbance & WQ effect sizes (as per MacNeil et al. in prep)
bleaching.mn.sd <- c(-0.19, 0.01)
COTS.mn.sd <- c(-0.54, 0.04)
disease.mn.sd <- c(-0.13, 0.01)
storms.mn.sd <- c(-0.64, 0.01)
unknown.mn.sd <- c(-0.16, 0.01)
WQ.mn.sd <- c(-0.68, 0.03)

####################!
# 2 INITIALIZE MODEL ----
####################!

initializeModel = function(PopData,COTSabund,CoralCover){

# Probably Change Storage to an array
Results = data.frame(sapply(PopData, rep.int, times=NYEARS*NSEASONS),
                     Year=rep(1996:2015,each=2*npops), Season=rep(c("summer", "winter"),each=npops), 
                     COTSJ1=NA, COTSJ2=NA, COTSA=NA, CoralCover=NA, DistCOTS=NA, DistCYCL=NA, DistBLCH=NA)

for(year in 1996){                  # loop through years
  for(season in SEASONS){               # loop through seasons
    COTSabund = doCOTSDispersal(season,COTSabund,SexRatio,ConnMat, avgPCF, sdPCF)
    COTSabund = doCOTSDemography(season, COTSabund, COTSmort, COTSremain)
    #COTSabund = doPredPreyDynamics(season, year, COTSabund,Reults)
    CoralCover = doCoralConsumption(year, season,  COTSabund, CoralCover, ConsRateS, ConsRateW)
    CoralCover = doCoralGrowth(CoralCover, B0, WQ)
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

NREPS = 10

masterDF = MakeLHSSamples(NREPS = 10)

####################!
# 3 RUN MODEL ----
####################!

runModel = function(masterDF, PopData, rep) {
    SexRatio = masterDF[rep, "SexRatio"]
    ConsRateW = masterDF[rep, "ConsRateW"]
    ConsRateS = masterDF[rep, "ConsRateS"]
    avgPCF = masterDF[rep, "avgPCF"]
    sdPCF = masterDF[rep, "sdPCF"]
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
    COTSabund = initializeCOTSabund(PopData, COTS.data, 1996, stagenames, COTS_StableStage)   # initialize the COTS abundance object (for year 0) 
    Results = initializeModel(PopData, COTSabund, CoralCover)
       
    # year Loop
    for(year in 1997:2015){                  # loop through years
      for(season in SEASONS){               # loop through seasons
        COTSabund = doCOTSDispersal(season,COTSabund,SexRatio,ConnMat, avgPCF, sdPCF)
        COTSabund = doCOTSDemography(season, COTSabund, COTSmort, COTSremain)
        COTSabund = doPredPreyDynamics(season, year, COTSabund, Results)
        CoralCover = doCoralConsumption(year, season,  COTSabund, CoralCover, ConsRateS, ConsRateW)
        CoralCover = doCoralGrowth(CoralCover, B0, WQ)
        #CoralCover = doCoralDisturbance(season,CoralCover,...)           # coral disturbance processes, including from COTS
        Results[(Results$Year==year) & (Results$Season==season),
                c("COTSJ1", "COTSJ2", "COTSA", "CoralCover")] = cbind(COTSabund, CoralCover)
      }
    }
    setwd(RESULTS_DIRECTORY)
    name <- sprintf("Sample_%s.Rdata",rep)
    save(Results, file=name)
    
}

# for (i in 1:10){
#   runModel(masterDF, PopData,i)
# }
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

df = as.data.frame(CoralMat) %>% dplyr::filter(PIXEL_ID==1) %>% tidyr::gather(LHS, CoralCover, 4:(NREPS+3))
df$rownum = rep((1:(NYEARS*NSEASONS)), NREPS)

Years = as.character(1996:2015)
# plot
ggplot(data = df, aes(x=rownum, y=CoralCover)) + geom_line(aes(colour=LHS))+
  scale_x_continuous("Year", breaks = seq(1,40,4), labels = Years[c(TRUE,FALSE)]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####################!
# 3 ANALYSE BRT ----
####################!