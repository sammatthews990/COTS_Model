
####################################
###  Utility functions for the COTS population model
#################################### 


###################
# loadPackages: GENERIC FUNCTION FOR INSTALLING/LOADING PACKAGES FROM CRAN
########################################
##########
# OBJECTIVE:
#    loads a package from CRAN into the global environment.  
# PARAMS:
#    - none
# RETURNS:
#    - nothing. results in libraries being loaded into the global environment 
###################

loadPackage <- function(pkg){

  if(pkg %in% rownames(installed.packages()) == FALSE) {suppressMessages(suppressWarnings(install.packages(pkg)))}
  eval(parse(text=sprintf("suppressMessages(suppressWarnings(require(%s)))",pkg)), envir= .GlobalEnv)

}

###################
# loadPackages
##########
# OBJECTIVE:
#    loads all necessary packages into the global environment. Makes use of the 
#    "loadPackage" utility function.  
# PARAMS:
#    - none
# RETURNS:
#    - nothing. results in libraries being loaded into the global environment 
###################

loadPackages <- function(){
  loadPackage("lhs")            # for latin hypercube sampling
  loadPackage("RCurl")          # for loading source code from github
  loadPackage("raster")         # for managing raster data
  loadPackage("leaflet")
  loadPackage("rgdal")          # for reading and writing all sorts of spatial data   
  loadPackage("popbio")         # for doing basic matrix population modeling
  loadPackage("tidyverse")      # data manipulation
  loadPackage("rgeos")          # geometry applications
  loadPackage("plyr")           # data wrangling
  loadPackage("dplyr")          # data wrangling
  loadPackage("tidyr")
  loadPackage("reshape2")       # data wrangling
  loadPackage("gstat")          # performing interpolation
  loadPackage("ggplot2")        # plotting
  loadPackage("emdbook")        # support for ecological models with data book
  loadPackage("SpatialTools")   # compute pairwise distances
  loadPackage("Matrix")         # Sparse matrices
  loadPackage("maptools")
  loadPackage("PBSmapping")
  loadPackage("RColorBrewer")
  loadPackage("sampling")
  loadPackage("psych")
  loadPackage("dismo")
  loadPackage("foreach")
  loadPackage("doParallel")
  #loadPackage("mbcv")
  #loadPackage("gbm")
  loadPackage("MASS")
  #loadPackage("nlme")
  #loadPackage("lme4")
}



###################
# source_github
##########
# OBJECTIVE:
#    Read source code from a GIThub repository 
# PARAMS:
#    - baseurl: URL for the repository
#    - scriptname: name of the script/library to load
# RETURNS:
#    - nothing. loads functions from GitHub 
###################

source_github <- function(baseurl,scriptname) {
  # load package
  loadPackage(RCurl)
 
  # read script lines from website
  url <- sprintf("%s%s",baseurl,scriptname)
  script <- getURL(url, ssl.verifypeer = FALSE)
  
  script <- gsub("\r\n", "\n", script)     # get rid of carriage returns (not sure why this is necessary...)
 
  # parse lines and evaluate in the global environement
  eval(parse(text = script), envir= .GlobalEnv)
}


###################
# saveWorkspace
##########
# OBJECTIVE:
#    save the current workspace or a set of R objects on the workspace 
# PARAMS:
#    - objNames: (optional) list of names of objects to save.
# RETURNS:
#    - nothing. saves an RData file to the R_Workspaces directory in Dropbox 
###################

 saveWorkspace <- function(objNames=NULL,filename="R_workspace", dir){
  if(is.null(objNames)) objNames = ls(all = TRUE,name=.GlobalEnv)
  setwd(dir)
  filename <- sprintf("%s_%s.RData",filename,Sys.Date())
  save(list = objNames, file = filename)
}
 

###########
## FUNCTION "specifyLHSParam"
##
## Information necessary to translate standard uniform LHS sample into parameters of interest for paleo project 
###########

specifyLHSParam <- function(paramslist,name,type,lb,ub){
  newlist <- paramslist
  eval(parse(text=sprintf("newlist$%s <- list()",name)))
  eval(parse(text=sprintf("newlist$%s$type <- \"%s\"",name,type)))
  eval(parse(text=sprintf("newlist$%s$lb <- %s",name,lb)))
  eval(parse(text=sprintf("newlist$%s$ub <- %s",name,ub))) 	
  return(newlist)
}

###############################
#        SPECIFY PARAMETER RANGES
################################


MakeLHSSamples <- function(NREPS){
  
  specifyLHSParam <- function(paramslist,name,type,lb,ub){
    newlist <- paramslist
    eval(parse(text=sprintf("newlist$%s <- list()",name)))
    eval(parse(text=sprintf("newlist$%s$type <- \"%s\"",name,type)))
    eval(parse(text=sprintf("newlist$%s$lb <- %s",name,lb)))
    eval(parse(text=sprintf("newlist$%s$ub <- %s",name,ub))) 	
    return(newlist)
  }
  
  LHSParms <- list()    # initialize the container for parameter bounds
  
  ### SEXRATIO : 1 = 0.1M:0.9F
  LHSParms <- specifyLHSParam(paramslist=LHSParms,name="SexRatio",type="CAT",lb=0,ub=9)
  
  ####  WINTER CONSUMPTION RATE
  LHSParms <- specifyLHSParam(LHSParms,"ConsRateW",type="CONT",lb=100,ub=200)
  
  #### SUMMER CONSUMPTION RATE
  LHSParms <- specifyLHSParam(LHSParms,"ConsRateS",type="CONT",lb=150,ub=400)
  
  ### AVERAGE PER CAPITA FECUNDITY 
  LHSParms <- specifyLHSParam(LHSParms,"avgPCF",type="CONT",lb=42e5,ub=21e6)       # KTS: changed to 500
  
  ### STD DEV PER CAPITA FECUNDITY
  LHSParms <- specifyLHSParam(LHSParms,"sdPCF",type="CONT",lb=6e4,ub=14.6e4)     # was 1000 to 40000  
  
  ### % JUVENILE1 MORTALITY PER TIME STEP
  LHSParms <- specifyLHSParam(LHSParms,"mortJ1",type="CONT",lb=0.5,ub=0.99)
  
  ### % JUVENILE2 MORTALITY PER TIME STEP
  LHSParms <- specifyLHSParam(LHSParms,"mortJ2",type="CONT",lb=0.4,ub=0.99) 
  
  ### % ADULT MORTALITY PER TIME STEP
  LHSParms <- specifyLHSParam(LHSParms,"mortA",type="CONT",lb=0.1,ub=0.6)
  
  ### % JUVENILE1 TO REMIAIN during transition
  LHSParms <- specifyLHSParam(LHSParms,"remJ1",type="CONT",lb=0.01,ub=0.5)
  
  ### % JUVENILE2 TO REMIAIN during transition 
  LHSParms <- specifyLHSParam(LHSParms,"remJ2",type="CONT",lb=0.1,ub=0.5)       # KTS: changed to 500
  
  ### % ADULT TO REMIAIN during transition
  LHSParms <- specifyLHSParam(LHSParms,"remA",type="CONT",lb=1,ub=1)     # was 1000 to 40000  
  
  ### PROPORTIONAL STABLE STAGE DISTRIBUTION J1
  LHSParms <- specifyLHSParam(LHSParms,"cssJ1",type="CONT",lb=0.9803,ub=0.9803)
  
  ### PROPORTIONAL STABLE STAGE DISTRIBUTION J2
  LHSParms <- specifyLHSParam(LHSParms,"cssJ2",type="CONT",lb=0.0171,ub=0.0171) 
  
  ### PROPORTIONAL STABLE STAGE DISTRIBUTION a
  LHSParms <- specifyLHSParam(LHSParms,"cssA",type="CONT",lb=0.0026,ub=0.0026)
  
  
  ### DENSITY DEPENDENCE ON HARVEST (y intercept of the harvest rate/abundance relationship)
  
  #  # HUNTDD_posvals <- c(seq(-0.05,0,length=10),seq(0.1,1,length=10),seq(1.1,2,length=10))
  # LHSParms <- specifyLHSParam(LHSParms,"HUNTDD",type="CONT",lb=-0.05,ub=2)
  
  #LHSParms <- specifyLHSParam(LHSParms,"HARVZ",type="CONT",lb=1,ub=2)
  
  #### HUMAN ARRIVAL (0 represents lower bound on a per-population basis, 1 represents upper bound)
  #LHSParms <- specifyLHSParam(LHSParms,"HUMAN",type="CONT",lb=0,ub=1)    # NOTE: was "CAT" not sure why
  
  ##################
  ##### GENERATE LATIN HYPERCUBE SAMPLE
  
  nVars <- length(names(LHSParms))  
  
  LHS <- lhs::randomLHS(NREPS, nVars )   # generate multiple samples from parameter space according to a LHS sampling scheme
  
  masterDF <- as.data.frame(LHS)    #  storage container (data frame) to record relevant details for each MP file. Rows:MP file/LHS samples. Cols: relevant variables
  
  
  ### translate raw lhs samples into desired parameter space
  colnames(masterDF) <- names(LHSParms)
  parm=1
  for(parm in 1:nVars){
    if(LHSParms[[parm]]$type=="CONT"){
      masterDF[,parm] <- LHSParms[[parm]]$lb + LHS[,parm]*(LHSParms[[parm]]$ub-LHSParms[[parm]]$lb)
    }
    if(LHSParms[[parm]]$type=="CAT"){
      masterDF[,parm] <- ceiling(LHSParms[[parm]]$lb + LHS[,parm]*(LHSParms[[parm]]$ub-LHSParms[[parm]]$lb))
    }
  }
  
  
  setwd(RESULTS_DIRECTORY)
  ## name file for LHS parameters 
  write.csv(masterDF,"masterDF_prelimCOTS.csv",row.names=F)
  
  return(masterDF)
}

################!
# MakeWorker ----
################!

MakeWorker = function (masterDF, npops, seasons){
  
  force(masterDF)
  force(npops)
  force(seasons)
  
  source("C://Users//jc312264//Documents//GitHub//COTS_Model//COTSModel_PrepareWorkspace.R", local=T)
  # dirs = list(BASE=BASE_DIRECTORY, CODE=CODE_DIRECTORY, RESULTS=RESULTS_DIRECTORY, 
  # DATA=DATA_DIRECTORY, ENVDATA=ENVDATA_DIRECTORY,SPATIALDATA=SPATIALDATA_DIRECTORY)
  setwd(BASE_DIRECTORY)
  load(file = "R_Objects/ConnMat.Rdata")
  setwd(CODE_DIRECTORY)
  source("COTSModel_LoadSmallObjectsForModelling.R", local=T)
  setwd(CODE_DIRECTORY)
  source("COTSModel_CoralFunctions.R", local=T)
  source("COTSModel_COTSFunctions.R", local=T)
  # ConnMat=ConnMat
  # SEASONS=seasons
  # PopData = PopData[1:npops,]
  # COTS.data = COTS.data[1:npops,]
  # data.grid = data.grid[1:npops,]
  # FvDParams=FvDParams
  #force(PopData);force(COTS.data); force(data.grid)
  
  browser()
  Worker = function(i){
    runModel(masterDF=masterDF, PopData=PopData,COTS.data = COTS.data, 
             data.grid = data.grid, ConnMat = ConnMat, npops=npops, rep=i, 
             seasons = seasons, FvDParams = FvDParams)
  }
  return(Worker)
}
