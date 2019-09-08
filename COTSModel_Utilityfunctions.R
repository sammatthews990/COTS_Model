
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
  # loadPackage("leaflet")
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
  loadPackage("Matrix")  # Sparse matrices
  loadPackage("pscl")
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
  LHSParms <- specifyLHSParam(paramslist=LHSParms,name="SexRatio",type="CAT",lb=5,ub=5)
  
  ####  WINTER CONSUMPTION RATE
  LHSParms <- specifyLHSParam(LHSParms,"ConsRateW",type="CONT",lb=0,ub=0)
  
  #### SUMMER CONSUMPTION RATE
  LHSParms <- specifyLHSParam(LHSParms,"ConsRateS",type="CONT",lb=0,ub=0)
  
  ### AVERAGE PER CAPITA FECUNDITY 
  LHSParms <- specifyLHSParam(LHSParms,"avgPCF",type="CONT",lb=2.82e7,ub=3.81e7)       # KTS: changed to 500
  
  ### STD DEV PER CAPITA FECUNDITY
  LHSParms <- specifyLHSParam(LHSParms,"sdPCF",type="CONT",lb=1.15e7,ub=0.85e7)     # was 1000 to 40000  
  
  ### % JUVENILE1 MORTALITY PER TIME STEP
  LHSParms <- specifyLHSParam(LHSParms,"mortJ1",type="CONT",lb=0.98,ub=0.98)
  
  ### % JUVENILE2 MORTALITY PER TIME STEP
  LHSParms <- specifyLHSParam(LHSParms,"mortJ2",type="CONT",lb=0.7,ub=0.7) 
  
  ### % ADULT MORTALITY PER TIME STEP
  LHSParms <- specifyLHSParam(LHSParms,"mortA",type="CONT",lb=0.43,ub=0.58)
  
  ### % JUVENILE1 TO REMIAIN during transition
  LHSParms <- specifyLHSParam(LHSParms,"remJ1",type="CONT",lb=0,ub=0)
  
  ### % JUVENILE2 TO REMIAIN during transition 
  LHSParms <- specifyLHSParam(LHSParms,"remJ2",type="CONT",lb=0,ub=0)       # KTS: changed to 500
  
  ### % ADULT TO REMIAIN during transition
  LHSParms <- specifyLHSParam(LHSParms,"remA",type="CONT",lb=1,ub=1)     # was 1000 to 40000  
  
  ### PROPORTIONAL STABLE STAGE DISTRIBUTION J1
  LHSParms <- specifyLHSParam(LHSParms,"cssJ1",type="CONT",lb=0.9803,ub=0.9803)
  
  ### PROPORTIONAL STABLE STAGE DISTRIBUTION J2
  LHSParms <- specifyLHSParam(LHSParms,"cssJ2",type="CONT",lb=0.0171,ub=0.0171) 
  
  ### PROPORTIONAL STABLE STAGE DISTRIBUTION a
  LHSParms <- specifyLHSParam(LHSParms,"cssA",type="CONT",lb=0.0026,ub=0.0026)
  
  ### PROPORTIONAL Larval Mortality FIXED
  LHSParms <- specifyLHSParam(LHSParms,"Pred",type="CONT",lb=0.995,ub=0.95)
  
  ### % Coral Cover at which to crash a COTS population entirely  
  LHSParms <- specifyLHSParam(LHSParms,"Crash",type="CONT",lb=0,ub=0)

  ### Number of years of high densities COTS after which to force crash population --CRUDE
  LHSParms <- specifyLHSParam(LHSParms,"OutbreakCrash",type="CONT",lb=Inf,ub=Inf)
  
  ### Lowest possible fecundity as a proportion of PCF (Per Capita Fecundity)  
  LHSParms <- specifyLHSParam(LHSParms,"Fbase",type="CONT",lb=0,ub=0)
  
  ### Coral-COTS Ratio at which COral becomes lmiting resource ~25
  LHSParms <- specifyLHSParam(LHSParms,"CCRatioThresh",type="CONT",lb=40,ub=30)
  
  ### Coral-COTS Ratio below which COTS mortality reaches its maximum
  LHSParms <- specifyLHSParam(LHSParms,"CCRatioThresh2",type="CONT",lb=4.13,ub=5.58)
  
  ### Maximum mortality experienced by adult COTS under severely resource limited conditions
  LHSParms <- specifyLHSParam(LHSParms,"maxmort",type="CONT",lb=1,ub=1)
  
  ### Proportion of estimated selfseeding larvae to remiain (assuming that KArlo's model overestimates)
  LHSParms <- specifyLHSParam(LHSParms,"selfseed",type="CONT",lb=0.85,ub=1)
  
  ### Intercept for the chlorophyll model, increasing it allows more larve to survive in oligotrophic cond's
  LHSParms <- specifyLHSParam(LHSParms,"chl.int",type="CONT",lb=1.7,ub=2.4)
  
  ### Proportion of CMax to be consumed at lowest Coral-COTS Ratio
  LHSParms <- specifyLHSParam(LHSParms,"Cbase",type="CONT",lb=0.115,ub=0.085)
  
  ### Max coral consumed (cm2) per day by COTS
  LHSParms <- specifyLHSParam(LHSParms,"CMax",type="CONT",lb=170,ub=230)
  
  ### Logistic parameter for J2 mortality
  LHSParms <- specifyLHSParam(LHSParms,"J2M",type="CONT",lb=-0.24e5,ub=-0.33e5)
  
  ### Logistic parameter for J1 mortality
  LHSParms <- specifyLHSParam(LHSParms,"J1M",type="CONT",lb=-1.7e7,ub=-2.3e7)
  
  ### Logistic parameter for J2 mortality
  LHSParms <- specifyLHSParam(LHSParms,"J2R",type="CONT",lb=0.0000081,ub=0.000011)
  
  ### Logistic parameter for J1 mortality
  LHSParms <- specifyLHSParam(LHSParms,"J1R",type="CONT",lb=1.69e-07,ub=2.29e-07)
  
  ### Logistic parameter for A mortality
  LHSParms <- specifyLHSParam(LHSParms,"AM",type="CONT",lb=115,ub=85)
  
  ### Logistic parameter for A mortality
  LHSParms <- specifyLHSParam(LHSParms,"AR",type="CONT",lb=0.039,ub=0.052)
  
  ### Von Bertanlanffy Growth of Fertilisation by density - MAX
  LHSParms <- specifyLHSParam(LHSParms,"Linf",type="CONT",lb=0.805,ub=0.595)
  
  ### on Bertanlanffy Growth of Fertilisation by density - Growth Rate
  LHSParms <- specifyLHSParam(LHSParms,"K",type="CONT",lb=3.7e-04,ub=5e-04)
  
  ### on Bertanlanffy Growth of Fertilisation by density - Value at 0
  LHSParms <- specifyLHSParam(LHSParms,"t0",type="CONT",lb=0,ub=0)
  
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
  
  
  #setwd(RESULTS_DIRECTORY)
  ## name file for LHS parameters 
  #write.csv(masterDF,"masterDF_prelimCOTS.csv",row.names=F)
  
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

inv.logit = function(x) {
  exp(x)/(1+exp(x))
}


#' @title R2 for models with zero-inflation
#' @name r2_zeroinflated
#'
#' @description Calculates R2 for models with zero-inflation component, including mixed effects models.
#'
#' @param model A model.
#' @param method Indicates the method to calculate R2. See 'Details'. May be abbreviated.
#'
#' @return For the default-method, a list with the R2 and adjusted R2 values.
#'   For \code{method = "correlation"}, a named numeric vector with the
#'   correlation-based R2 value.
#'
#' @details The default-method calculates an R2 value based on the residual
#'   sums of squares (using Pearson residuals), divided by the total sum of
#'   squares. For \code{method = "correlation"}, R2 is a correlation-based measure,
#'   which is rather crude. It simply computes the squared correlation between the
#'   model's actual and predicted reponse.
#'
#' @examples
#' library(pscl)
#' data(bioChemists)
#' model <- zeroinfl(
#'   art ~ fem + mar + kid5 + ment | kid5 + phd,
#'   data = bioChemists
#' )
#'
#' r2_zeroinflated(model)
#'
#' @importFrom stats cor predict residuals fitted
#' @importFrom insight model_info get_response find_parameters n_obs
#' @export
r2_zeroinflated <- function(model, method = c("default", "correlation")) {
  method <- match.arg(method)
  
  mi <- insight::model_info(model)
  if (!mi$is_zero_inflated) {
    warning("Model has no zero-inflation component.")
  }
  
  if (method == "default")
    .r2_zi_default(model)
  else
    .r2_zi_correlation(model)
}


.r2_zi_correlation <- function(model) {
  r2_zi <- stats::cor(insight::get_response(model), stats::predict(model, type = "response"))^2
  names(r2_zi) <- "R2 for ZI-models"
  r2_zi
}


.r2_zi_default <- function(model) {
  n <- insight::n_obs(model)
  p <- length(insight::find_parameters(model)[["conditional"]])
  
  r2_zi <- 1 - (sum(stats::fitted(model)^2) /
                  (sum(stats::fitted(model)^2) + sum(stats::residuals(model, type = "pearson")^2)))
  r2_zi_adj <- 1 - (1 - r2_zi) * (n - 1) / (n - p - 1)
  
  out <- list(R2 = r2_zi, R2_adjusted = r2_zi_adj)
  
  names(out$R2) <- "R2"
  names(out$R2_adjusted) <- "adjusted R2"
  
  attr(out, "model_type") <- "Zero-Inflated and Hurdle"
  structure(class = "r2_generic", out)
}


df1 = data.frame(MT = c(0,0.1,0.22,1), DENS = c(0,3500,4900,11000))

MTCalib.gam = lm(DENS~sqrt(MT), data=df1)
MTCalib.gaminv = lm(sqrt(MT)~DENS, data=df1)
