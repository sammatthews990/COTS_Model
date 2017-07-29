
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
  loadPackage("rgdal")          # for reading and writing all sorts of spatial data   
  loadPackage("popbio")         # for doing basic matrix population modeling
  loadPackage("tidyverse")      # data manipulation
  loadPackage("rgeos")          # geometry applications
  loadPackage("plyr")           # data wrangling
  loadPackage("dplyr")          # data wrangling
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
  loadPackage("mbcv")
  loadPackage("gbm")
  loadPackage("MASS")
  loadPackage("nlme")
  loadPackage("lme4")
  
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

saveWorkspace <- function(objNames=NULL,filename="R_workspace"){
  if(is.null(objNames)) objNames = ls(all = TRUE,name=.GlobalEnv)
  setwd(RDATA_DIRECTORY)
  filename <- sprintf("%s_%s.RData",filename,Sys.Date())
  save(list = objNames, file = filename)
}
 

#function to ignore NA's when summing
plus <- function(x) {
  if(all(is.na(x))){
    c(x[0],NA)} else {
      sum(x,na.rm = TRUE)}
}


