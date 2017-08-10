##########################
#  R script: prepare workspace for modeling COTS outbreaks in the Great Barrier Reef 
#
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
#  
#  ***IMPORTANT*** This is the ONLY script where you should need to reference local file trees. All other
#                  subsequent scripts should be able to be run without referencing user-specific directories
# 
#  08 April 2015 -- started scripting

##########################


#######################
#   CLEAR THE WORKSPACE 
#######################

#rm(list=ls())


#######################
#   SET USER 
#######################

USER = "SAM"
# USER = "SAM_RENO"


####################
#  SET GLOBAL VARIABLES  (USER SPECIFIED PARAMS)
####################

# NREPS <- 10      
NYEARS <- 20
NSEASONS <- 2
SEASONS <- c("summer","winter")
npops = 50 #number of reefs we want to test

VERBOSE <- TRUE        # flag whether functions should return detailed information
DEBUG <- TRUE          # flag whether to output debug files etc. 

projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   #"+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

stagenames <- c('J_1', 'J_2', 'A')
COTSmort <- c(0.8,0.7,0.2)
names(COTSmort) <- stagenames
COTSremain <- c(0.02,0.2,1) # proportion remianing in each life stage --> we can make this a function of resource
names(COTSremain) <- stagenames
#VBG.Params = VBG.Models[[2]]

COTS_StableStage <- c(0.9803, 0.0171, 0.0026)   # very approximate stable stage distribution (J1, J2, Adult: see below for back-of-the-envelope calculation)



#########################
# SET PROJECT DIRECTORIES (this should be the only place where local directories should be referenced)
#########################

if(USER=="KEVIN") BASE_DIRECTORY <- "C:\\Users\\Kevin\\Dropbox\\CoTS_Model"             # NOTE: this should link to the Dropbox folder with shared project resources	                                                                        
if(USER=="KEVIN") CODE_DIRECTORY <- "C:\\Users\\Kevin\\GIT\\COTS_Model"              # NOTE: code directory should be your local copy of the GitHub repository

if(USER=="SAM") BASE_DIRECTORY <- "C:\\Users\\jc312264\\Dropbox\\CoTS_Model"
if(USER=="SAM") CODE_DIRECTORY <- "C:\\Users\\jc312264\\Documents\\GitHub\\COTS_Model"

if(USER=="SAM_RENO") BASE_DIRECTORY <- "C:\\Users\\kshoemaker\\Dropbox\\CoTS_Model"
if(USER=="SAM_RENO") CODE_DIRECTORY <- "C:\\Users\\kshoemaker\\Documents\\GitHub\\COTS_Model"

SPATIALDATA_DIRECTORY <- paste(BASE_DIRECTORY,"\\Spatial Layers",sep="")                          # directory for storing relevant spatial data (ASC, SHP files)
if(is.na(file.info(SPATIALDATA_DIRECTORY)[1,"isdir"])) dir.create(SPATIALDATA_DIRECTORY)

DATA_DIRECTORY <- paste(BASE_DIRECTORY,"\\Data",sep="")                                 # directory for storing data (CSV files)
if(is.na(file.info(DATA_DIRECTORY)[1,"isdir"])) dir.create(DATA_DIRECTORY)

ENVDATA_DIRECTORY <- paste(DATA_DIRECTORY,"\\Environmental",sep="")                                 # directory for storing data (CSV files)
if(is.na(file.info(ENVDATA_DIRECTORY)[1,"isdir"])) dir.create(ENVDATA_DIRECTORY)

FIGURES_DIRECTORY <- paste(BASE_DIRECTORY,"\\Figures\\RawFigures",sep="")               # directory for storing raw figures generated from R
if(is.na(file.info(FIGURES_DIRECTORY)[1,"isdir"])) dir.create(FIGURES_DIRECTORY)

RESULTS_DIRECTORY <- paste(BASE_DIRECTORY,"\\results",sep="")                           # directory for storing relevant results
if(is.na(file.info(RESULTS_DIRECTORY)[1,"isdir"])) dir.create(RESULTS_DIRECTORY)

RDATA_DIRECTORY <- paste(BASE_DIRECTORY,"\\R_Workspaces",sep="")                        # directory for storing .RData files (R workspaces and data objects)
if(is.na(file.info(RDATA_DIRECTORY)[1,"isdir"])) dir.create(RDATA_DIRECTORY)

cat(sprintf("The current user is %s",USER))

####################
#  LOAD FUNCTIONS AND SCRIPTS FROM SOURCE CODE (your local GitHub repository)
#################### 

setwd(CODE_DIRECTORY)
source("COTSModel_Utilityfunctions.R")   # load utility functions, e.g., for loading packages etc. 
#source("COTSModel_COTSfunctions.R")      # load functions for implementing COTS demography and dispersal
#source("COTSModel_Coralfunctions.R")     # load functions for implemeting coral growth/recovery and dispersal
source("COTSModel_GISfunctions.R")     # load functions for implemeting coral growth/recovery and dispersal

#############################
#  LOAD PACKAGES
#############################
# note: 'loadPackage' should install the package from CRAN automatically if it is not already installed

loadPackages()   # load all packages into the global environment 
# THis only installs packages now.. everything else will be defined inline
#library(dplyr)
#library(ggplot2)
`%>%` <- magrittr::`%>%`

# save global params
saveWorkspace(filename="GlobalParams.RData", dir=RDATA_DIRECTORY)



#############
#  END SCRIPT
#############




