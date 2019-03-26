# COTSMod: Spatially Explicit COTS-Coral Metapopulation Model ----

# Load Pre Defined Data ----

rm(list=ls())

load(...)
source("COTSModel_LoadObjectsForModelling.R")

# Load Functions ----

source(...)

setwd(CODE_DIRECTORY)
source("COTSModel_Utilityfunctions.R")   # load utility functions, e.g., for loading packages etc.
source("COTSModel_COTSfunctions.R")      # load functions for implementing COTS demography and dispersal
source("COTSModel_Coralfunctions.R")     # load functions for implemeting coral growth/recovery and dispersal
source("COTSModel_GISfunctions.R")     # load functions for implemeting coral growth/recovery and dispersal

library(...)

`%>%` <- magrittr::`%>%` # bringt the pipe operator into the global environment

# Set Global Parameters ----

NREPS <- 10      
NYEARS <- 22
NSEASONS <- 2
seasons <- c("summer","winter")
npops = 15802 #number of reefs we want to test
nsimul <- 100

VERBOSE <- TRUE        # flag whether functions should return detailed information
DEBUG <- TRUE          # flag whether to output debug files etc. 

projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   

stagenames <- c('J_1', 'J_2', 'A')
COTSmort <- c(0.8,0.7,0.2)
names(COTSmort) <- stagenames
COTSremain <- c(0.02,0.2,1) 
names(COTSremain) <- stagenames
COTS_StableStage <- c(0.9803, 0.0171, 0.0026)   

# Resample Disturbance ----


# Save Workspace for import within each thread ----

saveWorkspace(filename="GlobalParams.RData", dir=RDATA_DIRECTORY)

# PARALLEL


load(...)










