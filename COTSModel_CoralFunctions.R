##########################
#  Functions for modeling coral dynamics in the great barrier reef as part of a framework for modeling COTS outbreaks  
#
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
# 
#  05 April 2015 -- started scripting
#  28 September 2016 -- intialize coral cover function
#  

  ## for coral dispersal: see kinlan and Gaines paper


###################
# initializeCoralCover
##########
# OBJECTIVE:
#    generate an matrix for storing the coral cover in each pixel. 
# PARAMS:
#    - reefmap: raster template for the study region: NA outside of reefs, reef ID value within reefs 
#    - initCoralCover: for every pixel in the study area, a vector of initial coral cover
#
# RETURNS:
#    - CoralCover: spatially-structured Coral Cover
#           CoralCover[,'F']: vector representing spatially structured cover for fast growing corals
#           CoralCover[,'S']: vector representing spatially structured abundance of Juvenile stage 2 individuals
#           COTSabund$A: vector representing spatially structured abundance of reproductive adult individuals
#           COTSabund$S: vector representing spatially structured abundance of senile adult individuals
#           NOTE: larvae are not considered explicitly here. 
###################



#########################