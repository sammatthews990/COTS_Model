#################!
# CONNECTIVITY MATRICES ----
#################!


source("C://Users//jc312264//Documents//GitHub//COTS_Model//COTSModel_PrepareWorkspace.R")
loadPackage("R.matlab")

setwd(DATA_DIRECTORY)
ConnMat = R.matlab::readMat("cotsconn/cotsconn.mat")
