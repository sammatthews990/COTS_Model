##########################
#  Functions and Code to Visualise Outputs  
#
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
# 
#  14 Novemeber 2018 -- started scripting
#  15 November-- view timestep with leaflet --- NEEEDS COVRTE TO FUNCTION
#  



#################!
# ViewLeaflet ----
#################!
# OBJECTIVE:
#    View coral cover + COTS numbers at a given time step using the leaflet package 
# PARAMS:
#    - Year: Which year to view data for
#    - reefmap: GBRMPA reef outlines 
#    - CoralCover: Estimtated Coral Cover for the time step
#    - COTS: Estimated COTS numbers for the time step
#    - Disturbance: Data frame containing disturbance data for the timestep
#
# RETURNS:
#    - Interactive map of coral cover and COTS numbers
###################

library(leaflet)

pal = colorNumeric(palette = "RdYlBu",
                   domain=CoralCoverParams[[2]])


leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group= "ESRI World Imagery") %>%
  addPolygons(data=shp, group="Reefs", label=~LOC_NAME_S,
              weight = 1,
              highlightOptions = highlightOptions(color = "white", weight = 2, bringToFront = FALSE)) %>%
  addCircles(data = data.grid, lng=~data.grid$lon, lat=~data.grid$lat, group="HCINI",
           radius=5, label = ~paste(round(CoralCoverParams[[3]][,1]),1),
           color =~pal(CoralCoverParams[[2]]),
           fillOpacity = ~0.5) %>%
  addMiniMap(tiles = providers$Esri.WorldImagery, width =150, height=150)

CAvailable = (30*data.grid$PercentReef/10000)*1e6*1e4 # in cm2
CConsumed = 300*1000*365
CRemaining=((CAvailable-CConsumed)/1e10)*(10000/data.grid$PercentReef)
CRemaining[CRemaining < 0.5] <- 0.5

doCoralGrowth = function(CoralCover, B0, WQ, HCMAX) {
  b0.wq <- B0 + WQ * rnorm(length(WQ), mean=WQ.mn.sd[1], sd=WQ.mn.sd[2])
  b1.wq <- b0.wq / log(HCMAX)
  CoralCover <- log(CoralCover)
  CoralCover <- b0.wq + (1 - b1.wq)*CoralCover
  return(exp(CoralCover))
}
doCoralGrowth(CRemaining, B0 = CoralCoverParams[[4]][,1], WQ = CoralCoverParams[[1]]$WQ, HCMAX = CoralCoverParams[[3]][,1])
