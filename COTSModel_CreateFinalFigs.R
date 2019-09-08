# Make FInal PLots For Publication

# Load Data & Functions ----
DIRECTORY = "C:/Users/jc312264/OneDrive - James Cook University/GitHub/COTS_Model"
setwd(DIRECTORY)
source("COTSMod_PlotSpatial.R")
source("COTSModel_UtilityFunctions.R")
loadPackages()
load("Results/Archive/Validation_2019-08-29_5.RData")

# Make fig of  SECTOR COTS/Manta Tow

data.manta.plot = data.manta.valid.all %>%
  group_by(A_SECTOR, Year) %>% 
  dplyr::summarise(MEAN = mean(MEAN_COTS, na.rm=T),
                   SE=sd(MEAN_COTS, na.rm=T)/sqrt(n()),
                   n = n()) %>% 
  dplyr::filter(A_SECTOR %in% c("CL", "CA", "TO", "WH", "SW", "CB"))

data.manta.plot$SECTOR = factor(data.manta.plot$A_SECTOR, levels = c("CL", "CA", "TO", "WH", "SW", "CB"),
                                labels = c("Cooktown/Lizard Island", "Cairns", "Townsville","Whitsundays", "Swains", "Capricorn Bunker"))
GG.SECTOR = ggplot(data.manta.plot, aes(x=Year, y=MEAN, fill=SECTOR)) +
  geom_point(aes(y=MEAN)) + 
  # geom_line(aes(x=REPORT_YEAR-1, y=MEAN_LIVE.corr)) +
  geom_line() + geom_point() + 
  geom_ribbon(aes(ymin = MEAN-SE, ymax=MEAN+SE,colour=SECTOR),alpha=0.4) +
  facet_wrap(~SECTOR, scales = "free_y", ncol = 2) + ylab("Mean COTS/Manta Tow (+/- SE)") +
  theme_bw(base_size=14) + theme(legend.position = "none") + xlim(c(1996,2017))
GG.SECTOR

export::graph2tif(file = "Results/Figures/SectorPlot", width =6, height=9, dpi=300)


# Filter Validation Table for Best Replicate ----
BESTREP = MPE.overall %>% 
  dplyr::arrange(Overall) %>%
  dplyr::top_n(1, Overall) %>% dplyr::pull(REP)

BESTVALIDATION = MPE %>% dplyr::filter(REP %in% BESTREP)
write.csv(BESTVALIDATION, file = "Results/Archive/Validation_2019-08-29_5.csv", row.names = F)
# Calculate % Coral Loss Per Pixel ----
normalized = function(x){(x-min(x))/(max(x)-min(x))}

res.plot = res.plot %>% dplyr::filter(REP %in% BESTREP)
coral.delta = res.plot %>% ungroup() %>%
  dplyr::select(REEF_NAME, Year, CC.mn) %>% 
  dplyr::filter(Year %in% c(1996, 2017)) %>%
  distinct() %>%
  tidyr::spread(Year, CC.mn) %>%
  mutate(delta.CC = (`2017`-`1996`),
         delta.CC.pa = delta.CC/22,
         delta.CC.norm = normalized(delta.CC.pa))



# Calculate Max and Mean COTS Abundance per Pixel----
cots.delta =  res.plot %>% ungroup() %>%
  dplyr::select(REEF_NAME, Year, COTS.mn) %>% dplyr::group_by(REEF_NAME) %>%
  dplyr::summarise(MEAN_COTS = mean(COTS.mn),
                   MAX_COTS = max(COTS.mn))

# Make Final Figures for 6 Reefs ----

# CORAL ----

resreefs = c("Martin Reef (14-123)",
             "Linnet Reef (14-126)",
             "Hastings Reef (16-057)", 
             "Moore Reef (16-071)",
             "Bowden Reef (19-019)",
             "East Cay Reef (21-305)")

res.fig = dplyr::select(data.grid, REEF_NAME, SECTOR, CROSS_SHELF, lat, lon) %>% 
  inner_join(cots.delta) %>% inner_join(coral.delta)
res.fig.GBR = res.plot %>% group_by(Year) %>%
  dplyr::summarise(CC.md = median(CC.mn),
                   CC.Q05 = quantile(CC.mn, probs=0.05),
                   CC.Q25 = quantile(CC.mn, probs=0.25),
                   CC.Q75 = quantile(CC.mn, probs=0.75),
                   CC.Q95 = quantile(CC.mn, probs=0.95),
                   COTS.md = mean(COTS.mn),
                   COTS.se = sd(COTS.mn)/sqrt(n()),
                   COTS.Q05 = quantile(COTS.mn, probs=0.05),
                   COTS.Q25 = quantile(COTS.mn, probs=0.25),
                   COTS.Q75 = quantile(COTS.mn, probs=0.75),
                   COTS.Q95 = quantile(COTS.mn, probs=0.95))


dev.off()
f = colorRamp(c("red", "yellow", "green"))
layout(rbind(c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(3,3,3,3)))

# PLOT CORAL COVER
color = rgb(f(res.fig$delta.CC.norm)/255)
par(mai=c(0,0,0,0))
plotPolys(shape45[1:1000,], col="gray95", border="gray70", xlim=c(148,158), ylim=c(-36, -21),
          bg="white", cex=1.2, axes=F, ylab=NA, plt=c(0,1,0,1), colHoles=NA, projection = 1)
lines(gridln.elide, col="lightgrey", xlim=c(148,158), ylim=c(-36, -21))
addPolys(shape45, col="gray95", border="gray70", bg="white", xlim=c(148,158), ylim=c(-36, -20),
         colHoles=NA)
points(grid45$x, grid45$y, pch=19, cex=0.5, col=addalpha(color, .5))
northarrow(c(150,-23.5), 1, bearing =5.5, cex = 1.8)

# text(149.9,-22.5,i, cex=4, font=2)

text(155.3, -22.5, "150E", col="lightgrey", srt=40, cex=1)
text(155.3, -29.4, "155E", col="lightgrey", srt=40, cex=1)
text(155.3, -26, "15S", col="lightgrey", srt=-40, cex=1)
text(155.3, -33.2, "20S", col="lightgrey", srt=-40, cex=1)

# pnts = cbind(x =c(149.3,150.2, 150.2,149.3), y =c(-35.6, -35.6,-34,-34))
pnts = cbind(x =c(153.2,154.2, 154.2,153.2), y =c(-28.8, -28.8,-27,-27))
# rect(149.2, -35.8, 151.1, -33.2, col="white", border = "black")
# rect(153, -29.1, 155.3, -26.2, col="white", border = "black")
color.grad = rgb(f(seq(0,1,0.01))/255)
SDMTools::legend.gradient(pnts,color.grad,
                          c(round(min(res.fig$delta.CC.pa),1),
                            round(max(res.fig$delta.CC.pa),1)), 
                          title = "% Coral Loss p.a.", cex=1.8)

# Make GGPLOT FIGS
res.plot.reef = data.manta.valid %>%
  filter(REP %in% BESTREP &
           REEF_NAME %in% resreefs)
res.plot.reef$REEF_NAME = factor(res.plot.reef$REEF_NAME, 
                                    levels = c("Martin Reef (14-123)",
                                    "Linnet Reef (14-126)",
                                    "Hastings Reef (16-057)", 
                                    "Moore Reef (16-071)",
                                    "Bowden Reef (19-019)",
                                    "East Cay Reef (21-305)"))

GG.CORAL = ggplot(res.plot.reef, aes(x=Year, y=CC.mn, fill=REEF_NAME)) +
  geom_point(aes(y=MEAN_LIVE.corr)) + geom_line(aes(x=REPORT_YEAR-1, y=MEAN_LIVE.corr)) +
  geom_line(aes(colour=as.factor(REEF_NAME))) + geom_ribbon(aes(ymin = CC.25, ymax=CC.75,colour=as.factor(REEF_NAME)),alpha=0.4) +
  facet_wrap(~REEF_NAME, scales = "free_y", ncol = 2) + ylim(c(0,60)) + ylab("Mean % Coral Cover") +
  theme_bw(base_size=14) + theme(legend.position = "none") + xlim(c(1996,2017))

GG.GBR = ggplot(res.fig.GBR, aes(x=Year, y=CC.md)) +
  geom_line(aes(x=Year, y=CC.md)) +
  geom_ribbon(aes(ymin = CC.Q05, ymax=CC.Q95),alpha=0.2) +
  geom_ribbon(aes(ymin = CC.Q25, ymax=CC.Q75),alpha=0.5) +
  ylab("Mean % Coral Cover") +
  theme_bw(base_size=14) + theme(legend.position = "none",
                                 axis.title.x = element_blank()) + xlim(c(1996,2017))


# GG.CORAL
## PLot GGPLOTs using viewport
vp<-grid::viewport(x=0.7,y=0.6,width=0.6,height = 0.8)
print(GG.CORAL, vp = vp)
vp<-grid::viewport(x=0.5,y=0.1,width=1,height = 0.2)
print(GG.GBR, vp = vp)


export::graph2tif(file = "Results/Figures/FinalPlots_2019-09-02", width =11, height=9, dpi=300)

# COTS ----

dev.off()
f = colorRamp(c("green", "yellow", "red"))
layout(rbind(c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(3,3,3,3)))

# PLOT CORAL COVER
color = ifelse(res.fig$MAX_COTS<0.01,"lightblue",
               ifelse(res.fig$MAX_COTS <0.11,"green",
                      ifelse(res.fig$MAX_COTS <0.22, "yellow",
                             ifelse(res.fig$MAX_COTS <1, "orange", "red"))))
# color = rgb(f(res.fig$MEAN_COTS/max(res.fig$MEAN_COTS))/255)
par(mai=c(0,0,0,0))
plotPolys(shape45[1:1000,], col="gray95", border="gray70", xlim=c(148,158), ylim=c(-36, -21),
          bg="white", cex=1.2, axes=F, ylab=NA, plt=c(0,1,0,1), colHoles=NA, projection = 1)
lines(gridln.elide, col="lightgrey", xlim=c(148,158), ylim=c(-36, -21))
addPolys(shape45, col="gray95", border="gray70", bg="white", xlim=c(148,158), ylim=c(-36, -20),
         colHoles=NA)
points(grid45$x, grid45$y, pch=19, cex=0.5, col=addalpha(color, .5))
northarrow(c(150,-23.5), 1, bearing =5.5, cex=1.8)

# text(149.9,-22.5,i, cex=4, font=2)

text(155.3, -22.5, "150E", col="lightgrey", srt=40, cex=1)
text(155.3, -29.4, "155E", col="lightgrey", srt=40, cex=1)
text(155.3, -26, "15S", col="lightgrey", srt=-40, cex=1)
text(155.3, -33.2, "20S", col="lightgrey", srt=-40, cex=1)

# pnts = cbind(x =c(149.3,150.2, 150.2,149.3), y =c(-35.6, -35.6,-34,-34))
pnts = cbind(x =c(153.2,154.2, 154.2,153.2), y =c(-28.8, -28.8,-27,-27))
# rect(149.2, -35.8, 151.1, -33.2, col="white", border = "black")
# rect(153, -29.1, 155.3, -26.2, col="white", border = "black")
color.grad = rgb(f(seq(0,1,0.01))/255)
legend(x='bottomleft', legend = c("NC", "NO", "PO", "EO", "SO"), 
       fill = c("lightblue", "green", "yellow", "orange", "red"), cex = 1)

# Make GGPLOT FIGS
res.plot.reef = data.manta.valid %>%
  filter(REP %in% BESTREP &
           REEF_NAME %in% resreefs)
res.plot.reef$REEF_NAME = factor(res.plot.reef$REEF_NAME, 
                                 levels = c("Martin Reef (14-123)",
                                            "Linnet Reef (14-126)",
                                            "Hastings Reef (16-057)", 
                                            "Moore Reef (16-071)",
                                            "Bowden Reef (19-019)",
                                            "East Cay Reef (21-305)"))

GG.COTS = ggplot(res.plot.reef, aes(x=Year, y=COTS.mn, fill=as.factor(REEF_NAME))) +
  geom_point(aes(y=MEAN_COTS)) + geom_line(aes(x=REPORT_YEAR-1, y=MEAN_COTS)) +
  geom_line(aes(colour=as.factor(REEF_NAME))) + geom_ribbon(aes(ymin = COTS.25, ymax=COTS.75,colour=as.factor(REEF_NAME)),alpha=0.4) +
  facet_wrap(~REEF_NAME, scales = "free_y", ncol = 2)  + ylab("Mean COTS/Manta Tow") +
  theme_bw(base_size=14) + theme(legend.position = "none") + xlim(c(1996,2017))

GG.GBR.COTS = ggplot(res.fig.GBR, aes(x=Year, y=COTS.md)) +
  geom_line(aes(x=Year, y=COTS.md)) +
  # geom_ribbon(aes(ymin = COTS.Q05, ymax=COTS.Q95),alpha=0.2) +
  geom_ribbon(aes(ymin = COTS.md-COTS.se, ymax=COTS.md+COTS.se),alpha=0.3) +
  ylab("Mean COTS/Manta Tow") +
  theme_bw(base_size=14) + theme(legend.position = "none", 
                                 axis.title.x = element_blank()) + xlim(c(1996,2017))


# GG.CORAL
## PLot GGPLOTs using viewport
vp<-grid::viewport(x=0.7,y=0.6,width=0.6,height = 0.8)
print(GG.COTS, vp = vp)
vp<-grid::viewport(x=0.5,y=0.1,width=1,height = 0.2)
print(GG.GBR.COTS, vp = vp)


export::graph2tif(file = "Results/Figures/FinalPlotsCOTS_2019-09-02", width =11, height=9, dpi=300)

