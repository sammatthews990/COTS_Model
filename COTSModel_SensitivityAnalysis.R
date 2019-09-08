# Sensitivity Anlaysis

# LOAD Validation Data ----
DIRECTORY = "C:/Users/jc312264/OneDrive - James Cook University/GitHub/COTS_Model"
setwd(DIRECTORY)
source("COTSModel_UtilityFunctions.R")
loadPackages()
load("Results/Archive/SensitivityAnalysis_2019-09-03_1.RData")

BESTREP = MPE.overall %>% 
  dplyr::arrange(Overall) %>%
  dplyr::top_n(1, Overall) %>% dplyr::pull(REP)

BESTPARAMS = masterDF[8,]

PARAMRANGE = rbind(BESTPARAMS*1.15, BESTPARAMS*0.85)

# Run Parameter Range on HPC

# Load Sensitivity Results ----
masterDF$avgPCF = signif(masterDF$avgPCF, 2)

df.sens = bind_cols(dplyr::select(masterDF, avgPCF, sdPCF, Pred, CCRatioThresh, CCRatioThresh2,
                                  selfseed, chl.int, CMax, J2M, J1M, J2R, J1R, AM, AR, Linf, K),
                    dplyr::select(MPE.overall, MPE.CC.ALL, ACC.Pres.ALL, KAP.Pres.ALL, ACC.Out.ALL, KAP.Out.ALL, Overall, ACC.CLASS))
# df.sens = df.sens %>% mutate_at(.vars = 1:16, .funs=scale)
df.sens$Overall = car::logit(df.sens$Overall)
vars.gbm = colnames(df.sens[1:16])#[-which(colnames(df.sens[1:16])%in% c("sdPCF","Pred","CCRatioThresh","CMax","AM","Linf"))]
gbm.sens = gbm.step(data = df.sens, gbm.x = vars.gbm, gbm.y = "Overall",tree.complexity = 3,
         family = "gaussian")
# options(scipen = 2)
# Optimizing for Coral Cover
gbm.sens.CC = gbm.step(data = df.sens, gbm.x = vars.gbm, gbm.y = "MPE.CC.ALL",tree.complexity = 3,
                    family = "gaussian")
gbm.sens.CC$contributions
# Optimizing for COTS Classification
gbm.sens.COTS = gbm.step(data = df.sens, gbm.x = vars.gbm, gbm.y = "ACC.CLASS",tree.complexity = 3,
                    family = "gaussian")
gbm.sens.COTS$contributions
gbm.sens$contributions
gbm.plot(gbm.sens.COTS, n.plots = 3)

# Plot Partial Dependency for Chapter ----
library(scales)
gbm.sens$contributions$var = factor(gbm.sens$contributions$var,as.character(gbm.sens$contributions$var))
gg_part = gg_partial_plot(gbm.sens,gbm.sens$contributions$var[1:6], 
                          var.labs = c("Pred_Larv", "MortJ1_k", "MortJ1_x0", "CCRatio_2","Fert_Linf", "Fec_Max")) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  theme_bw(base_size=14)

gg_imp = importance_plot(gbm.sens) + theme_light(base_size = 14) + 
  # scale_size_manual(values=3)+
  scale_x_discrete("Model Parameters", 
                   labels = c(expression(MortA[x[0]]),expression(CCRatio[1]),expression(MortJ2[x[0]]),
                              expression(MortA[k]),expression(Fec[sd]),expression(Fert[K]),  
                              expression(SelfSeed),expression(Chl[intercept]), expression(MortJ2[k]),
                              expression(Cons[max]),  expression(Fec[max]), expression(Fert[Linf]),
                              expression(CCRatio[2]),expression(MortJ1[x[0]]),expression(MortJ1[k]), expression(Pred[Larv])))

gridExtra::grid.arrange(gg_imp,gg_part,
                        widths = c(1, 2))
export::graph2tif(file = "Results/Figures/Sensitivity/SensitivityResults", width=8, height=8, dpi=300)
export::graph2svg(file = "Results/Figures/Sensitivity/SensitivityResults", width=8, height=8)

# Look at interactions in 
gbm.sens.int = gbm.interactions(gbm.sens)
gbm.sens.int$interactions
gbm.sens.int$rank.list

# Find mean oubreak Density
COTS.dens = data.manta.valid %>% dplyr::filter(MEAN_COTS >0.22) %>%
  dplyr::group_by(SECTOR) %>% dplyr::summarise(mn = mean(MEAN_COTS))
COTS.dens.mod = data.manta.valid %>% dplyr::filter(COTS.mn >0.22) %>%
  dplyr::group_by(SECTOR) %>% dplyr::summarise(mn.mod = mean(COTS.mn)) %>% 
  left_join(COTS.dens) %>% dplyr::mutate(delta = mn.mod-mn)
# GBR level
COTS.dens.GBR = data.manta.valid %>% dplyr::filter(MEAN_COTS >0.22) %>%
  dplyr::summarise(mn = mean(MEAN_COTS)) %>% mutate(GBR = "GBR")
COTS.dens.mod.GBR = data.manta.valid %>% dplyr::filter(COTS.mn >0.22) %>%
  dplyr::summarise(mn.mod = mean(COTS.mn)) %>% mutate(GBR = "GBR") %>%
  left_join(COTS.dens.GBR) %>% dplyr::mutate(delta = mn.mod-mn)
