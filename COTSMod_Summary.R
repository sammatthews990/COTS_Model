
load("RData/ModelWorkspace_FULL.RData")
DIRECTORY = getwd()
setwd("Results")
load("Sample_1.Rdata")
DIRECTORY

ForGraph = ResultsDash
ForGraph$REP=1
files = list.files()
files = files[grep("Sample",files)]
NREPS = nrow(masterDF)

for (i in 2:NREPS) {
  load(sprintf("Sample_%s.Rdata", i))
  ResultsDash$REP=i
  # ForDashboard = cbind(ForDashboard, ResultsDash[7:14])
  ForGraph = rbind(ForGraph, ResultsDash)
}

setwd(DIRECTORY)

ForGraph = ForGraph %>% 
  mutate(COTS.mn = round(predict(MTCalib.gaminv, newdata=data.frame(DENS=COTS.mn))^2,3),
         COTS.Q50 = round(predict(MTCalib.gaminv, newdata=data.frame(DENS=COTS.Q50))^2,3),
         COTS.Q25 = round(predict(MTCalib.gaminv, newdata=data.frame(DENS=COTS.Q25))^2,3),
         COTS.Q75 = round(predict(MTCalib.gaminv, newdata=data.frame(DENS=COTS.Q75))^2,3))


ForGraph$YearDec = as.numeric(ifelse(ForGraph$Season=="summer", ForGraph$Year, c(ForGraph$Year, ".5")))
# colnames(ForDashboard)[7:length(colnames(ForDashboard))] = paste0(rep(1:NREPS, each=8),"_", colnames(ResultsDash[7:14]))



#### Validation ----
manta.reefs = data.manta %>% 
  dplyr::filter(REEF_NAME.y %in% unique(data.grid$REEF_NAME) & REPORT_YEAR > 1995) %>% 
  mutate(REEF_NAME = as.factor(REEF_NAME.y)) %>%
  group_by(A_SECTOR, REEF_NAME) %>% 
  dplyr::summarise(n=length(REPORT_YEAR)) %>% filter(n>8) %>% top_n(3, n) %>%
  pull(REEF_NAME) %>% as.character()

manta.reefs.all = data.manta %>% 
  dplyr::filter(REEF_NAME.y %in% unique(data.grid$REEF_NAME) & REPORT_YEAR > 1995) %>% 
  mutate(REEF_NAME = as.factor(REEF_NAME.y)) %>%
  group_by(A_SECTOR, REEF_NAME) %>% 
  dplyr::summarise(n=length(REPORT_YEAR)) %>% filter(n>3) %>%
  pull(REEF_NAME) %>% as.character()


res.plot = ForGraph %>% filter(REEF_NAME %in% data.grid$REEF_NAME) %>% 
  dplyr::group_by(REEF_NAME,REP,Year,SECTOR, CROSS_SHELF) %>%
  dplyr::summarise(COTS.mn = mean(COTS.mn), # add in detection
                   COTS.25 = mean(COTS.Q25),
                   COTS.75 = mean(COTS.Q75),
                   CC.mn = mean(CC.mn),
                   CC.25 = mean(CC.Q25),
                   CC.75 = mean(CC.Q75))
res.plot.manta = res.plot[which(res.plot$REEF_NAME %in% manta.reefs),]
res.plot.manta.all = res.plot[which(res.plot$REEF_NAME %in% manta.reefs.all),]

data.manta.valid = data.manta %>% filter(REEF_NAME.y %in% manta.reefs & REPORT_YEAR > 1995) %>% 
  mutate(REEF_NAME = REEF_NAME.y,
         Year = REPORT_YEAR - 1) %>% 
  right_join(res.plot.manta, by=c("REEF_NAME", "Year"))

data.manta.valid.all = data.manta %>% filter(REEF_NAME.y %in% manta.reefs.all & REPORT_YEAR > 1995) %>% 
  mutate(REEF_NAME = REEF_NAME.y,
         Year = REPORT_YEAR - 1) %>% 
  right_join(res.plot.manta.all, by=c("REEF_NAME", "Year"))


SECTORS = unique(as.character(data.manta.valid$SECTOR))
MPE = data.frame(SECTOR= rep(SECTORS, each=NREPS), 
                 REP = 1:NREPS,
                 MPE.CC = vector(mode = "numeric", NREPS*length(SECTORS)),
                 MPE.CC.c = vector(mode = "numeric", NREPS*length(SECTORS)),
                 MPE.CC.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 MPE.CC.ALL.c = vector(mode = "numeric", NREPS*length(SECTORS)),
                 MPE.COTS = vector(mode = "numeric", NREPS*length(SECTORS)),
                 MPE.COTS.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 ACC.Pres = vector(mode = "numeric", NREPS*length(SECTORS)),
                 ACC.Pres.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 KAP.Pres = vector(mode = "numeric", NREPS*length(SECTORS)),
                 KAP.Pres.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 P.Pres = vector(mode = "numeric", NREPS*length(SECTORS)),
                 P.Pres.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 KAP.Out = vector(mode = "numeric", NREPS*length(SECTORS)),
                 KAP.Out.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 P.Out = vector(mode = "numeric", NREPS*length(SECTORS)),
                 P.Out.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 ACC.Out = vector(mode = "numeric", NREPS*length(SECTORS)),
                 ACC.Out.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 pR2.COTS = vector(mode = "numeric", NREPS*length(SECTORS)),
                 pR2.COTS.ALL = vector(mode = "numeric", NREPS*length(SECTORS)),
                 pR2sqrt.COTS = vector(mode = "numeric", NREPS*length(SECTORS)),
                 pR2sqrt.COTS.ALL = vector(mode = "numeric", NREPS*length(SECTORS)))



#Compute validation table at sector level
for (i in 1:NREPS) {
  # browser()
  print(i)
  manta.i = dplyr::filter(data.manta.valid, REP==i)
  lm = lm(CC.mn~MEAN_LIVE.corr, data=manta.i)
  summary(lm)$r.squared  ## R2 = 0.637
  MPE[MPE$REP==i,"MPE.CC.ALL"] = mean(abs(summary(lm)$residuals)) ## Mean prediction error = 7.7 %
  MPE[MPE$REP==i,"MPE.CC.ALL.c"] = lm$coefficients[2]
  lm = lm(COTS.mn~MEAN_COTS, data=manta.i)
  summary(lm)$r.squared  ## R2 = 0.637
  MPE[MPE$REP==i,"MPE.COTS.ALL"] = mean(abs(summary(lm)$residuals))
  # make confusion matrix --> get KAPPA and Accuracy
  obs.bin = factor(ifelse(manta.i$MEAN_COTS>0.01,1,0),levels = c(0,1)) ; pred.bin = factor(ifelse(manta.i$COTS.mn > 0.01, 1,0),levels = c(0,1))
  obs.binOUT = factor(ifelse(manta.i$MEAN_COTS >0.22,1,0), levels = c(0,1)) ; pred.binOUT = factor(ifelse(manta.i$COTS.mn > 0.22, 1,0),levels = c(0,1))
  pres.conf = caret::confusionMatrix(pred.bin, obs.bin)
  out.conf = caret::confusionMatrix(pred.binOUT, obs.binOUT) #### Confusion Matrix
  MPE[MPE$REP==i,c("ACC.Pres.ALL", "KAP.Pres.ALL", "P.Pres.ALL")] = round(rep(pres.conf$overall[c(1:2, 6)], each=length(SECTORS)),3)
  MPE[MPE$REP==i,c("ACC.Out.ALL", "KAP.Out.ALL", "P.Out.ALL")] = round(rep(out.conf$overall[c(1:2, 6)], each=length(SECTORS)),3)
  df = manta.i %>%  
    dplyr::select(REEF_NAME, Year, REP, COTS.mn, MEAN_COTS) %>% 
    mutate(COTS.mn = predict(MTCalib.gam, newdata=data.frame(MT=COTS.mn)),
           MEAN_COTS = predict(MTCalib.gam, newdata=data.frame(MT=MEAN_COTS)),
           pred.bin = pred.bin) %>%
    mutate(COTS.mn = ifelse(COTS.mn <=0.01, 0, COTS.mn),
           MEAN_COTS = round(ifelse(MEAN_COTS <0, 0, MEAN_COTS),0))
  # run Zeron-infl neg bin on results using pscl to get pseudo R2 --> convert to integer value
  try({
    m1 = pscl::zeroinfl(MEAN_COTS ~ COTS.mn | pred.bin,
                        data = df, dist = "negbin", EM = TRUE)
    m2 = pscl::zeroinfl(round(sqrt(MEAN_COTS),0) ~ round(sqrt(COTS.mn),0) | pred.bin,
                        data = df, dist = "negbin", EM = TRUE)
    MPE[MPE$REP==i,"pR2.COTS.ALL"] = round(r2_zeroinflated(m1)$R2,3)
    MPE[MPE$REP==i,"pR2sqrt.COTS.ALL"] = round(r2_zeroinflated(m2)$R2,3)
  }, silent=T)
  
  
  for (j in 1:length(SECTORS)) {
    # browser()
    print(j)
    manta.i = dplyr::filter(data.manta.valid, REP==i & SECTOR == SECTORS[j])
    lm = lm(CC.mn~MEAN_LIVE.corr, data=manta.i)
    summary(lm)$r.squared  ## R2 = 0.637
    MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,"MPE.CC"] = mean(abs(summary(lm)$residuals)) ## Mean prediction error = 7.7 %
    MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,"MPE.CC.c"] = lm$coefficients[2]
    lm = lm(COTS.mn~MEAN_COTS, data=manta.i)
    summary(lm)$r.squared  ## R2 = 0.637
    MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,"MPE.COTS"] = mean(abs(summary(lm)$residuals)) ##
    # make confusion matrix --> get KAPPA and Accuracy
    obs.bin = factor(ifelse(manta.i$MEAN_COTS>0.01,1,0),levels = c(0,1)) ; pred.bin = factor(ifelse(manta.i$COTS.mn > 0.01, 1,0),levels = c(0,1))
    obs.binOUT = factor(ifelse(manta.i$MEAN_COTS >0.22,1,0), levels = c(0,1)) ; pred.binOUT = factor(ifelse(manta.i$COTS.mn > 0.22, 1,0),levels = c(0,1))
    pres.conf = caret::confusionMatrix(pred.bin, obs.bin)
    out.conf = caret::confusionMatrix(pred.binOUT, obs.binOUT)#### Confusion Matrix
    MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,c("ACC.Pres", "KAP.Pres", "P.Pres")] = pres.conf$overall[c(1:2, 6)]
    MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,c("ACC.Out", "KAP.Out", "P.Out")] = out.conf$overall[c(1:2, 6)]
    df = manta.i %>%  
      dplyr::select(REEF_NAME, Year, REP, COTS.mn, MEAN_COTS) %>% 
      mutate(COTS.mn = predict(MTCalib.gam, newdata=data.frame(MT=COTS.mn)),
             MEAN_COTS = predict(MTCalib.gam, newdata=data.frame(MT=MEAN_COTS)),
             pred.bin = pred.bin) %>%
      mutate(COTS.mn = ifelse(COTS.mn <=0.01, 0, COTS.mn),
             MEAN_COTS = round(ifelse(MEAN_COTS <0, 0, MEAN_COTS),0))
    
    # run Zeron-infl neg bin on results using pscl to get pseudo R2 --> convert to integer value
    try({
      m1 = pscl::zeroinfl(MEAN_COTS ~ COTS.mn | pred.bin,
                          data = df, dist = "negbin", EM = TRUE)
      m2 = pscl::zeroinfl(round(sqrt(MEAN_COTS),0) ~ round(sqrt(COTS.mn),0) | pred.bin,
                          data = df, dist = "negbin", EM = TRUE)
      MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,"pR2.COTS"] = round(r2_zeroinflated(m1)$R2,3)
      MPE[MPE$SECTOR==SECTORS[j] & MPE$REP==i,"pR2sqrt.COTS"] = round(r2_zeroinflated(m2)$R2,3)
    }, silent = T)
    
  }
}  

range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
MPE.overall = MPE %>% group_by(REP) %>% 
  dplyr::select(contains("ALL")) %>%
  # dplyr::select(MPE.CC.ALL, MPE.COTS.ALL, ACC.Pres.ALL, KAP.Pres.ALL, KAP.Out.ALL, 
  #               ACC.Out.ALL, pR2sqrt.COTS.ALL,P.Pres.ALL, P.Out.ALL) %>%
  summarise_all(mean) %>%
  mutate(Overall = (range01(ACC.Pres.ALL) + range01(KAP.Pres.ALL) +range01(ACC.Out.ALL) +range01(KAP.Out.ALL) + (1-range01(MPE.CC.ALL)))/5)


save(res.plot, MPE, MPE.overall,  data.manta.valid, data.manta.valid.all, masterDF,
     file = paste0("Results/Archive/Validation_", Sys.Date(),"_", nruns, ".RData"))
