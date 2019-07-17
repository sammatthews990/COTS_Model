### COTS Manta Tow Times Series Anlaysis ---

reefs.manta = unique(data.manta$REEF_NAME.y)

data.manta$DATE = as.Date(substring(data.manta$SAMPLE_DATE,1,
                            nchar(as.character(data.manta$SAMPLE_DATE))-5), "%d/%m/%y")
data.manta2 = data.manta %>%
  arrange(DATE) %>%
  group_by(REEF_NAME.y) %>%
  mutate(Ratio = MEAN_LIVE.corr/MEAN_COTS,
         Max.COTS = max(MEAN_COTS, na.rm = T),
         Loss.COTS = lead(MEAN_COTS)-MEAN_COTS,
         Loss.COTS.Prop = (lead(MEAN_COTS)-MEAN_COTS)/MEAN_COTS, 
         Loss.COTS.Prop2 = ifelse(Loss.COTS.Prop >1,1, Loss.COTS.Prop),
         Loss.COTS.Prop3 = ifelse(Loss.COTS.Prop >0.5,1, ifelse(Loss.COTS.Prop < (-0.5), 0, NA)),
         Loss.CC = lead(MEAN_LIVE.corr)-MEAN_LIVE.corr,
         Loss.CC.Prop = (lead(MEAN_LIVE.corr)-MEAN_LIVE.corr)/MEAN_LIVE.corr,
         NextSurvey = lead(DATE) - DATE,
         Ratio.sqrt = sqrt(Ratio))
data.manta2$Loss.COTS.Binom = BBmisc::normalize(data.manta2$Loss.COTS.Prop2, method="range")
mantasub = dplyr::filter(data.manta2, REEF_NAME.y %in% reefs.manta[i]) %>%
  dplyr::select(REEF_ID, REEF_NAME.y, A_SECTOR, SHELF, REEF_LAT, REEF_LONG, DATE, 
                REPORT_YEAR, MEAN_LIVE.corr, MEAN_COTS, SE_COTS, Ratio, Loss.COTS:NextSurvey)

library(tidyverse)

# OUtbreak Must Have been established (i.e > 1 COTS/Manta Tow)

ggplot(data.manta2 %>% filter(NextSurvey < 480 & MEAN_COTS >5 & MEAN_COTS <20), aes(x=MEAN_COTS, y=Loss.COTS)) + 
  geom_point() + 
  geom_smooth(method = "lm", fill="green", colour="black") +
  labs(x="COTS/Manta Tow", y="Change in COTS Population", title = "COTS Population Crashes") +
  theme_classic()
ggplot(data.manta2 %>% filter(NextSurvey < 480 & MEAN_COTS >0.22 & MEAN_COTS < 5), aes(x=MEAN_COTS, y=Loss.COTS.Prop2)) + 
  geom_point() + 
  geom_smooth(method = "lm", fill="green", colour="black") +
  labs(x="COTS/Manta Tow", y="Change in COTS Population", title = "COTS Population Crashes") +
  theme_classic() 
ggplot(data.manta2 %>% filter(NextSurvey < 480 & Ratio < 5), aes(x=MEAN_COTS, y=Loss.COTS.Prop)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=4), fill="green", colour="black") +
  labs(x="COTS/Manta Tow", y="Change in COTS Population", title = "COTS Population Crashes") +
  theme_classic()
ggplot(data.manta2 %>% filter(NextSurvey < 480 & MEAN_COTS > 1), aes(x=Ratio, y=Loss.CC)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=4), fill="green", colour="black") +
  labs(x="COTS/Manta Tow", y="Change in COTS Population", title = "COTS Population Crashes") +
  theme_classic()

library(interplot)

## Filter only by reefs that have had an outbreak

dat = data.manta2 %>% filter(NextSurvey<420 & MEAN_COTS > 5)

# dat = data.manta2 %>% filter(NextSurvey<420 & Loss.COTS.Prop3 == 1)
# summary(dat$Ratio)
model1 = lm(Loss.COTS~MEAN_COTS, data = dat)
model2 = lm(Loss.COTS~MEAN_COTS*MEAN_LIVE.corr, data = dat)
model3 = glm(Loss.COTS.Prop3~Ratio, data = dat, family = "binomial")
interplot(m = model2, var1 = "MEAN_COTS", var2 = "MEAN_LIVE.corr") +
  xlab("SQRT-Ratio") + ylab("Estimated Coefficient for MEAN_COTS")
summary(model1)
plot(model2)

ggplot(dat %>% filter(Ratio < 100), aes(x=Ratio, y=Loss.COTS.Prop3)) + 
  geom_point() + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), fill="green", colour="black") +
  labs(x="Coral-Cover-COTS/Manta Tow", y="Loss of COTS", title = "COTS Population Crashes") +
  theme_classic() 

ggplot(dat %>% filter(Ratio < 100), aes(x=Ratio, y=Loss.COTS.Prop3)) + 
  geom_point(aes(colour=A_SECTOR)) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), aes(fill=A_SECTOR), colour="black") +
  labs(x="Coral-Cover-COTS/Manta Tow", y="Loss of COTS", title = "COTS Population Crashes") +
  theme_classic() 

ggplot(mantasub, aes(x=DATE)) + 
  geom_ribbon(aes(ymin = MEAN_COTS*5 - SE_COTS*5, ymax = MEAN_COTS*5 + SE_COTS*5), fill="purple", alpha=0.3) +
  geom_line(aes(y=MEAN_LIVE.corr, colour = "Mean Live Coral")) + 
  geom_line(aes(y=MEAN_COTS*5, colour = "Mean COTS")) + 
  scale_y_continuous(sec.axis = sec_axis(~./5, name = "Mean COTS/Manta Tow")) +
  scale_colour_manual(values = c("purple", "coral")) +
  labs(y="Mean Live Coral") +
  ggtitle(reefs.manta[i]) +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.17, 0.12),
        legend.title = element_blank()) 

hist(data.manta2$Ratio[which(data.manta2$Ratio <25)])
