# COTS Chlorophyll Functions


# install.packages("DHARMa")
# install.packages("effects")
# install.packages("sjPlot")
library(ggplot2)
# library(DHARMa)
# library(sjPlot)
library(lme4)
library(car)
library(effects)
library(MASS)

inv.logit = function(x) {
  exp(x)/(1+exp(x))
}



fab = read.csv("Data/ChlorophyllModel/Fabricius.csv")
fab$weights = 100
fab=fab[which(fab$keep),]
car::logit(fab$surv/100)
# chl.lm$coefficients = c(1,2.2)
chl.lm = lm(car::logit(surv/100, adjust = 0.001)~log(chl,base = 2), data=fab)

xlevels = as.list(data.frame(chl = seq(0, 6,len = 100)))
newdata = as.data.frame(Effect("chl", chl.lm, xlevels = xlevels))
newdata$fit = inv.logit(newdata$fit)
newdata$lower = inv.logit(newdata$lower)
newdata$upper = inv.logit(newdata$upper)

ggplot(data = newdata, aes(y = fit, x = chl)) + 
  geom_point(data = fab,aes(y = surv/100)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper),fill = "blue", alpha = 0.3) + 
  geom_line() + scale_y_continuous("P(Survival") +
  theme_classic() + scale_x_continuous(trans="log2")

coefs <- MASS::mvrnorm(n = 10000, mu = coefficients(chl.lm), Sigma = vcov(chl.lm))
plot(coefs[,1]~coefs[,2])
mean(coefs[,1]); mean(coefs[,2])

coefs.quants = rbind(quantile(coefs[,2], probs = c(0.005,0.01, 0.015)),
                     quantile(coefs[,2], probs = c(0.045,0.05, 0.055)),
                     quantile(coefs[,2], probs = c(0.245,0.25, 0.255)),
                     quantile(coefs[,2], probs = c(0.495,0.5, 0.505)),
                     quantile(coefs[,2], probs = c(0.745,0.75, 0.755)),
                     quantile(coefs[,2], probs = c(0.945,0.95, 0.955)),
                     quantile(coefs[,2], probs = c(0.985,0.99, 0.995)))

quants = data.frame(Int=NA, chl=NA, Quant=NA)[-1,]
for (i in 1:nrow(coefs.quants)){
  quants[i,2] = coefs.quants[i,2]
  quants[i,1] = median(coefs[which(coefs[,2] >= coefs.quants[i,1] & 
                                     coefs[,2] <= coefs.quants[i,3]),1])  
}
#Quantile Coefficient

quants$Quant = c("Q0.01","Q0.05", "Q0.25", "Q0.5", "Q0.75", "Q0.95", "Q0.99")
quants = rbind(quants, c(coef(chl.lm), "model"))

# inv.logit(predict(chl.lm, newdata = data.frame(chl=10),type = "response"))

# save(chl.lm, coefs.quants, file="Data/ChlorophyllModel/ChlModelParams.RData")

modelfits = data.frame(chl=NA, fit=NA, lower=NA, upper=NA, model=NA)[-1,]
# Check to see what it looks like
for (i in 1:nrow(quants)) {
  chl.lm$coefficients = as.numeric(quants[i,1:2])
  newdata = as.data.frame(Effect("chl", chl.lm, xlevels = xlevels))
  newdata$fit = inv.logit(newdata$fit)
  newdata$lower = inv.logit(newdata$lower)
  newdata$upper = inv.logit(newdata$upper)
  newdata$model = quants[i,3]
  modelfits = rbind(modelfits, newdata)
}

modelfits$model = factor(modelfits$model)
#newdata$chl = log(newdata$chl, base = 2)
ggplot(data = modelfits, aes(y = fit, x = chl)) + 
  geom_line(size=1) +
  geom_ribbon(aes(ymin = lower, ymax = upper,  fill=model), alpha = 0.3) +
  scale_y_continuous("P(Survival)") +
  geom_point(data = fab,aes(y = surv/100, x=chl)) + 
  theme_classic() + 
  scale_x_continuous(trans="log2", 
                     labels = scales::number_format(accuracy = 0.01)) +
  facet_wrap(~model, ncol = 2)

save(quants, chl.lm, file =  "Data/ChlorophyllModel/ChlModels.RData")
# Sven Model

data=read.csv("Data/ChlorophyllModel/COTS contour.csv", header=T)
attach(data)
library(lattice)

x<-data$Temp
y<-data$Algae
z<-data$First_LB

ggplot(data[-c(2,5),], aes(x=Temp, y=First_LB, color=as.factor(Algae))) + 
  geom_point(position=position_dodge(width = 0.5))
ggplot(data[-c(2,5),], aes(x=Algae, y=First_LB, color=as.factor(Temp))) + 
  geom_point(position=position_dodge(width = 500))
## First late brachs

First_LB.loess = loess(z ~ x*y, degree = 1, span = 1,control=loess.control(surface = "direct"))
fit.glm = lm(z~x*y, data = data)
summary(fit.glm)
First_LB.fit = expand.grid(list(x = seq(26, 31, 0.01), y = seq(1100, 9800, 10)))

z = predict(First_LB.loess, newdata = First_LB.fit)
z2 = predict(fit.glm, newdata = First_LB.fit)
First_LB.fit$LB = as.numeric(z)
First_LB.fit$LB2 = as.numeric(z2)


#par(mfrow=c(2,1), mar=c(4,4,1.5,0), oma=c(4,5,1,5), las=1)
col.l <- colorRampPalette(c('green', 'yellow', 'orange','red'))
levelplot(LB ~ x*y, data = First_LB.fit,
          xlab = "Temperature", ylab = "Algae (ml-1)",
          main = "Time to first Late Brachiolaria (d)", cuts = 10, pretty=T, 
          colourkey="bottom",
          col.regions = col.l
)

# ## Just to get individual values out
# First_LB2.fit = expand.grid(list(x = seq(26, 30, 1), y = seq(1100, 9800, 100)))
# z2 = predict(First_LB.loess, newdata = First_LB2.fit)
# First_LB2.fit$LB = as.numeric(z2)
# First_LB2.fit$Survival<-exp(-0.16*(First_LB2.fit$LB+2))

## test survival prabability: assume: fist brach + 2d can settle

## A) M= 0.05 (asteroid average Rumrill) / B = 0.16 (change in script)

First_LB.fit$Survival<-exp(-0.16*(First_LB.fit$LB+2))
First_LB.fit$Survival2<-exp(-0.16*(First_LB.fit$LB2+2))
levelplot(Survival2 ~ x*y, data = First_LB.fit,
          xlab = "Temperature", ylab = "Algae (ml-1)",
          main = "Survival probability (M = 0.16)", cuts = 15, pretty=T, 
          colourkey="bottom",
          col.regions = col.l
)
