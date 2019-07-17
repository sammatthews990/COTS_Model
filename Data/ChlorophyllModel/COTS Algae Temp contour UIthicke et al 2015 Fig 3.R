## set the working directory


setwd("C:/Users/suthicke/Documents/R/COTS experiments 2013")


##### clear workspace, read in file define variables
## after first checks 4 extreme outliers (all in Blue 30 removed)


rm(list=ls())

data<-read.csv("COTS contour.csv", header=T)
attach(data)


##library(ggplot2)

library(lattice)


x<-data$Temp
y<-data$Algae
z<-data$First_LB


## First late brachs

First_LB.loess = loess(z ~ x*y, degree = 1, span = 1)

First_LB.fit = expand.grid(list(x = seq(28, 30, 0.01), y = seq(1100, 9800, 10)))
z = predict(First_LB.loess, newdata = First_LB.fit)
First_LB.fit$LB = as.numeric(z)


#par(mfrow=c(2,1), mar=c(4,4,1.5,0), oma=c(4,5,1,5), las=1)
col.l <- colorRampPalette(c('green', 'yellow', 'orange','red'))
levelplot(LB ~ x*y, data = First_LB.fit,
          xlab = "Temperature", ylab = "Algae (ml-1)",
          main = "Time to first Late Brachiolaria (d)", cuts = 10, pretty=T, 
          colourkey="bottom",
          col.regions = col.l
)

## Just to get individual values out
First_LB2.fit = expand.grid(list(x = seq(28, 30, 1), y = seq(1100, 9800, 100)))
z2 = predict(First_LB.loess, newdata = First_LB2.fit)
First_LB2.fit$LB = as.numeric(z2)
First_LB2.fit$Survival<-exp(-0.16*(First_LB2.fit$LB+2))

## test survival prabability: assume: fist brach + 2d can settle

## A) M= 0.05 (asteroid average Rumrill) / B = 0.16 (change in script)

First_LB.fit$Survival<-exp(-0.16*(First_LB.fit$LB+2))

levelplot(Survival ~ x*y, data = First_LB.fit,
          xlab = "Temperature", ylab = "Algae (ml-1)",
          main = "Survival probability (M = 0.16)", cuts = 15, pretty=T, 
          colourkey="bottom",
          col.regions = col.l
)
