######################!
# COTSFertilisation By Distance
######################!

###################!
# Fertilisation by distance function ----


Fert.data <- data.frame(Dist = c(0,2,4,8,16,32,64,100), PercFert = c(90,86.5, 71.8,71.9,41.5,26.8,20.5,5.8))

m1 <- as.formula(PercFert ~ p * exp(k * Dist)) #standard 
m2 <- as.formula(PercFert ~ p * exp(k * Dist) + q) #standard 
m3 <- as.formula(PercFert ~ p * exp(k * Dist) + (92-p)) #fixed intercept of 92%

em <-function(x,p,k,q) {(p*exp(k*x)) + q}
em.fixed <-function(x,p,k,f) {(p*exp(k*x)) + (f-p)}
em2<-function(x,p,k) {(p*exp(k*x))}


nls1 <- nls(m1,start=list(p=80,k=-0.05), data = Fert.data)
nls2 <- nls(m2,start=list(p=90,k=-0.05, q=10), data = Fert.data, control = list(maxiter=500))
nls3 <- nls(m3,start=list(p=80,k=-0.05), data = Fert.data)


###################!
###################!
# Retrieve best fit model parameters ----
BestFitPars <- function(nls.object){
  confints <- confint(nls.object)
  bestpars <- nls.object$m$getPars()
  upperpars <- confints[,2] 
  lowerpars <- confints[,1]
  pars <- data.frame(bestpars,upperpars,lowerpars)
  return(pars)
}

BestFitPars(nls1)
BestFitPars(nls2)
BestFitPars(nls3)


###################!
###################!
# Plot Fertilisation Function ----
Data <- Fert.data
nls.object <- nls1

nlsCIplot <- function(nls.object, Data) {
  
  (confints <- confint(nls.object))
  (bestpars <- nls.object$m$getPars())
  upperpars <- confints[,2] 
  lowerpars <- confints[,1]
  (pars <- data.frame(bestpars,upperpars,lowerpars))
  
  fit.data <- data.frame(x=seq(0,200,len=100), 
                         best=NA, upper=NA, lower=NA)
  fit.data$best <- em(fit.data$x, bestpars[1], bestpars[2], bestpars[3])
  
  fit.data$upper <- em(fit.data$x, upperpars[1], upperpars[2], upperpars[3])
  fit.data$lower <- em(fit.data$x, lowerpars[1], lowerpars[2], lowerpars[3])
  ggplot(fit.data, aes(y=best, x=x)) +
    geom_line() + theme_classic() +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) +
    geom_point(data=Data, aes(x=Dist,y=PercFert)) +
    xlab("Distance") +
    ylab("Percentage Eggs Fertilised") +
    ggtitle("Fertilisation Function")
}

nlsCIplot2 <- function(nls.object, Data) {
  
  (confints <- confint(nls.object))
  (bestpars <- nls.object$m$getPars())
  upperpars <- confints[,2] 
  lowerpars <- confints[,1]
  (pars <- data.frame(bestpars,upperpars,lowerpars))
  
  fit.data <- data.frame(x=seq(0,200,len=100), 
                         best=NA, upper=NA, lower=NA)
  fit.data$best <- em2(fit.data$x, bestpars[1], bestpars[2])
  
  fit.data$upper <- em2(fit.data$x, upperpars[1], upperpars[2])
  fit.data$lower <- em2(fit.data$x, lowerpars[1], lowerpars[2])
  ggplot(fit.data, aes(y=best, x=x)) +
    geom_line() + theme_classic() +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) +
    geom_point(data=Data, aes(x=Dist,y=PercFert)) +
    xlab("Distance") +
    ylab("Percentage Eggs Fertilised") +
    ggtitle("Fertilisation Function")
}

nlsCIplot(nls2, Fert.data)
nlsCIplot2(nls1, Fert.data)

(bestpars <- BestFitPars(nls1))

###################!
###################!
# Fertilisation by Density ---- 

# This section will determine the Von Bertalanffy Growth Parameters for 10 different sex ratios and plot the results

geomSeries <- function(base, max) {
  base^(0:floor(log(max, base)))
}


Densities <- geomSeries(2, 10000)
Densities <- c(Densities, 1536, 2560, 3072, 3584, 5120, 6144,7168,9216)
Densities <- sort(Densities)
SexRatio <- c('M' = 0.2, 'F' = 0.8)
SexRatios <- list(c('M' = 0.1, 'F' = 0.9), c('M' = 0.2, 'F' = 0.8), c('M' = 0.3, 'F' = 0.7),
                  c('M' = 0.4, 'F' = 0.6), c('M' = 0.5, 'F' = 0.5), c('M' = 0.6, 'F' = 0.4),
                  c('M' = 0.7, 'F' = 0.3), c('M' = 0.8, 'F' = 0.2), c('M' = 0.9, 'F' = 0.1))

expParams <- bestpars[,1]

sexes <- seq(0.1,0.9, 0.1)
sexes <- paste("M", sexes, sep="")  
# I will have to convert fertilisation based on reefpercent

FertVsDensity <- function (SexRatio, Densities, expParams) {
  # Place all adults amongst landscape
  props <- data.frame(Dens=Densities, Female = NA, Male = NA)
  fertbydens <- list()
  
  for (i in 1:length(Densities)){
    coordsx <- runif(Densities[i], min = 0, max = 1000)
    coordsy <- runif(Densities[i], min = 0, max = 1000)
    coords <- data.frame(x = coordsx, y = coordsy)
    # randomly assign each coordinate pair as male or female
    coords$sex <- sample(c('M', 'F'), length(coordsx), replace = T, prob = SexRatio)
    #determine number of males and females
    props[i,2:3] <- prop.table(table(coords$sex))[1:2]*Densities[i]
    
    if(sum(is.na(props[i,2:3]))==1) {
      fertbydens[[i]] <- rep(NA, Densities[i])
      next
    } else {
      
      # now for each female, we calculate the distance and the probability that her eggs were NOT
      # fertilised by that male
      
      Females <- as.matrix(coords[coords[,'sex']=='F',][,1:2])
      Males <- as.matrix(coords[coords[,'sex']=='M',][,1:2])
      # Each row has the distance between each female and male
      DistMat <- SpatialTools::dist2(Females, Males)
      # convert distance into probability of not fertalising
      em2 <- function(x,p,k) {(1- (p*exp(k*x))/100)}
      DistMat.NotFert <- matrix(sapply(DistMat, FUN = em2, expParams[1], expParams[2]), 
                                nrow = length(Females[,1]), ncol = length(Males[,1]))
      #Multiply across each row to find out the probability that egss from that female were not fertilised
      # by any male
      NotFert <- apply(DistMat.NotFert, MARGIN = 1, FUN = prod)
      fertbydens[[i]] <- 1-NotFert
    }
  }
  # name the list by the desnity of the population
  names(fertbydens) <- Densities
  #convert into a dataframe for plotting
  dens <- rep(Densities, each = 1, times = props[,'Female'])
  fert <- unlist(fertbydens)
  fert.df <- data.frame(dens,fert)
  return(list(Props = props, FertbyDens.df = fert.df, FertbyDens = fertbydens))
}

f1 <- FertVsDensity(SexRatio, Densities, expParams) 

l1 <- lapply(SexRatios, FUN = FertVsDensity, expParams=expParams, Densities=Densities)

names(l1) <- sexes 



# Von Bertalanffy Growth
#
# This function takes the data created by the FertVsDensity function to define the parameters of the 
# Von Bertalanffy Growth for each of our potential densities
#

FvD <- l1

SSQ <- function(theta, dens, fert) {
  Linf <- theta[1]
  K <- theta[2]
  t0 <- theta[3]
  epsilon <- rep(0, length(dens))
  lpred <- rep(0, length(dens))
  for (i in 1:length(dens)) {
    lpred[i] <- Linf * (1 - exp(-K * (dens[i] - t0)))
    epsilon[i] <- (fert[i] - lpred[i])^2
  }
  ssq <- sum(na.omit(epsilon))
  return(ssq)
}

VonBertalannfyGrowth <- function(FvD){
  Model <- list()
  Parameters <- data.frame(SexRatio=names(FvD), Linf = NA, K = NA, t0 = NA)
  for (i in 1:length(FvD)) {
    fert <- FvD[[i]][[2]][,2]
    dens <- FvD[[i]][[2]][,1]
    
    theta <- c(1, 0.001, 0.1)
    
    out <- optim(theta, fn = SSQ, method = "BFGS", 
                 dens = na.omit(dens), fert=fert,  hessian = TRUE)
    #out$V <- solve(out$hessian)  #solve the hessian
    #out$S <- sqrt(diag(out$V))  #Standard Error
    #out$R <- out$V/(out$S %o% out$S)  #Correlation
    Model[[i]] <- out
    Parameters[i,2:4] <- out$par
  }
  return(list(Models=Model, Parameters=Parameters))
}


VBG.Models <- VonBertalannfyGrowth(FvD)




###################!
###################!
# VBGPlot ----
VBGPlot <- function(SR, FvD, Parameters) {
  #Sex Ratio is integer from 1 to 9 relating to proportion of Males --> i.e 1 = 0.1M, 0.9F Sex Ratio
  data <- FvD[[SR]][[2]]
  
  fit <- Parameters[SR,2] * (1 - exp(-Parameters[SR, 3] * (data$dens - Parameters[SR,4])))
  fit.data <- data.frame(dens=data$dens, fit=fit)
  
  ggplot(data=data, aes(x=dens, y=fert)) +
    geom_point() +
    geom_smooth() +
    geom_line(data=fit.data, aes(x=dens, y=fit), col="green", size=1) +
    labs(x=expression(CoTS~Density~(km^{-2})),
         y="% of Eggs Fertilised") +
    ggtitle(paste(SR,"M:",10-SR, "F", sep="")) +
    theme_classic() 
  
}

VBGPlots = list()
for (i in 1:9){ 
  VBGPlots[[i]] = VBGPlot(i, FvD, Parameters = VBG.Models[["Parameters"]])
}

library(ggpubr)
ggarrange(VBGPlots[[1]], VBGPlots[[2]],VBGPlots[[3]], VBGPlots[[4]],
          VBGPlots[[5]], VBGPlots[[6]], VBGPlots[[7]], VBGPlots[[8]],
          VBGPlots[[9]], ncol=3, nrow=3)

VBG = function (Linf, K, dens, t0) {Linf * (1 - exp(-K * (dens - t0)))}

plot(1:10000,VBG(0.7,2e-4, 1:10000, 0))
points(1:10000,VBG(0.7,6e-4, 1:10000, 0))
