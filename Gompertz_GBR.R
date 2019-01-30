
##############################################################
# DYNAMIC MODEL OF CORAL COVER - CM 05/06/18                 #
##############################################################

# Load libraries ------------
rm(list = ls())

library(maptools)
library(PBSmapping)
library(RColorBrewer)
library(sampling)
library(psych)
library(mgcv)
library(lhs)
library(gbm)
library(MASS)
library(nlme)
library(lme4)
library(extrafont)
library(FactoMineR)
library(sm)
library(png)

# Functions ------------

vioplot <- function(x,...,range=1.5,h=NULL,ylim=NULL,names=NULL, horizontal=FALSE,
                    col="magenta", border="black", lty=1, lwd=1, rectCol="black", colMed="white", pchMed=19, at, add=FALSE, wex=1,
                    drawRect=TRUE)
{
  # process multiple datas
  datas <- list(x,...)
  n <- length(datas)
  
  if(missing(at)) at <- 1:n
  
  # pass 1
  #
  # - calculate base range
  # - estimate density
  #
  
  # setup parameters for density estimation
  upper  <- vector(mode="numeric",length=n)
  lower  <- vector(mode="numeric",length=n)
  q1     <- vector(mode="numeric",length=n)
  q3     <- vector(mode="numeric",length=n)
  med    <- vector(mode="numeric",length=n)
  base   <- vector(mode="list",length=n)
  height <- vector(mode="list",length=n)
  baserange <- c(Inf,-Inf)
  
  # global args for sm.density function-call
  args <- list(display="none")
  
  if (!(is.null(h)))
    args <- c(args, h=h)
  
  for(i in 1:n) {
    data<-datas[[i]]
    
    # calculate plot parameters
    #   1- and 3-quantile, median, IQR, upper- and lower-adjacent
    data.min <- min(data)
    data.max <- max(data)
    q1[i]<-quantile(data,0.25, na.rm=T)
    q3[i]<-quantile(data,0.75, na.rm=T)
    med[i]<-median(data)
    iqd <- q3[i]-q1[i]
    upper[i] <- min( q3[i] + range*iqd, data.max )
    lower[i] <- max( q1[i] - range*iqd, data.min )
    
    #   strategy:
    #       xmin = min(lower, data.min))
    #       ymax = max(upper, data.max))
    #
    
    est.xlim <- c( min(lower[i], data.min), max(upper[i], data.max) )
    
    # estimate density curve
    smout <- do.call("sm.density", c( list(data, xlim=est.xlim), args ) )
    
    # calculate stretch factor
    #
    #  the plots density heights is defined in range 0.0 ... 0.5
    #  we scale maximum estimated point to 0.4 per data
    #
    hscale <- 0.4/max(smout$estimate) * wex
    
    # add density curve x,y pair to lists
    base[[i]]   <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    
    # calculate min,max base ranges
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1],t[1])
    baserange[2] <- max(baserange[2],t[2])
    
  }
  
  # pass 2
  #
  # - plot graphics
  
  # setup parameters for plot
  if(!add){
    xlim <- if(n==1)
      at + c(-.5, .5)
    else
      range(at) + min(diff(at))/2 * c(-1,1)
    
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  } else {
    label <- names
  }
  
  boxwidth <- 0.05 * wex
  
  # setup plot
  if(!add)
    plot.new()
  if(!horizontal) {
    if(!add){
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1,at = at, label=label )
    }
    
    box()
    for(i in 1:n) {
      # plot left/right density curve
      polygon( c(at[i]-height[[i]], rev(at[i]+height[[i]])),
               c(base[[i]], rev(base[[i]])),
               col = col[i %% length(col) + 1], border=border, lty=lty, lwd=lwd)
      
      if(drawRect){
        # plot IQR
        lines( at[c( i, i)], c(lower[i], upper[i]) ,lwd=lwd, lty=lty)
        
        # plot 50% KI box
        rect( at[i]-boxwidth/2, q1[i], at[i]+boxwidth/2, q3[i], col=rectCol)
        
        # plot median point
        points( at[i], med[i], pch=pchMed, col=colMed )
      }
    }
    
  }
  else {
    if(!add){
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2,at = at, label=label )
    }
    
    box()
    for(i in 1:n) {
      # plot left/right density curve
      polygon( c(base[[i]], rev(base[[i]])),
               c(at[i]-height[[i]], rev(at[i]+height[[i]])),
               col = col[i %% length(col) + 1], border=border, lty=lty, lwd=lwd)
      
      if(drawRect){
        # plot IQR
        lines( c(lower[i], upper[i]), at[c(i,i)] ,lwd=lwd, lty=lty)
        
        # plot 50% KI box
        rect( q1[i], at[i]-boxwidth/2, q3[i], at[i]+boxwidth/2,  col=rectCol)
        
        # plot median point
        points( med[i], at[i], pch=pchMed, col=colMed )
      }
    }
  }
  invisible (list( upper=upper, lower=lower, median=med, q1=q1, q3=q3))
}  

addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

colorRampPaletteAlpha <- function(colors, n=32, interpolate='linear') {
  # Create the color ramp normally
  cr <- colorRampPalette(colors, interpolate=interpolate)(n)
  # Find the alpha channel
  a <- col2rgb(colors, alpha=T)[4,]
  # Interpolate
  if (interpolate=='linear') {
    l <- approx(a, n=n)
  } else {
    l <- spline(a, n=n)
  }
  l$y[l$y > 255] <- 255 # Clamp if spline is > 255
  cr <- addalpha(cr, l$y/255.0)
  return(cr)
}

scale <- function(x) {
  x <- (x- min(na.omit(x))) / (max(na.omit(x))- min(na.omit(x)))
  x
}

scale.global <- function(x,y) {
  x <- (x- min(na.omit(y))) / (max(na.omit(y))- min(na.omit(y)))
  x
}

northarrow <- function(loc,size,bearing=0,cols,cex=1,...) {
  # checking arguments
  if(missing(loc)) stop("loc is missing")
  if(missing(size)) stop("size is missing")
  # default colors are white and black
  if(missing(cols)) cols <- rep(c("white","black"),8)
  # calculating coordinates of polygons
  radii <- rep(size/c(1,4,2,4),4)
  x <- radii[(0:15)+1]*cos((0:15)*pi/8+bearing)+loc[1]
  y <- radii[(0:15)+1]*sin((0:15)*pi/8+bearing)+loc[2]
  # drawing polygons
  for (i in 1:15) {
    x1 <- c(x[i],x[i+1],loc[1])
    y1 <- c(y[i],y[i+1],loc[2])
    polygon(x1,y1,col=cols[i])
  }
  # drawing the last polygon
  polygon(c(x[16],x[1],loc[1]),c(y[16],y[1],loc[2]),col=cols[16])
  # drawing letters
  b <- c("E","N","W","S")
  for (i in 0:3) text((size+par("cxy")[1])*cos(bearing+i*pi/2)+loc[1],
                      (size+par("cxy")[2])*sin(bearing+i*pi/2)+loc[2],b[i+1],
                      cex=cex)
}

scalebar <- function(loc,length,unit="km",division.cex=.8,...) {
  if(missing(loc)) stop("loc is missing")
  if(missing(length)) stop("length is missing")
  x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
  y <- c(0,length/(10*3:1))+loc[2]
  cols <- rep(c("black","white"),2)
  for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
  for (i in 1:5) segments(x[i],y[2],x[i],y[3])
  labels <- (x[c(1,3)]-loc[1])*100
  labels <- append(labels,paste((x[5]-loc[1])*100,unit))
  text(x[c(1,3,5)],y[4],labels=labels,adj=.5,cex=division.cex)
}

map.pred <- function(k) {
  layout(rbind(c(1, 2), c(1, 3), c(1, 4)))
  plotMap(shape, col="lightgrey", border="darkgrey", xlim=c(143,155), ylim=c(-25,-12), bg="white", cex=1.2)
  points(data.grid$LONG, data.grid$LAT, pch=19, cex=.3, col=rgb(yellow(scale.global(res.med[,k], res.med)),max=255))
  title(main=years[k], cex.main=3)
  
  for (i in c(1,8,16)) {
    x.mn <- mean(obs.manta$REEF_LONG[obs.manta$REEF_ID==res.manta.fid[i]])
    y.mn <- mean(obs.manta$REEF_LAT[obs.manta$REEF_ID==res.manta.fid[i]])
    rect(x.mn-.25, y.mn-.25, x.mn+.25, y.mn+.25, lwd=2, border=rgb(palette[[i]](.7),max=255))
  }
  
  for (i in c(1,8,16)) {
    par(mai=c(.5,.5,.5,.5))  
    plot(NA, xlim=c(1996,2015), ylim=c(10,80), xlab="Time (yr)", ylab="Coral cover (%)", xaxs="i", yaxs="i", las=1)
    
    # Add individual grid cells
    data.reef.i <- res.med[data.grid$REEF_ID==res.manta.fid[i],]
    if (!is.vector(data.reef.i)) {
      for (j in 1: dim(data.reef.i)[1]) {lines(1996:years[k], data.reef.i[j,1:k], col=rgb(palette[[i]](.5),max=255), lty=2, lwd=.5)}
    } else {lines(1996:years[k], data.reef.i[1:k], col=rgb(palette[[i]](.5),max=255), lty=2, lwd=.5)}
    
    # Add confidence interval and median at reef level
    polygon (c(1996:years[k], years[k]:1996), c(res.manta.min[i,1:k], res.manta.max[i,k:1]), border = NA,col=rgb(palette[[i]](.5),max=255), fillOddEven=T)
    polygon (c(1996:years[k], years[k]:1996), c(res.manta.25[i,1:k], res.manta.75[i,k:1]), border = NA,col=rgb(palette[[i]](.7),max=255), fillOddEven=T)
    lines(1996:years[k], res.manta.med[i,1:k], col=rgb(palette[[i]](.9),max=255), lwd=3)
    title(main= unique(obs.manta$REEF_NAME[obs.manta$REEF_ID == res.manta.fid[i]]))
    
    # Add Manta observations
    manta.i <- subset(obs.manta, REEF_ID == res.manta.fid[i] & REPORT_YEAR <= years[k], select=c(REPORT_YEAR, MEAN_LIVE.corr))
    manta.i <- manta.i[order(manta.i$REPORT_YEAR),]
    lines(manta.i$REPORT_YEAR, manta.i$MEAN_LIVE.corr)
    points(manta.i$REPORT_YEAR, manta.i$MEAN_LIVE.corr, pch=19, cex=1.5, col="black")
    }
}

HC.model <- function(bleaching.freq, cots.freq, storms.freq) {
for (j in 1:nsimul) {

  # Sample initial parameters from normal distributions
  HC.asym <- HCMAX[,j]
  HC.1996 <- HCINI[,j]
  b0 <- B0[,j]
  b0.wq <- b0 + WQ * rnorm(length(WQ), mean=WQ.mn.sd[1], sd=WQ.mn.sd[2])
  b1.wq <- b0.wq / log(HC.asym)
  res[,1,j] <- as.numeric(HC.1996)
  HC.tmp <- log(HC.1996)

  # Re-initialize disturbance data
  data.bleaching <- data.bleaching.bckp
  data.COTS <- data.COTS.bckp
  data.storms <-  data.storms.bckp

  # Resample disturbance data in each year
  data.unknown <- data.disease <- data.COTS
  data.unknown[,3:22] <- data.disease[,3:22] <- 0
  for (i in 1:20) {
    ## Simulate disease and unknown disturbance based on observed frequencies
    data.unknown[,i+2] <- srswor(round(length(data.unknown[,i+2])*0.01),length(data.unknown[,i+2]))
    data.disease[,i+2] <- srswor(round(length(data.disease[,i+2])*0.01),length(data.disease[,i+2]))
    ## Resample other disturbance based on P(Impact|Disturbance) and as a function of distance to observations
    count.cots <- length(data.COTS[,i+2][data.COTS[,i+2]>0])
    count.storms <- length(data.storms[,i+2][data.storms[,i+2]>0])
    if (count.cots>0)  data.COTS[,i+2][data.COTS[,i+2]>0][sample(count.cots, count.cots*cots.freq)] <- 0
    if (count.storms>0)  data.storms[,i+2][data.storms[,i+2]>0][sample(count.storms, round(count.storms*storms.freq))] <- 0
    if (i < 3) {
      count.bleaching <- length(data.bleaching[,i+2][data.bleaching[,i+2]>0])
      if (count.bleaching>0)  data.bleaching[,i+2][data.bleaching[,i+2]>0][sample(count.bleaching, count.bleaching*bleaching.freq)] <- 0
    }
  }

  # Add known disturbance for LTMP reefs
  data.bleaching[,-(1:2)][!is.na(data.ltmp.bleaching[,6]),] <- data.ltmp.bleaching[,c(7,11)][!is.na(data.ltmp.bleaching[,6]),]
  data.COTS[,4:22][!is.na(data.ltmp.COTS[,6]),] <- data.ltmp.COTS[,6:24][!is.na(data.ltmp.COTS[,6]),]
  data.storms[,4:22][!is.na(data.ltmp.storms[,6]),] <- data.ltmp.storms[,6:24][!is.na(data.ltmp.storms[,6]),]
  data.disease[,4:22][!is.na(data.ltmp.disease[,6]),] <- data.ltmp.disease[,6:24][!is.na(data.ltmp.disease[,6]),]
  data.unknown[,4:22][!is.na(data.ltmp.unknown[,6]),] <- data.ltmp.unknown[,6:24][!is.na(data.ltmp.unknown[,6]),]

  data.bleaching[,-(1:2)][data.bleaching[,-(1:2)]>60] <- 60
  data.COTS[,3:22][data.COTS[,3:22]>86] <- 86
  data.storms[,3:22][data.storms[,3:22]>75] <- 75
  data.disease[,3:22][data.disease[,3:22]>1] <- 1
  data.unknown[,3:22][data.unknown[,3:22]>1] <- 1

  data.COTS[data.grid$REEF_ID==2577,] <- data.COTS.bckp[data.grid$REEF_ID==2577,]
  data.COTS[data.grid$REEF_ID==179,] <- data.COTS.bckp[data.grid$REEF_ID==179,]
  data.COTS[data.grid$REEF_ID==1056, c("X2009","X2010")]<- c(2,3)
  data.storms[data.grid$REEF_ID==1056,"X1997"] <- 0
  data.storms[data.grid$REEF_ID==1056,"X2009"] <- 0
  data.bleaching[data.grid$REEF_ID==1056,3:4] <- c(0,0)
  data.bleaching[data.grid$REEF_ID==1121,3:4] <- c(5,10)
  #data.storms[data.grid$REEF_ID==1121,] <- data.storms.bckp[data.grid$REEF_ID==1121,]
  data.storms[data.grid$REEF_ID==1121,c("X2009","X2010","X2013","X2014")] <- c(1,0,1,1)
  data.unknown[data.grid$REEF_ID==1056,-(1:2)] <- 0
  data.disease[data.grid$REEF_ID==1056,-(1:2)] <- 0


  # Year loop
  for (i in 2:nyears) {

    ### Apply disturbances (bleaching, CoTS, disease, storms) in year 1994+i (starting from 1996)
    HC.tmp[data.COTS[,i+1]>0] <- HC.tmp[data.COTS[,i+1]>0] + data.COTS[,i+1][data.COTS[,i+1]>0] * rnorm(length(HC.tmp[data.COTS[,i+1]>0]), mean=COTS.mn.sd[1], sd=COTS.mn.sd[2])
    HC.tmp[data.storms[,i+1]>0] <- HC.tmp[data.storms[,i+1]>0] + data.storms[,i+1][data.storms[,i+1]>0] * rnorm(length(HC.tmp[data.storms[,i+1]>0]), mean=storms.mn.sd[1], sd=storms.mn.sd[2])
    if (i==3) HC.tmp[data.bleaching[,"X1998"] > 5] <- HC.tmp[data.bleaching[,"X1998"] > 5] + data.bleaching[,"X1998"][data.bleaching[,"X1998"] > 5] * rnorm(length(HC.tmp[data.bleaching[,"X1998"] > 5]), mean=bleaching.mn.sd[1], sd=bleaching.mn.sd[2])
    if (i==7) HC.tmp[data.bleaching[,"X2002"] > 5] <- HC.tmp[data.bleaching[,"X2002"] > 5] + data.bleaching[,"X2002"][data.bleaching[,"X2002"] > 5] * rnorm(length(HC.tmp[data.bleaching[,"X2002"] > 5]), mean=bleaching.mn.sd[1], sd=bleaching.mn.sd[2])
    HC.tmp[data.disease[,i+1]>0] <- HC.tmp[data.disease[,i+1]>0] + rnorm(length(HC.tmp[data.disease[,i+1]>0]), mean=disease.mn.sd[1], sd=disease.mn.sd[2])
    HC.tmp[data.unknown[,i+1]>0] <- HC.tmp[data.unknown[,i+1]>0] + rnorm(length(HC.tmp[data.unknown[,i+1]>0]), mean=unknown.mn.sd[1], sd=unknown.mn.sd[2])
    HC.tmp[HC.tmp < log(0.01)] <- log(0.01) # sets minimal value to 0.5% (as 0% does not allow for recovery. 0.5% is the minimum HC cover observed in the LTMP data)

    ### Make coral grow/recover
    HC.tmp <- b0.wq + (1 - b1.wq)* HC.tmp

    ### Store results
    res[,i,j] <- exp(HC.tmp)

  }

  bleaching.rsmpl[,3,j] <- data.bleaching[,3]
  bleaching.rsmpl[,7,j] <- data.bleaching[,4]
  COTS.rsmpl[,,j] <- as.matrix(data.COTS[,-(1:2)])
  disease.rsmpl[,,j] <- as.matrix(data.disease[,-(1:2)])
  storms.rsmpl[,,j] <- as.matrix(data.storms[,-(1:2)])
  unknown.rsmpl[,,j] <- as.matrix(data.unknown[,-(1:2)])
  #print(j)

}


  ### Compute HC stats for each grid cell and reef/cluster/bioregion and in each year (e.g. median, quartiles, 95% CI) and store them each in a matrix
  res.mn <- apply(res, c(1,2), mean, na.rm=T)
  delta.SA <- res.mn[,20] - res.mn[,1]

  ### Compute mean disturbance impact across simulations
  for (i in 1:nsimul) bleaching.rsmpl[,,i][is.na(bleaching.rsmpl[,,i])] <- 0

  bleaching.mn <- apply(bleaching.rsmpl, c(1,2), mean, na.rm=T)
  COTS.mn <- apply(COTS.rsmpl, c(1,2), mean, na.rm=T)
  disease.mn <- apply(disease.rsmpl, c(1,2), mean, na.rm=T)
  storms.mn <- apply(storms.rsmpl, c(1,2), mean, na.rm=T)
  unknown.mn <- apply(unknown.rsmpl, c(1,2), mean, na.rm=T)

  # Calculate severity of disturbance events (no of events x average impact on CC)
  cnt.bleaching <- bleaching.mn
  cnt.bleaching <- data.frame(cnt.bleaching, cnt = rowSums(cnt.bleaching, na.rm=T)*18.8)

  cnt.disease <- disease.mn
  cnt.disease <- data.frame(cnt.disease, cnt = rowSums(cnt.disease, na.rm=T)*18.6)

  cnt.COTS <- COTS.mn
  cnt.COTS <- data.frame(cnt.COTS, cnt = rowSums(cnt.COTS)*35)

  cnt.storms <- storms.mn
  cnt.storms <- data.frame(cnt.storms, cnt = rowSums(cnt.storms)*67)

  cnt.unknown <- unknown.mn
  cnt.unknown <- data.frame(cnt.unknown, cnt = rowSums(cnt.unknown)*17.4)

  cnt.TOT.SA <- cnt.COTS$cnt + cnt.bleaching$cnt + cnt.storms$cnt + cnt.unknown$cnt + cnt.disease$cnt

  ### Export results
  HC.model.out <- data.frame(delta.SA, cnt.TOT.SA, cnt.bleaching$cnt, cnt.COTS$cnt, cnt.storms$cnt)
  HC.model.out
}

lmer.predict <- function(model, terms, newdat) {
  X <- model.matrix(terms, data = newdat)
  newdat$pred <- X %*% fixef(model)
  pvar1 <- diag(X %*% tcrossprod(vcov(model), X))
  newdat <- data.frame(newdat,
                       plo = newdat$pred - (1.96 * sqrt(pvar1)),
                       phi = newdat$pred + (1.96 * sqrt(pvar1)))
  return(newdat)
  
}

#display.brewer.all() 
yellow=colorRamp(brewer.pal(9, "YlOrRd"), interpolate="linear", space="rgb", bias=.8)
blue=colorRamp(brewer.pal(9, "Blues")[3:8], interpolate="linear", space="rgb", bias=.8)
green=colorRamp(brewer.pal(9, "Greens")[5:9], interpolate="linear", space="rgb", bias=.8)
orange=colorRamp(brewer.pal(9, "YlOrBr")[2:6], interpolate="linear", space="rgb", bias=.8)
red=colorRamp(brewer.pal(9, "Reds"), interpolate="linear", space="rgb", bias=.8)
purple=colorRamp(brewer.pal(9, "Purples"), interpolate="linear", space="rgb", bias=.8)
brown=colorRamp(brewer.pal(9, "BrBG")[6:4], interpolate="linear", space="rgb", bias=.8)
greys=colorRamp(brewer.pal(9, "Greys")[2:6], interpolate="linear", space="rgb", bias=.8)

green2 <- colorRamp(brewer.pal(9, "PiYG")[6:9], interpolate="linear", space="rgb", bias=.8)
green2.nopale <- colorRamp(brewer.pal(9, "PiYG")[7:9], interpolate="linear", space="rgb", bias=.8)
magenta <- colorRamp(brewer.pal(9, "PiYG")[4:1], interpolate="linear", space="rgb", bias=.8)
magenta.nopale <- colorRamp(brewer.pal(9, "PiYG")[3:1], interpolate="linear", space="rgb", bias=.8)

# Load & filter data ------------

## Mac 
setwd("/Users/Camille/Desktop/Coral Cover model/Gompertz_GBR")

#load("Gompertz_update2.RData")
source("input/brt.functions.R")
source("input/ModelFun.R")


data.manta <- read.table("input/Manta.csv", header = TRUE, sep = ",")
data.manta.env <- read.table("input/Manta_ENV.txt", header = TRUE, sep = "\t")
data.bleaching <- read.table("input/Disturb_bleaching.txt", header = TRUE, sep = "\t")
data.COTS <- read.table("input/Disturb_COTS.txt", header = TRUE, sep = "\t")
data.disease <- read.table("input/Disturb_disease.txt", header = TRUE, sep = "\t")
data.storms <-  read.table("input/Disturb_cyclones.txt", header = TRUE, sep = "\t")
data.grid <- read.table("input/XYZ_BRTpred_MRTpred_GG.csv", header = TRUE, sep = ",")
data.grid <- data.grid[order(data.grid$LONG, data.grid$LAT),]

data.ltmp.bleaching <- read.table("input/LTMP_B_XYZ_updated.txt", header = TRUE, sep = "\t")
data.ltmp.COTS <- read.table("input/LTMP_C_XYZ_updated.txt", header = TRUE, sep = "\t")
data.ltmp.disease <- read.table("input/LTMP_D_XYZ_updated.txt", header = TRUE, sep = "\t")
data.ltmp.storms <-  read.table("input/LTMP_S_XYZ_updated.txt", header = TRUE, sep = "\t")
data.ltmp.unknown <-  read.table("input/LTMP_U_XYZ_updated.txt", header = TRUE, sep = "\t")

#data.bioreg <- read.table("bioreg.deviation_mpa.txt", header = TRUE, sep = "\t")
data.WQ <- read.table("input/pst_grid.txt", header = TRUE, sep = "\t")

data.reef <- read.table("input/MRT_reef_clust_Gompertz.txt", header = TRUE, sep = "\t") 
names(data.reef)[5:7] <- c("clust.bent", "b0.mu", "b0.sd")

data.rap <- read.table("input/KERCOD_GP_Bent_Env3.txt", header = TRUE, sep = "\t")
rap.reefs <- data.rap$REEF[data.rap$P_CODE=="RAP"]


# REMOVE NORTHERNMOST REEFS FROM SPATIAL TABLES
data.ltmp.bleaching <- data.ltmp.bleaching[data.grid$LAT < (-14),]
data.ltmp.COTS <- data.ltmp.COTS[data.grid$LAT < (-14),]
data.ltmp.disease <- data.ltmp.disease[data.grid$LAT < (-14),]
data.ltmp.storms <- data.ltmp.storms[data.grid$LAT < (-14),]
data.ltmp.unknown <- data.ltmp.unknown[data.grid$LAT < (-14),]

data.bleaching <- data.bleaching[data.grid$LAT < (-14),]
data.COTS <- data.COTS[data.grid$LAT < (-14),]
data.disease <- data.disease[data.grid$LAT < (-14),]
data.storms <- data.storms[data.grid$LAT < (-14),]
data.grid <- data.grid[data.grid$LAT < (-14),]


WQ <- data.grid$Primary + data.grid$Secondary + data.grid$Tertiary

shape <- importShapefile("/Users/Camille/Dropbox/My documents/Projects/CoTS/Modelling/Coral Cover Model_update/shapefiles/SDE_OWNER_crcgis_land250", readDBF=FALSE)
shape.poly <- readShapePoly("/Users/Camille/Dropbox/My documents/Projects/CoTS/Modelling/Coral Cover Model_update/shapefiles/SDE_OWNER_crcgis_land250", proj4string = CRS("+proj=longlat +datum=WGS84"))
shape.elide <- elide(shape.poly, rotate=45)
shape45 <- SpatialPolygons2PolySet(shape.elide)

grid <- SpatialPoints(cbind(data.grid$LONG, data.grid$LAT), proj4string = CRS("+proj=longlat +datum=WGS84"))
grid.elide <- elide(grid, bb=bbox(shape.poly), rotate=45)
grid45 <- as.data.frame(grid.elide)

pst.grid <- SpatialPoints(cbind(data.WQ$X, data.WQ$Y), proj4string = CRS("+proj=longlat +datum=WGS84"))
pst.elide45 <- elide(pst.grid, bb=bbox(shape.poly), rotate=45)
pst45 <- as.data.frame(pst.elide45)
 
gridln <- gridlines(shape.poly)
gridln.elide <- elide(gridln, bb=bbox(shape.poly), rotate=45)
gridln45 <- SpatialLines2PolySet(gridln.elide)

extrapol <- importShapefile("/Users/Camille/Dropbox/My documents/Projects/CoTS/Modelling/Coral Cover Model_update/shapefiles/Extrapolated_buffer_20km", readDBF=FALSE)
extrapol.poly <- readShapePoly("/Users/Camille/Dropbox/My documents/Projects/CoTS/Modelling/Coral Cover Model_update/shapefiles/Extrapolated_buffer_20km", proj4string = CRS("+proj=longlat +datum=WGS84"))
extrapol.elide <- elide(extrapol.poly, bb=bbox(shape.poly), rotate=45)
extrapol45 <- SpatialPolygons2PolySet(extrapol.elide)


# Convert bleaching scores to mid points ------

#data.bleaching[,c("X1998","X2002","X2016")] <- round(data.bleaching[,c("X1998","X2002","X2016")])

data.bleaching.bckp <- data.bleaching
data.bleaching[,-(1:2)][data.bleaching.bckp[,-(1:2)]==1] <- 0.05
data.bleaching[,-(1:2)][data.bleaching.bckp[,-(1:2)]==2] <- 0.15
data.bleaching[,-(1:2)][data.bleaching.bckp[,-(1:2)]==3] <- 0.45
data.bleaching[,-(1:2)][data.bleaching.bckp[,-(1:2)]==4] <- 0.80


# Create back up files of disturbance tables ----
data.bleaching[is.na(data.bleaching)] <- 0
data.COTS[is.na(data.COTS)] <- 0
data.storms[is.na(data.storms)] <- 0

data.bleaching.bckp <- data.bleaching
data.COTS.bckp <- data.COTS
data.storms.bckp <- data.storms

# Sample model parameters from theirs distributions and store them in arrays --------

# N <- 100
# HCINI <- HCMAX <- B0 <- matrix(NA, ncol = N, nrow = dim(data.grid)[1])
# 
# for (i in 1:dim(data.grid)[1]) {
# #Define sigma for ith grid cell based on covariances 
# sigma <- matrix(c(I(data.grid$pred.HCini.sd[i])^2, 0.54*data.grid$pred.HCini.sd[i]*data.grid$pred.HCmax.sd[i], 0.14*data.grid$pred.HCini.sd[i]*data.grid$pred.b0.sd[i],
#                   0.54*data.grid$pred.HCini.sd[i]*data.grid$pred.HCmax.sd[i], I(data.grid$pred.HCmax.sd[i])^2, -0.13*data.grid$pred.HCmax.sd[i]*data.grid$pred.b0.sd[i], 
#                   0.14*data.grid$pred.HCini.sd[i]*data.grid$pred.b0.sd[i], -0.13*data.grid$pred.HCmax.sd[i]*data.grid$pred.b0.sd[i], I(data.grid$pred.b0.sd[i])^2),
#                   3,3)
# 
# #Pick N random parameters for ith grid cell
# pick <- mvrnorm(n=N, mu=c(data.grid$pred.HCini.mean[i], data.grid$pred.HCmax.mean[i], data.grid$pred.b0.mean[i]), Sigma = sigma)
# HCINI[i,] <- pick[,1]
# HCMAX[i,] <- pick[,2]
# B0[i,] <- pick[,3]
# }


# Gompertz parameters: disturbance & WQ effect sizes (as per MacNeil et al. in prep)
bleaching.mn.sd <- c(-1.0182, 0.1028)
COTS.mn.sd <- c(-0.1713, 0.0097)
disease.mn.sd <- c(-0.1681, 0.0303)
storms.mn.sd <- c(-0.5232, 0.01173)
unknown.mn.sd <- c(-0.2681, 0.0331)

WQ_bleach <- c(0.851062345,0.096005095)
WQ_CoTS <- c(-2.034101163,0.134657287)
WQ_Cyclone <- c(0.649030968,0.056093201)
WQ_Disease <- c(-0.369621141,0.15480862)
#WQ_Unknown <- c(-0.242144071,0.3328742) # Overlaps zero

nsimul <- 100
LHS <- randomLHS(nsimul,9)
B.BLEACHING <- qnorm(LHS[,1], mean=bleaching.mn.sd[1], sd=bleaching.mn.sd[2])
B.COTS <- qnorm(LHS[,2], mean=COTS.mn.sd[1], sd=COTS.mn.sd[2])
B.DISEASE <- qnorm(LHS[,3], mean=disease.mn.sd[1], sd=disease.mn.sd[2])
B.STORMS <- qnorm(LHS[,4], mean=storms.mn.sd[1], sd=storms.mn.sd[2])
B.UNKNOWN <- qnorm(LHS[,5], mean=unknown.mn.sd[1], sd=unknown.mn.sd[2])
B.WQ <- qnorm(LHS[,6], mean=WQ.mn.sd[1], sd=WQ.mn.sd[2])

A <- B0 <- HCINI <- HCMAX <- matrix(NA, ncol = nsimul, nrow = dim(data.grid)[1])

for (i in 1:nsimul) {
A[,i] <- qnorm(LHS[i,7], mean=data.grid$pred.a.mean, sd=data.grid$pred.a.sd)
B0[,i] <- qnorm(LHS[i,8], mean=data.grid$pred.b0.mean, sd=data.grid$pred.b0.sd)
HCINI[,i] <- qnorm(LHS[i,9], mean=data.grid$pred.HCini.mean, sd=data.grid$pred.HCini.sd)  
HCMAX[,i] <- qnorm(LHS[i,9], mean=data.grid$pred.HCmax.mean, sd=data.grid$pred.HCmax.sd) 
}
HCMAX[HCINI > HCMAX] <- HCINI[HCINI > HCMAX]

# Loop through years and make coral grow (100 stochastic simulations) ---------------
  
  nyears <- length(1996:2017)
  nsimul <- 100
  res <- array(NA, dim=c(dim(data.grid)[1], nyears, nsimul)) ## Stores the coral cover values for each grid cell (rows), year (columns) and simulation (third dimension)
  bleaching.rsmpl <- COTS.rsmpl <- disease.rsmpl <- storms.rsmpl <- unknown.rsmpl <- array(NA, dim=c(dim(data.grid)[1], nyears, nsimul)) ## Stores the actual (i.e. resampled) disturbance intensity values for each grid cell (rows), year (columns) and simulation (third dimension)
  
  
  for (j in 1:nsimul) {

  # Define initial parameters for jth simulation
  HC.1996 <- HCINI[,j]
  b0 <- B0[,j]
  #b1 <- A[,j]
  b1 <- B0[,j]/log(HCMAX[,j])
  res[,1,j] <- as.numeric(HC.1996)
  HC.tmp <- log(HC.1996)
  
  # Re-initialize disturbance data
  data.bleaching <- data.bleaching.bckp
  data.COTS <- data.COTS.bckp
  data.storms <-  data.storms.bckp
  
  # Resample disturbance data in each year
  data.unknown <- data.disease <- data.COTS
  data.unknown[,3:24] <- data.disease[,3:24] <- 0
  for (i in 1:nyears) {
    ## Simulate disease and unknown disturbance based on observed frequencies
    data.unknown[,i+2] <- srswor(round(length(data.unknown[,i+2])*0.01),length(data.unknown[,i+2]))
    data.disease[,i+2] <- srswor(round(length(data.disease[,i+2])*0.01),length(data.disease[,i+2]))

    ## Resample other disturbance based on P(Impact|Disturbance)
    count.cots <- length(data.COTS[,i+2][data.COTS[,i+2]>0])  
    count.storms <- length(data.storms[,i+2][data.storms[,i+2]>0])
    count.bleaching <- length(data.bleaching[,i+2][data.bleaching[,i+2]>0])
    if (count.cots>0)  data.COTS[,i+2][data.COTS[,i+2]>0][sample(count.cots, count.cots*.5)] <- 0
    if (count.storms>0)  data.storms[,i+2][data.storms[,i+2]>0][sample(count.storms, round(count.storms*.5))] <- 0     
    if (count.bleaching>0)  data.bleaching[,i+2][data.bleaching[,i+2]>0][sample(count.bleaching, count.bleaching*.1)] <- 0
  }
  data.disease[,-(1:2)][data.disease[,-(1:2)]>1] <- 1
  data.unknown[,-(1:2)][data.unknown[,-(1:2)]>1] <- 1
  
  data.storms[,"X2009"] <- data.storms[,"X2009"]*.5
  data.storms[,"X2013"] <- data.storms[,"X2013"]*.5
  data.storms[,"X2014"] <- data.storms[,"X2014"]*.5
  data.storms[,"X2015"] <- data.storms[,"X2015"]*.5
  
  # Add known disturbance for LTMP reefs
  data.bleaching[,-(1:2)][!is.na(data.ltmp.bleaching[,-(1:4)])] <- data.ltmp.bleaching[,-(1:4)][!is.na(data.ltmp.bleaching[,-(1:4)])]
  data.COTS[,-(1:2)][!is.na(data.ltmp.COTS[,-(1:4)])] <- data.ltmp.COTS[,-(1:4)][!is.na(data.ltmp.COTS[,-(1:4)])]
  data.storms[,-(1:2)][!is.na(data.ltmp.storms[,-(1:4)])] <- data.ltmp.storms[,-(1:4)][!is.na(data.ltmp.storms[,-(1:4)])]
  data.disease[,-(1:2)][!is.na(data.ltmp.disease[,-(1:4)])] <- data.ltmp.disease[,-(1:4)][!is.na(data.ltmp.disease[,-(1:4)])]
  data.unknown[,-(1:2)][!is.na(data.ltmp.unknown[,-(1:4)])] <- data.ltmp.unknown[,-(1:4)][!is.na(data.ltmp.unknown[,-(1:4)])]
  

  
  
  
  # Year loop
  for (i in 2:nyears) {
    
    ### Apply disturbances (bleaching, CoTS, disease, storms) in year 1994+i (starting from 1996)
    data.storms[,i+1][WQ > (-1)*B.STORMS[j]/WQ_Cyclone[1]] <- 0
    HC.tmp[data.storms[,i+1]>0] <- HC.tmp[data.storms[,i+1]>0] + data.storms[,i+1][data.storms[,i+1]>0] * (B.STORMS[j] + WQ[data.storms[,i+1]>0] * WQ_Cyclone[1])
    HC.tmp[data.COTS[,i+1]>0] <- HC.tmp[data.COTS[,i+1]>0] + data.COTS[,i+1][data.COTS[,i+1]>0] * (B.COTS[j] + WQ[data.COTS[,i+1]>0] * WQ_CoTS[1])
    #HC.tmp[data.COTS[,i+1]>0] <- HC.tmp[data.COTS[,i+1]>0] + data.COTS[,i+1][data.COTS[,i+1]>0] * B.COTS[j]
    HC.tmp[data.bleaching[,i+1]>0] <- HC.tmp[data.bleaching[,i+1]>0] + data.bleaching[,i+1][data.bleaching[,i+1]>0] * (B.BLEACHING[j] + WQ[data.bleaching[,i+1]>0] * WQ_bleach[1])
    HC.tmp[data.disease[,i+1]>0] <- HC.tmp[data.disease[,i+1]>0] + B.DISEASE[j] 
    HC.tmp[data.unknown[,i+1]>0] <- HC.tmp[data.unknown[,i+1]>0] + B.UNKNOWN[j] 
    HC.tmp[HC.tmp < log(0.1)] <- log(0.1) # sets minimal value to 10% (as 0% does not allow for recovery. 10% is the minimum HC cover observed in the LTMP data)
 
    ### Make coral grow/recover
    HC.tmp <- b0 + (1 - b1)* HC.tmp
    
    ### Store results
    res[,i,j] <- exp(HC.tmp)
    
  }
  
  bleaching.rsmpl[,,j] <- as.matrix(data.bleaching[,-(1:2)])
  COTS.rsmpl[,,j] <- as.matrix(data.COTS[,-(1:2)])
  disease.rsmpl[,,j] <- as.matrix(data.disease[,-(1:2)])
  storms.rsmpl[,,j] <- as.matrix(data.storms[,-(1:2)])
  unknown.rsmpl[,,j] <- as.matrix(data.unknown[,-(1:2)])
  print(j)

  }
  

# Compute summary stats ---------------------  
  
  ### Compute HC stats at reef, cluster and bioregion levels in each year {
  res.reef <- array(0, dim=c(length(unique(data.grid$REEF_ID)), nyears, nsimul))  
  res.cluster <- array(0, dim=c(length(unique(data.grid$bent.clust)), nyears, nsimul))
  # res.bioregion <- array(0, dim=c(length(unique(data.grid$BIOREGION)), nyears, nsimul))
  
  for (i in 1:dim(res)[3]) {
    res.reef[,,i] <- as.matrix(aggregate(res[,,i], by=list(data.grid$REEF_ID), FUN=mean, na.rm=T)[-1])
    res.cluster[,,i] <- as.matrix(aggregate(res[,,i], by=list(data.grid$bent.clust), FUN=mean, na.rm=T)[-1])
  #   res.bioregion[,,i] <- as.matrix(aggregate(res[,,i], by=list(data.grid$BIOREGION), FUN=mean, na.rm=T)[-1])
   }
  
  ### Compute HC stats for each grid cell and reef/cluster/bioregion and in each year (e.g. median, quartiles, 95% CI) and store them each in a matrix
  res.mn <- apply(res, c(1,2), mean, na.rm=T)
  res.med <- apply(res, c(1,2), median, na.rm=T)
  res.min <- apply(res, c(1,2), quantile, probs=0.05, na.rm=T)
  res.max <- apply(res, c(1,2), quantile, probs=0.95, na.rm=T)
  res.25 <- apply(res, c(1,2), quantile, probs=0.25, na.rm=T)
  res.75 <- apply(res, c(1,2), quantile, probs=0.75, na.rm=T)
  
  res.gbr.mn <- colMeans(res.mn, na.rm=T)
  res.gbr.med <- colMeans(res.med, na.rm=T)
  res.gbr.min <- colMeans(res.min, na.rm=T)
  res.gbr.max <- colMeans(res.max, na.rm=T)
  res.gbr.25 <- colMeans(res.25, na.rm=T)
  res.gbr.75 <- colMeans(res.75, na.rm=T)

  res.reef.med <- apply(res.reef, c(1,2), median, na.rm=T)
  res.reef.min <- apply(res.reef, c(1,2), quantile, probs=0.05, na.rm=T)
  res.reef.max <- apply(res.reef, c(1,2), quantile, probs=0.95, na.rm=T)
  res.reef.25 <- apply(res.reef, c(1,2), quantile, probs=0.25, na.rm=T)
  res.reef.75 <- apply(res.reef, c(1,2), quantile, probs=0.75, na.rm=T)
  
  res.cluster.med <- apply(res.cluster, c(1,2), median, na.rm=T)
  res.cluster.mn <- apply(res.cluster, c(1,2), mean, na.rm=T)
  res.cluster.min <- apply(res.cluster, c(1,2), quantile, probs=0.05, na.rm=T)
  res.cluster.max <- apply(res.cluster, c(1,2), quantile, probs=0.95, na.rm=T)
  res.cluster.25 <- apply(res.cluster, c(1,2), quantile, probs=0.25, na.rm=T)
  res.cluster.75 <- apply(res.cluster, c(1,2), quantile, probs=0.75, na.rm=T)
  # 
  # res.bioregion.med <- apply(res.bioregion, c(1,2), median)
  # res.bioregion.min <- apply(res.bioregion, c(1,2), quantile, probs=0.05, na.rm=T)
  # res.bioregion.max <- apply(res.bioregion, c(1,2), quantile, probs=0.95, na.rm=T)
  # res.bioregion.25 <- apply(res.bioregion, c(1,2), quantile, probs=0.25, na.rm=T)
  # res.bioregion.75 <- apply(res.bioregion, c(1,2), quantile, probs=0.75, na.rm=T)
  
  # Compute mean disturbance impact across simulations
  for (i in 1:nsimul) bleaching.rsmpl[,,i][is.na(bleaching.rsmpl[,,i])] <- 0
  
  bleaching.mn <- apply(bleaching.rsmpl, c(1,2), mean, na.rm=T)
  COTS.mn <- apply(COTS.rsmpl, c(1,2), mean, na.rm=T)
  disease.mn <- apply(disease.rsmpl, c(1,2), mean, na.rm=T)
  storms.mn <- apply(storms.rsmpl, c(1,2), mean, na.rm=T)
  unknown.mn <- apply(unknown.rsmpl, c(1,2), mean, na.rm=T)
  
  
# Exploratory plots: HC trajectories for the first 10 reefs -----------------
  
  palette=list(eval(bquote(purple)), eval(bquote(orange)), eval(bquote(green)), eval(bquote(blue)), eval(bquote(red)), eval(bquote(purple)), eval(bquote(magenta)),
               eval(bquote(blue)), eval(bquote(purple)), eval(bquote(blue)), eval(bquote(green)), eval(bquote(orange)), eval(bquote(red)), eval(bquote(purple)), eval(bquote(magenta)),
               eval(bquote(green)), eval(bquote(blue)), eval(bquote(green)), eval(bquote(orange)), eval(bquote(red)), eval(bquote(purple)), eval(bquote(magenta)))
               
  par(mfcol=c(5,2), mai=c(0,0,.2,0))

  for (i in 1:10) {
  plot(NA, xlim=c(1996,2017), ylim=c(0,100), xlab="Time (yr)", ylab="Coral cover (%)", xaxs="i", yaxs="i", las=1)
  polygon (c(1996:2017, 2017:1996), c(res.reef.min[i,], res.reef.max[i,][nyears:1]), border = NA,col=rgb(palette[[i]](.5),max=255), fillOddEven=T)
  polygon (c(1996:2017, 2017:1996), c(res.reef.25[i,], res.reef.75[i,][nyears:1]), border = NA,col=rgb(palette[[i]](.7),max=255), fillOddEven=T)
  lines(1996:2017, res.reef.med[i,], col=rgb(palette[[i]](.9),max=255), lwd=3)
  
  # Add individual grid cells (dashed lines)
  reef.i <- aggregate(res[,,i], by=list(data.grid$REEF_ID), FUN=mean)[i,1]
  title(main= paste("Reef", reef.i))
  data.reef.i <- res.med[data.grid$REEF_ID==reef.i,]
  if (!is.vector(data.reef.i)) {
    for (j in 1: dim(data.reef.i)[1]) {lines(1996:2017, data.reef.i[j,], col=rgb(palette[[i]](.5),max=255), lwd=2, lty=2)}
  } else {lines(1996:2017, data.reef.i, col=rgb(palette[[i]](.5),max=255), lwd=2, lty=2)}
  }
  

# Plot time series of observed vs. predicted coral cover for validation manta-tow reefs -----------
  
  ## Select manta reefs by FULLREEF_ID
  data.manta <- subset(data.manta, REEF_LAT < (-14)) # Remove norternmost reefs (off the grid)
  manta.recent <- data.manta[data.manta$REPORT_YEAR>1995,]
  reef10.ls <- names(table(manta.recent$FULLREEF_ID)[table(manta.recent$FULLREEF_ID)>9]) ## n=54 reefs with at least 10 years of data post 1995
  reef.calib.ls <- reef10.ls[reef10.ls %in% data.reef$FULLREEF_ID]
  reef.valid.ls <- reef10.ls[!(reef10.ls %in% data.reef$FULLREEF_ID)]

  ## Link to REEF_ID through data.manta.env and create lists of FID (REEF_ID)
  obs.manta.tmp <- subset(data.manta, FULLREEF_ID %in% reef10.ls)
  obs.manta <- merge.data.frame(obs.manta.tmp, data.manta.env[,c("FULLREEF_ID", "REEF_ID")])
  obs.manta <- obs.manta[order(obs.manta$REEF_ID, obs.manta$REPORT_YEAR),]
  obs.fid <- unique(obs.manta$REEF_ID)[order(unique(obs.manta$REEF_ID))]

  obs.calib <- obs.manta[obs.manta$FULLREEF_ID %in% reef.calib.ls & obs.manta$REEF_LAT < (-14),]
  obs.valid <- obs.manta[obs.manta$FULLREEF_ID %in% reef.valid.ls & obs.manta$REEF_LAT < (-14),]
  obs.calib.fid <- unique(obs.calib$REEF_ID)
  obs.valid.fid <- unique(obs.valid$REEF_ID)

  ## Select model predictions corresponding to CALIBRATION and VALIDATION reef sets
  res.reef.fid <- aggregate(res[,,i], by=list(data.grid$REEF_ID), FUN=mean, na.rm=T)[,1]
  
  res.calib.med <- res.reef.med[res.reef.fid %in% obs.calib.fid,]
  res.calib.min <- res.reef.min[res.reef.fid %in% obs.calib.fid,]
  res.calib.max <- res.reef.max[res.reef.fid %in% obs.calib.fid,]
  res.calib.25 <- res.reef.25[res.reef.fid %in% obs.calib.fid,]
  res.calib.75 <- res.reef.75[res.reef.fid %in% obs.calib.fid,]
  res.calib.fid <- res.reef.fid[res.reef.fid %in% obs.calib.fid]

  res.valid.med <- res.reef.med[res.reef.fid %in% obs.valid.fid,]
  res.valid.min <- res.reef.min[res.reef.fid %in% obs.valid.fid,]
  res.valid.max <- res.reef.max[res.reef.fid %in% obs.valid.fid,]
  res.valid.25 <- res.reef.25[res.reef.fid %in% obs.valid.fid,]
  res.valid.75 <- res.reef.75[res.reef.fid %in% obs.valid.fid,]
  res.valid.fid <- res.reef.fid[res.reef.fid %in% obs.valid.fid]
  
  res.manta.med <- res.reef.med[res.reef.fid %in% obs.fid,]
  res.manta.min <- res.reef.min[res.reef.fid %in% obs.fid,]
  res.manta.max <- res.reef.max[res.reef.fid %in% obs.fid,]
  res.manta.25 <- res.reef.25[res.reef.fid %in% obs.fid,]
  res.manta.75 <- res.reef.75[res.reef.fid %in% obs.fid,]
  res.manta.fid <- res.reef.fid[res.reef.fid %in% obs.fid]
  

  ## Plot 10 validation reefs and compute validation metrics
  
  k <- 22 # Change to plot incomplete time series
  years <- 1996:2017
  pred.obs.all <- ""
  
  par(mfcol=c(5,2), mai=c(.3,.4,.3,.3), family="Calibri")
  
  for (i in 1:10) {
    
    print(unique(obs.valid$REEF_NAME[obs.valid$REEF_ID == obs.valid.fid[i]]))
    print(data.grid$bent.clust[data.grid$REEF_ID==obs.valid.fid[i]][1])
    
    plot(NA, xlim=c(1996,2017), ylim=c(0,80), xlab="Time (yr)", ylab="Coral cover (%)", xaxs="i", yaxs="i", las=1)
    
    # Add disturbance events
    abline(v = years[colMeans(data.COTS[data.grid$REEF_ID==obs.valid.fid[i],-(1:2)]) > .2], col="darkgoldenrod1", lwd=2)
    abline(v = years[colSums(data.bleaching[data.grid$REEF_ID==obs.valid.fid[i],-(1:2)]) > 0], col="powderblue", lwd=2)
    abline(v = years[colSums(data.ltmp.disease[data.grid$REEF_ID==obs.valid.fid[i],-(1:4)]) > 0], col="lightgreen", lwd=2)
    abline(v = years[colSums(data.ltmp.unknown[data.grid$REEF_ID==obs.valid.fid[i],-(1:4)]) > 0], col="lightgrey", lwd=2)
    abline(v = years[colSums(data.storms[data.grid$REEF_ID==obs.valid.fid[i],-(1:2)]) > 0], col="lightcoral", lwd=2)
    
    # Add individual grid cells
    # data.reef.i <- res.med[data.grid$REEF_ID==obs.valid.fid[i],]
    # if (!is.vector(data.reef.i)) {
    #   for (j in 1: dim(data.reef.i)[1]) {lines(1996:years[k], data.reef.i[j,1:k], col=rgb(blue(.5),max=255), lty=2, lwd=.5)}
    # } else {lines(1996:years[k], data.reef.i[1:k], col=rgb(blue(.5),max=255), lty=2, lwd=.5)}
    
    # Add confidence interval and median at reef level
    polygon (c(1996:years[k], years[k]:1996), c(res.valid.min[i,1:k], res.valid.max[i,k:1]), border = NA,col=rgb(blue(.4),max=255), fillOddEven=T)
    polygon (c(1996:years[k], years[k]:1996), c(res.valid.25[i,1:k], res.valid.75[i,k:1]), border = NA,col=rgb(blue(.6),max=255), fillOddEven=T)
    lines(1996:years[k], res.valid.med[i,1:k], col=rgb(blue(.9),max=255), lwd=3)
    text(2016.5, 65, pos=2, unique(obs.valid$REEF_NAME[obs.valid$REEF_ID == obs.valid.fid[i]]), cex=1.3, font=2)
    
    # Add Manta observations
    manta.i <- subset(obs.manta, REEF_ID == obs.valid.fid[i] & REPORT_YEAR <= years[k], select=c(REEF_NAME, REEF_LONG, REEF_LAT, REPORT_YEAR, MEAN_LIVE.corr))
    manta.i <- manta.i[order(manta.i$REPORT_YEAR),]
    lines(manta.i$REPORT_YEAR, manta.i$MEAN_LIVE.corr)
    points(manta.i$REPORT_YEAR, manta.i$MEAN_LIVE.corr, pch=19, cex=1.5, col="black")
    
    #Compute validation table
    pred.obs.i <- data.frame(manta.i[,1:3][manta.i$REPORT_YEAR %in% 1996:2017,], year = as.numeric(manta.i$REPORT_YEAR[manta.i$REPORT_YEAR %in% 1996:2017]), pred = as.numeric(res.valid.med[i,][1996:2017 %in% manta.i$REPORT_YEAR]), obs = as.numeric(manta.i$MEAN_LIVE.corr[manta.i$REPORT_YEAR %in% 1996:2017]))
    pred.obs.all <- rbind(pred.obs.all, pred.obs.i)
    #print(paste(i, obs.valid.fid[i], mean(data.grid$LAT[data.grid$REEF_ID==obs.valid.fid[i]])), sep="          ")
    
  }
  
  pred.obs.all <- pred.obs.all[-1,]
  pred.obs.all$pred <- as.numeric(pred.obs.all$pred)
  pred.obs.all$obs <- as.numeric(pred.obs.all$obs)
  #plot(pred.obs.all$pred, pred.obs.all$obs, xlim=c(10,70), ylim=c(10,70))
  lm <- lm(pred ~ obs, data = pred.obs.all)
  summary(lm)$r.squared  ## R2 = 0.637
  mean(abs(summary(lm)$residuals)) ## Mean prediction error = 7.7 %
  
  


# Plot time series of observed vs. predicted coral cover for calibration manta-tow reefs -----------

k <- 22 # Change to plot partial time series
years <- 1996:2017

par(mfcol=c(5, 4), mai=c(.2,.3,.2,.2))

for (i in c(1:5, 8:14,16:23)) {
  plot(NA, xlim=c(1996,2017), ylim=c(0,80), xlab="Time (yr)", ylab="Coral cover (%)", xaxs="i", yaxs="i", las=1)
  
  # Add disturbance events
  abline(v = years[colMeans(data.COTS[data.grid$REEF_ID==obs.calib.fid[i],-(1:2)]) > .2], col="darkgoldenrod1", lwd=2)
  abline(v = years[colSums(data.bleaching[data.grid$REEF_ID==obs.calib.fid[i],-(1:2)]) > 0], col="powderblue", lwd=2)
  abline(v = years[colSums(data.ltmp.disease[data.grid$REEF_ID==obs.calib.fid[i],-(1:4)]) > 0], col="lightgreen", lwd=2)
  abline(v = years[colSums(data.ltmp.unknown[data.grid$REEF_ID==obs.calib.fid[i],-(1:4)]) > 0], col="lightgrey", lwd=2)
  abline(v = years[colSums(data.storms[data.grid$REEF_ID==obs.calib.fid[i],-(1:2)]) > 0], col="lightcoral", lwd=2)
  
  # Add individual grid cells
  data.reef.i <- res.med[data.grid$REEF_ID==obs.calib.fid[i],]
  if (!is.vector(data.reef.i)) {
    for (j in 1: dim(data.reef.i)[1]) {lines(1996:years[k], data.reef.i[j,1:k], col=rgb(blue(.5),max=255), lty=2, lwd=.5)}
  } else {lines(1996:years[k], data.reef.i[1:k], col=rgb(blue(.5),max=255), lty=2, lwd=.5)}
  
  # Add confidence interval and median at reef level
  polygon (c(1996:years[k], years[k]:1996), c(res.calib.min[i,1:k], res.calib.max[i,k:1]), border = NA,col=rgb(blue(.5),max=255), fillOddEven=T)
  polygon (c(1996:years[k], years[k]:1996), c(res.calib.25[i,1:k], res.calib.75[i,k:1]), border = NA,col=rgb(blue(.7),max=255), fillOddEven=T)
  lines(1996:years[k], res.calib.med[i,1:k], col=rgb(blue(.9),max=255), lwd=3)
  text(2016.5, 65, pos=2, unique(obs.calib$REEF_NAME[obs.calib$REEF_ID == obs.calib.fid[i]]), cex=1.3)
  
  # Add Manta observations
  manta.i <- subset(obs.manta, REEF_ID == obs.calib.fid[i] & REPORT_YEAR <= years[k], select=c(REEF_NAME, REEF_LONG, REEF_LAT, REPORT_YEAR, MEAN_LIVE.corr))
  manta.i <- manta.i[order(manta.i$REPORT_YEAR),]
  lines(manta.i$REPORT_YEAR, manta.i$MEAN_LIVE.corr)
  points(manta.i$REPORT_YEAR, manta.i$MEAN_LIVE.corr, pch=19, cex=1.5, col="black")
  
  print(paste(unique(obs.calib$REEF_NAME[obs.calib$REEF_ID == obs.calib.fid[i]]), obs.calib.fid[i], mean(data.grid$LAT[data.grid$REEF_ID==obs.calib.fid[i]])), sep="          ")

  
  #Compute validation table
  pred.obs.i <- data.frame(manta.i[,1:3][manta.i$REPORT_YEAR %in% 1996:2017,], year = as.numeric(manta.i$REPORT_YEAR[manta.i$REPORT_YEAR %in% 1996:2017]), pred = as.numeric(res.calib.med[i,][1996:2017 %in% manta.i$REPORT_YEAR]), obs = as.numeric(manta.i$MEAN_LIVE.corr[manta.i$REPORT_YEAR %in% 1996:2017]))
  pred.obs.all <- rbind(pred.obs.all, pred.obs.i)
  
  }

  #write.table(pred.obs.all, file= "Residuals by reefNyear.txt", quote= FALSE, sep="\t", row.names=FALSE)

# Plot time series of observed vs. predicted coral cover for RAP reefs -----------
  
  ## Select RAP reefs and link to REEF_ID
  obs.rap.tmp <- subset(data.manta, REEF_NAME %in% rap.reefs & REPORT_YEAR>1995)
  obs.rap <- merge.data.frame(obs.rap.tmp, data.manta.env[,c("FULLREEF_ID", "REEF_ID")])
  obs.rap <- obs.rap[order(obs.rap$REEF_ID, obs.rap$REPORT_YEAR),]
  obs.rap.fid <- unique(obs.rap$REEF_ID)[order(unique(obs.rap$REEF_ID))]
  obs.rap$REEF_NAME[obs.rap$REEF_ID == obs.rap.fid[40]] <- "FORK REEF"
  
  ## Select model predictions corresponding to CALIBRATION and VALIDATION reef sets
  res.reef.fid <- aggregate(res[,,i], by=list(data.grid$REEF_ID), FUN=mean, na.rm=T)[,1]

  res.rap.med <- res.reef.med[res.reef.fid %in% obs.rap.fid,]
  res.rap.min <- res.reef.min[res.reef.fid %in% obs.rap.fid,]
  res.rap.max <- res.reef.max[res.reef.fid %in% obs.rap.fid,]
  res.rap.25 <- res.reef.25[res.reef.fid %in% obs.rap.fid,]
  res.rap.75 <- res.reef.75[res.reef.fid %in% obs.rap.fid,]
  res.rap.fid <- res.reef.fid[res.reef.fid %in% obs.rap.fid]
  
  
  ## Plot 10 validation reefs and compute validation metrics
  
  k <- 22 # Change to plot incomplete time series
  years <- 1996:2017
  
  par(mfcol=c(5, 4), mai=c(.2,.3,.2,.2), family="Calibri")
  
  
  for (i in 1:20) {
    
    print(unique(obs.rap$REEF_NAME[obs.rap$REEF_ID == obs.rap.fid[i]]))
    print(obs.rap.fid[i])
    print(data.grid$bent.clust[data.grid$REEF_ID==obs.rap.fid[i]][1])
    
    plot(NA, xlim=c(1996,2017), ylim=c(0,80), xlab="Time (yr)", ylab="Coral cover (%)", xaxs="i", yaxs="i", las=1)
    
    # Add disturbance events
    abline(v = years[colMeans(data.COTS[data.grid$REEF_ID==obs.rap.fid[i],-(1:2)]) > .2], col="darkgoldenrod1", lwd=2)
    abline(v = years[colSums(data.bleaching[data.grid$REEF_ID==obs.rap.fid[i],-(1:2)]) > 0], col="powderblue", lwd=2)
    abline(v = years[colSums(data.ltmp.disease[data.grid$REEF_ID==obs.rap.fid[i],-(1:4)]) > 0], col="lightgreen", lwd=2)
    abline(v = years[colSums(data.ltmp.unknown[data.grid$REEF_ID==obs.rap.fid[i],-(1:4)]) > 0], col="lightgrey", lwd=2)
    abline(v = years[colSums(data.storms[data.grid$REEF_ID==obs.rap.fid[i],-(1:2)]) > 0], col="lightcoral", lwd=2)
    
    # Add individual grid cells
    # data.reef.i <- res.med[data.grid$REEF_ID==obs.rap.fid[i],]
    # if (!is.vector(data.reef.i)) {
    #   for (j in 1: dim(data.reef.i)[1]) {lines(1996:years[k], data.reef.i[j,1:k], col=rgb(blue(.5),max=255), lty=2, lwd=.5)}
    # } else {lines(1996:years[k], data.reef.i[1:k], col=rgb(blue(.5),max=255), lty=2, lwd=.5)}
    
    # Add confidence interval and median at reef level
    polygon (c(1996:years[k], years[k]:1996), c(res.rap.min[i,1:k], res.rap.max[i,k:1]), border = NA,col=rgb(blue(.4),max=255), fillOddEven=T)
    polygon (c(1996:years[k], years[k]:1996), c(res.rap.25[i,1:k], res.rap.75[i,k:1]), border = NA,col=rgb(blue(.6),max=255), fillOddEven=T)
    lines(1996:years[k], res.rap.med[i,1:k], col=rgb(blue(.9),max=255), lwd=3)
    text(2016.5, 65, pos=2, unique(obs.rap$REEF_NAME[obs.rap$REEF_ID == obs.rap.fid[i]]), cex=1.3, font=2)
    
    # Add Manta observations
    manta.i <- subset(obs.rap, REEF_ID == obs.rap.fid[i] & REPORT_YEAR <= years[k], select=c(REEF_NAME, REEF_LONG, REEF_LAT, REPORT_YEAR, MEAN_LIVE.corr))
    manta.i <- manta.i[order(manta.i$REPORT_YEAR),]
    lines(manta.i$REPORT_YEAR, manta.i$MEAN_LIVE.corr)
    points(manta.i$REPORT_YEAR, manta.i$MEAN_LIVE.corr, pch=19, cex=1.5, col="black")
    
  }
  
  
  
  
  
# Map coral decline and plot time series for selected reefs ---------------------

delta <- res.mn[,22] - res.mn[,1]
data.map <- data.frame(x = data.grid$LONG, y = data.grid$LAT, z = delta)
data.map.neg <- data.map[data.map$z<=0,]
data.map.pos <- data.map[data.map$z>0,]
  
  
fid <- c(2295, 3028, 1121)
yellow=colorRamp(brewer.pal(9, "YlOrRd"), interpolate="linear", space="rgb", bias=.8)
turquoise=colorRamp(brewer.pal(9, "PRGn")[5:10], interpolate="linear", space="rgb", bias=.8)
orange=colorRamp(brewer.pal(9, "Oranges")[2:7], interpolate="linear", space="rgb", bias=.8)

layout(rbind(c(1, 2), c(1, 3), c(1,4)))

# GBR map
par(mai=c(.8,.8,.2,.2), family="Calibri")
plotPolys(shape, plt=NULL, col="gray95", border="gray85", xlim=c(143,155), ylim=c(-25,-14), bg="white", cex=1.2, xlab="Longitude", ylab="Latitude", lwd=.5)
points(data.map.neg$x, data.map.neg$y, pch=19, cex=.2, col=rgb(orange(scale(-1*data.map.neg$z)),max=255))
points(data.map.pos$x, data.map.pos$y, pch=19, cex=.2, col=rgb(turquoise(scale(data.map.pos$z)),max=255))
#title(main="Annual change in coral cover (% y-1)", cex.main=2)
#Add legend
for (i in 0:100) {segments(x0 = 151, y0 = -16 - i/100, x1 = 152, col=rgb(orange(i/100),max=255))}
for (i in 0:50) {segments(x0 = 151, y0 = -16 + i/100, x1 = 152, col=rgb(turquoise(i/50),max=255))}
text(152, -15.5, pos=4, labels = paste("+",round(max(delta, na.rm=T)/20,1),sep=""), cex=1.1)
text(152, -17, pos=4, labels = round(min(delta, na.rm=T)/20,1), cex=1.1)
text(152, -16, pos=4, labels="0")
#text(150, -15, labels = "Annual rate of change (% y-1)", font=2, cex=1.3)

for (i in 1:3) {
  x.mn <- mean(obs.manta$REEF_LONG[obs.manta$REEF_ID==fid[i]])
  y.mn <- mean(obs.manta$REEF_LAT[obs.manta$REEF_ID==fid[i]])
  rect(x.mn-.7, y.mn-.5, x.mn+.7, y.mn+.5, lwd=2, border=rgb(blue(.9),max=255))
  text(x.mn+.4, y.mn+.3, labels=as.character(c("A","B","C")[i]), cex=1.2, font=2)
}

# Add extrapolation contour
addPolys(extrapol, col="lightgrey", border="lightgrey", lwd=2, density=8, colHoles=NA)

points(145.78, -16.92, pch=17, cex=1.5)
points(146.82, -19.26, pch=17, cex=1.5)
points(149.18, -21.14, pch=17, cex=1.5)
points(150.51, -23.28, pch=17, cex=1.5)

text(145.78, -16.92, "Cairns", pos=2)
text(146.82, -19.26, "Townsville", pos=2)
text(149.18, -21.14, "Mackay", pos=2)
text(150.51, -23.28, "Rockampton", pos=2)

text(145, -24, expression(italic("AUSTRALIA")), cex=1.5, col="dimgrey")
text(151.5, -18, expression(italic("CORAL SEA")), cex=1.5, col="dimgrey")

text(154.5, -15, expression(bold("Annual change in coral cover (% y-1)")), pos= 2, cex=1.3)

scalebar(loc=c(147.5, -24.1), length=2, unit="km", cex=1.2)
northarrow(loc=c(145,-22), size=.8)
box()

for (i in 1:3) {
  
  par(mai=c(.5,.5,.5,.5), family="Calibri")  
  if (i == 3)   {plot(NA, xlim=c(1996,2017), ylim=c(0,60), xlab="Time (yr)", ylab="Coral cover (%)", xaxs="i", yaxs="i", las=1)
  } else 
  { plot(NA, xlim=c(1996,2017), ylim=c(0,60), xlab=NA, ylab="Coral cover (%)", xaxt="n", xaxs="i", yaxs="i", las=1)
  }
  
  # Add disturbance events
  abline(v = years[colMeans(data.COTS[data.grid$REEF_ID==fid[i],-(1:2)]) > .2], col="darkgoldenrod1", lwd=2)
  abline(v = years[colSums(data.bleaching[data.grid$REEF_ID==fid[i],-(1:2)]) > 0], col="powderblue", lwd=2)
  abline(v = years[colSums(data.ltmp.disease[data.grid$REEF_ID==fid[i],-(1:4)]) > 0], col="lightgreen", lwd=2)
  abline(v = years[colSums(data.ltmp.unknown[data.grid$REEF_ID==fid[i],-(1:4)]) > 0], col="lightgrey", lwd=2)
  abline(v = years[colSums(data.storms[data.grid$REEF_ID==fid[i],-(1:2)]) > 0], col="lightcoral", lwd=2)
  
  # Add individual grid cells
  # data.reef.i <- res.med[data.grid$REEF_ID==fid[i],]
  # if (!is.vector(data.reef.i)) {
  #   for (j in 1: dim(data.reef.i)[1]) {lines(1996:years[k], data.reef.i[j,1:k], col=rgb(palette[[i]](.5),max=255), lty=2, lwd=.5)}
  # } else {lines(1996:years[k], data.reef.i[1:k], col=rgb(palette[[i]](.5),max=255), lty=2, lwd=.5)}
  
  # Add confidence interval and median at reef level
  polygon (c(1996:years[k], years[k]:1996), c(res.manta.min[res.manta.fid==fid[i],1:k], res.manta.max[res.manta.fid==fid[i],k:1]), border = NA,col=rgb(blue(.4),max=255), fillOddEven=T)
  polygon (c(1996:years[k], years[k]:1996), c(res.manta.25[res.manta.fid==fid[i],1:k], res.manta.75[res.manta.fid==fid[i],k:1]), border = NA,col=rgb(blue(.6),max=255), fillOddEven=T)
  lines(1996:years[k], res.manta.med[res.manta.fid==fid[i],1:k], col=rgb(blue(.9),max=255), lwd=3)
  title(main= paste(as.character(c("A","B","C")[i]), ") ", unique(obs.manta$REEF_NAME[obs.manta$REEF_ID == fid[i]]), sep=""))
  
  # Add Manta observations
  manta.i <- subset(obs.manta, REEF_ID == fid[i] & REPORT_YEAR <= years[k], select=c(REPORT_YEAR, MEAN_LIVE.corr))
  manta.i <- manta.i[order(manta.i$REPORT_YEAR),]
  lines(manta.i$REPORT_YEAR, manta.i$MEAN_LIVE.corr)
  points(manta.i$REPORT_YEAR, manta.i$MEAN_LIVE.corr, pch=19, cex=1.5, col="black")
}


# GBR plot
KolT=brewer.pal(11, "Spectral")[c(11,10,8,5,3,1)]

dev.new(width=10, height=5)
par(mai=c(.5,1,.5,1)) 
k=22
plot(NA, xlim=c(1996,2017), ylim=c(0,60), xlab="Time (yr)", ylab="Coral cover (%)", xaxs="i", yaxs="i", las=1)
polygon (c(1996:years[k], years[k]:1996), c(res.gbr.min[1:k], res.gbr.max[k:1]), border = NA, col=addalpha("lightcyan3",.5), fillOddEven=T)
polygon (c(1996:years[k], years[k]:1996), c(res.gbr.25[1:k], res.gbr.75[k:1]), border = NA, col=addalpha("lightcyan4",.5), fillOddEven=T)
lines(1996:years[k], res.gbr.med[1:k], col="gray25", lwd=4, lty=3)
lines(1996:years[k], res.gbr.mn[1:k], col="gray25", lwd=4)
for (j in 1: dim(res.cluster.med)[1]) {lines(1996:years[k], res.cluster.med[j,1:k], col=KolT[j], lwd=2)}
title(main= "GREAT BARRIER REEF", cex.main=1.5)
#dev.off()

#for (j in 1: dim(res.cluster.med)[1]) {print((res.cluster.mn[j,20]-res.cluster.mn[j,1])/20)}


# Map uncertainty based on CV(%) in model predictions across simulations -------------------

res.mn <- apply(res, c(1,2), mean, na.rm=T) # mean HC over time in each grid cell
res.sd <- apply(res, c(1,2), sd, na.rm=T) # std dev in HC over time in each grid cell (among simulations)
res.cv <- res.sd / res.mn
res.cv.avg <- rowMeans(res.cv)
res.sd.avg <- rowMeans(res.sd)

dev.new()
par(mai=c(1,1,.2,.2), family="Calibri")
hist(res.cv.avg, breaks=100, col = "black", cex.lab=4, xlab="CV (%)", main = NULL, xaxt="n", yaxt="n")
axis(1, at=c(0,.8), cex.axis=3)
axis(2, at=c(0,400), cex.axis=3)

insert <- readPNG("Input/Uncertainty_Inset.png")


#Reef level
res.reef.mn <- apply(res.reef, c(1,2), mean, na.rm=T)
res.reef.sd <- apply(res.reef, c(1,2), sd, na.rm=T)
res.reef.cv <- res.reef.sd/res.reef.mn
res.reef.cv.avg <- rowMeans(res.reef.cv)
res.reef.sd.avg <- rowMeans(res.reef.sd)
hist(exp(res.reef.cv.avg+1), breaks=100)

# res.reef.fid <- aggregate(data.grid$LONG, by=list(data.grid$REEF_ID), FUN="mean")[,1]
# res.reef.cv.gd <- rep(0, dim(data.grid)[1])
# for (i in 1:dim(data.grid)[1])  res.reef.cv.gd[i] <- res.reef.cv.avg[res.reef.fid==data.grid$REEF_ID[i]]

par(mai=c(.8,.8,.2,.2), family="Calibri")
plotPolys(shape, plt=NULL, col="gray95", border="gray85", xlim=c(143,155), ylim=c(-25,-14), bg="white", cex=1.2, xlab="Longitude", ylab="Latitude", lwd=.5)
#points(data.grid$LONG, data.grid$LAT, pch=19, cex=.2, col=rgb(yellow(scale(res.sd.avg)),max=255))
points(data.grid$LONG, data.grid$LAT, pch=19, cex=.2, col=rgb(yellow(scale(exp(res.cv.avg))),max=255))

addPolys(extrapol, col="lightgrey", border="lightgrey", lwd=2, density=8, colHoles=NA)
points(obs.calib$REEF_LONG, obs.calib$REEF_LAT, pch=21, cex=1)

rasterImage(insert, 144, -24, 148, -21)
box()

#Legend
x=3; y=1
for (i in 1:100) {segments(x0 = 148+x, y0 = -16.5 +y - i/100, x1 = 148.5+x, col=rgb(yellow(1- i/100),max=255))}
text(148.7+x, -16.5+y, labels = round(max(res.cv.avg),2)*100, cex=1.5, pos=4, font = 1.5)
text(148.7+x, -17.5+y, labels = round(min(res.cv.avg),2)*100, cex=1.5, pos=4, font = 1.5)
text(146.5+x, -15.9+y, labels = "Uncertainty (CV; %)", cex=1.5, pos=4)


# Map disturbance vs. decline + uncertainty -------------------------

# Calculate severity of disturbance events (cumul severity * effect size)
cnt.bleaching <- bleaching.mn * abs(bleaching.mn.sd[1] + WQ*WQ_bleach[1]) 
cnt.bleaching[is.na(cnt.bleaching)] <- 0
cnt.bleaching <- data.frame(cnt.bleaching, cnt = rowSums(cnt.bleaching, na.rm=T))

cnt.disease <- disease.mn * abs(disease.mn.sd[1] + WQ*WQ_Disease[1]) 
cnt.disease <- data.frame(cnt.disease, cnt = rowSums(cnt.disease, na.rm=T))

cnt.COTS <- COTS.mn * abs(COTS.mn.sd[1] + WQ*WQ_CoTS[1])
cnt.COTS[is.na(cnt.COTS)] <- 0
cnt.COTS <- data.frame(cnt.COTS, cnt = rowSums(cnt.COTS))

cnt.storms <- storms.mn * ifelse((storms.mn.sd[1] + WQ*WQ_Cyclone[1]) > 0, 0, abs(storms.mn.sd[1] + WQ*WQ_Cyclone[1]))
cnt.storms[is.na(cnt.storms)] <- 0
cnt.storms <- data.frame(cnt.storms, cnt = rowSums(cnt.storms))

cnt.unknown <- unknown.mn * abs(unknown.mn.sd[1])
cnt.unknown <- data.frame(cnt.unknown, cnt = rowSums(cnt.unknown))

disturb <- cnt.COTS$cnt + cnt.bleaching$cnt + cnt.storms$cnt + cnt.unknown$cnt + cnt.disease$cnt
disturb.reef <- aggregate(disturb, by=list(data.grid$REEF_ID), FUN="mean")[,2]
decline.reef <- res.reef.mn[,1] - res.reef.mn[,22]

res.reef.fid <- aggregate(data.grid$LONG, by=list(data.grid$REEF_ID), FUN="mean")[,1]
disturb.gd <- decline.gd <- rep(0, dim(data.grid)[1])
for (i in 1:dim(data.grid)[1])  {
  disturb.gd[i] <- disturb.reef[res.reef.fid==data.grid$REEF_ID[i]]
  decline.gd[i] <- decline.reef[res.reef.fid==data.grid$REEF_ID[i]]
}


decline.bin <- cut(decline.gd, quantile(decline.gd, probs=(0:2)/2), labels=F, include.lowest = T, na.rm=T)
disturb.bin <- cut(disturb.gd, quantile(disturb.gd, probs=(0:2)/2, na.rm=T), labels=F, include.lowest = T, na.rm=T)

decdis <- decdis.kol <- ""
decdis[decline.bin ==1 & disturb.bin ==1] <- "Lo_disturb_Lo_decl" 
decdis[decline.bin ==1 & disturb.bin ==2] <- "Hi_disturb_Lo_decl" 
decdis[decline.bin ==2 & disturb.bin ==1] <- "Lo_disturb_Hi_decl" 
decdis[decline.bin ==2 & disturb.bin ==2] <- "Hi_disturb_Hi_decl" 

decdis.kol[decline.bin ==1 & disturb.bin ==1] <-  "mediumseagreen"
decdis.kol[decline.bin ==1 & disturb.bin ==2] <- "darkseagreen1" 
decdis.kol[decline.bin ==2 & disturb.bin ==1] <-  "steelblue3"
decdis.kol[decline.bin ==2 & disturb.bin ==2] <- "lightskyblue2" 


extrapol_X5 <- extrapol
extrapol_X5$X <- extrapol$X + 5

dev.off()
par(mai=c(.8,.8,.2,.2), family="Calibri", cex.lab=2, cex.sub=2, cex.axis=1.5)
plotPolys(shape, col="lightgrey", border="darkgrey", xlim=c(144,158), ylim=c(-25,-14), bg="white", xlab=list("Longitude", cex=2, font=2), ylab = list("Latitude", cex=2, font=2), xaxt = "n", yaxt = "n", colHoles=NA)
axis(1, at=seq(142,152, by=2))
axis(2, at=seq(-24, -14, by=2))

points(data.grid$LONG, data.grid$LAT, pch=19, cex=.3, col=decdis.kol)
addPolys(extrapol, col=addalpha("white",.4), border="white", lwd=2, colHoles=NA)
addPolys(extrapol, col="darkgrey", border="darkgrey", lwd=2, density=8, colHoles=NA)

points(data.grid$LONG+5, data.grid$LAT, pch=19, cex=.3, col=rgb(yellow(scale(exp(res.cv.avg))),max=255))
addPolys(extrapol_X5, col=addalpha("white",.4), border="white", lwd=2, colHoles=NA)
addPolys(extrapol_X5, col="darkgrey", border="darkgrey", lwd=2, density=8, colHoles=NA)


#scalebar(loc=c(145, -24), length=2, unit="km", cex=2)
#northarrow(loc=c(146,-22), size=.5)

#Legend
rect(145,-24,145.5,-23.5,col="darkseagreen1")
rect(145,-23.4,145.5,-22.9,col="steelblue3")
rect(145.6,-24,146.1,-23.5,col="mediumseagreen")
rect(145.6,-23.4,146.1,-22.9,col="lightskyblue2")

text(145.25,-24.2, "L")
text(145.85,-24.2, "H")
text(144.8,-23.75, "L")
text(144.8,-23.15, "H")
text(145.55,-24.5, "Disturbance", cex=1.2)
text(144.4,-23.55, "Decline",srt=90, cex=1.2)

for (i in 1:100) {segments(x0 = 154.8, y0 = -16 - i/100, x1 = 155.5, col=rgb(yellow(1- i/100),max=255))}
text(156.5, -16, labels = round(max(res.cv.avg)*100), cex=1.2, pos=2, font = 1.5)
text(156.5, -17, labels = "0", cex=1.2, pos=2, font = 1.5)
text(155.8, -15.5, labels = "Uncertainty (CV%)", cex=1.2, font=2)

# Resilience ordination and correlates --------------

# Remove grid cells in extrapolated areas
data.grid.EID <- as.EventData(data.frame(EID=data.grid$PIXEL_ID, X = data.grid$LONG, Y= data.grid$LAT))
data.grid.EID.extrapol <- findPolys(data.grid.EID, extrapol)

# PCA and definition of resilience
pca <- PCA(data.frame(decdis, disturb.gd, decline.gd)[!(data.grid$PIXEL_ID %in% data.grid.EID.extrapol[,1]) & decdis %in% c("Lo_disturb_Hi_decl", "Hi_disturb_Lo_decl"),], quali.sup=1)
plot.PCA(pca, choix="ind", habillage=1, label="none", col.hab=c("lightskyblue2","mediumseagreen","steelblue3","darkseagreen1"))
resilience <- -1*pca$ind$coord[,1]         
resilience[resilience > quantile(resilience, probs=.995)] <- quantile(resilience, probs=.995)

# Relationship with correlates
#resilience.sub <- resilience[decdis %in% c("Lo_disturb_Hi_decl", "Hi_disturb_Lo_decl")]
resilience.sub <- resilience
WQ.sub <- WQ[!(data.grid$PIXEL_ID %in% data.grid.EID.extrapol[,1]) & decdis %in% c("Lo_disturb_Hi_decl", "Hi_disturb_Lo_decl")]
green.sub <- data.grid$Green[!(data.grid$PIXEL_ID %in% data.grid.EID.extrapol[,1]) & decdis %in% c("Lo_disturb_Hi_decl", "Hi_disturb_Lo_decl")]
gravity.sub <- data.grid$Gravity[!(data.grid$PIXEL_ID %in% data.grid.EID.extrapol[,1]) & decdis %in% c("Lo_disturb_Hi_decl", "Hi_disturb_Lo_decl")]
gravity.sub[gravity.sub > 0.015] <- NA # remove outliers
gravity.bin <- ifelse(gravity.sub < 0.005, 0, cut(gravity.sub[gravity.sub>0.005], quantile(gravity.sub[gravity.sub>0.005], probs=(0:5)/5, na.rm=T), labels=F, include.lowest = T))


boxplot(resilience.sub ~ green.sub)
vioplot(resilience.sub[green.sub==0], resilience.sub[green.sub==1], col=c("darkseagreen","darkseagreen1"), names=c("open","closed"))

boxplot(resilience.sub ~ gravity.bin)
vioplot(resilience.sub[gravity.bin==0],
        resilience.sub[gravity.bin==1],
        resilience.sub[gravity.bin==2],
        resilience.sub[gravity.bin==3],
        resilience.sub[gravity.bin==4],
        resilience.sub[gravity.bin==5],
        col=brewer.pal(9,"YlGnBu")[c(2,7,6,5,4,3)],
        names=c("<50","[50;55[","[55;60[","[60-80[","[80-125[","[126-148["))
   
                                
plot(resilience.sub ~ WQ.sub, pch=19)
gam <- gam(resilience.sub ~ s(WQ.sub, k=3, bs="cs"))
plot(gam)

new.data <- data.frame(WQ.sub=seq(0,1,by=0.01))
pred <- predict.gam(gam, new.data, se.fit=T)
pred.mn <- pred$fit
pred.hi <- pred$fit + 1.96*pred$se.fit
pred.lo <- pred$fit - 1.96*pred$se.fit

par(mai=c(.8,.8,.2,.2), family="Calibri")
plot(resilience.sub ~ WQ.sub, type="n", xlab="PFc", ylab="Resilience", xaxs="i", yaxs="i", xlim=c(0,1), ylim=c(-3,3))
polygon(c(new.data[,1], rev(new.data[,1])), c(pred.hi, rev(pred.lo)), col="lightgrey", border=NA)
lines(new.data$WQ.sub, pred$fit, col="darkgrey", lwd=2)
points(WQ.sub, resilience.sub, pch=16, col=rgb(0, 0, 1, 0.25))
rug(WQ.sub)

gravity.log <- log10(gravity.sub+1)*100
plot(resilience.sub ~ gravity.log, pch=19)
gam <- gam(resilience.sub ~ s(gravity.log, k=3, bs="cs"))
plot(gam)

new.data <- data.frame(gravity.log=seq(0,0.6, by=0.01))
pred <- predict.gam(gam, new.data, se.fit=T)
pred.mn <- pred$fit
pred.hi <- pred$fit + 1.96*pred$se.fit
pred.lo <- pred$fit - 1.96*pred$se.fit

par(mai=c(.8,.8,.2,.2), family="Calibri")
plot(resilience.sub ~ gravity.log, type="n", xlab="Human density", ylab="Resilience", xaxs="i", yaxs="i", xlim=c(0,0.5), ylim=c(-3,3))
polygon(c(new.data[,1], rev(new.data[,1])), c(pred.hi, rev(pred.lo)), col="lightgrey", border=NA)
lines(new.data$gravity.log, pred$fit, col="darkgrey", lwd=2)
points(gravity.log, resilience.sub, pch=16, col=rgb(0, 0, 1, 0.25))
rug(gravity.log)


## MPA -- suppl violin plots
par(mfrow=c(2,2), mai=c(.5,.5,.5,.5))
vioplot(COTS.sp[data.grid$Green==0], COTS.sp[data.grid$Green==1], 
        col=c("seashell1","seashell3"), names=c("Open","Closed"))
title(main="COTS", xlab="Reef zoning")
means.COTS <- tapply(COTS.sp,data.grid$Green,mean)
points(means.COTS,col="red", pch=19, cex=3)
grid()

vioplot(STORMS.sp[data.grid$Green==0], STORMS.sp[data.grid$Green==1], 
        col=c("seashell1","seashell3"), names=c("Open","Closed"))
title(main="Cyclones", xlab="Reef zoning")
means.storms <- tapply(STORMS.sp,data.grid$Green,mean)
points(means.storms,col="red", pch=19, cex=3)
grid()

vioplot(BLEACH.sp[data.grid$Green==0], BLEACH.sp[data.grid$Green==1], 
        col=c("seashell1","seashell3"), names=c("Open","Closed"))
title(main="Bleaching", xlab="Reef zoning")
means.bleach <- tapply(BLEACH.sp,data.grid$Green,mean)
points(means.bleach,col="red", pch=19, cex=3)
grid()

vioplot(WQ.sub[green.sub==0], WQ.sub[green.sub==1], 
        col=c("seashell1","seashell3"), names=c("Open","Closed"))
title(main="PFc", xlab="Reef zoning")
means.WQ <- tapply(WQ.sub,green.sub,mean)
points(means.WQ,col="red", pch=19, cex=3)
grid()


## Boxplots
boxplot(BLEACH.sp ~ data.grid$Green, col=c("seashell1","seashell3"), names=c("Open","Closed"), main="Bleaching")
kruskal.test(BLEACH.sp ~ data.grid$Green)
means.bleach <- tapply(BLEACH.sp,data.grid$Green,mean)
points(means.bleach,col="red", pch=19, cex=3)

boxplot(COTS.sp ~ data.grid$Green, col = "lightgray", main="CoTS")
kruskal.test(COTS.sp ~ data.grid$Green)
means.COTS <- tapply(COTS.sp,data.grid$Green,mean)
points(means.COTS,col="red", pch=19, cex=3)

boxplot(STORMS.sp ~ data.grid$Green, col = "lightgray", main="Cyclones")
kruskal.test(STORMS.sp ~ data.grid$Green)
means.STORMS <- tapply(STORMS.sp,data.grid$Green,mean)
points(means.STORMS,col="red", pch=19, cex=3)

boxplot(WQ.sub ~ green.sub, col = "lightgray", main="PFc")
kruskal.test(WQ.sub ~ green.sub)
means.WQ <- tapply(WQ.sub,green.sub,mean)
points(means.WQ,col="red", pch=19, cex=3)







# Sensitivity analysis ----------------

# GBR-wide trajectories across simulations
res.gbr <- t(apply(res, c(2,3), mean))

plot(1996:2017, colMeans(res.gbr), type="l", ylim=c(10,40))
for (i in 1:100) lines(1996:2017, res.gbr[i,], col="grey")
lines(1996:2017, colMeans(res.gbr))

decline <- res.gbr[,1] - res.gbr[,22]
hist(decline, breaks=20)


#Collate data for BRT
BLEACH <- bleaching.rsmpl * abs(bleaching.mn.sd[1] + WQ*WQ_bleach[1]) 
BLEACH[is.na(BLEACH)] <- 0
BLEACH.sum <- apply(BLEACH, 3, mean, na.rm=T)

COTS <- COTS.rsmpl * abs(COTS.mn.sd[1] + WQ*WQ_CoTS[1])
COTS[is.na(COTS)] <- 0
COTS.sum <- apply(COTS, 3, mean, na.rm=T)

STORMS <- storms.rsmpl * ifelse((storms.mn.sd[1] + WQ*WQ_Cyclone[1]) > 0, 0, abs(storms.mn.sd[1] + WQ*WQ_Cyclone[1]))
STORMS[is.na(STORMS)] <- 0
STORMS.sum <- apply(STORMS, 3, mean, na.rm=T)

data.brt <- data.frame(decline, HCMAX = LHS[,7], B0 = LHS[,8], HCINI = LHS[,9],
                    bleach = scale(BLEACH.sum), COTS = scale(COTS.sum), storms=scale(STORMS.sum))

gbm.SA <- gbm.step(data = data.brt,
                   gbm.x = 2:7,
                   gbm.y = 1,
                   family = "gaussian",
                   tree.complexity = 3,
                   learning.rate = 0.001,
                   bag.fraction = 0.5)

gbm.plot(gbm.SA, n.plots=6, write.title=F, plot.layout = c(3,3))

plot(colMeans(B0), decline, pch=19, col="black", ylab="Coral decline (%)", xlab="B0")

summary(gbm.SA)

#find.int <- gbm.interactions(gbm.SA)
#find.int$rank.list

#Interaction between HCini and B0, size=1.70
gbm.perspec(gbm.SA,2,3, z.range=c(12,16))


# FIG 1 -------------------

# Calculate combined effect of cyclones, bleahcing, COTS
BLEACH <- bleaching.rsmpl * abs(bleaching.mn.sd[1] + WQ*WQ_bleach[1]) 
BLEACH[is.na(BLEACH)] <- 0
BLEACH.sp <- rowSums(apply(BLEACH, c(1,2), mean, na.rm=T))/3
BLEACH.sp[BLEACH.sp > quantile(BLEACH.sp,.99)] <- quantile(BLEACH.sp,.99)

COTS <- COTS.rsmpl * abs(COTS.mn.sd[1] + WQ*WQ_CoTS[1])
COTS[is.na(COTS)] <- 0
COTS.sp <- apply(COTS, 1, mean, na.rm=T)
COTS.sp[COTS.sp > quantile(COTS.sp,.99)] <- quantile(COTS.sp,.99)

storms.mn.sd <- c(-0.2832, 0.01173)
STORMS <- storms.rsmpl * ifelse((storms.mn.sd[1] + WQ*WQ_Cyclone[1]) > 0, 0, abs(storms.mn.sd[1] + WQ*WQ_Cyclone[1]))
STORMS[is.na(STORMS)] <- 0
STORMS.sp <- apply(STORMS, 1, mean, na.rm=T)
STORMS.sp[STORMS.sp > quantile(STORMS.sp,.99)] <- quantile(STORMS.sp,.99)

# Tabulate combined effect per year
COTS.tp <- apply(COTS, 2, mean, na.rm=T)
COTS.tp[COTS.tp > quantile(COTS.tp,.99)] <- quantile(COTS.tp,.99)

STORMS.tp <- apply(STORMS, 2, mean, na.rm=T)
STORMS.tp[STORMS.tp > quantile(STORMS.tp,.99)] <- quantile(STORMS.tp,.99)

BLEACH.tp <- apply(BLEACH, 2, mean, na.rm=T)
BLEACH.tp[BLEACH.tp > quantile(BLEACH.tp,.99)] <- quantile(BLEACH.tp,.99)


# Map combined effects
orange=colorRamp(brewer.pal(9, "YlOrBr")[2:8], interpolate="linear", space="rgb", bias=.8)
purple=colorRamp(brewer.pal(9, "BuPu")[2:7], interpolate="linear", space="rgb", bias=.8)
azures=colorRamp(c("lightcyan1","lightcyan2","lightcyan3","lightcyan4","lightskyblue4"), interpolate="linear", space="rgb", bias=.8)
turquoise=colorRamp(brewer.pal(9, "YlGn")[3:5], interpolate="linear", space="rgb", bias=.8)

data.map.pos$z[data.map.pos$z > quantile(data.map.pos$z,.9)] <- quantile(data.map.pos$z,.9)
data.map.neg$z[data.map.neg$z < quantile(data.map.neg$z,.05)] <- quantile(data.map.neg$z,.05)

extrapol_X17 <- extrapol
extrapol_X17$X <- extrapol$X + 17

insert <- readPNG("Input/Australia_Inset.png")

#FIG1
{
dev.off()
par(mai=c(.5,.6,.3,.3), omi=c(0,0,0,0), family="Calibri", cex=1, cex.lab=1.5, cex.sub=1.5, cex.axis=1.5)
layout.show(layout(rbind(1,2,3), widths= c(1,1), heights = c(2,1,1)))

plotPolys(shape, col="lightgrey", border="darkgrey", xlim=c(144,170), ylim=c(-25,-14), bg="white",  xlab=NA, ylab=NA, xaxt = "n", yaxt = "n", las=1, colHoles=NA)
axis(1, at=seq(142,152, by=2))
axis(2, at=seq(-24, -14, by=2))
title(xlab="Longitude", line=4, cex.lab=1.5)
title(ylab="Latitude", line=4, cex.lab=1.5)

points(data.grid$LONG, data.grid$LAT, pch=19, cex=.3, col=rgb(orange(scale(STORMS.sp)),max=255))
points(data.grid$LONG+5, data.grid$LAT, pch=19, cex=.3, col=rgb(orange(scale(COTS.sp)),max=255))
points(data.grid$LONG+10, data.grid$LAT, pch=19, cex=.3, col=rgb(orange(scale(BLEACH.sp)),max=255))

points(data.map.neg$x+17, data.map.neg$y, pch=19, cex=.2, col=rgb(azures(scale(-1*data.map.neg$z)),max=255))
points(data.map.pos$x+17, data.map.pos$y, pch=19, cex=.2, col=rgb(turquoise(scale(data.map.pos$z)),max=255))

addPolys(extrapol_X17, col=addalpha("white",.4), border="white", lwd=2, colHoles=NA)
addPolys(extrapol_X17, col="lightgrey", border="lightgrey", lwd=2, density=8, colHoles=NA)

#insert
rasterImage(insert, 144.1, -24.8, 147.8, -21.8)
box()


#scalebar(loc=c(145, -24), length=2, unit="km", cex=2)
#northarrow(loc=c(146,-22), size=.5)

#Legend
text(146.2, -15, "A) Cyclones", cex=1.5, font=2, pos=4)
text(146.2 +5, -15, "B) CoTS", cex=1.5, font=2, pos=4)
text(146.2 +10, -15, "C) Bleaching", cex=1.5, font=2, pos=4)
text(146.2 +17, -15, "D) Change in coral cover", cex=1.5, font=2, pos=4)

# x=0
# for (i in 1:100) {segments(x0 = 147.75+x, y0 = -16.5 - i/100, x1 = 148.25+x, col=rgb(orange(1- i/100),max=255))}
# text(148.2+x, -16.5, labels = round(max(STORMS.sp),1), cex=1.5, pos=4, font = 1.5)
# text(148.2+x, -17.5, labels = "0", cex=1.5, pos=4, font = 1.5)
# text(148+x, -15.8, labels = "Effect size", cex=1.5, font=2)

# x=5
# for (i in 1:100) {segments(x0 = 147.75+x, y0 = -16.5 - i/100, x1 = 148.25+x, col=rgb(orange(1- i/100),max=255))}
# text(148.2+x, -16.5, labels = round(max(COTS.sp),1), cex=1.5, pos=4, font = 1.5)
# text(148.2+x, -17.5, labels = "0", cex=1.5, pos=4, font = 1.5)
# text(148+x, -15.8, labels = "Effect size", cex=1.5, font=2)
# 
x=12
for (i in 1:100) {segments(x0 = 147.75+x, y0 = -16.5 - i/100, x1 = 148.25+x, col=rgb(orange(1- i/100),max=255))}
text(148.2+x, -16.5, labels = round(max(BLEACH.sp),1), cex=1.5, pos=4, font = 1.5)
text(148.2+x, -17.5, labels = "0", cex=1.5, pos=4, font = 1.5)
text(148+x, -15.8, labels = "Impact", cex=1.5, font=2)

x=19
for (i in 0:100) {segments(x0 = 147.75+x, y0 = -17 - i/100, x1 = 148.25+x, col=rgb(azures(i/100),max=255))}
for (i in 0:50) {segments(x0 = 147.75+x, y0 = -17 + i/100, x1 = 148.25+x, col=rgb(turquoise(i/50),max=255))}
text(148.2+x, -16.5, pos=4, labels = paste("+",round(max(delta, na.rm=T)/22,1),sep=""), cex=1.5)
text(148.2+x, -18, pos=4, labels = round(min(delta, na.rm=T)/22,1), cex=1.5)
text(148.2+x, -17, pos=4, labels="0", cex=1.5)
text(148+x, -15.8, labels = "(% y-1)", cex=1.5, font=2)

points(145.78, -16.92, pch=17, cex=1.5)
points(146.82, -19.26, pch=17, cex=1.5)
points(149.18, -21.14, pch=17, cex=1.5)
points(150.51, -23.28, pch=17, cex=1.5)

text(145.78, -16.92, "Cairns", cex=1.2, pos=2)
text(146.82, -19.26, "Townsville", cex=1.2, pos=2)
text(149.18, -21.14, "Mackay", cex=1.2, pos=2)
text(150.51, -23.28, "Rockhampton", cex=1.2, pos=2)



# Bar plot and time series
data.plot <- data.frame(STORMS.tp, COTS.tp, BLEACH.tp)
data.plot[22,] <- 0
kol <- c(rgb(orange(.75),max=255), rgb(orange(.5),max=255), "khaki1")

par(mai=c(.3,.7,.5,.3))
barplot(t(data.plot), ylab=NA, 
        col=kol, border="darkgrey", space=0.5, las=1, xaxs="i", yaxs="i", ylim=c(0,.5))
text(32.5,.1,"NA", cex=1.5)
title(ylab="Impact", line=6, cex.lab=1.5)

#legend(4, .85, c("Cyclones","CoTS","Bleaching"), cex=1.2, fill=kol)
legend(13, .45, c("Cyclones","CoTS","Bleaching"), cex=1.2, fill=kol)
text(1.25,.45,"E)", cex=1.5, font=2)
box()

k=22
par(mai=c(.6,.7,.2,.3))
plot(NA, xlim=c(1996,2017), ylim=c(0,50), xlab="Time (yr)", ylab=NA, xaxs="i", yaxs="i", las=1)
polygon (c(1996:years[k], years[k]:1996), c(res.gbr.min[1:k], res.gbr.max[k:1]), border = NA, col=addalpha("lightcyan3",.5), fillOddEven=T)
polygon (c(1996:years[k], years[k]:1996), c(res.gbr.25[1:k], res.gbr.75[k:1]), border = NA, col=addalpha("lightcyan4",.5), fillOddEven=T)
lines(1996:years[k], res.gbr.mn[1:k], col="gray25", lwd=4)
title(ylab="Coral cover (%)", line=6, cex.lab=1.5)
text(1996.5,45,"F)", cex=1.5, font=2)
}

# FIG 3 ----------

dev.off()
par(mai=c(.5,.6,.3,.3), omi=c(0,0,0,0), family="Calibri", cex=1, cex.lab=1.5, cex.sub=1.5, cex.axis=1.5)
layout.show(layout(rbind(c(1,2),c(1,3),c(1,4)), widths= c(1,1), heights = c(1,1,1)))

par(mai=c(0,0,0,0))
plotPolys(shape45, col="gray95", border="gray70", xlim=c(151,156), ylim=c(-36, -22), bg="white", cex=1.2, axes=F, ylab=NA, plt=c(0,1,0,1), colHoles=NA)
points(pst45$x, pst45$y, pch=19, cex=1, col=rgb(greys(scale(data.WQ$PST)),max=255))
lines(gridln.elide, col="lightgrey")
addPolys(shape45, col="gray95", border="gray70", bg="white", colHoles=NA)
points(grid45$x, grid45$y, pch=19, cex=1, col=decdis.kol)
addPolys(extrapol45, col=addalpha("white",.4), border="white", lwd=2, colHoles=NA)
addPolys(extrapol45, col="lightgrey", border="lightgrey", lwd=2, density=8, colHoles=NA)

northarrow(loc=c(152,-24), size=.6, cex=2, bearing=-pi/4)
text(155.7,-22.5,"A)", cex=2, font=2)

text(155.3, -22.5, "150E", col="lightgrey", srt=40, cex=1.5)
text(155.3, -29.4, "155E", col="lightgrey", srt=40, cex=1.5)
text(155.3, -26, "15S", col="lightgrey", srt=-40, cex=1.5)
text(155.3, -33.2, "20S", col="lightgrey", srt=-40, cex=1.5)

#Legend
x0 <- 155.7
y0 <- -35.8

arrows(x0-1.2, y0+1.7, x0+.1, y0+.4, length=.1, lwd=2)

rect(x0-1.1, y0+.5, x0-.6, y0+1, col="darkseagreen1")
rect(x0-1.1, y0+1.1, x0-.6, y0+1.6, col="steelblue3")
rect(x0-.5, y0+.5, x0, y0+1, col="mediumseagreen")
rect(x0-.5, y0+1.1, x0, y0+1.6, col="lightskyblue2")

text(x0-.85, y0+.3, "L", cex=1.5)
text(x0-.25, y0+.3, "H", cex=1.5)
text(x0-1.3, y0+.75, "L", cex=1.5)
text(x0-1.3, y0+1.35, "H", cex=1.5)
text(x0-.55, y0, "Disturbance", cex=1.5, font=2)
text(x0-1.5, y0+.95, "Decline",srt=90, cex=1.5, font=2)
text(x0+.2, y0+.3, "R", cex=2, font=4)


# WQ plot
gam <- gam(resilience.sub ~ s(WQ.sub, k=3, bs="cs"))
new.data <- data.frame(WQ.sub=seq(0,1,by=0.01))
pred <- predict.gam(gam, new.data, se.fit=T)
pred.mn <- pred$fit
pred.hi <- pred$fit + 1.96*pred$se.fit
pred.lo <- pred$fit - 1.96*pred$se.fit

par(mai=c(.6,.6,.1,.1))
plot(resilience.sub ~ WQ.sub, type="n", xlab="PFc", ylab="Resilience (R)", xaxs="i", yaxs="i", xlim=c(0,1), ylim=c(-2,3))
polygon(c(new.data[,1], rev(new.data[,1])), c(pred.hi, rev(pred.lo)), col="lightgrey", border=NA)
points(WQ.sub, resilience.sub, pch=16, col=addalpha("seashell2",.5))
lines(new.data$WQ.sub, pred$fit, col="darkgrey", lwd=2)
rug(WQ.sub)
text(.9,2.5, "B)", cex=2, font=2)
box()

# Gravity plot
gam <- gam(resilience.sub ~ s(gravity.log, k=3, bs="cs"))
new.data <- data.frame(gravity.log=seq(0,0.6, by=0.01))
pred <- predict.gam(gam, new.data, se.fit=T)
pred.mn <- pred$fit
pred.hi <- pred$fit + 1.96*pred$se.fit
pred.lo <- pred$fit - 1.96*pred$se.fit

par(mai=c(.6,.6,.1,.1))
plot(resilience.sub ~ gravity.log, type="n", xlab="Human density", ylab="Resilience (R)", xaxs="i", yaxs="i", xlim=c(0,0.5), ylim=c(-2,3))
polygon(c(new.data[,1], rev(new.data[,1])), c(pred.hi, rev(pred.lo)), col="lightgrey", border=NA)
points(gravity.log, resilience.sub, pch=16, col=addalpha("seashell2",.5))
lines(new.data$gravity.log, pred$fit, col="darkgrey", lwd=2)
rug(gravity.log)
text(.45,2.5, "C)", cex=2, font=2)
box()

# Green zone plot
par(mai=c(.6,.6,.1,.1))
vioplot(resilience.sub[green.sub==0], resilience.sub[green.sub==1], 
        col=c("seashell1","seashell3"), names=c("Open","Closed"))
title(ylab="Resilience (R)", xlab="Reef zoning")
text(2.4,2.5, "D)", cex=2, font=2)



# Gridded plots of disturbance severity (FIG S1) -------------------

storm.mat <- as.matrix(data.storms.bckp[,-(1:2)])
storm.mat[storm.mat>quantile(storm.mat,.98)] <- quantile(storm.mat,.98)

cots.mat <- as.matrix(data.COTS.bckp[,-(1:2)])
cots.mat[cots.mat>quantile(cots.mat,.95)] <- quantile(cots.mat,.95)

bleach.mat <- as.matrix(data.bleaching.bckp[,-(1:2)])
#bleach.mat[bleach.mat>quantile(bleach.mat,.95)] <- quantile(bleach.mat,.95)

dev.off()
par(mai=c(0,0,0,0), omi=c(0,0,0,0), family="Calibri", cex=1, cex.lab=1.5, cex.sub=1.5, cex.axis=1.5)
layout.show(layout(rbind(1,2,3)))

## STORMS
plotPolys(shape, col="lightgrey", border="darkgrey", xlim=c(144,258), ylim=c(-25,-14), bg="white",  xlab=NA, ylab=NA, xaxt = "n", yaxt = "n", las=1, colHoles=NA)
points(data.grid$LONG, data.grid$LAT, pch=19, cex=.3, col=rgb(orange(scale.global(log10(storm.mat[,1]+2),log10(storm.mat+2))),max=255))
text(152,-24.5,"1996", cex=1.5, font=2)
text(145, -21, "Cyclones", cex=2, font=2, srt=90)

for (i in 2:21) {
  points(data.grid$LONG+5*(i-1), data.grid$LAT, pch=19, cex=.3, col=rgb(orange(scale.global(log10(storm.mat[,i]+2),log10(storm.mat+2))),max=255))
  text(152+5*(i-1),-24.5,(1996:2017)[i], cex=1.5, font=2)
}

x=105
for (i in 1:100) {segments(x0 = 146+x, y0 = -16.5 - i/100, x1 = 148+x, col=rgb(orange(1- i/100),max=255))}
text(148.2+x, -16.5, labels = ">17", cex=1.5, pos=4, font = 1.5)
text(148.2+x, -17.5, labels = "0", cex=1.5, pos=4, font = 1.5)
text(148+x, -15.8, labels = "Severity", cex=1.5, font=2)


## COTS
plotPolys(shape, col="lightgrey", border="darkgrey", xlim=c(144,258), ylim=c(-25,-14), bg="white",  xlab=NA, ylab=NA, xaxt = "n", yaxt = "n", las=1, colHoles=NA)
points(data.grid$LONG, data.grid$LAT, pch=19, cex=.3, col=rgb(orange(scale.global(log10(cots.mat[,1]+2),log10(cots.mat+2))),max=255))
text(152,-24.5,"1996", cex=1.5, font=2)
text(145, -21, "COTS", cex=2, font=2, srt=90)

for (i in 2:21) {
  points(data.grid$LONG+5*(i-1), data.grid$LAT, pch=19, cex=.3, col=rgb(orange(scale.global(log10(cots.mat[,i]+2),log10(cots.mat+2))),max=255))
  text(152+5*(i-1),-24.5,(1996:2017)[i], cex=1.5, font=2)
}

x=105
for (i in 1:100) {segments(x0 = 146+x, y0 = -16.5 - i/100, x1 = 148+x, col=rgb(orange(1- i/100),max=255))}
text(148.2+x, -16.5, labels = ">1", cex=1.5, pos=4, font = 1.5)
text(148.2+x, -17.5, labels = "0", cex=1.5, pos=4, font = 1.5)
text(148+x, -15.8, labels = "Severity", cex=1.5, font=2)

## BLEACHING
plotPolys(shape, col="lightgrey", border="darkgrey", xlim=c(144,258), ylim=c(-25,-14), bg="white",  xlab=NA, ylab=NA, xaxt = "n", yaxt = "n", las=1, colHoles=NA)
points(data.grid$LONG, data.grid$LAT, pch=19, cex=.3, col=rgb(orange(scale.global(bleach.mat[,1],bleach.mat)),max=255))
text(152,-24.5,"1996", cex=1.5, font=2)
text(145, -21, "Bleaching", cex=2, font=2, srt=90)

for (i in 2:21) {
  points(data.grid$LONG+5*(i-1), data.grid$LAT, pch=19, cex=.3, col=rgb(orange(scale.global(bleach.mat[,i],bleach.mat)),max=255))
  text(152+5*(i-1),-24.5,(1996:2017)[i], cex=1.5, font=2)
}

x=105
for (i in 1:100) {segments(x0 = 146+x, y0 = -16.5 - i/100, x1 = 148+x, col=rgb(orange(1- i/100),max=255))}
text(148.2+x, -16.5, labels = "80", cex=1.5, pos=4, font = 1.5)
text(148.2+x, -17.5, labels = "0", cex=1.5, pos=4, font = 1.5)
text(148+x, -15.8, labels = "Severity", cex=1.5, font=2)

