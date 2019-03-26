rm(list=ls())

## External Packages:
library("glmulti")    # glm: subset selection 
library("scales")     # colorscales
library("lubridate")  # working with times and dates
library("stargazer")  # LaTeX tables
library("rcompanion") # stats for glm
library("Hmisc")      # further stats for glm
library("devtools")   # install packages

## Set 'your_path' to 'GFLMPOI_blinded':
my_path <- "your_path/GFLMPOI_blinded/"

## Install the a ccompaning R-package GFLMPOI
devtools::install(pkg=paste0(my_path,"GFLMPOI"))

## Load the accompaning R-package GFLMPOI
library("GFLMPOI")    

## Attach the data
data("VideoRatingData")

## Variables: #################
Y     <- VideoRatingData$Y
X.mat <- t(as.matrix(VideoRatingData[,-1]))
##
(p     <- nrow(X.mat))
(N     <- ncol(X.mat))
##
a      <- 0
b      <- 1
t.grid <- (1:p-1)/(p-1)
## #############################


## Figure of rating trajectories
par(mar=c(4,5,4,1)+.1, cex.lab=1.3, cex=1.3)
matplot(t.grid, X.mat, type = "l", lty=1, lwd=1, col=alpha("black", 0.55), xlab="", ylab="")
mtext(side = 1, text = "Time (Standardized)", line = 3, at = 0.5, cex=1.3)
mtext(side = 3, text = "very\npositive", line = -1.25, at = -.1, cex=1.3)
mtext(side = 1, text = "very\nnegative", line = -.5,  at = -.1, cex=1.3)
par(mar=c(5,4,4,2)+.1)
dev.off()


## Estimation
result.obj  <- FUN_PoI_BIC(Y             = Y, 
                           X.mat         = X.mat, 
                           S.max         = 3)

## 
summary(result.obj$glm.obj)


## Video time-points (min, sec) of the selected PoIs
PoI_seconds <- TauIndHat_to_VideoSec(TauIndHat = result.obj$tau.ind.hat)
seconds_to_period(PoI_seconds)

## Link to the video: https://youtu.be/9F6UpuJIFaY

## 1st: after 'even genetiles' 
## 2nd: after 'Selling his brother's body parts'

## Peak-End Rule Predictors
X.peak.abs   <- apply(X.mat[-p,],2,function(x)max(abs(x)))
X.peak.pos   <- apply(X.mat[-p,],2,function(x)max(x))
X.peak.neg   <- apply(X.mat[-p,],2,function(x)min(x))
##
X.peak.t.abs <- apply(X.mat[-p,],2,function(x)t.grid[which.max(abs(x))])
X.peak.t.pos <- apply(X.mat[-p,],2,function(x)t.grid[which.max(x)])
X.peak.t.neg <- apply(X.mat[-p,],2,function(x)t.grid[which.min(x)])
##
X.end        <- X.mat[p,]

## Selected PoI-Predictors
X.PoI.1      <- X.mat[result.obj$tau.ind.hat[1], ]
X.PoI.2      <- X.mat[result.obj$tau.ind.hat[2], ]

## Comparision of three Models (PER-1, PER-2, and POI) 
PER_1_glm    <- glm(Y ~ X.peak.abs                + X.end, family = "binomial")
PER_2_glm    <- glm(Y ~ X.peak.pos  + X.peak.neg  + X.end, family = "binomial")
PoI_glm      <- glm(Y ~ X.PoI.1     + X.PoI.2            , family = "binomial")


## LaTeX Table
stargazer(PoI_glm, PER_1_glm, PER_2_glm, title="Results", align=TRUE)

## Comparision stats:
nagelkerke( PoI_glm )$Pseudo.R.squared.for.model.vs.null
nagelkerke(PER_1_glm)$Pseudo.R.squared.for.model.vs.null
nagelkerke(PER_2_glm)$Pseudo.R.squared.for.model.vs.null
##
somers2(fitted(PoI_glm), Y)
somers2(fitted(PER_1_glm), Y)
somers2(fitted(PER_2_glm), Y)

## Figure of estimation results:
par(mar=c(4,5,4,1)+.1, cex=1.3, cex.lab=1.3)
matplot(t.grid, X.mat, type = "l", lty=1, lwd=1, col=alpha("black", 0.35), xlab="", ylab="")
mtext(side = 1, text = "Time (Standardized)", line = 3,     at = 0.5, cex=1.3)
mtext(side = 3, text = "very\npositive", line = -1.0, at = -.11, cex=1.3)
mtext(side = 1, text = "very\nnegative", line = -.25,   at = -.11, cex=1.3)
##
points(x = X.peak.t.pos, y = X.peak.pos, pch="p")
points(x = X.peak.t.neg, y = X.peak.neg, pch="n")
abline(v=t.grid[result.obj$tau.ind.hat])
axis(side = 3, at = t.grid[result.obj$tau.ind.hat[1]], 
     labels = "''Even genitals''", 
     line = 0)
axis(side = 3, at = t.grid[result.obj$tau.ind.hat[2]], 
     labels = "''Selling his brother's\nbody parts''", 
     line = 0)
par(mar=c(5,4,4,2)+.1)
dev.off()

