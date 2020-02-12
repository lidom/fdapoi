## External Packages:
library("glmulti")    # glm: subset selection 
library("scales")     # colorscales
library("lubridate")  # working with times and dates
library("stargazer")  # LaTeX tables
library("rcompanion") # stats for glm
library("Hmisc")      # further stats for glm
library("devtools")   # to install the fdapoi package from GitHub

## Install the a ccompaning R-package fdapoi
# devtools::install("fdapoi") # install from local pkg
devtools::install_github("lidom/fdapoi/fdapoi") # install from Gitub

## Load the accompaning R-package fdapoi
library("fdapoi")    

## Attach the data
data("emotion_rating")

## Variables: #################
Y     <- emotion_rating$Y
X.mat <- t(as.matrix(emotion_rating[,-1]))
##
(p     <- nrow(X.mat))
(N     <- ncol(X.mat))
##
a      <- 0
b      <- 1
t.grid <- (1:p-1)/(p-1)
## #############################

## Estimation
result.obj  <- FUN_PoI_BIC(Y             = Y, 
                           X.mat         = X.mat)

## 
summary(result.obj$glm.obj)

## Video time-points (min, sec) of the selected PoIs
PoI_seconds <- TauIndHat_to_VideoSec(TauIndHat = result.obj$tau.ind.hat)
seconds_to_period(PoI_seconds)

## Link to the video: https://youtu.be/9F6UpuJIFaY

## 1st point of impact (movie scene):   Portrait shot of African albino nervously moving eyes 
## 2nd point of impact (spoken workds): [...] selling his brother's body parts

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
stargazer(PoI_glm, PER_1_glm, PER_2_glm, title="Results", align=TRUE, digits=2)

## Comparision stats:
round(nagelkerke( PoI_glm )$Pseudo.R.squared.for.model.vs.null, 2)
round(nagelkerke(PER_1_glm)$Pseudo.R.squared.for.model.vs.null, 2)
round(nagelkerke(PER_2_glm)$Pseudo.R.squared.for.model.vs.null, 2)
##
round(somers2(fitted(PoI_glm), Y), 2)
round(somers2(fitted(PER_1_glm), Y), 2)
round(somers2(fitted(PER_2_glm), Y), 2)


## Figure of estimation results:
cex_value <- 1.6
par(mar=c(4,5.1,5,1)+.1, cex=cex_value, cex.lab=cex_value, family="serif")
matplot(t.grid, X.mat, type = "l", lty=1, lwd=.25, col=gray(.5), xlab="", ylab="")
mtext(side = 1, text = "Time (Standardized)", line = 3,     at = 0.5,  cex=cex_value)
mtext(side = 3, text = "very\npositive",      line =  -1.8, at = -.19, cex=cex_value)
mtext(side = 1, text = "very\nnegative",      line = -.5,   at = -.19, cex=cex_value)
##
points(x = X.peak.t.pos, y = X.peak.pos, pch="p")
points(x = X.peak.t.neg, y = X.peak.neg, pch="n")
abline(v=t.grid[result.obj$tau.ind.hat])
axis(side = 3, at = t.grid[result.obj$tau.ind.hat[1]], labels = "", line = 0)
axis(side = 3, at = t.grid[result.obj$tau.ind.hat[1]],
     labels = expression(hat(tau)[1]),
     line = -.5, tick = F)
axis(side = 3, at = t.grid[result.obj$tau.ind.hat[2]], labels = "", line = 0)
axis(side = 3, at = t.grid[result.obj$tau.ind.hat[2]],
     labels = expression(hat(tau)[2]),
     line = -.5, tick = F)
mtext(side = 3, at = t.grid[result.obj$tau.ind.hat[1]], 
      text = "Movie scene:\nPortrait shot of African albino\nnervously moving eyes", 
      line = 2, cex=cex_value)
mtext(side = 3, at = t.grid[result.obj$tau.ind.hat[2]],
      text = "Spoken works:\n''Selling his body parts''",
      line = 2, cex=cex_value)
par(mar=c(5,4,4,2)+.1)
dev.off()
