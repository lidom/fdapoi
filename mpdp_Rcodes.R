#############################################################################
#############################################################################
#############################################################################
######## F. FERRATY, P. HALL, P. VIEU (2010). Most predictive       #########
######## design points for functional data predictors. Biometrika,  #########
######## 97, 807-824                                                #########
#############################################################################
#############################################################################
#############################################################################
############################            #####################################
############################ R ROUTINES #####################################
############################            #####################################
#############################################################################
#############################################################################
#############################################################################

#############################################################################
#############################################################################
#############################################################################
############################                    #############################
############################ STEPWISE ALGORITHM #############################
############################                    #############################
#############################################################################
#############################################################################
#############################################################################

#####################################################################
#####################################################################
############### IDENTIFYING INFLUENCE POINTS THROUGH ################
############### LOCAL LINEAR MULTIVARIATE REGRESSION ################
############### WITH SIMPLIFIED BANDWIDTHS:          ################
############### STEPWISE ALGORITHM                   ################
#####################################################################
#####################################################################
mpdp = function(CURVES, Response, nbmax=5, nbbw=5, Grid=(1:ncol(CURVES)), pcvpar=1, kind.of.kernel="quadratic2")
##############################################################################
# Select most-predictive design points ("mpdp") when regressing a saclar response 
# "Response" on functional data predictors "CURVES"
# "CURVES" matrix (n rows and p columns) containing the functional data predictors 
#          X1, X2,...,Xn (n curves) stored row by row: 
#          row 1: X1(t1),...,X1(tp) 
#          row 2: X2(t1),...,X2(tp) 
#          .......................
#          row p: Xn(t1),...,Xn(tp) 
#          it is worth noting that this sample of curves must be observed at the 
#          same p design points t1,...,tp
# "Response" vector of length n containing the corresponding responses Y1,...,Yn 
# "nbmax" maximum number of selected design points (default value: 5)
# "nbbw" size of the grid of bandwidths (h1,h2,..,hnbbw) used to choose the optimal 
#          one. From "nbbw", a grid of bandwidths is built in terms of k-nearest
#          neighbours: 
#          h1 takes into account 2% of the sample (for building the kernel estimator)
#          ......................................  
#          hk takes into account [2+{(50-2)/(nbbw-1)}*(k-1)]% of the sample
#          ......................................  
#          hnbbw takes into account 50% of the sample
#          (default value 5 ==> 2%, 14%, 26%, 38%, 50%)
# "Grid" vector containg a subset of design points in which one looks for the 
#          "nbmax" most-predicitve ones (default value: the whole set of 
#          design points is considered). 
# "pcvpar" parameter used in the penalized corssvalidation criterion (delta0 in 
#          the corresponding paper). The default value is 1. 
# "kind.of.kernel" gives the kernel used for computing the modes (default kernel 
#          is the quadratic one)          
##############################################################################
# Returns an object of class "mpdp": 
# "$Bwdselection"  vector identifying the most-predicitve design points retained 
#                    definitely after the backward deletion stage
# "$Fwdselection"  vector identifying the most-predicitve design points retained 
#                    just after the forward addition stage
# "$bandwidth"     optimal bandwidth (see paper for more details)
# "$Bwd.estimated.values" vector containing the leave-one-out local linear 
#                    estimations of the responses obtained by using the  
#                    selected most-predictive design points after backward stage 
# "$Bwdpcv"        vector containing the values of the penalized crossvalidation 
#                    criterion at each step of the backward deletion stage 
# "$Fwd.estimated.values" vector containing the leave-one-out local linear 
#                    estimations of the responses obtained by using the  
#                    selected most-predictive design points just after forward stage 
# "$Fwdpcv"        vector containing the values of the penalized crossvalidation 
#                    criterion at each step of the forward addition stage 
#                          final penalized crossvalidation criterion
# "$fpcv.status"   "yes" means pcv reached its minimum at the forward addition stage;
#                  "no" means that the number of selected design points to reach the 
#                    minimum must be greater than "nbmax" (i.e. start again "mpdp" 
#                    with a larger "nbmax")
# REMARK : when only one design point is selected in the forward stage, the procedure 
#          stops without launching backward stage (consequently all backward stage's 
#          outputs are not available)
##############################################################################
{
        Selected.points <- NULL
        XX <- NULL
        kernel <- get(kind.of.kernel)
        nind <- nrow(CURVES)
        nind2 <- nind^2
        Ind11 <- rep(1:nind,nind)
        Ind12 <- rep(1:nind,rep(nind,nind))
        Gridnew <- Grid
        Gridpoints <- NULL
        poszero <- ((0:(nind-1))*nind)+(1:nind)
        Nbknn <- ceiling(seq(0.02,0.5,length=nbbw)*nind)
        varesp <- var(Response)
        Leaveoneout.estimated.values <- list()
        Meancurve <- apply(CURVES,2,mean)
        CCURVES <- t(t(CURVES)-Meancurve)
        Stdev <- sqrt(apply(CCURVES^2,2,mean))
        STDCURVES <- t(t(CCURVES)/Stdev)
        Bandwidths <- 0
        Fcv <- 0
        Fpcv <- 0
        minimum.reached.for.pcv <- "no"
        fpcv.ini <- sum((Response-mean(Response))^2)/nind
##########################################################
# FORWARD ADDITION
##########################################################
        for(kk in 1:nbmax){
               cat(paste("Forward addition in progress"),sep="\n")
               if(length(Selected.points)!=0){Gridnew <- Gridnew[-selectedpoint]}
               lgridnew <- length(Gridnew)
               INVSUMOFSQUARES <- matrix(0,nrow=nbbw,ncol=lgridnew)
               ESTIMATIONS <- array(0,c(nbbw,lgridnew,nind))
               jj=0
               for(tt in Gridnew){
                      jj=jj+1
                      Xxtt <- STDCURVES[,tt]
                      XXNEW <- cbind(XX, Xxtt)
                      XXNEW.ALT.ROW <- XXNEW[Ind11,]
                      XXNEW.REP.ROW <- XXNEW[Ind12,]
                      D2 <- (XXNEW.ALT.ROW-XXNEW.REP.ROW)^2
                      Prox <- sqrt(apply(as.matrix(D2),1,sum))
                      PROX <- matrix(Prox, nrow=nind)
                      Max <- apply(PROX,1,max) 
                      diag(PROX)=Max 
                      PROXSORTED <- apply(PROX,1,sort)
                      Bw <- apply(PROXSORTED[Nbknn,],1,max)
                      ll=0
                      for(ll in 1:nbbw){
                             Dist <- Prox/Bw[ll]
                             Ker <- rep(0,nind2)
                             Ker[-poszero] <- kernel(Dist[-poszero])
                             KERNEL <- matrix(Ker, nrow=nind)
                             Denom <- apply(KERNEL,1,sum)
                             if(sum(Denom==0)>0){
                                     INVSUMOFSQUARES[ll,jj] <- 0
                                     next
                             }else{
                                     NUM <- KERNEL %*% XXNEW
                                     XBAR <- NUM/Denom
                                     Ynum <- as.vector(KERNEL %*% Response) 
                                     Ybar <- Ynum/Denom
                                     T11 <- crossprod(t(KERNEL)*Response,XXNEW)
                                     T12 <- XBAR*Ynum
                                     T1 <- T11-T12
                                     for(ii in 1:nind){
                                           CXXNEW <- t(t(XXNEW)-XBAR[ii,])
                                           T2 <- crossprod(CXXNEW,(CXXNEW*KERNEL[ii,]))
                                           if(abs(det(T2))<1e-8){
                                                 ESTIMATIONS[ll,jj,ii] <- Ybar[ii]
                                           }else{
                                                 Bb <- solve(T2,T1[ii,])
                                                 ESTIMATIONS[ll,jj,ii] <- Ybar[ii]+crossprod(Bb,CXXNEW[ii,])
                                           } 
                                     }
                                     INVSUMOFSQUARES[ll,jj] <- 1/sum((Response-ESTIMATIONS[ll,jj,])^2)
                             }
                      }
               }
               nr <- nrow(INVSUMOFSQUARES)
               Numcolofmax <- max.col(INVSUMOFSQUARES)
               numrowofmax <- which.max(INVSUMOFSQUARES[cbind(1:nr,Numcolofmax)])
               selectedpoint <- Numcolofmax[numrowofmax]
               Selected.points <- c(Selected.points,selectedpoint)
               Gridpoints <- c(Gridpoints,Gridnew[selectedpoint])
               Xxtt <- STDCURVES[,Gridnew[selectedpoint]]
               XX <- cbind(XX, Xxtt)
               XX.ALT.ROW <- XX[Ind11,]
               XX.REP.ROW <- XX[Ind12,]
               D2 <- (XX.ALT.ROW-XX.REP.ROW)^2
               Prox <- sqrt(apply(as.matrix(D2),1,sum))
               PROX <- matrix(Prox, nrow=nind)
               Max <- apply(PROX,1,max) 
               diag(PROX) <- Max 
               PROXSORTED <- apply(PROX,1,sort)
               Bw <- apply(PROXSORTED[Nbknn,],1,max)
               Bandwidths[kk] <- Bw[numrowofmax]
               Leaveoneout.estimated.values[[kk]] <- ESTIMATIONS[numrowofmax,Numcolofmax[numrowofmax],]
               Fcv[kk] <- sum((Response-Leaveoneout.estimated.values[[kk]])^2)/nind
               Fpcv[kk] <- Fcv[kk]*(1+kk*pcvpar/log(nind))
               if(kk==1) {
                    if(Fpcv[kk]>fpcv.ini) {
                         stop("No design point selected: stepwise algorithm over")
                    } else {
                         cat(paste("Design point",Gridnew[selectedpoint],"selected"),sep="\n")
                    }
               } else {
                    if(Fpcv[kk]>Fpcv[kk-1]) {
                         minimum.reached.for.pcv <- "yes"
                         break      
                    } else {
                         cat(paste("Design point",Gridnew[selectedpoint],"selected"),sep="\n")
                    }
               }
        }
        if(minimum.reached.for.pcv=="yes") {
                  Forwardselection <- Gridpoints[-kk]
                  Fbwstd <- Bandwidths[kk-1]
                  Fleaveoneout.estimated.values <- Leaveoneout.estimated.values[[kk-1]]
                  Fpcv <- Fpcv[-kk]
        } else {
             cat("Minimum for PCV not confirmed; increase nbmax", sep="\n")
             Forwardselection <- Gridpoints
             Fbwstd <- Bandwidths[kk]
             Fleaveoneout.estimated.values <- Leaveoneout.estimated.values[[nbmax]]
        }
##########################################################
# BACKWARD DELETION
##########################################################
        lforwardselection <- length(Forwardselection)
        if(lforwardselection>=2){
             Selected.points <- NULL
             Bandwidths <- 0
             XX <- STDCURVES[,Forwardselection]
             Gridpoints <- NULL
             Gridnew <- Forwardselection
             Bcv <- 0
             Bpcv <- Fpcv[lforwardselection]
             exitloop <- F
             for(kk in 1:(lforwardselection-1)){
                    if(kk==1) {abbrev <- "st"}
                    if(kk==2) {abbrev <- "nd"}
                    if(kk==3) {abbrev <- "rd"}
                    if(kk>3) {abbrev <- "th"}
                    cat(paste("Backward deletion: is the deletion of a ",kk,abbrev," design point necessary?", sep=""),sep="\n")
                    if(length(Selected.points)!=0){Gridnew <- Gridnew[-selectedpoint]}
                    lgridnew <- length(Gridnew)
                    SUMOFSQUARES <- matrix(0,nrow=lgridnew,ncol=nbbw)
                    ESTIMATIONS <- array(0,c(nbbw,lgridnew,nind))
                    jj=0
                    for(tt in 1:lgridnew){
                           jj=jj+1
                           XXNEW <- as.matrix(XX[,-tt])
                           XXNEW.ALT.ROW <- XXNEW[Ind11,]
                           XXNEW.REP.ROW <- XXNEW[Ind12,]
                           D2 <- (XXNEW.ALT.ROW-XXNEW.REP.ROW)^2
                           Prox <- sqrt(apply(as.matrix(D2),1,sum))
                           PROX <- matrix(Prox, nrow=nind)
                           Max <- apply(PROX,1,max) 
                           diag(PROX)=Max 
                           PROXSORTED <- apply(PROX,1,sort)
                           Bw <- apply(PROXSORTED[Nbknn,],1,max)
                           ll=0
                           for(ll in 1:nbbw){
                                  Dist <- Prox/Bw[ll]
                                  Ker <- rep(0,nind2)
                                  Ker[-poszero] <- kernel(Dist[-poszero])
                                  KERNEL <- matrix(Ker, nrow=nind)
                                  Denom <- apply(KERNEL,1,sum)
                                  if(sum(Denom==0)>0){
                                          SUMOFSQUARES[jj,ll] <- sum((Response-mean(Response))^2)
                                          next
                                  }else{
                                          NUM <- KERNEL %*% XXNEW
                                          XBAR <- NUM/Denom
                                          Ynum <- as.vector(KERNEL %*% Response) 
                                          Ybar <- Ynum/Denom
                                          T11 <- crossprod(t(KERNEL)*Response,XXNEW)
                                          T12 <- XBAR*Ynum
                                          T1 <- T11-T12
                                          for(ii in 1:nind){
                                                #cat(paste("ii ",ii),sep="\n")
                                                CXXNEW <- t(t(XXNEW)-XBAR[ii,])
                                                T2 <- crossprod(CXXNEW,(CXXNEW*KERNEL[ii,]))
                                                if(abs(det(T2))<1e-8){
                                                      ESTIMATIONS[ll,jj,ii] <- Ybar[ii]
                                                }else{
                                                      Bb <- solve(T2,T1[ii,])
                                                      ESTIMATIONS[ll,jj,ii] <- Ybar[ii]+crossprod(Bb,CXXNEW[ii,])
                                                } 
                                          }
                                          SUMOFSQUARES[jj,ll] <- sum((Response-ESTIMATIONS[ll,jj,])^2)
                                  }
                           }
                    }
                    nr <- nrow(SUMOFSQUARES)
                    Numcolofmax <- max.col(1/SUMOFSQUARES)
                    selectedpoint <- which.min(SUMOFSQUARES[cbind(1:nr,Numcolofmax)])
                    Selected.points <- c(Selected.points,selectedpoint)
                    Gridpoints <- Gridnew[-selectedpoint]
                    XX <- as.matrix(STDCURVES[,Gridpoints]) 
                    XX.ALT.ROW <- XX[Ind11,]
                    XX.REP.ROW <- XX[Ind12,]
                    D2 <- (XX.ALT.ROW-XX.REP.ROW)^2
                    Prox <- sqrt(apply(as.matrix(D2),1,sum))
                    PROX <- matrix(Prox, nrow=nind)
                    Max <- apply(PROX,1,max) 
                    diag(PROX) <- Max 
                    PROXSORTED <- apply(PROX,1,sort)
                    Bw <- apply(PROXSORTED[Nbknn,],1,max)
                    Bandwidths[kk] <- Bw[selectedpoint]
                    Leaveoneout.estimated.values[[kk]] <- ESTIMATIONS[Numcolofmax[selectedpoint],selectedpoint,]
                    Bcv[kk] <- sum((Response-Leaveoneout.estimated.values[[kk]])^2)/nind
                    Bpcv[kk+1] <- Bcv[kk]*(1+(lforwardselection-kk)*pcvpar/log(nind))
                    if((kk<lforwardselection) && (Bpcv[kk]<Bpcv[kk+1])) {
                         exitloop <- T
                         cat(paste("No!"),sep="\n")
                         break      
                    } else {
                         cat(paste("Yes: design point",Gridnew[selectedpoint], "deleted"),sep="\n")
                    }
             }
             if(exitloop) {
                  Backwardselection <- Gridnew
                  if(kk==1) {
                     Bbwstd <- Fbwstd
                     Bleaveoneout.estimated.values <- Fleaveoneout.estimated.values
                  } else {
                     Bbwstd <- Bandwidths[kk-1]
                     Bleaveoneout.estimated.values <- Leaveoneout.estimated.values[[kk-1]]
                  }
             } else {
                  Backwardselection <- Gridpoints
                  Bbwstd <- Bandwidths[kk]
                  Bleaveoneout.estimated.values <- Leaveoneout.estimated.values[[kk]]
             }
             cat("Selected most-predicitve design points:", sep="\n")
             cat(Backwardselection, sep="\n")
             object <- list()
             class(object) <- "mpdp"
             object$Bwdselection <- Backwardselection
             object$Fwdselection <- Forwardselection
             object$bandwidth <- Bbwstd
             object$Bwd.estimated.values <- Bleaveoneout.estimated.values
             object$Bwdpcv <- Bpcv
             object$Fwd.estimated.values <- Fleaveoneout.estimated.values
             object$Fwdpcv <- Fpcv
             object$fpcv.status <- minimum.reached.for.pcv
             object$predictors <- CURVES
             object$responses <- Response
             object$kind.of.kernel <- kind.of.kernel
             return(object) 
         }else{
             cat("Selected most-predicitve design points:", sep="\n")
             cat(Forwardselection, sep="\n")
             object <- list()
             class(object) <- "mpdp"
             object$Bwd.estimated.values <- "No backward estimations (no backward stage)"
             object$Bwdselection <- "No backward selection (no backward stage)"
             object$Bwdpcv <- "No backward penalized cross-validation (no backward stage)"
             object$Fwdselection <- Forwardselection
             object$bandwidth <- Fbwstd
             object$Fwd.estimated.values <- Fleaveoneout.estimated.values
             object$Fwdpcv <- Fpcv
             object$fpcv.status <- minimum.reached.for.pcv
             object$predictors <- CURVES
             object$responses <- Response
             object$kind.of.kernel <- kind.of.kernel
         return(object) 
         }
}
#####################################################################
#####################################################################
############### LOCAL LINEAR MULTIVARIATE PREDICTION ################
############### USING AN OBJECT OF CLASS "mpdp"      ################
#####################################################################
#####################################################################
predict.mpdp <- function(object, CURVPRED)
{
      if(is.vector(CURVPRED)) CURVPRED <- as.matrix(CURVPRED) 
      VECEST <- object$predictors[,object$Bwdselection]
      VECPRED <- CURVPRED[,object$Bwdselection]
      Meancurve <- apply(VECEST,2,mean)
      CCOVARIATES <- t(t(VECEST)-Meancurve)
      Stdev <- sqrt(apply(CCOVARIATES^2,2,mean))
      COVARIATE <- t(t(CCOVARIATES)/Stdev)
      CVECPRED <- t(t(VECPRED)-Meancurve)
      PRED <- t(t(CVECPRED)/Stdev)
      nind <- nrow(COVARIATE)
      kernel <- get(object$kind.of.kernel)
      nind2 <- nrow(PRED)
      Ind1 <- rep(1:nind2,nind)
      Ind2 <- rep(1:nind,rep(nind2,nind))
      COVAR.REP.ROW <- COVARIATE[Ind2,]
      PRED.ALT.ROW <- PRED[Ind1,]
      D2 <- (PRED.ALT.ROW-COVAR.REP.ROW)^2
      D2OVERBW <- D2/(object$bandwidth^2)
      Dist <- sqrt(apply(D2OVERBW,1,sum))
      Ker <- rep(0,nind2*nind)
      Ker[Dist>0] <- kernel(Dist[Dist>0])
      KERNEL <- matrix(Ker, nrow=nind2)
      Denom <- apply(KERNEL,1,sum)
      NUM <- KERNEL %*% COVARIATE
      XBAR <- NUM/Denom
      Ynum <- as.vector(KERNEL %*% object$responses) 
      Ybar <- Ynum/Denom
      TP1 <- crossprod(t(KERNEL)*object$responses,COVARIATE)
      TP2 <- XBAR*Ynum
      TP1M2 <- TP1-TP2
      Predictions <- 0
      for(ii in 1:nind2){
             CCOVAR <- t(t(COVARIATE)-XBAR[ii,])
             TP3 <- crossprod(CCOVAR,(CCOVAR*KERNEL[ii,]))
             if(abs(det(TP3))<1e-8){
                   Predictions[ii] <- Ybar[ii]
             }else{
                   Bb <- solve(TP3,TP1M2[ii,])
                   Cpred <- t(t(PRED)-XBAR[ii,])
                   Predictions[ii] <- Ybar[ii]+crossprod(Bb,Cpred[ii,]) 
             } 
      }
      return(Predictions) 
}    
####################################
### ASYMMETRICAL QUADRATIC KERNEL  #
####################################
quadratic2 <- function(u) 
{
        u[u<0] <- 1
        u[u>1] <- 1 
        return(1 - (u)^2) 
}

