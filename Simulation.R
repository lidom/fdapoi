## Load packages:
library("glmulti")   # glm: subset selection 
library("MASS")      # multivariate normals
library("np")        # for nonparametric regression
options(np.messages = FALSE) # surpress run-time messages from np (like "Multistart x of y")

## Install and load the accompaning R-package 
## devtools::install(pkg="fdapoi_pkg")
library("fdapoi")   



## ###############################
## Simulation parameters #########
## ###############################
## Number of MC-repetitions
B               <- 1000
## Domain
a               <- 0
b               <- 1
##
DGP.seq         <- c(1, 2, 3, 4, 5)

N.seq           <- c(100,200,500,1000,3000)
p.seq           <- c(100,500,1000)
##
##
(Start.Time <- Sys.time())
##
for (DGP in DGP.seq) { 
  ## ####################################################
  ## set seed
  set.seed(DGP)
  ## ####################################################
  ## Feedback
  cat("\nDGP =", DGP)
  ## ####################################################
  for (p in p.seq) { # p <- 100
    ## ####################################################
    ## Feedback
    cat("\np =", p,"\n")
    ## ####################################################
    ##
    if(DGP==1){
      beta0         <- 1
      beta          <- c(4)
      S             <- length(beta)
      tau.true      <- c(1/2) 
      t.grid        <- (1:p-1)/(p-1)
      tau.ind.true  <- rep(NA,S)
      for(s in 1:S){
        tau.ind.true[s]  <- which.min(abs(t.grid - tau.true[s]))
      }
    }
    if(any(DGP==c(2,4,5))){
      beta0         <- 1
      beta          <- c(-6, 5)
      S             <- length(beta)
      tau.true      <- c(1/3,2/3) 
      t.grid        <- (1:p-1)/(p-1)
      tau.ind.true  <- rep(NA,S)
      for(s in 1:S){
        tau.ind.true[s]  <- which.min(abs(t.grid - tau.true[s]))
      }
      if (DGP == 4) {
        ## covariance matrix
        calcSigma <- function(X1, X2, dd = 1 / 10) {
          Sigma <- matrix(rep(0, length(X1) * length(X2)), nrow = length(X1))
          for (i in 1:nrow(Sigma)) {
            for (j in 1:ncol(Sigma)) {
              Sigma[i, j] <- exp(-(abs(X1[i] - X2[j]) / dd) ^ 2)
            }
          }
          return(Sigma)
        }
        CoVGP <- calcSigma(t.grid, t.grid, dd = 1 / 10)  
      }
    }
    if(DGP==3){
      beta0         <- 1
      beta          <- c(-6, 6, -5, 5)
      S             <- length(beta)
      tau.true      <- c(1/6, 2/6, 4/6, 5/6) 
      t.grid        <- (1:p-1)/(p-1)
      tau.ind.true  <- rep(NA,S)
      for(s in 1:S){
        tau.ind.true[s]  <- which.min(abs(t.grid - tau.true[s]))
      }
    }
    ##
    for(N in N.seq){ # N <- 100
      ## ####################################################
      ## Feedback
      cat("N =", N," ")
      ## ####################################################
      if (DGP == 1) {
        sim.results <- matrix(NA, nrow = B, ncol = 6)
        colnames(sim.results) <- c("beta0.hat.PoI",  "beta1.hat.PoI",  "tau1.ind.hat.PoI", 
                                   "beta0.hat.LMcK", "beta1.hat.LMcK", "tau1.ind.hat.LMcK")
      } else {
        sim.results <- matrix(NA, nrow = B, ncol = 20)
        colnames(sim.results) <- c("beta0.hat.PoI", "beta1.hat.PoI",   "beta2.hat.PoI",   "beta3.hat.PoI",   "beta4.hat.PoI",
                                   "tau1.ind.hat.PoI","tau2.ind.hat.PoI","tau3.ind.hat.PoI","tau4.ind.hat.PoI", "S.hat.PoI", 
                                   "beta0.hat.TRH", "beta1.hat.TRH",   "beta2.hat.TRH",   "beta3.hat.TRH",   "beta4.hat.TRH",
                                   "tau1.ind.hat.TRH","tau2.ind.hat.TRH","tau3.ind.hat.TRH","tau4.ind.hat.TRH", "S.hat.TRH")
      }
      ## NP-Regression
      sim.np.results           <- matrix(NA, nrow = B, ncol = 2)
      colnames(sim.np.results) <- c("np.mase.PoI", "np.mase.TRH")
      ## ######################################################
      for(repet in 1:B){
        ## ####################################################
        ## Generate Data:
        ## ####################################################
        ## Discretization points
        t.vec  <- seq(a,b,length.out=p)
        ## Random functions:
        if(any(DGP==c(1,2,3))){
          X.mat  <- OU.SIM(a=a,b=b,p=p, N=N)
        }
        if(DGP==4){
          X.mat <- t(MASS::mvrnorm(N, rep(0, p), CoVGP))
        }
        if(DGP==5){
          X.mat  <- matrix(NA,nrow = p, ncol = N)	
          for(i in 1:N){X.mat[,i] <- exp(BM.SIM(a=a,b=b,p=p))}
        }
        ##
        X.mat  <- X.mat - rowMeans(X.mat)
        ## PoIs:
        X.tau  <- t(X.mat[tau.ind.true, 1:N, drop=F])
        ## Generate dependent variable:
        pi.x   <- exp(beta0 + X.tau %*% beta) / (1 + exp(beta0 + X.tau %*% beta))
        Y      <- rbinom(n=length(1:N), size=1, prob=pi.x)
        
        ## Estimations
        if (DGP == 1) {
          ## ####################################################
          ## Estimation procedure of Lindquist & McKeague
          ## ####################################################
          LMcK.estim        <- FUN_LMcK(Y=Y, X.mat=X.mat) 
          beta0.hat.LMcK    <- LMcK.estim$beta0.hat
          beta.hat.LMcK     <- LMcK.estim$beta.hat
          tau.ind.hat.LMcK  <- LMcK.estim$tau.ind.slct
        }
        ## ####################################################
        ## Our PoI estimation procedure
        ## ####################################################
        logit.PoI.estim  <- try(
          FUN_PoI_BIC(
            Y          = Y, 
            X.mat      = X.mat, 
            S.max      = ifelse(DGP==1, S, S+2))
        )
        
        if(!Error_Checker(logit.PoI.estim)) {
          ## Allocation of results
          alloc.obj.PoI    <- Allocator(tau.ind.hat  = logit.PoI.estim$tau.ind.hat, 
                                        beta.hat     = logit.PoI.estim$beta.hat, 
                                        tau.ind.true = tau.ind.true, 
                                        p            = p)
          ## Results
          beta0.hat.PoI    <- logit.PoI.estim$beta0.hat
          beta.hat.PoI     <- alloc.obj.PoI$beta.hat.alloc
          tau.ind.hat.PoI  <- alloc.obj.PoI$tau.ind.hat.alloc
          S.hat.PoI        <- logit.PoI.estim$S.hat
        } else {
          beta0.hat.PoI    <- NA
          beta.hat.PoI     <- NA
          tau.ind.hat.PoI  <- NA
          S.hat.PoI        <- NA
        }
        
        ## ####################################################
        ## Estimation procedure PoI with Lambda-threshold
        ## ####################################################
        logit.TRH.estim <- FUN_PoI_TRH(Y = Y, X.mat = X.mat, scale = ifelse(DGP==5,3,1.5)) 
        ## Allocation of results
        
        alloc.obj.TRH   <- Allocator(tau.ind.hat  = logit.TRH.estim$tau.ind.hat, 
                                     beta.hat     = logit.TRH.estim$beta.hat, 
                                     tau.ind.true = tau.ind.true, 
                                     p            = p)
        ## Results
        beta0.hat.TRH    <- logit.TRH.estim$beta0.hat
        beta.hat.TRH     <- alloc.obj.TRH$beta.hat.alloc
        tau.ind.hat.TRH  <- alloc.obj.TRH$tau.ind.hat.alloc
        S.hat.TRH        <- logit.TRH.estim$S.hat
        ##
        ## Equalize lengths of vectors (4 max. length of DGP3)
        if(DGP != 1) {
          length(beta.hat.PoI)     <- 4
          length(beta.hat.TRH)     <- 4
          length(tau.ind.hat.PoI)  <- 4
          length(tau.ind.hat.TRH)  <- 4
        }
        ## ####################################################
        ## Feedback
        ## if(repet %% 200 ==0){cat("repet/B =", repet,"/", B,"\n")}
        ## ####################################################
        ## Collect simulation results:
        if(DGP==1){
          sim.results[repet,] <- c(
            ## PoI
            beta0.hat.PoI,
            beta.hat.PoI,
            tau.ind.hat.PoI,
            ## LMcK
            beta0.hat.LMcK,
            beta.hat.LMcK,
            tau.ind.hat.LMcK
          )
        }else{
          sim.results[repet,] <- c(
            ## PoI
            beta0.hat.PoI,
            beta.hat.PoI,
            tau.ind.hat.PoI,
            S.hat.PoI,
            ## THR
            beta0.hat.TRH,
            beta.hat.TRH,
            tau.ind.hat.TRH,
            S.hat.TRH
          )
        }
        ## ####################################################
        ## NP-regression
        ## ####################################################
        if (Error_Checker(logit.PoI.estim) | N >= N.seq[4] | length(c(na.omit(tau.ind.hat.PoI))) == 0) {
          np.mase.PoI <- NA
        } else {
          ## Create data frame consisting of Y and X evaluated at the estimated points of impact:
          np.dat           <- data.frame(Y, t(X.mat[c(na.omit(tau.ind.hat.PoI)),, drop = F])) 
          colnames(np.dat) <- c(paste(c("Y", paste("X", rep(1:length(logit.PoI.estim$tau.ind.hat)), sep = ""))))
          ## Derive optimal bandwidths by aic.CV:
          bw.all           <- npregbw(formula = as.formula(paste("Y~", paste(names(np.dat)[-1], collapse = "+"))),
                                      regtype = "lc", bwmethod = "cv.aic", data = np.dat)
          ## Execute non parametric regression:
          model.np         <- npreg(bws = bw.all)
          ## Save the mean average squared error
          np.mase.PoI      <- mean((pi.x - fitted(model.np))^2)
        }
        ## ####################################################
        ## Plotting: let's look at a slice of the nonparametric regression along one direction
        ## x.eval <- data.frame(X1 = seq(-2, 2, length.out = 100))
        ## model.np.pred <-predict(model.np, newdata = x.eval)
        ## plot(as.matrix(x.eval),model.np.pred,type="l")
        ## lines(as.matrix(x.eval), exp(beta0 + as.matrix(x.eval) * beta) / (1 + exp(beta0 + as.matrix(x.eval) * beta)), type = "l", col = "red")
        #######################################################################
        ##
        if (Error_Checker(logit.TRH.estim) | N >= N.seq[4] | length(c(na.omit(tau.ind.hat.THR))) == 0) {
          np.mase.TRH <- NA
        } else {
          ## Create data frame consisting of Y and X evaluated at the estimated points of impact:
          np.dat.TRH    <- data.frame(Y, t(X.mat[c(na.omit(tau.ind.hat.THR)),, drop = F])) 
          ## Derive optimal bandwidths by aic.CV:
          bw.all.TRH    <- npregbw(formula = as.formula(paste("Y~", paste(names(np.dat.TRH)[-1], collapse = "+"))),
                                   regtype = "lc", bwmethod = "cv.aic", data = np.dat.TRH)
          ## Execute non parametric regression:
          model.np.TRH  <- npreg(bws = bw.all.TRH)
          ## Save the mean average squared error
          np.mase.TRH   <- mean((pi.x - fitted(model.np.TRH))^2)
        }
        ## ####################################################
        ## Plotting: let's look at a slice of the nonparametric regression along one direction
        ## x.eval <- data.frame(X1 = seq(-2, 2, length.out = 100))
        ## model.np.pred<-predict(model.np.TRH, newdata = x.eval)
        ## plot(as.matrix(x.eval),model.np.pred,type="l")
        ## lines(as.matrix(x.eval), exp(beta0 + as.matrix(x.eval) * beta) / (1 + exp(beta0 + as.matrix(x.eval) * beta)), type = "l", col = "red")
        #######################################################################
        ## ####################################################
        ## Collect simulation results for np-reg:
        sim.np.results[repet,] <- c(np.mase.PoI, np.mase.TRH)
      }## End of foreach-loop over 1:B
      ## ###################
      ## Save results
      ## ###################
      save(sim.results,    file = paste0("Simulation_Results/Sim_Results_DGP=",       DGP, "_N=", N, "_p=", p, ".RData"))
      save(sim.np.results, file = paste0("Simulation_Results/Sim_npReg_Results_DGP=", DGP, "_N=", N, "_p=", p, ".RData"))
    }## End of N.seq loop
  }## End of p.seq loop
}## End of DGP.seq loop
## ######################################################
## ######################################################
End.Time <- Sys.time()
##
## Run-time:
round(End.Time - Start.Time, 2)
##