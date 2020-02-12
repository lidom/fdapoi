#' FUN_PoI_BIC
#'
#' This function is an implementation of the fully data driven estimator POI.
#' @param Y dependent variable (n-vector)
#' @param X.mat pxn matrix of discretized functional predictors
#' @param a lower interval border (default a=0)
#' @param b upper interval border (default b=1)
#' @param CoVars additinal covariables (default CoVars=NULL)
#' @param S.max maximum number of points of impact (default S.mat=5)
#' @param show.pb show a progress bar (default show.pb = FALSE)
#' @param k.seq define k.seq (i.e., the delta values - in terms of grid-points - to be considered). The default (k.seq=NULL) means: k.seq = ceiling(p * seq(0.02, 0.40, by = 0.02))
#' @param IC information criterion ("bic" (default), or "aic")
#' @param standardize.X standardize (to unit variance) the functional predictors (default standardize.X  = FALSE)
#' @param center.X center the functional predictors (default center.X = FALSE)
#' @param family any family of ?stats::family (default family = "binomial")
#' @export FUN_PoI_BIC
FUN_PoI_BIC <- function(Y = Y, X.mat = X.mat, a = 0, b = 1,
                        CoVars         = NULL, 
                        S.max          = 10, 
                        show.pb        = FALSE, 
                        k.seq          = NULL,
                        IC             = c("bic", "aic")[1],
                        standardize.X  = FALSE,
                        center.X       = FALSE,
                        family         = "binomial"){
  
  ## Default for glmuli:
  glmulti.method = c("h", "g", "l", "d")[1]
  
  ## For future use:
  nbasis.max     = 0    # if nbasis.max==0, no FLR is estimated
  
  ##
  if (!is.null(CoVars)) { if (is.null(nrow(CoVars))) { stop("CoVars must a matrix with ncol(CoVars)==N.") }}
  if (S.max <= 0 & nbasis.max <=0) {stop("Either S.max or nbasis.max must be strictly geater than 0.") }
  ##
  
  ## Dimensions
  N <- ncol(X.mat)
  p <- nrow(X.mat)
  ##
  if(!is.null(CoVars)){p.CoVars <- nrow(CoVars)}else{p.CoVars <- 0}
  ## t-grid
  t.vec <- seq(a, b, length.out = p) # grid in the domain [a,b]
  
  ## center X?
  if (center.X) {X.mat <- X.mat - rowMeans(X.mat)}

  ## Standardize X?
  if (standardize.X) {
    ## Standardization of X is only for PoI-selection.
    # X.st         <- t(apply(X.mat, MARGIN = 1, FUN = function(x){x/sd(x)}))# old standardization
    X.st <- apply(X.mat, MARGIN = 2, FUN = function(x) { x / stats::sd(x) }) # As suggested in Lindquist & McKeague
  }else{
    X.st <- X.mat
  }
  ##
  if(nbasis.max > 0) {
    ## Generate Fourier basis functions with nbasis.max basis functions:
    Phi <- matrix(0, nrow = p, ncol = nbasis.max)
    for (jj in 1:nbasis.max) {
      Phi[, jj] <- sqrt(2) * sin(jj * pi * t.vec)
    }
    ## Compute Fourier basis coefficients
    FBscores <- matrix(c((1 / p) * t(X.mat) %*% Phi), nrow=N, ncol=nbasis.max)
  }
  ## k (i.e., delta in terms of grid points) sequence for optimizing
  if(is.null(k.seq)){
    ## k.seq <- seq(max(floor(p * 0.01), 1), floor(p * 0.3), by = 4)
    k.seq <- ceiling(p * seq(0.02, 0.40, by = 0.02))
  } 
  if(S.max==0){k.seq <- 1}
  
  ## Progress bar
  if(show.pb){ pb <- utils::txtProgressBar(min = 0, max = length(k.seq), style = 3) }
  
  ############################################## #NEUNEUNEUNEU
  # need list of  k.seq x p matrix!
  if(length(k.seq)>1){
    ind.selected.list<- lapply(seq(nbasis.max+1), function(x) replicate(nbasis.max+1,matrix(0, nrow = length(k.seq), ncol = p))[ , , x])
  }else{
    ind.selected.list<- lapply(seq(nbasis.max+1), function(x) t(as.matrix(replicate(nbasis.max+1,matrix(0, nrow = length(k.seq), ncol = p))[ , , x])))
  } 
  ############################################# #NEUNEUNEUNEU
  ## Container for BIC-values
  IC.mat <- matrix(NA, nrow = length(k.seq), ncol = nbasis.max + 1)
  ##
  for(k in k.seq) { # k <- k.seq[1]
    ## Idea: 
    ## 1. For each k (i.e., delta), compute tau-candidates
    ## 2. For each k force 0:nbasis.max components into the model (approx. functional linear regression part) and save the bic
    ## 3. Choose model with minimal bic, estimate, and save results
    ##
    k.c <- which(k == k.seq) # just a counter - only important if k.seq != 1:k.max
    ## 
    
    for(jj in 0:nbasis.max) { # jj <- 0
      ###
      if(S.max>0){
        ## Compute PoI-candidates (their indices):
        tau.ind.cand <- PoIMaker(k = k, xc = X.st, y = Y, a = a, b = b, plotting = FALSE)
        ## Add the true PoIs to the plot of PoIMaker: axis(side=3,at=t.vec[tau.ind.true])
        ## At most S.max PoI-candidates
        tau.ind.cand <- tau.ind.cand[1:min(length(tau.ind.cand), S.max)]
        ##
        ## Compute X(tau.candidates)
        XX <- X.st[tau.ind.cand,, drop = FALSE] # standardized since only selection - not beta estimation    
        ##
        if (jj == 0) {
          if (is.null(CoVars)){
            Xy <- as.data.frame(cbind(t(XX), Y)) # all as a data-frame 
            colnames(Xy) <- c(paste("TAU", 1:length(tau.ind.cand), sep = ""), "Y") # give standardized 
          } else {
            Xy <- as.data.frame(cbind(t(XX), t(CoVars), Y))
            colnames(Xy) <- c(paste0("TAU", 1:length(tau.ind.cand)), paste0("CoVar", 1:p.CoVars), "Y") # give standardized colnames
          }
        }else{
          if(is.null(CoVars)){
            Xy <- as.data.frame(cbind(t(XX), FBscores[, 1:jj, drop=FALSE], Y)) # all as a data-frame 
            colnames(Xy) <- c(paste0("TAU", 1:length(tau.ind.cand)), paste0("B.fct", 1:jj), "Y") # give standardized colnames
          }else{
            Xy <- as.data.frame(cbind(t(XX), FBscores[, 1:jj, drop=FALSE], t(CoVars), Y))
            colnames(Xy) <- c(paste0("TAU", 1:length(tau.ind.cand)), paste0("B.fct", 1:jj), paste0("CoVar", 1:p.CoVars), "Y") # give standardized colnames
          }
        }
        ## Subset selection only with respect to tau-candidates
        formula.obj <- stats::as.formula(paste("Y~", paste(colnames(Xy)[1:length(tau.ind.cand)], collapse = "+"))) 
        ##
        ## Define 'always.arg', i.e., the variables that are forced in
        if(is.null(CoVars)){
          always.arg <- ifelse(jj > 0, paste(paste0("+B.fct", 1:jj), collapse = ""), "")
        }else{
          always.arg <- ifelse(jj > 0, paste(c(paste0("+CoVar", 1:p.CoVars), paste0("+B.fct", 1:jj)), collapse = ""), paste(paste0("+CoVar", 1:p.CoVars), collapse = ""))
        }
        ##
        glm.force.in.out <- do.call("glmulti",
                                    list(
                                      y           = formula.obj,
                                      data        = Xy,
                                      confsetsize = 50,
                                      maxsize     = S.max + p.CoVars + nbasis.max + 1,
                                      maxK        = S.max + p.CoVars + nbasis.max + 1,
                                      intercept   = TRUE,
                                      fitfunc     = glm.redefined,
                                      always      = always.arg,
                                      crit        = IC,
                                      family      = family,
                                      method      = glmulti.method,
                                      level       = 1, # only main effects (no interactions)
                                      plotty = FALSE, report = FALSE
                                    ))
        ## Save value of information criterion
        if(IC=="bic"){IC.mat[k.c, jj + 1] <- glmulti::weightable(glm.force.in.out)$bic[1]}
        if(IC=="aic"){IC.mat[k.c, jj + 1] <- glmulti::weightable(glm.force.in.out)$aic[1]}
        ##
        tau.ind.cand.slct <- which(paste0("TAU", 1:length(tau.ind.cand)) %in% names(stats::coef(glm.force.in.out@objects[[1]])))
        ind.selected.list[[jj+1]][k.c, tau.ind.cand[tau.ind.cand.slct]]<-1 ##NEUNEUNEU
        ######ind.selected.mat[k.c, tau.ind.cand[tau.ind.cand.slct]] <- 1 ##ALTALT ACHTUNG FEHLER DIE WAHL HIER MUSS AUCH VON JJ ABHAENGEN!!!
        ##
      }# closing 'if(S.max>0)'
      ####
      if(S.max==0){
        if (jj == 0) {
          if(is.null(CoVars)){
            Xy <- as.data.frame(Y) # all as a data-frame 
            colnames(Xy) <- "Y"    # give standardized colnames
          }else{
            Xy <- as.data.frame(cbind(t(CoVars), Y))
            colnames(Xy) <- c(paste0("CoVar", 1:p.CoVars), "Y") # give standardized colnames
          }
        }else{
          if(is.null(CoVars)){
            Xy <- as.data.frame(cbind(FBscores[, 1:jj, drop=FALSE], Y)) # all as a data-frame 
            colnames(Xy) <-  c(paste0("B.fct", 1:jj), "Y")   # give standardized colnames
          }else{
            Xy <- as.data.frame(cbind(t(CoVars), FBscores[, 1:jj, drop=FALSE], Y))
            colnames(Xy) <- c(paste0("CoVar", 1:p.CoVars), paste0("B.fct", 1:jj), "Y") # give standardized colnames
          } 
        }
        ##
        if(IC=="bic"){
          IC.mat[k.c, jj + 1] <- glmulti::bic(stats::glm(Y ~ ., family = family, data = Xy))
        }
        if(IC=="aic"){
          IC.mat[k.c, jj + 1] <- glmulti::aic(stats::glm(Y ~ ., family = family, data = Xy))
        }
        ##
      }# closing 'if(S.max==0)'
    }# closing loop over jj
  }# closing the loop over k
  
  
  ## Select best k and best number of basis functions
  opt.k_j <- which(IC.mat == min(IC.mat), arr.ind = TRUE)
  if(is.matrix(opt.k_j)){opt.k_j <- opt.k_j[1,]}
  k.best  <- opt.k_j[1]     # Index of optimal k (i.e., delta in terms of grid points)
  j.best  <- opt.k_j[2] - 1 # optimal number of Fourier basis functions
  ##
  tau.ind.slct <- c(1:p)[ind.selected.list[[j.best+1]][k.best,,drop=F ] == 1]
  #######################################
  ## Re-fitting the selected glm model 
  ## ####################################
  ## Only Intercept 
  if(length(tau.ind.slct)==0 & is.null(CoVars) & j.best == 0){
    DF              <- data.frame("Y"=Y)
    glm.obj         <- stats::glm(Y ~ 1, family=family, data = DF)
    beta0.hat       <- stats::coef(glm.obj)
    ##
    beta.hat        <- NA
    scores.hat      <- NA
    beta.CoVar.hat  <- NA
  }
  ## PoIs
  if(length(tau.ind.slct)>0 & is.null(CoVars) & j.best == 0){
    DF              <- as.data.frame(cbind(t(X.mat[tau.ind.slct,,drop=FALSE]),Y))
    #print(tau.ind.slct)
    #   colnames(DF)    <- c(paste0("TAU", 1:length(tau.ind.cand)),"Y")
    colnames(DF)    <- c(paste0("TAU", 1:length(tau.ind.slct)),"Y")
    glm.obj         <- stats::glm(Y ~ . , family=family, data = DF)
    beta0.hat       <- stats::coef(glm.obj)[1]
    beta.hat        <- stats::coef(glm.obj)[2:(length(tau.ind.slct)+1)]
    ##
    scores.hat      <- NA
    beta.CoVar.hat  <- NA
  }
  ## PoIs and CoVars
  if(length(tau.ind.slct)>0 & !is.null(CoVars) & j.best == 0){
    DF              <- as.data.frame(cbind(t(X.mat[tau.ind.slct,,drop=FALSE]),t(CoVars),Y))
    colnames(DF)    <- c(paste0("TAU", 1:length(tau.ind.slct)), paste0("CoVar", 1:p.CoVars), "Y")
    glm.obj         <- stats::glm(Y ~ . , family=family, data = DF)
    beta0.hat       <- stats::coef(glm.obj)[1]
    beta.hat        <- stats::coef(glm.obj)[2:(length(tau.ind.slct)+1)]
    beta.CoVar.hat  <- stats::coef(glm.obj)[  (length(tau.ind.slct)+1+1):(length(tau.ind.slct)+1+p.CoVars)]
    ##
    scores.hat      <- NA
  }
  ## FLR
  if(length(tau.ind.slct)==0 & is.null(CoVars) & j.best > 0){
    DF              <- as.data.frame(cbind(FBscores[, 1:j.best, drop=FALSE],Y)) 
    colnames(DF)    <- c(paste0("B.fct", 1:j.best), "Y")
    glm.obj         <- stats::glm(Y ~ . , family=family, data = DF)
    beta0.hat       <- stats::coef(glm.obj)[1]
    scores.hat      <- stats::coef(glm.obj)[2:(j.best+1)]
    ##
    beta.hat        <- NA
    beta.CoVar.hat  <- NA
  }
  ## FLR and CoVars
  if(length(tau.ind.slct)==0 & !is.null(CoVars) & j.best > 0){
    DF              <- as.data.frame(cbind(FBscores[, 1:j.best, drop=FALSE], t(CoVars),Y))
    colnames(DF)    <- c(paste0("B.fct", 1:j.best), paste0("CoVar", 1:p.CoVars), "Y")
    glm.obj         <- stats::glm(Y ~ . , family=family, data = DF)
    beta0.hat       <- stats::coef(glm.obj)[1]
    scores.hat      <- stats::coef(glm.obj)[2:(j.best+1)]
    beta.CoVar.hat  <- stats::coef(glm.obj)[  (j.best+1+1):(j.best+1+p.CoVars)]
    ##
    beta.hat        <- NA
  }
  ## CoVars
  if(length(tau.ind.slct)==0 & !is.null(CoVars) & j.best == 0){
    DF              <- as.data.frame(cbind(t(CoVars),Y))
    colnames(DF)    <- c(paste0("CoVar", 1:p.CoVars), "Y")
    glm.obj         <- stats::glm(Y ~ . , family=family, data = DF)
    beta0.hat       <- stats::coef(glm.obj)[1]
    scores.hat      <- NA
    beta.CoVar.hat  <- stats::coef(glm.obj)[-1]
    ##
    beta.hat        <- NA
  }
  ## PoIs and FLR
  if(length(tau.ind.slct)>0 & is.null(CoVars) & j.best > 0){
    DF              <- as.data.frame(cbind(t(X.mat[tau.ind.slct,,drop=FALSE]), FBscores[, 1:j.best, drop=FALSE],Y))
    colnames(DF)    <- c(paste0("TAU", 1:length(tau.ind.slct)), paste0("B.fct", 1:j.best), "Y")
    glm.obj         <- stats::glm(Y ~ . , family=family, data = DF)
    beta0.hat       <- stats::coef(glm.obj)[1]
    beta.hat        <- stats::coef(glm.obj)[2:(length(tau.ind.slct)+1)]
    scores.hat      <- stats::coef(glm.obj)[  (length(tau.ind.slct)+1+1):(length(tau.ind.slct)+1+j.best)]
    ##
    beta.CoVar.hat <- NA
  }
  ## PoIs, CoVars, and FLR
  if(length(tau.ind.slct)>0 & !is.null(CoVars) & j.best > 0){
    DF              <- as.data.frame(cbind(t(X.mat[tau.ind.slct,,drop=FALSE]), t(CoVars), FBscores[, 1:j.best, drop=FALSE],Y))
    colnames(DF)    <- c(paste0("TAU", 1:length(tau.ind.slct)), paste0("CoVar", 1:p.CoVars), paste0("B.fct", 1:j.best), "Y")
    glm.obj         <- stats::glm(Y ~ . , family=family, data = DF)
    beta0.hat       <- stats::coef(glm.obj)[1]
    beta.hat        <- stats::coef(glm.obj)[2:(length(tau.ind.slct)+1)]
    scores.hat      <- stats::coef(glm.obj)[  (length(tau.ind.slct)+1+1):(length(tau.ind.slct)+1+j.best)]
    beta.CoVar.hat  <- stats::coef(glm.obj)[(length(tau.ind.slct)+1+j.best+1):(length(tau.ind.slct)+1+j.best+p.CoVars)]
  }
  
  ## Beta.fct (FLR)
  if (j.best > 0) {
    beta.fct.hat <- Phi[, 1:j.best, drop = FALSE] %*% scores.hat
    #coef(glm.force.in.out@objects[[1]])[names(coef(glm.force.in.out@objects[[1]])) %in% paste0("B.fct", 1:j.best)]
    #plot(beta.fct.hat)
  } else {
    beta.fct.hat <- rep(0,p)
  }
  ## 
  if (S.max != 0) {
    S.hat <- length(tau.ind.slct)
    k.opt <- k.seq[k.best]
  } else {
    S.hat <- NA
    k.opt <- NA
  }
  ##
  if (length(tau.ind.slct)==0) {
    tau.ind.slct <- NA
  }
  ##
  if(j.best == 0) {
    scores.out <- rep(0,nbasis.max)
  } else {
    scores.out <- rep(0, nbasis.max)
    scores.out[1:j.best]<- scores.hat
  }
  ## ###########################
  ## Return Estimation Result
  ## ###########################
  return(list("tau.ind.hat"    = tau.ind.slct,   # Estimated PoI-locations (indices on grid)
              "beta0.hat"      = beta0.hat,      # Estimated intercept coefficient
              "beta.hat"       = beta.hat,       # Estimated PoI-slope coefs
              "beta.CoVar.hat" = beta.CoVar.hat, # Estimated PoI-slope coefs
              "beta.fct.hat"   = beta.fct.hat,   # Estimated regression function
              "scores.fct.hat" = scores.out,     # Estimated scores 
              "S.hat"          = S.hat,          # Number of slected PoIs
              "k.opt"          = k.opt,          # Optimal k (i.e., optimal delta in terms of grid points)
              "nbasis.opt"     = ifelse(nbasis.max==0, NA, j.best),        # Optimal number of Fourier basis functions
              "glm.obj" = glm.obj # glm object
  ))
}

PoIMaker <- function(k, xc, y, a = 0, b = 1, plotting = TRUE) {
  
  N      <- ncol(xc) # number of functions 
  p      <- nrow(xc) # number of discretization points 
  t.vec  <- seq(a, b, length.out = p) # grid in the domain [a,b]
  xc     <- xc - rowMeans(xc) # making sure that X is centered
  
  delta  <- k * (b - a) / (p - 1) # translate k to delta
  Jdelta <- (k + 1):(p - k) # only 'inner' grid-point in [a,b]
  
  CorXY  <- as.numeric(1 / N * ((xc %*% y))) 
  CorZY  <- CorXY[Jdelta] - 1 / 2 * (CorXY[Jdelta - k] + CorXY[Jdelta + k])
  
  COR     <- abs(CorZY) #Criteria
  COR.aux <- COR
  Jdelta.aux <- Jdelta
  t.delta <- t.vec[Jdelta]
  t.aux   <- t.delta
  taucand <- vector() #container for hat(tau)
  ## Estimate tau
  while (sum(COR.aux > .Machine$double.eps^(1/2)) > 1) {   # COR.aux > 0 (should also work)
    tau.aux <- which.max(COR.aux)
    taucand <- c(taucand, t.aux[tau.aux])
    indi    <- (abs(t.aux - t.aux[tau.aux]) >= sqrt(delta) / 2) # theoretical threshold
    t.aux   <- t.aux[indi]
    COR.aux <- COR[t.delta %in% t.aux]
  }
  if (plotting) {
    graphics::plot(t.delta, COR, xlab = "", ylab = "", type = "l") #main=paste("k =", k , " delta= ", round(delta,3)))
    #abline(v=t.grid[tau.ind.true],col="black",lwd=1)
    graphics::abline(v = taucand, col = "red")
  }
  ## Get the indexes corresponding to the tau-candidates and return them:
  ind.x <- which(t.vec %in% taucand)[rank(taucand)]
  return(ind.x)
}

glm.redefined = function(formula, data, always = "", ...) {
  ## wrapper for glmulti in order to force 
  ## in the first 1:jj basis functions:
  stats::glm(stats::as.formula(paste(deparse(formula), always)), data = data, ...) 
}


#' FUN_PoI_THR
#'
#' This function is an implementation of the threshold estimator THR.
#' @param Y dependent variable (n-vector)
#' @param X.mat pxn matrix of discretized functional predictors
#' @param a lower interval border (default a=0)
#' @param b upper interval border (default b=1)
#' @param scale scaling parameter c_delta (default scale=1.5)
#' @param force.1st.PoI force selection of the first impact point (default force.1st.PoI = FALSE)
#' @param center.X center the functional predictors (default center.X = FALSE)
#' @export FUN_PoI_TRH
FUN_PoI_TRH <- function(Y, X.mat, a = 0, b = 1, 
                        scale          = 1.5, 
                        force.1st.PoI  = FALSE, 
                        center.X       = FALSE){
  ## center X?
  if (center.X) {X.mat <- X.mat - rowMeans(X.mat)}
  ## Dimensions
  N <- ncol(X.mat)
  p <- nrow(X.mat)
  ## 
  delta   <- scale / sqrt(N)
  k.delta <- ceiling(delta * (p - 1) / (b - a))
  out.aux <- PoIMaker(k = k.delta, xc = X.mat, y = Y, a = a, b = b, plotting = FALSE)
  Z       <- X.mat[out.aux,, drop = F] - (1 / 2) * (X.mat[out.aux - k.delta,, drop = F] + X.mat[out.aux + k.delta,, drop = F])
  f_XY    <- abs((1 / N * Z %*% Y) / sqrt(apply(Z, 1, function(x) mean(x ^ 2))))
  ## Compute lambda threshold
  lambda  <- sqrt(2 * sqrt(3)) * sqrt(sqrt(mean(Y ^ 4)) / N * log((b - a) / delta))
  
  # plot(f_XY,ylim=c(0,lambda+0.5))
  # abline(h = lambda)
  
  if (force.1st.PoI == TRUE) {
    X.mat.tau   <- (X.mat[out.aux[1],, drop = F])
    tau.ind.hat <- out.aux[1]
    S.hat       <- 1
  }
  if (force.1st.PoI == FALSE) { 
    ## Compute tau and X(tau) via threshold criterion
    if (any(f_XY > lambda)){# at least one candidate above threshold lambda
      if (all(f_XY > lambda)){# all candidates above threshold lambda
        poi.slct    <- length(f_XY)
        X.mat.tau   <- (X.mat[out.aux[1:poi.slct],, drop = F])
        tau.ind.hat <- out.aux[1:poi.slct]
        S.hat       <- length(1:poi.slct)
      } else {# not all, but some, candidates above threshold lambda
        poi.slct <- min(which(f_XY < lambda)) - 1
        if (poi.slct > 0){
          X.mat.tau   <- (X.mat[out.aux[1:poi.slct],, drop = F])
          tau.ind.hat <- out.aux[1:poi.slct]
          S.hat       <- length(1:poi.slct)
        } else {
          X.mat.tau   <- NULL
          tau.ind.hat <- NA
          S.hat       <- 0
        }
      }
    } else {# no candidate above threshold lambda
      X.mat.tau   <- NULL
      tau.ind.hat <- NA
      S.hat       <- 0
    }
  }
  if (!is.null(X.mat.tau)){
    DF              <- as.data.frame(cbind(t(X.mat.tau),Y))
    colnames(DF)    <- c(paste0("TAU", 1:length(tau.ind.hat)),"Y")
    glm.obj         <- stats::glm(Y ~ . , family=stats::binomial(link = "logit"), data = DF)
    beta0.hat       <- stats::coef(glm.obj)[1]
    beta.hat        <- stats::coef(glm.obj)[2:(length(tau.ind.hat)+1)]
  } else {
    beta0.hat <- NA
    beta.hat  <- NA
  }
  ## ###########################
  ## Return Estimation Result
  ## ###########################
  return(list("tau.ind.hat"    = tau.ind.hat, # Estimated PoI-locations (indices on grid)
              "beta0.hat"      = beta0.hat,   # Estimated intercept coefficient
              "beta.hat"       = beta.hat,    # Estimated PoI-slope coefs
              "S.hat"          = S.hat        # Number of slected PoIs
  ))
}


#' FUN_LMcK
#'
#' This function is an implementation of the estimation procedure described in 
#' 'Logistic Regression With Brownian-Like Predictors' 
#' by Martin Lindquist and Ian McKeague
#' @param Y dependent variable (n-vector)
#' @param X.mat pxn matrix of discretized functional predictors
#' @export FUN_LMcK
FUN_LMcK <- function(Y=Y, X.mat=X.mat){
  log.Lik.vec  <- rep(NA,nrow(X.mat))
  beta.hat.vec <- matrix(NA,ncol=2,nrow=nrow(X.mat))
  for(r in 1:nrow(X.mat)){# slct <- 1
    fit  <- stats::glm(Y~X.mat[r,],family=stats::binomial(link = "logit"))
    log.Lik.vec[r]   <- stats::logLik(fit)
    beta.hat.vec[r,] <- stats::coef(fit)
  }
  tau.ind.slct <- which.max(log.Lik.vec)
  beta.hat     <- beta.hat.vec[tau.ind.slct,]
  ## Return estim result:
  return(list("tau.ind.slct" = tau.ind.slct,
              "beta0.hat"    = beta.hat[1],
              "beta.hat"     = beta.hat[-1]))
}


#' TauIndHat_to_VideoSec
#'
#' This function allows to compute the video time points for each PoI grid-point. 
#' @param TauIndHat Selected PoI indexes from FUN_PoI_BIC() or FUN_PoI_TRH()
#' @export TauIndHat_to_VideoSec
TauIndHat_to_VideoSec <- function(TauIndHat){
  gap_len       <- 5     # orig data was interpolated piecwiese local constant (knot-distance: 5)
  grid_len_orig <- 862   # grid length of original data
  video_len_sec <- 113   # length (in seconds) of the video
  sec_per_ind   <- video_len_sec/grid_len_orig # seconds per grid ('ind-ex') point
  slct_grid     <- seq(1,831,by=gap_len)       # used grid points from the original grid
  ##
  PoI_ind_orig  <- slct_grid[TauIndHat] + 30 # Selected PoI in term of the original grid
  ## start and end seconds for each PoI (interpolation introduced a small blurring)
  PoI_mat       <- cbind(PoI_ind_orig * sec_per_ind, (PoI_ind_orig + gap_len) * sec_per_ind)
  ## Return secs
  return(floor(PoI_mat[,1]))
}


#' BM.SIM 
#'
#' This function allows to simulate Brownian motion processes.
#' @param a lower interval boundary
#' @param b upper interval boundary
#' @param p length of the path
#' @export BM.SIM
BM.SIM <- function(a,b,p){
  dt <-(b-a)/(p-1)
  W  <- cumsum(stats::rnorm(p,0,sqrt(dt)))
  return(W)
}


#' OU.SIM
#'
#' This function allows to simulate Ornstein–Uhlenbeck processes.
#' @param a lower interval boundary
#' @param b upper interval boundary
#' @param p length of the path
#' @param N number of paths
#' @param theta parameter theta
#' @param sigmaou parameter sigma
#' @export OU.SIM
OU.SIM  <- function(a, b, p, N, theta=5, sigmaou=3.5){
  dt <- (b-a)/(p-1)
  mu <- exp(-theta*dt)
  sigmasquared <- sigmaou^2/(2*theta)*(1-mu^2)
  XX <- matrix(0,p,N)
  ZZ <- matrix(stats::rnorm(N*(p-1)),p-1,N)
  for(j in 2:p){
    XX[j,] <- XX[j-1,] * mu + sqrt(sigmasquared) * ZZ[j-1,]
  }
  return(XX)
}


#' Error_Checker
#'
#' This function allows to check for error messages from try().
#' @param x value from try()
#' @export Error_Checker
Error_Checker <- function(x){
  bool.result <- inherits(x, "try-error")
  bool.result <- unname(bool.result)
  return(bool.result)
}

## result <- try(log("a"), silent=TRUE); Error_Checker(result)


#' Allocator
#'
#' Allocates the estimated PoIs to the true PoIs 
#' @param tau.ind.hat indices of estimated PoIs
#' @param beta.hat empirical slope coefficients
#' @param tau.ind.true indices of true PoIs
#' @param p number of discretization points
#' @export Allocator
## ##########################################
## Allocating estimated PoIs to true PoIs  ##
## ##########################################
Allocator <- function(tau.ind.hat, beta.hat, tau.ind.true = tau.ind.true, p = p){
  ##
  S.hat <- length(tau.ind.hat)
  S     <- length(tau.ind.true)
  ## 
  if(any(diff(tau.ind.true)<=0)){stop("tau.ind.true must be strictly increasing.")}
  if(length(tau.ind.hat)!=length(beta.hat)){stop("tau.ind.hat and beta.hat must have equal lengths")}
  ##
  tau.ind.hat.alloc <- rep(NA, S)
  beta.hat.alloc    <- rep(NA, S)
  ##
  if ( all(!is.na(c(tau.ind.hat,beta.hat))) ) {# Only if there are no NAs
    ##
    if(S==1){
      ind.aux           <- which.min(abs(tau.ind.hat - tau.ind.true))
      tau.ind.hat.alloc <- tau.ind.hat[ind.aux]
      beta.hat.alloc    <- beta.hat[ind.aux]
    }
    ##
    if(S > 1){
      ## define intervals for allocating the estimated PoI candidates (one interval for each true PoI): 
      vbin     <- c(0, sapply(2:length(tau.ind.true), function(x) 0.5 * (tau.ind.true[x] + tau.ind.true[x - 1])), p)
      ## 1st column: empirical PoI indices, 2nd column corresponding allocation intervals
      aux.tau  <- cbind(tau.ind.hat, findInterval(tau.ind.hat, vbin, all.inside = T))
      ## 1st column: slope-estimates of PoIs, 2nd column corresponding allocation intervals
      aux.beta <- cbind(beta.hat,    findInterval(tau.ind.hat, vbin, all.inside = T))
      ## Allocating the empirical PoIs and their slope-estimates to the 1st, 2nd, ... allocation interval
      for (vv in 1:S)
        if (any(aux.tau[, 2] == vv)) {#vv=1
          ind.aux               <- which.min(abs(aux.tau[which(aux.tau[, 2] == vv), 1] - tau.ind.true[vv]))
          tau.ind.hat.alloc[vv] <- aux.tau[ which(aux.tau[ , 2] == vv), 1][ind.aux]
          beta.hat.alloc[vv]    <- aux.beta[which(aux.beta[, 2] == vv), 1][ind.aux]
        }
    }
  }
  return(list("tau.ind.hat.alloc" = tau.ind.hat.alloc,
              "beta.hat.alloc"    = beta.hat.alloc))
}



Allocator2 <- function(tau.ind.hat, tau.ind.true = tau.ind.true, p = p){
  ##
  S.hat <- length(tau.ind.hat)
  S     <- length(tau.ind.true)
  ## 
  if(any(diff(tau.ind.true)<=0)){stop("tau.ind.true must be strictly increasing.")}
  ##
  tau.ind.hat.alloc <- rep(NA, S)
  ##
  if ( all(!is.na(tau.ind.hat)) ) {# Only if there are no NAs
    ##
    if(S==1){
      ind.aux           <- which.min(abs(tau.ind.hat - tau.ind.true))
      tau.ind.hat.alloc <- tau.ind.hat[ind.aux]
    }
    ##
    if(S > 1){
      ## define intervals for allocating the estimated PoI candidates (one interval for each true PoI): 
      vbin     <- c(0, sapply(2:length(tau.ind.true), function(x) 0.5 * (tau.ind.true[x] + tau.ind.true[x - 1])), p)
      ## 1st column: empirical PoI indices, 2nd column corresponding allocation intervals
      aux.tau  <- cbind(tau.ind.hat, findInterval(tau.ind.hat, vbin, all.inside = T))
      ## Allocating the empirical PoIs and their slope-estimates to the 1st, 2nd, ... allocation interval
      for (vv in 1:S)
        if (any(aux.tau[, 2] == vv)) {#vv=1
          ind.aux               <- which.min(abs(aux.tau[which(aux.tau[, 2] == vv), 1] - tau.ind.true[vv]))
          tau.ind.hat.alloc[vv] <- aux.tau[ which(aux.tau[ , 2] == vv), 1][ind.aux]
        }
    }
  }
  return(list("tau.ind.hat.alloc" = tau.ind.hat.alloc))
}



#' Video Rating Data
#'
#' A dataset containing the emotion ratings reported from n=65 participants 
#' while watching an affective online video on the persecution of African albinos. 
#' A version of the video can be found online at YouTube (www.youtube.com). 
#' Link to the video: https://youtu.be/9F6UpuJIFaY 
#' The video clip used in the experiment corresponds to the first 112 sec. 
#' The first six data points (<1 sec.) are removed from the trajectories as they 
#' contain some obviously erratic components.
#'
#' \itemize{
#'   \item Y: Final overall rating (Y=0 denotes 'I feel negative' and Y=1 denotes 'I feel not negative')
#'   \item X_1 to X_167: Continuous emotion ratings (discretization points 1 to 167) 
#' }
#'
#' @docType data
#' @keywords datasets
#' @name emotion_rating
#' @usage data(emotion_rating)
#' @format A data frame with 65 rows and 168 variables
"emotion_rating"




