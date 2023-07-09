#########################################################
#
################ START SUB-ROUTINES #####################
#
#########################################################
# N.B. we should avoid to use "gam" insetad of "gamma"
# for the gamma paratemer since that gamma is also a
# function
### A) CENSORED LOG-LIKELIHOOD BASED ON THE GEV DISTRIBUTION
### A.1) DEFINE THE GEV CDF
gevcdf.root <- function(x, mu, sigma, gam)
{
  # define a zero neighborhood
  eps <- .Machine$double.eps^.3
  # set the three domains
  # GUMBEL CASE
  if(abs(gam)<=eps) return (exp(-exp(-(x-mu)/sigma)))
  else{
    # FRECHET CASE
    if(gam>0){
      if(x<=(mu-sigma/gam)) return(0)
      else return(exp(-(1+gam*(x-mu)/sigma)^(-1/gam)))
    }
    # REVERSE-WEIBULL CASE
    if(gam<0){
      if(x<=(mu-sigma/gam)) return(exp(-(1+gam*(x-mu)/sigma)^(-1/gam)))
      else return(1)
    }
  }
}
### A.2) DEFINE THE GEV PDF
gevpdf.root <- function(x, mu, sigma, gam)
{
  # define a zero neighborhood
  eps <- .Machine$double.eps^.3
  # set the three domains
  # GUMBEL CASE
  if(abs(gam)<=eps){
    z <- (x-mu)/sigma
    return (exp(-exp(-z)-z)/sigma)
  }
  else{
    # FRECHET CASE
    if(gam>0){
      if(x<=(mu-sigma/gam)) return(0)
      else{
        z <- 1 + gam*(x-mu)/sigma
        return(1/sigma * exp(-z^(-1/gam)) * z^(-1/gam-1))
      }
    }
    # REVERSE-WEIBULL CASE
    if(gam<0){
      if(x<=(mu-sigma/gam)){
        z <- 1 + gam*(x-mu)/sigma
        return(1/sigma * exp(-z^(-1/gam)) * z^(-1/gam-1))
      }
      else return(0)
    }
  }
}
gevcdf <- Vectorize(gevcdf.root, "x")
gevpdf <- Vectorize(gevpdf.root, "x")
# A.3) DEFINE THE CENSORED LOG-LIKELIHOOD GEV FUNCTION
gev.cens.llik <- function(param, data, th){
  # return value when something goes wrong
  LOW <- -1e300
  # parameters definition
  mu <- param[1]
  sigma <- param[2]
  gam <- param[3]
  # scale parameter check
  if(sigma<0){return(LOW)}
  # derive the exceedances based on the threhold "t"
  I <- data>th
  exceed <- data[I]
  # set sample size and effective sample sizes:
  n <- length(data)
  k <- sum(I)
  # new parametrization
  mu <- mu - sigma * (1 - (k/n)^gam) / gam
  sigma <- sigma * (k/n)^(gam)
  # log-likelihood function concerning the exceedances
  up <- sum(log(gevpdf(exceed, mu, sigma, gam)))
  # log-likelihood concerning the small values
  low <- (n-k) * log(gevcdf(th, mu, sigma, gam))
  # define the ultimate censored log-likelihood
  res <- up+low
  # utilmate check in the case the log-likelihood in -Inf
  if(is.infinite(res))return(-1e300)
  return(res)
}
###    MCMC algorithm to derive an approximate posterior distribution
###    for the parameters (location, scale and shape) of the GEV distribution
###    The MCMC algirithm is a Random Walk Metropolis-Hastings
RWMH.CenLik <- function(data, start, th, p, sig0, nsim, mle, prior){
  # A) INITIALIZATION
  #start <- par0
  #p <- etastar <- 0.234
  #alpha <- zeta0
  alpha  <- -qnorm(p/2)
  d <- length(start) # Dimension of the vector of parameters
  #if(!any(d==3, d==4)){stop("Wrong length of parameter vector")}
  if(!any(d==3)){stop("Wrong length of parameter vector")}
  #
  sig <- sig0 # Initial value of sigma
  sig.vec <- sig      #RETURN
  #
  #sigma <- sigMat
  sigMat <- diag(d)  # Initial proposal covarince matrix
  acc.vec <- rep(NaN, nsim) # Vector of acceptances
  accepted <- rep(0, nsim)
  straight.reject <- rep(0, nsim) # Monitor the number of proposed parameters that don't respect the constraints
  sig.start<- sig
  sig.restart<- sig       #RETURN
  #
  ####etastar!!!
  n0 <- round(5/(p*(1-p)))
  iMax <- 100 # iMax is the max number of iterations before the last restart
  Numbig <- 0
  Numsmall <- 0
  #
  param <- start
  msg <- "none"
  # sample size
  n <- length(data)
  # A) LIKELIHOOD DEFNITION, SET PARAMETERS CONSTRAINTS AND ACCEPTANCE PROBABILITY DEFINITION
  llikfun <- function(para) return(gev.cens.llik(para, data=data, th=th))
  parcheck <- function(para){
    res <- any(para[2]<=0,
              (para[3]>0 & (para[1]>min(data)+para[2]/para[3])),
              (para[3]<0 & any(para[1]<=data+para[2]/para[3])))
    return(res)
  }
  if(prior=="uniform"){
    accept_prob <- function(par, par_new){
      res <- min(exp(llikfun(par_new)-llikfun(par)+log(par[2])-log(par_new[2])), 1)
      return(res)
    }
  }
  if(prior=="empirical"){
    accept_prob <- function(par, par_new){
      scale <- mle[2]
      res <- min(exp(llikfun(par_new)-llikfun(par)+
                     dnorm(par_new[1], mle[1], scale^2, TRUE)-
                     dnorm(par[1], mle[1], scale^2, TRUE)+
                     dgamma(par_new[2], 1, scale = scale, log=TRUE)-
                     dgamma(par[2], 1, scale = scale, log=TRUE)+
                     dtt(par_new[3], df=1, left = -1, right=Inf, log=TRUE)-
                     dtt(par[3], df=1, left = -1, right=Inf, log=TRUE)), 1)
      return(res)
    }
  }
  # B) START ALGORITHM
  # set the chains
  sparam <- param
  # create a progress bar
  ###PROGRESS BAR!!!
  pb <- txtProgressBar(min = 0, max = nsim, initial = 0, style=3)
  # main algorithm loop
  for(i in 1:nsim){
    # simulation from the proposal (MVN):
    s <-  try(eigen(sig^2 * sigMat, symmetric=TRUE), silent=TRUE)
    if(class(s)=="try-error"){
      msg <- "unfortunately the algorihtm has failed"
      return(list(param_post=sparam,
                  accepted=accepted,
                  straight.reject=straight.reject,
                  nsim=nsim,
                  sig.vec=sig.vec,
                  sig.restart=sig.restart,
                  acc.vec=acc.vec,
                  msg=msg))
    }
    param_new <- as.vector(mvtnorm:::rmvnorm(1, mean = param, sigma = sig^2 * sigMat))
    if(any(is.na(param_new))){
      straight.reject[i] <- 1
      acc.vec[i] <- 0
    }
    else{
    # check basic condition
    if(parcheck(param_new)){
      straight.reject[i] <- 1
      acc.vec[i] <- 0
    }else{
      # compute the acceptance probability
      ratio <- accept_prob(param, param_new)
      if(is.na(ratio)){
        straight.reject[i] <- 1
        acc.vec[i] <- 0
      }
      else{
      acc.vec[i] <- ratio
      u <- runif(1)
      if(u<ratio){
        param <- param_new
        accepted[i] <- 1
      }}
    }}
    sparam <- rbind(sparam, param)
    # update covariance matrix with adaptive MCMC
    if(i>100) {
      if(i==101) {
        sigMat <- cov(sparam)
        thetaM <- apply(sparam, 2, mean)
      } else
      {
        tmp <- update.cov(sigMat = sigMat, i = i, thetaM = thetaM, theta = param, d = d)
        sigMat <- tmp$sigMat
        thetaM <- tmp$thetaM
      }
    }
    # Update sigma
    #n0 <- round(5/(etastar*(1-etastar)))
    if (i>n0) {
      #p <- etastar
      #alpha <- zeta0 <- -qnorm(etastar/2)
      sig <- update.sig(sig = sig, acc = ratio, d = d, p = p, alpha = alpha, i = i)
      #k.vec
      sig.vec <- c(sig.vec, sig)
      if ((i <= (iMax+n0)) && (Numbig<5 || Numsmall<5)) {
        Toobig <- (sig > (3*sig.start))
        Toosmall <-(sig < (sig.start/3))#
        if (Toobig || Toosmall) {
          #restart the algorithm
          message("restart the program at", i, "th iteration", "\n")
          #k.restart
          sig.restart <- c(sig.restart, sig)
          Numbig <- Numbig + Toobig
          Numsmall <- Numsmall + Toosmall
          i <- n0
          sig.start <- sig
        }
      } #end iMax
    }
    print.i <- seq(0, nsim, by=100)
    if(i %in% print.i) setTxtProgressBar(pb, i)
  }
  return(list(param_post=sparam,
              accepted=accepted,
              straight.reject=straight.reject,
              nsim=nsim,
              sig.vec=sig.vec,
              sig.restart=sig.restart,
              acc.vec=acc.vec,
              msg=msg))
}
#

# function to adjust the value of sigma (=log(theta)) in each iteration
update.sig <- function(sig, acc, d = d, p = p, alpha = alpha, i) {
  #alpha <- zeta0 <- -qnorm(etastar/2)
  #a <- ((2*pi)^(1/2))*exp((zeta0^2)/2)/(2*zeta0)
  c <- ((1-1/d) * sqrt(2*pi) * exp(alpha^2 / 2) / (2 * alpha) + 1 / (d*p*(1-p))) # Eq(7) of Garthwaite, Fan & Sisson (2016)
  Theta <- log(sig)
  # Theta = Theta + c * (acc-p) / max(200, i/d)
  Theta <- Theta + c * (acc-p) / i
  return(exp(Theta))
}
# function to update the covariance matrix
update.cov <- function(sigMat, i, thetaM, theta, d){
  epsilon <- 1/i
  thetaM2 <- ((thetaM*i)+theta)/(i+1)
  sigMat <- (i-1)/i*sigMat + thetaM%*%t(thetaM)-(i+1)/i*thetaM2%*%t(thetaM2)+1/i*theta%*%t(theta) + epsilon*diag(d)
  return(list(sigMat=sigMat, thetaM=thetaM2))
}


#########################################################
#
################ END SUB-ROUTINES #######################
#
#########################################################

#########################################################
#
################ START MAIN CODE ########################
#
#########################################################
#
# A.4) MAIN CODE FITTING PROCEDURE
FitCensLLik<- function(data, th=NULL, param0=NULL, llik.type="Gev-Cens", p=1/length(data),
                       T=50, t.prob=0.9, method="frequentist", optim.meth="Nelder-Mead",
                       control=NULL, sig0=NULL, nsim=NULL, burn=NULL, p.acc=0.234,
                       prior="empirical", val.show=FALSE, ...){
  # check whether starting values are provided
  if(is.null(param0)){stop("Need to specify param0 parameter")}
  # set the threshold if not provided
  if(is.null(th)){
    th <- quantile(data, probs=t.prob, type=3)
    message("threshold set as an empirical quantile of the data by default \n")
  }
  # DEFINE THE PROPORTION OF EXCEEDANCES
  exc.prop <- sum(data>th)/length(data)
  # BAYESIAN INFERENCE
  if(method=="bayesian"){
    # check the minimal required parameters
   if(is.null(sig0) || is.null(nsim) || is.null(p.acc)){stop("Missing initial values")}
   if(is.null(nsim)){stop("Missing number of replimessagees in algorithm")}
   # set the number of quantiles to extrapolate
   nQuant <- length(p)
   # ML ESTIMATION OF THE GEV PARAMETERS
   if(llik.type=="Gev-Cens"){
     optimfun <- function(par){
       return(gev.cens.llik(par, data=data, th=th))
   }}
   mle <- optim(param0, optimfun, control=list(fnscale=-1, maxit = 8e+5))$par
   # compute the expected information matrix at the ML estimates
   # MCMC ALGORITHM
   mcmc <- RWMH.CenLik(data=data, start=param0, p=p.acc, sig0=sig0, nsim=nsim, mle=mle,
                       prior=prior, th=th)

   meanacc <- rep(NA, length(mcmc$acc.vec))
   for(j in c(1:length(mcmc$acc.vec))){
     meanacc[j] <- mean(mcmc$acc.vec[round(j/2) :j])
   }
   index <- c(1:2000, seq(2001,length(mcmc$sig.vec), by=100))
   if(val.show){
     oldpar <- par(mar=c(4.4,4.8,0.5,0.5))
     on.exit(par(oldpar))
     par(mfrow=c(1,2))
     plot( cbind(index, mcmc$sig.vec[index]^2) , type="l", col=3, ylim=c(0, max(mcmc$sig.vec^2)), ylab=expression(tau^2), xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
     abline(h=0, lwd=2)
     plot(cbind(index, meanacc[index]), type="l", col=2, ylim=c(0,1), ylab="Acceptance Probability", xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
     abline(h=p.acc, lwd=2)
     # deduce approximate samples from the parameter's posterior distributions
   }
   post_sample <- mcmc$param_post[(burn+1):nsim,]
   mu.ps <-  post_sample[,1]
   sig.ps <- post_sample[,2]
   gam.ps <- post_sample[,3]
   # DEFINE THE EXTREME QUANTILE
   Q.ext <- matrix(NaN, nrow=nrow(post_sample), ncol=nQuant)
   for(j in 1:nQuant){
     Q.ext[,j] <- mu.ps + sig.ps * ((exc.prop/p[j])^gam.ps-1) / gam.ps
   }
   # DEFINE THE RETURN LEVEL
   R.lev <- mu.ps + sig.ps * ((-log(1-1/T))^(-gam.ps)-1) / gam.ps
   # DEFINE OUTPUT
   return(list(Q.extreme=Q.ext,
               R.level=R.lev,
               post_sample=post_sample,
               burn=burn,
               straight.reject=mcmc$straight.reject[-c(1:(burn-1))],
               sig.vec=mcmc$sig.vec,
               accept.prob=cbind(index, meanacc[index]),
               msg=mcmc$msg,
               mle=mle,
               exc.prop=exc.prop))
   # FREQUENTIST INFERENCE
   } else if(method=="frequentist"){
    # set the optimisation method if not provided
    if(is.null(optim.meth)){stop("Need to specify an optimisation method in frequentist setting")}
    # set further checks
    if(!("control" %in% names(sys.call())) ){
      control <- list(fnscale=-1, maxit = 8e+5)
    }else{
      if(!("fnscale" %in% names(control))){
        control[["fnscale"]] <- -1
        control[["maxit"]] <- 8e+5
      }
    }
    # outcome setting
    optim.msg <- "something went wrong"
    if(llik.type=="Gev-Cens"){
      optimfun <- function(par){
        return(gev.cens.llik(par, data=data, th=th))
      }
    }
    names(exc.prop) <- c("Excesses proportion")
    ### COMPUTE MAXIMUM CENSORED LOG-LIKELIHOOD ESTIMATES
    est <- optim(param0, optimfun, method=optim.meth, control=control)
    ### DEFINE THE ESTIMATED GEV PARAMETERS
    #scale <- est$par[2] * exc.prop^(-est$par[3])
    #loc <- est$par[1] + scale * (1 - exc.prop^est$par[3]) / est$par[3]
    #est$par[2] <- scale
    #est$par[1] <- loc
    loc <- est$par[1]
    scale <- est$par[2]
    names(est$par) <- c("location", "scale", "shape")
    optim.msg <- est$convergence
    names(optim.msg) <- c("optimization message")
    # DEFINE THE EXTREME QUANTILE
    Q.ext <- loc + scale * ((exc.prop/p)^est$par[3]-1) / est$par[3]
    names(Q.ext) <- c("Extreme-Quantile")
    # DEFINE THE RETURN LEVEL
    R.lev <- loc + scale * ((-log(1-1/T))^(-est$par[3])-1) / est$par[3]
    names(R.lev) <- c("Return-Level")
    return(list(estimate=est$par,
                Q.extreme=Q.ext,
                R.level=R.lev,
                optim.msg=optim.msg,
                exc.prop=exc.prop))
    }
}
