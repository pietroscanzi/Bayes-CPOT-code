#GEV distribution

library(evd)

#-------------------------------------------------------------------------------

#Probability density function

dgev1 <- function(x, par, log = FALSE){
  gam <- par[1]; mu <- par[2]; delta <- par[3]
  z <- (x - mu)/delta
  #Define 0 neighbourhood
  eps <- .Machine$double.eps^.3
  #the smallest positive floating-point number x such that 1 + x != 1
  pdf <- 0
  if(abs(gam) <= eps){
    pdf <- (1/delta)*exp(-(z + exp(-z)))
  }
  else{
    if((1 + gam*z) > 0){
      gev.k <- (1 + gam*z)^(-1/gam)
      pdf <- (1/delta)*(gev.k^(gam + 1))*exp(-gev.k)
    }
  }
  if(log == FALSE) return(pdf)
  if(log == TRUE) return(log(pdf))
}
dgev.V <- Vectorize(dgev1, "x")

#Comparison with evd functions
dgev1(1, c(0,2,4)); dgev1(1, c(0,2,4), log = TRUE)
dgev(1, loc = 2, scale = 4, shape = 0); dgev(1, loc = 2, scale = 4, shape = 0, log = TRUE)
dgev1(1, c(1,-4,6)); dgev1(1, c(1,-4,6), log = TRUE)
dgev(1, loc = -4, scale = 6, shape = 1); dgev(1, loc = -4, scale = 6, shape = 1, log = TRUE)
dgev1(-1, c(-1,0,1)); dgev1(-1, c(-1,0,1), log = TRUE)
dgev(-1, loc = 0, scale = 1, shape = -1); dgev(-1, loc = 0, scale = 1, shape = -1, log = TRUE)

#Comparison plots (domains of attraction: Gumbel, Frechet, Weibull)
#mu = 0, delta = 1, --> gam = c(-1, 0, 1)
plot(function(x) dgev.V(x, c(0, 0, 1)), -5, 5, xlab = "x", ylab = "pdf",
     main = bquote(mu == 0 ~","~ delta == 1), col = 2, lwd = 2,
     ylim = c(0, 1))
plot(function(x) dgev.V(x, c(1, 0, 1)), -5, 5, col = 3, lwd = 2, add = T)
plot(function(x) dgev.V(x, c(-1, 0, 1)), -5, 5, col = 4, lwd = 2, add = T)
legend("topleft", bty = "n", legend = c(expression(gamma  == 0),
                              expression(gamma == 1),
                              expression(gamma == -1)),
       col = c(2,3,4), lty = 1, lwd = 2)
#mu = 0, delta = 1, --> gamma = c(-1/2, 0, 1/2)
plot(function(x) dgev.V(x, c(0, 0, 1)), -5, 5, xlab = "x", ylab = "pdf",
     main = bquote(mu == 0 ~","~ delta == 1), col = 2, lwd = 2,
     ylim = c(0, 0.5))
plot(function(x) dgev.V(x, c(1/2, 0, 1)), -5, 5, col = 3, lwd = 2, add = T)
plot(function(x) dgev.V(x, c(-1/2, 0, 1)), -5, 5, col = 4, lwd = 2, add = T)
legend("topleft", bty = "n", legend = c(expression(gamma  == 0),
                                        expression(gamma == 1/2),
                                        expression(gamma == -1/2)),
       col = c(2,3,4), lty = 1, lwd = 2)

#-------------------------------------------------------------------------------

#Cumulative density function

pgev1 <- function(x, par, log = FALSE){
  gam <- par[1]; mu <- par[2]; delta <- par[3]
  z <- (x - mu)/delta
  #Define 0 neighbourhood
  eps <- .Machine$double.eps^.3
  #the smallest positive floating-point number x such that 1 + x != 1
  if(abs(gam) <= eps){
    cdf <- exp(-exp(-z))
  }
  else{
    gev.k <- (1 + gam*z)^(-1/gam)
    cdf <- exp(-gev.k)
    if(gam > 0 & ((1 + gam*z) <= 0)) cdf = 0
    if(gam < 0 & ((1 + gam*z) <= 0)) cdf = 1
  }
  if(log == FALSE) return(cdf)
  if(log == TRUE) return(log(cdf))
}
pgev.V <- Vectorize(pgev1, "x")

#Comparison with evd functions
pgev1(1, c(0,2,4)); pgev1(1, c(0,2,4), log = TRUE)
pgev(1, loc = 2, scale = 4, shape = 0)
pgev1(1, c(1,-4,6)); pgev1(1, c(1,-4,6), log = TRUE)
pgev(1, loc = -4, scale = 6, shape = 1)
pgev1(-1, c(-1,0,1)); pgev1(-1, c(-1,0,1), log = TRUE)
pgev(-1, loc = 0, scale = 1, shape = -1)

#Comparison plots (domains of attraction: Gumbel, Frechet, Weibull)
#mu = 1, delta = 1, --> gamma = c(-1, 0, 1)
plot(function(x) pgev.V(x, c(0, 1, 1)), -5, 5, xlab = "x", ylab = "cdf",
     main = bquote(mu == 1 ~","~ delta == 1), col = 2, lwd = 2,
     ylim = c(0, 1))
plot(function(x) pgev.V(x, c(1, 1, 1)), -5, 5, col = 3, lwd = 2, add = T)
plot(function(x) pgev.V(x, c(-1, 1, 1)), -5, 5, col = 4, lwd = 2, add = T)
legend("topleft", bty = "n", legend = c(expression(gamma  == 0),
                                        expression(gamma == 1),
                                        expression(gamma == -1)),
       col = c(2,3,4), lty = 1, lwd = 2)
#mu = 1, delta = 1, --> gamma = c(-1/2, 0, 1/2)
plot(function(x) pgev.V(x, c(0, 1, 1)), -5, 5, xlab = "x", ylab = "pdf",
     main = bquote(mu == 0 ~","~ delta == 1), col = 2, lwd = 2,
     ylim = c(0, 1))
plot(function(x) pgev.V(x, c(1/2, 1, 1)), -5, 5, col = 3, lwd = 2, add = T)
plot(function(x) pgev.V(x, c(-1/2, 1, 1)), -5, 5, col = 4, lwd = 2, add = T)
legend("topleft", bty = "n", legend = c(expression(gamma  == 0),
                                        expression(gamma == 1/2),
                                        expression(gamma == -1/2)),
       col = c(2,3,4), lty = 1, lwd = 2)

#-------------------------------------------------------------------------------

#Example

set.seed(1)
sam.exp <- rexp(100, rate = 1)
#Gumbel
dgev.V(sam.exp, c(0, 0, 1)); prod(dgev.V(sam.exp, c(0, 0, 1)))
dgev(sam.exp, loc = 0, scale = 1, shape = 0); prod(dgev(sam.exp, loc = 0, scale = 1, shape = 0))
dgev.V(sam.exp, c(0, 0, 1), log = T); sum(dgev.V(sam.exp, c(0, 0, 1), log = T))
dgev(sam.exp, loc = 0, scale = 1, shape = 0, log = T); sum(dgev(sam.exp, loc = 0, scale = 1, shape = 0, log = T))
#Frechet
dgev.V(sam.exp, c(1, 0, 1)); prod(dgev.V(sam.exp, c(1, 0, 1)))
dgev(sam.exp, loc = 0, scale = 1, shape = 1); prod(dgev(sam.exp, loc = 0, scale = 1, shape = 1))
dgev.V(sam.exp, c(1, 0, 1), log = T); sum(dgev.V(sam.exp, c(1, 0, 1), log = T))
dgev(sam.exp, loc = 0, scale = 1, shape = 1, log = T); sum(dgev(sam.exp, loc = 0, scale = 1, shape = 1, log = T))
#Weibull
dgev.V(sam.exp, c(-1, 0, 1)); prod(dgev.V(sam.exp, c(-1, 0, 1)))
dgev(sam.exp, loc = 0, scale = 1, shape = -1); prod(dgev(sam.exp, loc = 0, scale = 1, shape = -1))
dgev.V(sam.exp, c(-1, 0, 1), log = T); sum(dgev.V(sam.exp, c(-1, 0, 1), log = T))
dgev(sam.exp, loc = 0, scale = 1, shape = -1, log = T); sum(dgev(sam.exp, loc = 0, scale = 1, shape = -1, log = T))

#-------------------------------------------------------------------------------

#Censored GEV likelihood

gev.lik <- function(par, data, t = NULL, p = 0.95, log = FALSE){
  
  #Set the smallest value of the log-likelihood
  Low <- -1e300
  #Check the validity of the data
  if(!is.vector(data)){stop("Data must be a vector \n")}
  #Check the validity of the parameters
  if(length(par) != 3){stop("Wrong length of parameter vector \n")}
  #Check the support of the scale parameter delta
  if(par[3] <= 0){
    if(log == FALSE) return(0)
    if(log == TRUE) return(Low)
  }
  
  #Empirical threshold
  #Use discontinuous sample 95% quantile: type = 3
  if(is.null(t)){
    t <- as.numeric(quantile(data, probs = p, type = 3))
    #message("T set to 95% quantile by default \n")
  } 
  #Indicator of censoring
  cond <- as.numeric(data > t)
  
  #Likelihood
  if(log == FALSE){
    lik.cens <- (pgev.V(t, par, log = FALSE))^(sum(cond == 0))
    lik.obs <- prod(dgev.V(data[cond == 1], par, log = FALSE))
    lik <- prod(lik.cens, lik.obs)
    if(is.infinite(lik)) return(0)
    return(lik)
  }
  
  #Log-likelihood
  if(log == TRUE){
    llik.cens <- (sum(cond == 0))*(pgev.V(t, par, log = TRUE))
    llik.obs <- sum(dgev.V(data[cond == 1], par, log = TRUE))
    llik <- sum(llik.cens, llik.obs)
    if(is.infinite(llik)) return(Low)
    return(llik)
  }
}
#gev.lik <- Vectorize(gev.lik, "par")

#Using evd functions
gev.lik.evd <- function(par, data, log = FALSE){
  
  #Empirical threshold
  t <- as.numeric(quantile(data, 0.95, type = 3))
  #Indicator of censoring
  cond <- as.numeric(data > t)
  
  #Likelihood
  if(log == FALSE){
    lik.cens <- (pgev(t, loc = par[2], scale = par[3], shape = par[1]))^(sum(cond == 0))
    lik.obs <- prod(dgev(data[cond == 1], loc = par[2], scale = par[3], shape = par[1]))
    return(prod(lik.cens, lik.obs))
  }
  
  #Log-likelihood
  if(log == TRUE){
    lik.cens <- (sum(cond == 0))*log(pgev(t, loc = par[2], scale = par[3],
                                          shape = par[1]))
    lik.obs <- sum(dgev(data[cond == 1], loc = par[2], scale = par[3],
                        shape = par[1], log = TRUE))
    return(sum(lik.cens, lik.obs))
  }
}

#-------------------------------------------------------------------------------

#Data

n <- 1000
#Gumbel domain: gamma sample
set.seed(1)
sam.gum <- rgamma(n, shape = 2, rate = 2)
hist(sam.gum, nclass = 50, probability = T, main = "Gumbel domain")
#Frechet domain: half cauchy sample
set.seed(2)
library(LaplacesDemon)
sam.fre <- rhalfcauchy(n, scale = 1)
hist(sam.fre, nclass = 50, probability = T, main = "Frechet domain") 
summary(sam.fre)
#Weibull domain: reverse Weibull sample
set.seed(3)
#evd::rrweibull
sam.wei <- rrweibull(n, shape = 3)
hist(sam.wei, nclass = 50, probability = T, main = "Weibull domain")

#-------------------------------------------------------------------------------

#Optimization: MLE /Bayesian Mode (empirical bayes prior)

fit.gev.inference <- function(data, t = NULL, t.prob = 0.95, llik.type = "Gev-Cens",
                              T.ret = 50, p = 1/length(data),
                              par0 = NULL, optim.meth = "Nelder-Mead", control = NULL,
                              hessian = FALSE,
                              inf.type, 
                              #inf.type = c("Frequent", "Bayes")
                              ...){
  
  #Setting of the starting value of the optimizator
  if(is.null(par0)) {stop("Need to specify a starting value for the optimization")} 
  #Empirical threshold
  #Use discontinuous sample 95% quantile: type = 3
  if(is.null(t)){
    t <- as.numeric(quantile(data, probs = t.prob, type = 3))
    message(paste0("T set to the ", round(100*t.prob, 0), "% quantile by default"))
  }
  #Optimization method: 
  #c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")
  if(is.null(optim.meth)) {stop("Need to specify an existing optimisation method!")}
  #Flow control parameters
  if(!("control" %in% names(sys.call()))){
    #Maximization
    control <- list(fnscale = -1, maxit = 8e+5)
  }
  else{
    if(!("fnscale" %in% names(control))){
      #Maximization
      control[["fnscale"]] <- -1
      control[["maxit"]] <- 8e+5
    }
  }
  #Define the proportion of exceedances
  #Sample size
  n <- length(data)
  #Effective sample size
  k <- sum(data > t)
  exc.prop <- k/n
  names(exc.prop) <- c("Excesses proportion")
  #s = Reciprocal of proportion of exceedances
  s <- 1/exc.prop
  
  #Frequentist inference setting
  if(inf.type == "Frequent"){
    #Outcome setting
    optim.msg <- "Something went wrong"
    if(llik.type == "Gev-Cens"){
      optimfun <- function(par){
        return(gev.lik(par = par, data = data, t = t, p = t.prob, log = TRUE))
      }
    }
    #Compute the maximum consered likelihood estimates
    #Hessian = TRUE preferred for minimizations
    est <- optim(par0, optimfun, method = optim.meth, control = control,
                 hessian = hessian)
  }
  
  #Bayesian inference setting
  if(inf.type == "Bayes"){
    #Outcome setting
    optim.msg <- "Something went wrong"
    #Mle
    mle <- optim(par0, function(x) gev.lik(x, data = data, t = t, p = t.prob, log = TRUE),
                 method = optim.meth, control = control,
                 hessian = hessian)
    mle.est <- mle$par
    #Empirical Bayes
    #hatdelta
    hatdeltak <- mle.est[3] * s^(mle.est[1])
    #hatmu
    hatmuk <- mle.est[2] + hatdeltak*(1 - s^(-mle.est[1]))/mle.est[1]
    #Prior hyperparameters
    hyp <- c(hatmuk, hatdeltak)
    
    #Prior
    gev.prior <- function(par, hyp){
      
      #Parameters
      gam <- par[1]; mu <- par[2]; delta <- par[3]
      #Hyperparameters (hatmu, hatdelta)
      bk <- hyp[1]; ak <- hyp[2]
      
      #gamma prior
      pi.gam <- (1 - pt(-1, df = 1))^(-1)*dt(gam, df = 1)*I(gam > -1)
      #Mu prior
      pi.mu <- (1/ak)*dnorm((mu - bk)/ak)
      #Delta prior
      pi.del <- (1/ak)*dgamma(delta, shape = 1, scale = ak)
      
      #Prior
      return(sum(log(pi.gam), log(pi.mu), log(pi.del)))
    }
    
    gev.post <- function(par, data, hyp, t = t, p = t.prob){
      return(sum(gev.lik(par = par, data = data, t = t, p = t.prob, log = TRUE),
                 gev.prior(par, hyp)))
    }
    
    if(llik.type == "Gev-Cens"){
      optimfun <- function(par){
        return(gev.post(par = par, data = data, hyp = hyp, t = t, p = t.prob))
      }
    }
    #Compute the posterior mode estimates
    #Hessian = TRUE preferred for minimizations
    est <- optim(par0, optimfun, method = optim.meth, control = control,
                 hessian = hessian)
  }
  
  #Define the estimated GEV parameters: GEV fitting for the maximum asymptotic
  #distribution --> I have to refer to this estimates for return values.
  #I need these special estimates in order to compute extreme quantiles.
  scale <- est$par[3] * s^(est$par[1])
  loc <- est$par[2] + scale * (1 - s^(-est$par[1]))/est$par[1]
  est$par[3] <- scale
  est$par[2] <- loc
  names(est$par) <- c("Shape", "Location", "Scale")
  optim.msg <- est$convergence
  names(optim.msg) <- c("Optimization message")
  
  #Define the extreme quantiles
  Q.ext <- loc + scale * ((s*p)^(-est$par[1]) - 1) / est$par[1]
  #Taylor expansion from the following...?
  #Q.ext <- loc + scale * ((-log(1 - s*p))^(-est$par[1]) - 1) / est$par[1]
  names(Q.ext) <- c("Extreme-Quantile")
  #Define the return-level: G(R.lev) = 1 - 1/T.ex
  #Return level associated with the return period T.ex
  R.lev <- loc + scale * ((-log(1 - 1/T.ret))^(-est$par[1]) - 1) / est$par[1]
  names(R.lev) <- c("Return-Level")
  
  #Return object
  if(hessian == FALSE){
    return(list(inference = inf.type,
                estimate = est$par,
                value = est$value,
                hessian = est$hessian,
                Q.extreme = Q.ext,
                R.level = R.lev,
                optim.msg = optim.msg,
                exc.prop = exc.prop))
  }
  else{
    #Preferred for minimizations
    return(list(inference = inf.type,
                est = est$par,
                value = est$value,
                hessian = est$hessian,
                Q.extreme = Q.ext,
                R.level = R.lev,
                optim.msg = optim.msg,
                exc.prop = exc.prop))
  }
}

#fit.gev.inference trials

#Frequentist inference
mle.gum <- nlminb(c(1,0,1), function(x) -gev.lik(x, sam.gum, log = T))
mle.gum        #gamma --> 0
gum.freq <- fit.gev.inference(data = sam.gum, par0 = c(1,0,1), hessian = TRUE,
                  inf.type = "Frequent"); gum.freq

mle.fre <- nlminb(c(1,5,5), function(x) -gev.lik(x, sam.fre, log = T))
mle.fre        #gamma > 0
fre.freq <- fit.gev.inference(data = sam.fre, par0 = c(1,5,5), hessian = TRUE,
                  inf.type = "Frequent"); fre.freq

mle.wei <- nlminb(c(1,1,5), function(x) -gev.lik(x, sam.wei, log = T))
mle.wei        #gamma < 0
wei.freq <- fit.gev.inference(data = sam.wei, par0 = c(-0.5,0,1), hessian = TRUE,
                  inf.type = "Frequent"); wei.freq

#-------------------------------------------------------------------------------

#Using evd likelihood
mle.gum.evd <- nlminb(c(5,0,1), function(x) -gev.lik.evd(x, sam.gum, log = T), 
                  lower = c(-1+1e-8,-Inf,1e-5), upper = rep(Inf, 3))
mle.gum.evd      #Same as mle.gum

mle.fre.evd <- nlminb(c(1,5,5), function(x) -gev.lik.evd(x, sam.fre, log = T), 
                  lower = c(-1+1e-8,-Inf,1e-5), upper = rep(Inf, 3))
mle.fre.evd      #Same as mle.fre

mle.wei.evd <- nlminb(c(1,1,2), function(x) -gev.lik.evd(x, sam.wei, log = T), 
                  lower = c(-1+1e-8,-Inf,1e-5), upper = rep(Inf, 3))
mle.wei.evd     #Same as mle.wei

#-------------------------------------------------------------------------------

#Empirical Bayes procedure: needed quantities

#Sample size
n <- 1000
#Exceedance threshold 
p <- 0.95
#Effective sample size
k <- round((1 - p)*n); k
#Reciprocal of the proportion of exceedances
s <- n/k

#ML estimator function for hyperparameters in the prior
#Estimates that approach the true parameters with an appropriate rate

#hatdelta
hatdeltak <- function(par, s){
  par[3] * s^(par[1])
}
#Observed values for the three samples 
hatd.gum <- hatdeltak(mle.gum$par, s); hatd.gum
hatd.fre <- hatdeltak(mle.fre$par, s); hatd.fre
hatd.wei <- hatdeltak(mle.wei$par, s); hatd.wei
#hatmu
hatmuk <- function(par, s, hatdk){
  par[2] + hatdk*(1 - s^(-par[1]))/par[1]
}
#Observed values for the three samples
hatmu.gum <- hatmuk(mle.gum$par, s, hatd.gum); hatmu.gum
hatmu.fre <- hatmuk(mle.fre$par, s, hatd.fre); hatmu.fre
hatmu.wei <- hatmuk(mle.wei$par, s, hatd.wei); hatmu.wei

#Prior hyperparameters
hyp.gum <- c(hatmu.gum, hatd.gum)
hyp.fre <- c(hatmu.fre, hatd.fre)
hyp.wei <- c(hatmu.wei, hatd.wei)

#-------------------------------------------------------------------------------

#Prior distribution for (gamma, mu, delta)

gev.prior <- function(par, hyp, log = FALSE){
  
  #Parameters
  gam <- par[1]; mu <- par[2]; delta <- par[3]
  #Hyperparameters (hatmu, hatdelta)
  bk <- hyp[1]; ak <- hyp[2]
  
  #gamma prior
  pi.gam <- (1 - pt(-1, df = 1))^(-1)*dt(gam, df = 1)*I(gam > -1)
  #Mu prior
  pi.mu <- (1/ak)*dnorm((mu - bk)/ak)
  #Delta prior
  pi.del <- (1/ak)*dgamma(delta, shape = 1, scale = ak)
  
  #Prior
  if(log == FALSE) return(prod(pi.gam, pi.mu, pi.del))
  if(log == TRUE) return(sum(log(pi.gam), log(pi.mu), log(pi.del)))
}

#-------------------------------------------------------------------------------

#Posterior distribution for (hat, mu, delta)

gev.post <- function(par, data, hyp, log = FALSE){
  if(log == FALSE) return(prod(gev.lik(par, data, log = FALSE),
                               gev.prior(par, hyp, log = FALSE)))
  if(log == TRUE) return(sum(gev.lik(par, data, log = TRUE),
                               gev.prior(par, hyp, log = TRUE)))
}

#-------------------------------------------------------------------------------

#Bayesian inference
mode.gum <- nlminb(mle.gum$par, function(x) -gev.post(x, sam.gum, hyp.gum, log = T), 
                   lower = c(-1+1e-8,-Inf,1e-5), upper = rep(Inf, 3))
mode.gum
gum.bayes <- fit.gev.inference(data = sam.gum, par0 = mle.gum$par, hessian = TRUE,
                  inf.type = "Bayes"); gum.bayes;
mode.fre <- nlminb(mle.fre$par, function(x) -gev.post(x, sam.fre, hyp.fre, log = T), 
                   lower = c(-1+1e-8,-Inf,1e-5), upper = rep(Inf, 3))
mode.fre
fre.bayes <- fit.gev.inference(data = sam.fre, par0 = mle.fre$par, hessian = TRUE,
                  inf.type = "Bayes"); fre.bayes
mode.wei <- nlminb(mle.wei$par, function(x) -gev.post(x, sam.wei, hyp.wei, log = T), 
                   lower = c(-1+1e-8,-Inf,1e-5), upper = rep(Inf, 3))
mode.wei
wei.bayes <- fit.gev.inference(data = sam.wei, par0 = mle.wei$par, hessian = TRUE,
                  inf.type = "Bayes"); wei.bayes

#-------------------------------------------------------------------------------

#Adaptive Random Walk Metropolis-Hastings Markov Chain Monte Carlo

library(mvtnorm)
library(tictoc)
library(plyr)
library(coda)
library(TeachingDemos)

#Useful quantities
#Sampling size
R <- 10^4
#Deriserd overall sampler acceptance probability
etastar <- 0.234
#Quantities related to the adaptive update of k
zeta0 <- -qnorm(etastar/2)
#Steplength constant
a <- ((2*pi)^(1/2))*exp((zeta0^2)/2)/(2*zeta0)

#MH algorithm
adamh <- function(R, data, hyp, par0, sigma0, k0){
  
  p <- length(par0)
  out <- array(dim = c(R, p))
  #Acceptance vector
  accepted <- 0
  #Initialization
  par <- par0; sigma <- sigma0; k <- k0
  
  for(i in 1:R){
    #Proposal
    pars <- par + rmvnorm(1, sigma = k*sigma)
    #Acceptance proability
    eta <- min(1, gev.post(pars, data, hyp)/gev.post(par, data, hyp))
    #Accepatance/rejection step
    if(runif(1) < eta){
      par <- pars
      accepted <- accepted + 1
    }
    out[i,] <- par
    
    #Adaptive covariance matrix update (Haario et al. 2001)
    if(i <= 100) sigma = (1 + k^2/i)*diag(p)
    else{
      if(i == 1) theta.til <- out[i,]
      else theta.til <- (1/i)*colSums(out[1:i,])
      sigma.til <- matrix(0, nrow = p, ncol = p)
      for(j in 1:i){
        sigma.til <- sigma.til + ((out[j,] - theta.til) %*% t(out[j,] - theta.til))
      }
      sigma <- (1/(i - 1))*sigma.til + ((k^2/i)*diag(p))
    }
    
    #Adaptive k update via Robbins-Monro process (Garthwaite et al. 2016)
    kstar <- exp(log(k) + a*(eta - etastar))
    k <- kstar
    
    #Diagnostic
    if(i %% 500 == 0) cat(round((i/R)*100, 2), "%\n")
  }
  list(values = out, accepted = accepted/R)
}

#-------------------------------------------------------------------------------

#Gumbel domain of attraction

#Initialization
theta0.gum <- mode.gum$par; theta0.gum
sigma0.gum <- solve(-gum.bayes$hessian); sigma0.gum
k0.gum <- 2
#Posterior sample
set.seed(1)
tic()
post.gum <- adamh(R, sam.gum, hyp.gum, theta0.gum, sigma0.gum, k0.gum)
toc()
#193.03 sec elapsed
post.gum$accepted       #0.2343
#Trace plots
par(mfrow = c(3,1))
plot(post.gum$values[,1], type = "l")
abline(h = theta0.gum[1], col = 2, lwd = 2.5)
plot(post.gum$values[,2], type = "l")
abline(h = theta0.gum[2], col = 2, lwd = 2.5)
plot(post.gum$values[,3], type = "l")
abline(h = theta0.gum[3], col = 2, lwd = 2.5)
par(mfrow = c(1,1))

#Check of convergence through 5 different chains
set.seed(1)
starts.gum <- rmvnorm(5, mean = theta0.gum, sigma = 1.5*sigma0.gum); starts.gum
res.list.gum <- list()
tic()
for(i in 1:5){
  res.list.gum[[i]] <- adamh(R, sam.gum, hyp.gum, starts.gum[i,], sigma0.gum, k0.gum)$values
  cat(i, "")
}
toc()
#938.57 sec elapsed
res.mcmc.gum <- llply(res.list.gum, function(x) mcmc(window(x, start = floor(R/2) +1),
                                                     start = floor(R/2) +1))
res.mcmc.gum <- mcmc.list(res.mcmc.gum)
summary(res.mcmc.gum)
plot(res.mcmc.gum)
acfplot(res.mcmc.gum)
gelman.diag(res.mcmc.gum)
gelman.plot(res.mcmc.gum)
effectiveSize(res.mcmc.gum)
1 - rejectionRate(res.mcmc.gum)

#Final sample
resf.gum <- NULL
for(i in 1:length(res.mcmc.gum)){
  resf.gum <- rbind(resf.gum, res.mcmc.gum[[i]])
}
dim(resf.gum)
#Plots and statistics
theta0.gum
apply(resf.gum, 2, mean)
apply(resf.gum, 2, emp.hpd)
par(mfrow = c(1,3))
hist(resf.gum[,1], nclass = 50, probability = T, xlab = expression(gamma),
     main = "Posterior density")
hist(resf.gum[,2], nclass = 50, probability = T, xlab = expression(mu),
     main = "Posterior density")
hist(resf.gum[,3], nclass = 50, probability = T, xlab = expression(delta),
     main = "Posterior density")
par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#Fréchet domain of attraction

#Initialization
theta0.fre <- mode.fre$par; theta0.fre
sigma0.fre <- solve(-fre.bayes$hessian); sigma0.fre
k0.fre <- 2
#Posterior sample
set.seed(2)
tic()
post.fre <- adamh(R, sam.fre, hyp.fre, theta0.fre, sigma0.fre, k0.fre)
toc()
#214.68 sec elapsed
post.fre$accepted      #0.2369
#Trace plots
par(mfrow = c(3,1))
plot(post.fre$values[,1], type = "l")
abline(h = theta0.fre[1], col = 2, lwd = 2.5)
plot(post.fre$values[,2], type = "l")
abline(h = theta0.fre[2], col = 2, lwd = 2.5)
plot(post.fre$values[,3], type = "l")
abline(h = theta0.fre[3], col = 2, lwd = 2.5)
par(mfrow = c(1,1))

#Check of convergence through 5 different chains
set.seed(20)
starts.fre <- rmvnorm(5, mean = theta0.fre, sigma = sigma0.fre); starts.fre
res.list.fre <- list()
tic()
for(i in 1:5){
  res.list.fre[[i]] <- adamh(R, sam.fre, hyp.fre, starts.fre[i,], sigma0.fre, k0.fre)$values
  cat(i, "----------------------------------------------------------------", "\n")
}
toc()
#1584.3 sec elapsed
res.mcmc.fre <- llply(res.list.fre, function(x) mcmc(window(x, start = floor(R/2) +1),
                                                     start = floor(R/2) +1))
res.mcmc.fre <- mcmc.list(res.mcmc.fre)
summary(res.mcmc.fre)
plot(res.mcmc.fre)
acfplot(res.mcmc.fre)
gelman.diag(res.mcmc.fre)
gelman.plot(res.mcmc.fre)
effectiveSize(res.mcmc.fre)
1 - rejectionRate(res.mcmc.fre)

#Final sample
resf.fre <- NULL
for(i in 1:length(res.mcmc.fre)){
  resf.fre <- rbind(resf.fre, res.mcmc.fre[[i]])
}
dim(resf.fre)
#Plots and statistics
theta0.fre
apply(resf.fre, 2, mean)
apply(resf.fre, 2, emp.hpd)
par(mfrow = c(1,3))
hist(resf.fre[,1], nclass = 50, probability = T, xlab = expression(gamma),
     main = "Posterior density")
hist(resf.fre[,2], nclass = 50, probability = T, xlab = expression(mu),
     main = "Posterior density")
hist(resf.fre[,3], nclass = 50, probability = T, xlab = expression(delta),
     main = "Posterior density")
par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#Reverse Weibull domain of attraction

#Initialization
R <- 50000
theta0.wei <- mode.wei$par; theta0.wei
sigma0.wei <- solve(-wei.bayes$hessian); sigma0.wei
k0.wei <- 2
#Posterior sample
set.seed(3)
tic()
post.wei <- adamh(R, sam.wei, hyp.wei, theta0.wei, sigma0.wei, k0.wei)
toc()
#4631.1 sec elapsed (R = 50000)
post.wei$accepted     #0.2309
#Trace plots
par(mfrow = c(3,1))
plot(post.wei$values[,1], type = "l")
abline(h = theta0.wei[1], col = 2, lwd = 2.5)
plot(post.wei$values[,2], type = "l")
abline(h = theta0.wei[2], col = 2, lwd = 2.5)
plot(post.wei$values[,3], type = "l")
abline(h = theta0.wei[3], col = 2, lwd = 2.5)
par(mfrow = c(1,1))

#Check of convergence through 5 different chains
set.seed(80)
starts.wei <- rmvnorm(5, mean = theta0.wei, sigma = 0.5*sigma0.wei); starts.wei
res.list.wei <- list()
tic()
for(i in 1:5){
  res.list.wei[[i]] <- adamh(R, sam.wei, hyp.wei, starts.wei[i,], sigma0.wei, k0.wei)$values
  cat(i, "----------------------------------------------------------------", "\n")
}
toc()
#22870.86 sec elapsed
res.mcmc.wei <- llply(res.list.wei, function(x) mcmc(window(x, start = floor(R/2) +1),
                                                     start = floor(R/2) +1))
res.mcmc.wei <- mcmc.list(res.mcmc.wei)
summary(res.mcmc.wei)
plot(res.mcmc.wei)
acfplot(res.mcmc.wei)
gelman.diag(res.mcmc.wei)
gelman.plot(res.mcmc.wei)
effectiveSize(res.mcmc.wei)
1 - rejectionRate(res.mcmc.wei)

#Check for convergence of single chains via posterior
for(i in 1:length(res.mcmc.wei)) print(summary(res.mcmc.wei[[i]]))
post.dist.wei <- matrix(0, nrow = floor(R/2), ncol = length(res.mcmc.wei))
for(i in 1:length(res.mcmc.wei)){
  post.dist.wei[,i] <- apply(res.mcmc.wei[[i]], 1, function(x) gev.post(x, sam.wei, hyp.wei, log = T))
}
plot(post.dist.wei[,1], type = "l", ylab = "Log-posterior", main = "Convergence comparison")
for(i in 2:dim(post.dist.wei)[2]) lines(post.dist.wei[,i], col = i)

#Final sample
resf.wei <- NULL
for(i in 1:length(res.mcmc.wei)){
  resf.wei <- rbind(resf.wei, res.mcmc.wei[[i]])
}
dim(resf.wei)
#Plots and statistics
theta0.wei
apply(resf.wei, 2, mean)
apply(resf.wei, 2, emp.hpd)
par(mfrow = c(1,3))
hist(resf.wei[,1], nclass = 50, probability = T, xlab = expression(gamma),
     main = "Posterior density")
hist(resf.wei[,2], nclass = 50, probability = T, xlab = expression(mu),
     main = "Posterior density")
hist(resf.wei[,3], nclass = 50, probability = T, xlab = expression(delta),
     main = "Posterior density")
par(mfrow = c(1,1))

#-------------------------------------------------------------------------------