#Useful functions and libraries

library(evd)
library(mvtnorm)
library(plyr)
library(coda)
library(TeachingDemos)
library(tictoc)

#-------------------------------------------------------------------------------

#GEV distribution

#Probability density function

dgev1 <- function(x, par, log = FALSE){
  gam <- par[1]; mu <- par[2]; delta <- par[3]
  pdf <- 0
  if(delta <= 0){
    if(log == FALSE) return(pdf)
    if(log == TRUE) return(log(pdf))
  }
  else{
    z <- (x - mu)/delta
    #Define 0 neighbourhood
    eps <- .Machine$double.eps^.3
    #the smallest positive floating-point number x such that 1 + x != 1
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
}
dgev.V <- Vectorize(dgev1, "x")

#Cumulative density function

pgev1 <- function(x, par, log = FALSE){
  gam <- par[1]; mu <- par[2]; delta <- par[3]
  cdf <- 0
  if(delta <= 0){
    if(log == FALSE) return(cdf)
    if(log == TRUE) return(log(cdf))
  }
  else{
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
      if(gam > 0 & ((1 + gam*z) <= 0)) cdf <- 0
      if(gam < 0 & ((1 + gam*z) <= 0)) cdf <- 1
    }
    if(log == FALSE) return(cdf)
    if(log == TRUE) return(log(cdf))
  }
}
pgev.V <- Vectorize(pgev1, "x")

#-------------------------------------------------------------------------------

#GP distribution

#Probability density function

dgp1 <- function(x, par, log = FALSE){
  gam <- par[1]; delta <- par[2]
  pdf <- 0
  if(delta <= 0 | x <= 0){
    if(log == FALSE) return(pdf)
    if(log == TRUE) return(log(pdf))
  }
  else{
    z <- x/delta
    #Define 0 neighbourhood
    eps <- .Machine$double.eps^.3
    #the smallest positive floating-point number x such that 1 + x != 1
    if(abs(gam) <= eps){
      pdf <- (1/delta)*exp(-z)
    }
    else{
      if((1 + gam*z) > 0){
        pdf <- (1/delta)*(1 + gam*z)^(-(1/gam) - 1)
      }
    }
    if(log == FALSE) return(pdf)
    if(log == TRUE) return(log(pdf))
  }
}
dgp.V <- Vectorize(dgp1, "x")

#Cumulative density function

pgp1 <- function(x, par, log = FALSE){
  gam <- par[1]; delta <- par[2]
  cdf <- 0
  if(delta <= 0 | x <= 0){
    if(log == FALSE) return(cdf)
    if(log == TRUE) return(log(cdf))
  }
  else{
    z <- x/delta
    #Define 0 neighbourhood
    eps <- .Machine$double.eps^.3
    #the smallest positive floating-point number x such that 1 + x != 1
    if(abs(gam) <= eps){
      cdf <- 1 - exp(-z)
    }
    else{
      cdf <- 1 - (1 + gam*z)^(-1/gam)
      if(gam > 0 & ((1 + gam*z) <= 0)) cdf <- 0
      if(gam < 0 & ((1 + gam*z) <= 0)) cdf <- 1
    }
    if(log == FALSE) return(cdf)
    if(log == TRUE) return(log(cdf))
  }
}
pgp.V <- Vectorize(pgp1, "x")

#-------------------------------------------------------------------------------

#Censored GEV likelihood

gev.lik <- function(par, data, t = NULL, p = 0.95, log = FALSE,
                    #llik.type = Max-Gev-Cens/ Gev-Cens
                    llik.type = "Gev-Cens",
                    #param.type = tilde / bar
                    param.type = "tilde"){
  
  gam <- par[1]; mu <- par[2]; delta <- par[3]
  
  #Set the smallest value of the log-likelihood
  Low <- -1e300
  #Check the validity of the required likelihood
  if(!any(llik.type == "Max-Gev-Cens", llik.type == "Gev-Cens")){
    stop("Must choose the likelihood between Max-Gev-Cens and Gev-Cens \n")
  }
  #Check the validity of the parameterization
  if(!any(param.type == "tilde", param.type == "bar")){
    stop("Must choose the parameterization between tilde and bar \n")
  }
  #Check the validity of the data
  if(!is.vector(data)){stop("Data must be a vector \n")}
  #Check the validity of the parameters
  if(length(par) != 3){stop("Wrong length of parameter vector \n")}
  #Check the support of the scale parameter delta
  if(delta <= 0){
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
  cond <- data > t
  #Exceedances
  exceed <- data[cond]
  #Sample size
  n <- length(data)
  #Effective sample size
  k <- sum(cond)
  
  if(llik.type == "Max-Gev-Cens"){
    #Likelihood (tilde parameterization)
    if(log == FALSE){
      lik.cens <- (pgev.V(t, par, log = FALSE)^(k/n))^(n - k)
      #lik.cens <- (pgev.V(t, par.new, log = FALSE))^(n - k)
      lik.obs <- prod((k/n)*(pgev.V(exceed, par, log = FALSE)^(k/n - 1))*
                        dgev.V(exceed, par, log = FALSE))
      #lik.obs <- prod(dgev.V(exceed, par.new, log = FALSE))
      lik <- prod(lik.cens, lik.obs)
      if(is.infinite(lik)) return(0)
      return(lik)
    }
    
    #Log-likelihood (tilde parameterization)
    if(log == TRUE){
      llik.cens <- (n - k)*((k/n)*pgev.V(t, par, log = TRUE))
      #llik.cens <- (n - k)*(pgev.V(t, par.new, log = TRUE))
      llik.obs <- sum(log(k/n) + dgev.V(exceed, par, log = TRUE) +
                        (k/n - 1)*pgev.V(exceed, par, log = TRUE))
      #llik.obs <- sum(dgev.V(exceed, par.new, log = TRUE))
      llik <- sum(llik.cens, llik.obs)
      if(is.infinite(llik)) return(Low)
      return(llik)
    }
  }
  
  if(llik.type == "Gev-Cens"){
    #Define the new parameterization (GEV parameterization)
    if(param.type == "tilde"){
      delta.new <- delta*(k/n)^(gam)
      mu.new <- mu - delta * (1 - (k/n)^gam)/gam
      par.new <- c(gam, mu.new, delta.new)
      
      #Likelihood (tilde parameterization)
      if(log == FALSE){
        #lik.cens <- (pgev.V(t, par, log = FALSE))^(n - k)
        lik.cens <- (pgev.V(t, par.new, log = FALSE))^(n - k)
        #lik.obs <- prod(dgev.V(exceed, par, log = FALSE))
        lik.obs <- prod(dgev.V(exceed, par.new, log = FALSE))
        lik <- prod(lik.cens, lik.obs)
        if(is.infinite(lik)) return(0)
        return(lik)
      }
      
      #Log-likelihood (tilde parameterization)
      if(log == TRUE){
        #llik.cens <- (n - k)*(pgev.V(t, par, log = TRUE))
        llik.cens <- (n - k)*(pgev.V(t, par.new, log = TRUE))
        #llik.obs <- sum(dgev.V(exceed, par, log = TRUE))
        llik.obs <- sum(dgev.V(exceed, par.new, log = TRUE))
        llik <- sum(llik.cens, llik.obs)
        if(is.infinite(llik)) return(Low)
        return(llik)
      }
    }
    if(param.type == "bar"){
      
      #Likelihood (bar parameterization)
      if(log == FALSE){
        #lik.cens <- (pgev.V(t, par, log = FALSE))^(n - k)
        lik.cens <- (pgev.V(t, par, log = FALSE))^(n - k)
        #lik.obs <- prod(dgev.V(exceed, par, log = FALSE))
        lik.obs <- prod(dgev.V(exceed, par, log = FALSE))
        lik <- prod(lik.cens, lik.obs)
        if(is.infinite(lik)) return(0)
        return(lik)
      }
      
      #Log-likelihood (bar parameterization)
      if(log == TRUE){
        #llik.cens <- (n - k)*(pgev.V(t, par, log = TRUE))
        llik.cens <- (n - k)*(pgev.V(t, par, log = TRUE))
        #llik.obs <- sum(dgev.V(exceed, par, log = TRUE))
        llik.obs <- sum(dgev.V(exceed, par, log = TRUE))
        llik <- sum(llik.cens, llik.obs)
        if(is.infinite(llik)) return(Low)
        return(llik)
      }
    }
  }
}

#Profile likelihood
#(gamma, mu, \hat{delta}_{(gamma, mu)})
gev.lik.gammu <- function(par, data, t = NULL, p = 0.95, log = FALSE,
                          llik.type = "Gev-Cens", param.type = "tilde"){
  par0 <- 5
  prof.lik <- -nlminb(par0, function(x) -gev.lik(par = c(par[1], par[2], x),
                                                 data = data, t = t, p = t.prob,
                                                 log = log, llik.type = llik.type,
                                                 param.type = param.type),
                      lower = 1e-10, upper = Inf)$objective
  return(prof.lik)
}
#gev.lik.gammu <- Vectorize(gev.lik.gammu, "par")

#Profile likelihood
#(gamma, \hat{mu}_{(gamma, delta)}, delta)
gev.lik.gamdel <- function(par, data, t = NULL, p = 0.95, log = FALSE,
                           llik.type = "Gev-Cens", param.type = "tilde"){
  par0 <- 1
  prof.lik <- -nlminb(par0, function(x) -gev.lik(par = c(par[1], x, par[2]),
                                                 data = data, t = t, p = t.prob,
                                                 log = log, llik.type = llik.type,
                                                 param.type = param.type),
                      lower = -Inf, upper = Inf)$objective
  return(prof.lik)
}
#gev.lik.gamdel <- Vectorize(gev.lik.gamdel, "par")

#Profile likelihood
#(\hat{gamma}_{(mu, delta)}, mu, delta)
gev.lik.mudel <- function(par, data, t = NULL, p = 0.95, log = FALSE,
                          llik.type = "Gev-Cens", param.type = "tilde"){
  par0 <- 1
  prof.lik <- -nlminb(par0, function(x) -gev.lik(par = c(x, par[1], par[2]),
                                                 data = data, t = t, p = t.prob,
                                                 log = log, llik.type = llik.type,
                                                 param.type = param.type),
                      lower = -Inf, upper = Inf)$objective
  return(prof.lik)
}
#gev.lik.mudel <- Vectorize(gev.lik.mudel, "par")

#-------------------------------------------------------------------------------

#GEV distribution in the tilde parameterization

#Probability density function

dgev1.tilde <- function(x, par, s, log = FALSE){
  gam <- par[1]; mu <- par[2]; delta <- par[3]
  cdf <- 0
  pdf <- 0
  if(delta <= 0){
    if(log == FALSE) return(pdf)
    if(log == TRUE) return(log(pdf))
  }
  else{
    z <- (x - mu)/delta
    #Define 0 neighbourhood
    eps <- .Machine$double.eps^.3
    #the smallest positive floating-point number x such that 1 + x != 1
    if(abs(gam) <= eps){
      cdf <- (exp(-exp(-z)))^(1/s)
      pdf <- (1/(s*delta))*cdf*exp(-z)
    }
    else{
      gev.k <- (1 + gam*z)^(-1/gam)
      gev.k2 <- (1 + gam*z)^(-1/gam -1)
      cdf <- (exp(-gev.k))^(1/s)
      pdf <- (1/(s*delta))*cdf*gev.k2
      if(gam > 0 & ((1 + gam*z) <= 0)){
        cdf <- 0
        pdf <- 0
      }
      if(gam < 0 & ((1 + gam*z) <= 0)){
        cdf <-  1
        pdf <- 0
      } 
    }
    if(log == FALSE) return(pdf)
    if(log == TRUE) return(log(pdf))
  }
}
dgev.tilde.V <- Vectorize(dgev1.tilde, "x")

#Cumulative density function

pgev1.tilde <- function(x, par, s, log = FALSE){
  gam <- par[1]; mu <- par[2]; delta <- par[3]
  cdf <- 0
  if(delta <= 0){
    if(log == FALSE) return(cdf)
    if(log == TRUE) return(log(cdf))
  }
  else{
    z <- (x - mu)/delta
    #Define 0 neighbourhood
    eps <- .Machine$double.eps^.3
    #the smallest positive floating-point number x such that 1 + x != 1
    if(abs(gam) <= eps){
      cdf <- (exp(-exp(-z)))^(1/s)
    }
    else{
      gev.k <- (1 + gam*z)^(-1/gam)
      cdf <- (exp(-gev.k))^(1/s)
      if(gam > 0 & ((1 + gam*z) <= 0)) cdf = 0
      if(gam < 0 & ((1 + gam*z) <= 0)) cdf = 1
    }
    if(log == FALSE) return(cdf)
    if(log == TRUE) return(log(cdf))
  }
}
pgev.tilde.V <- Vectorize(pgev1.tilde, "x")

#-------------------------------------------------------------------------------

#Uniform prior for (gamma, mu, delta)

gev.unif.prior <- function(par, log = FALSE){
  
  #Same prior for tilde and bar parameterization
  #Jeffreys' prior in scale and location families
  #Jeffreys' prior is parameterization invariant
  delta <- par[3]
  pi.prior <- 1/delta
  if(log == FALSE) return(pi.prior)
  if(log == TRUE) return(log(pi.prior))
}

#Empirical Bayes prior distribution for (gamma, mu, delta)

gev.emp.prior <- function(par, hyp,
                          #param.type = tilde / bar
                          param.type = "tilde",
                          #s = n/k
                          s.frac = NULL,
                          log = FALSE){
  
  #Parameters
  gam <- par[1]; mu <- par[2]; delta <- par[3]
  #Hyperparameters (hatmu, hatdelta)
  bk <- hyp[1]; ak <- hyp[2]
  
  if(param.type == "tilde"){
    #gamma prior
    pi.gam <- (1 - pt(-1, df = 1))^(-1)*dt(gam, df = 1)*I(gam > -1)
    #Mu prior
    pi.mu <- (1/ak)*dnorm((mu - bk)/ak)
    #Delta prior
    pi.del <- (1/ak)*dgamma(delta, shape = 1, scale = ak)*I(delta > 0)
    
    #Prior
    if(log == FALSE) return(prod(pi.gam, pi.mu, pi.del))
    if(log == TRUE) return(sum(log(pi.gam), log(pi.mu), log(pi.del)))
  }
  
  if(param.type == "bar"){
    
    #Set s
    if(is.null(s.frac)){
      stop("Set the inverse of the proportion of exceedances s \n")
    }
    
    #Define 0 neighbourhood
    eps <- .Machine$double.eps^.3
    #the smallest positive floating-point number x such that 1 + x != 1
    if(abs(gam) <= eps){
      #Prior components
      pi.1 <- (1 - pt(-1, df = 1))^(-1)*dt(gam, df = 1)*I(gam > -1)
      pi.2 <- (1/ak)*dnorm((mu - bk + (delta*log(s.frac)))/ak)
      pi.3 <- (1/ak^2)*exp(-(1/ak)*delta)*I(delta > 0)
    }
    else{
      #Prior components
      pi.1 <- (s.frac^(gam))*(1 - pt(-1, df = 1))^(-1)*dt(gam, df = 1)*I(gam > -1)
      pi.2 <- (1/ak)*dnorm((mu - bk + ((s.frac^(gam))*delta*(1 - s.frac^(-gam)))/gam)/ak)
      pi.3 <- (1/ak^2)*exp(-(1/ak)*(s.frac^(gam))*delta)*I(delta > 0)
    }
    
    #Prior
    if(log == FALSE) return(prod(pi.1, pi.2, pi.3))
    if(log == TRUE) return(sum(log(pi.1), log(pi.2), log(pi.3)))
  }
}

#Posterior distribution for (hat, mu, delta)

gev.post <- function(par, data, hyp, t = NULL, p = 0.95, log = FALSE,
                     #llik.type = Max-Gev-Cens/ Gev-Cens
                     llik.type = "Gev-Cens",
                     #param.type = tilde / bar
                     param.type = "tilde",
                     #prior = "uniform" / "empirical"
                     prior = "empirical"){
  if(prior == "uniform"){
    if(log == FALSE) return(prod(gev.lik(par, data, t = t, p = p, log = FALSE,
                                         llik.type = llik.type,
                                         param.type = param.type),
                                 gev.unif.prior(par, log = FALSE)))
    if(log == TRUE) return(sum(gev.lik(par, data, t = t, p = p, log = TRUE,
                                       llik.type = llik.type,
                                       param.type = param.type),
                               gev.unif.prior(par, log = TRUE)))
  }
  if(prior == "empirical"){
    
    #Empirical threshold
    #Use discontinuous sample 95% quantile: type = 3
    if(is.null(t)){
      t <- as.numeric(quantile(data, probs = p, type = 3))
      #message("T set to 95% quantile by default \n")
    }
    #Indicator of censoring
    cond <- data > t
    #Sample size
    n <- length(data)
    #Effective sample size
    k <- sum(cond)
    #s = n/k
    s.frac <- n/k
    
    if(log == FALSE) return(prod(gev.lik(par, data, t = t, p = p, log = FALSE,
                                         llik.type = llik.type,
                                         param.type = param.type),
                                 gev.emp.prior(par, hyp,
                                               param.type = param.type,
                                               s.frac = s.frac,
                                               log = FALSE)))
    if(log == TRUE) return(sum(gev.lik(par, data, t = t, p = p, log = TRUE,
                                       llik.type = llik.type,
                                       param.type = param.type),
                               gev.emp.prior(par, hyp,
                                             param.type = param.type,
                                             s.frac = s.frac,
                                             log = TRUE)))
  }
}

#-------------------------------------------------------------------------------

#Adaptive Random Walk Metropolis-Hastings Markov Chain Monte Carlo

adamh <- function(R, data, hyp, par0, k0, t = NULL, p = p,
                  #llik.type = Max-Gev-Cens/ Gev-Cens
                  llik.type = "Gev-Cens",
                  #param.type = tilde / bar
                  param.type = "tilde",
                  #prior = "uniform" / "empirical"
                  prior = "empirical",
                  #etastar: desired overall sampler acceptance probability
                  etastar = 0.234){
  
  #Required packages
  require(mvtnorm)
  
  #Set parameters constraints
  parcheck <- function(para){
    res <- any(para[3] <= 0,
               (para[1] > 0 & (para[2] > min(data) + para[3]/para[1])),
               (para[1] < 0 & any(para[2] <= data + para[3]/para[1])))
    return(res)
  }
  n <- length(data)
  d <- length(par0)
  if(!any(d == 3)){stop("Wrong length of parameter vector")}
  #Correct specification of the prior
  if(!any(prior == "uniform", prior == "empirical")){
    stop("Specify a prior type within uniform and empirical")
  }
  
  #Quantities related to the adaptive update of k
  zeta0 <- -qnorm(etastar/2)
  #Steplength constant
  #a <- (sqrt(2*pi))*exp((zeta0^2)/2)/(2*zeta0) 
  #Suitable overstimate of steplength constant in multivariate framework
  a <- (1 - 1/d)*(sqrt(2*pi))*exp((zeta0^2)/2)/(2*zeta0) + 
    1/(d*etastar*(1 - etastar))
  
  out <- array(dim = c(R, d))
  #Acceptance vector
  accepted <- rep(0, R)
  #Acceptance probabilities vector
  accepted.prob <- rep(0, R)
  #Automatic rejection counter (proposed parameters that don't respect the constraints)
  straight.reject <- rep(0, R)
  #Vector of scaling parameters k
  k.vec <- k0
  #Monitor the adequate values of k
  k.start <- k0
  k.restart <- k0
  #Number of iterations before the beginning of RM search
  n0 <- round(5/(etastar*(1 - etastar)))
  #Max number of iterations before the last restart
  iMax <- 100
  Numbig <- 0
  Numsmall <- 0
  #Counter of the iteration number at each restart
  numRS <- 1
  #Message about the outcome of the algorithm
  msg <- "MCMC run without errors"
  #Initialization
  sigma0 <- diag(d)
  par <- par0; sigma <- sigma0; k <- k0
  
  #Progress bar
  pb <- txtProgressBar(min = 0, max = R, initial = 0, style = 3)
  
  for(i in 1:R){
    #Adequacy of the covariance matrix
    s <- try(eigen((k^2)*sigma, symmetric = TRUE), silent = TRUE)
    if(class(s) == "try-error"){
      msg <- "MCMC failed"
      return(list(values = out,
                  accepted.vec = accepted,
                  accepted = mean(accepted),
                  accepted.prob = accepted.prob,
                  straight.reject = straight.reject,
                  k.vec = k.vec,
                  k.restart = k.restart,
                  msg = msg))
    }
    
    #Proposal
    if(d == 1) pars <- rnorm(1, sigma = k)
    #pars <- par + rmvnorm(1, sigma = k*sigma)
    else pars <- par + rmvnorm(1, sigma = (k^2)*sigma)
    numRS <- numRS + 1
    #Check for NA proposed values
    if(any(is.na(pars))){
      straight.reject[i] <- 1
      accepted.prob[i] <- 0
    }
    #Check basis condition
    else{
      if(parcheck(pars)){
        straight.reject[i] <- 1
        accepted.prob[i] <- 0
      }
      #Compute the acceptance probability
      else{
        eta <- min(1, exp(gev.post(pars, data, hyp, t = t, p = p, log = TRUE,
                                   llik.type = llik.type, param.type = param.type,
                                   prior = prior) - 
                            gev.post(par, data, hyp,t = t, p = p,log = TRUE,
                                     llik.type = llik.type, param.type = param.type,
                                     prior = prior)))
        #Check for unsuccessful probability computation
        if(is.na(eta) || is.null(eta)){
          straight.reject[i] <- 1
          accepted.prob[i] <- 0
          eta <- 0
        }
        else{
          accepted.prob[i] <- eta
          #Accepatance/rejection step
          if(runif(1) < eta){
            par <- pars
            accepted[i] <- 1
          }
        }
      }
    }
    out[i,] <- par
    
    #Adaptive covariance matrix update (Haario et al. 2001)
    #if(i <= 100) sigma = (1 + k^2/i)*diag(d)
    if(i > 100 && d > 1){
      if(i == 101){
        #Update partial covariance matrix
        sigMat <- cov(out[1:i,])
        #Update partial mean vector
        thetaM <- apply(out[1:i,], 2, mean)
        sigma <- sigMat
        + diag(d)/i
        #+ ((k^2/i)*diag(d))
      }
      else{
        #Recursive update of partial mean vector
        thetaM2 <- (thetaM*(i - 1) + out[i,])/i
        #Recursive update of partial covariance matrix
        sigMat <- (i - 2)/(i - 1)*sigMat + thetaM%*%t(thetaM) -
          (i)/(i - 1)*thetaM2%*%t(thetaM2) + 1/(i - 1)*out[i,]%*%t(out[i,])
        #+ ((k^2/i)*diag(d))
        #sigMat <- ((i - 2)*sigMat)/(i - 1) +
        #  (thetaM - thetaM2)%*%t((thetaM - thetaM2)) +
        #  (1/(i - 1))*(out[i,] - thetaM2)%*%t((out[i,] - thetaM2)) 
        sigma <- sigMat
        + diag(d)/i
        #+ ((k^2/i)*diag(d))
        #Recursive update of partial mean vector
        thetaM <- thetaM2
      }
    }
    
    #Adaptive k update via Robbins-Monro process (Garthwaite et al. 2016)
    if(i > n0){
      #kstar <- exp(log(k) + a*(eta - etastar)/i)
      if(d == 1) kstar <- exp(log(k) + a*(eta - etastar)/i)
      else kstar <- exp(log(k) + a*(eta - etastar)/max(200, i/d))
      k <- kstar
      k.vec <- c(k.vec, k)
      if ((i <= (iMax + n0)) && (Numbig < 5 || Numsmall < 5)) {
        Toobig <- (k > (3*k.start))
        Toosmall <- (k < (k.start/3))
        if (Toobig || Toosmall) {
          #Restart the algorithm
          message("\n", "Restart the program at ", numRS, "th iteration")
          #k.restart
          k.restart <- c(k.restart, k)
          Numbig <- Numbig + Toobig
          Numsmall <- Numsmall + Toosmall
          i <- n0
          k.start <- k
        }
      }
    }
    
    #Progress bar
    print.i <- seq(0, R, by = 100)
    if(i %in% print.i) setTxtProgressBar(pb, i)
    
  }
  #Return object
  return(list(values = out,
              accepted.vec = accepted,
              accepted = mean(accepted),
              accepted.prob = accepted.prob,
              straight.reject = straight.reject,
              k.vec = k.vec,
              k.restart = k.restart,
              msg = msg))
}

#Adaptive k update via Robbins-Monro process (Garthwaite et al. 2016)
#a <- (sqrt(2*pi))*exp((zeta0^2)/2)/(2*zeta0) 
#kstar <- exp(log(k) + a*(eta - etastar))
#k <- kstar
#k.vec <- c(k.vec, k)

#-------------------------------------------------------------------------------

#DOUBTS ABOUT THE ITERATIVE UPDATE OF COVARIANCE MATRIX sigma AND STEPLENGHT k:

#- Why does the update of the covariance matrix start at the 100th iteration
#rather than the 1st?

#- Why do we start from an identity covariance matrix instead of (for example)
#the posterior hessian computed in the posterior mode?

#- Why do we compute the covariance matrix sigma adaptive update like this?
#Code for these computations is in lines 331-358, but in Padoan & Rizzelli (2022)
#expression (4.1) at page 12 would lead to the same code with uncommented lines
#332, 340, 354 (in place of line 353).
#Maybe the reference that could partially clarify my doubt is Garthwait et al. 
#(2016), expression (13) at page 10.

#- Why do we compute the steplenght constant k adaptive update like this?
#The code for these computations is in lines 360-381.
#If we just followed Padoan & Rizzelli (2022), equation (4.2) page 12
#the code would be:
#START CODE
#kstar <- exp(log(k) + a*(eta - etastar))
#k <- kstar
#k.vec <- c(k.vec, k)
#END CODE
#with constant a given by line 247 (Padoan & Rizzelli (2022)).
#Instead what we do is: 
# -- we compute constant a in lines 249-250 like in equation (8) at page six of 
#Garthwaite et al. (2016);
# -- we adatively update the steplength k in line 363 like in Garthwait et al. 
#(2016) expression (11), page 9. Our procedure seems to differ since we don't 
#consider different k updates for acceptance or rejection of the proposed value
#as it is showed in expression (11);
# -- we start updating k at the n0th (line 264) iteration;
# -- we restart updating k if k reaches more than 3*k.start or less than 
#k.start/3. Other quantities like iMax, Toobig, Toosmall, Numbig, Numsmall 
#control this process of updating.
#The reference for this algorithm should be the last two paragraphs of 
#subsection 4.1 at page 9 of Garthwaite et al. (2016).

#I'm not really sure to have understood how this works and I'm going to send 
#you via email a picture of the flow of iterations' numbers and restarts
#during one of the simulations (gumbel model).
#Unlike stated in the reference, in none of the used models the steplength 
#constant k shrinks to 0 as the number of iterations grows.

#-------------------------------------------------------------------------------

#Optimization: MLE /Bayesian Mode (uniform/empirical bayes prior)

fit.gev.inference <- function(data, t = NULL, t.prob = 0.95, llik.type = "Gev-Cens",
                              param.type = "tilde", T.ret = 50, p = 1/length(data),
                              par0 = NULL, optim.meth = "Nelder-Mead", control = NULL,
                              hessian = FALSE,
                              #inf.type = c("Frequent", "Bayes")
                              inf.type = NULL,
                              #prior = "uniform" / "empirical"
                              prior = "empirical",
                              #Posterior chain size
                              R = NULL,
                              #Burn-in
                              burn = NULL,
                              #Starting RW steplength k
                              k = 1,
                              #Desired overall sampler acceptance probability
                              etastar = 0.234,
                              #Plot about the MCMC algorithm
                              val.show = FALSE,
                              ...){
  
  #Choose one of the possible types of inference
  if(is.null(inf.type) || (inf.type != "Frequent" && inf.type != "Bayes")){
    stop("Need to specify the type of inference between: Frequent or Bayes")
  }
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
  
  #Sample size
  n <- length(data)
  #Effective sample size
  n.eff <- sum(data > t)
  #Define the proportion of exceedances
  #exc.prop <- k/n
  exc.prop <- n.eff/n
  names(exc.prop) <- c("Excesses proportion")
  #s = Reciprocal of proportion of exceedances
  #s = (n/k)
  s <- 1/exc.prop
  
  #Outcome setting
  optim.msg <- "Something went wrong"
  optimfun <- function(par){
    return(gev.lik(par = par, data = data, t = t, p = t.prob, log = TRUE,
                   llik.type = llik.type, param.type = param.type))
  }
  #Compute the maximum consered likelihood estimates
  #Hessian = TRUE preferred for minimizations
  mle <- optim(par0, optimfun, method = optim.meth, control = control,
               hessian = hessian)
  #Original estimate
  #mle.or <- mle$par
  
  #Define the estimated GEV parameters
  #scale <- mle$par[3] * exc.prop^(-mle$par[1])
  #loc <- mle$par[2] + scale * (1 - exc.prop^(mle$par[1]))/mle$par[1]
  #mle$par[3] <- scale
  #mle$par[2] <- loc
  names(mle$par) <- c("Shape", "Location", "Scale")
  
  #Frequentist inference setting
  if(inf.type == "Frequent"){
    
    optim.msg <- mle$convergence
    names(optim.msg) <- c("Optimization message")
    
    #Define the extreme quantiles
    #Q.ext <- loc + scale * ((s*p)^(-mle$par[1]) - 1) / mle$par[1]
    Q.ext <- mle$par[2] + mle$par[3] * ((s*p)^(-mle$par[1]) - 1) / mle$par[1]
    #Taylor expansion from the following...?
    #Q.ext <- loc + scale * ((-log(1 - s*p))^(-mle$par[1]) - 1) / mle$par[1]
    names(Q.ext) <- c("Extreme-Quantile")
    #Define the return-level: G(R.lev) = 1 - 1/T.ex
    #Return level associated with the return period T.ex
    #R.lev <- loc + scale * ((-log(1 - 1/T.ret))^(-mle$par[1]) - 1) / mle$par[1]
    R.lev <- mle$par[2] + mle$par[3] * ((-log(1 - 1/T.ret))^(-mle$par[1]) - 1) / mle$par[1]
    names(R.lev) <- c("Return-Level")
    
    #Return object
    if(hessian == FALSE){
      return(list(inference = inf.type,
                  #                  mle.or = mle.or,
                  mle = mle$par,
                  max.lik = mle$value,
                  hessian = mle$hessian,
                  Q.extreme = Q.ext,
                  R.level = R.lev,
                  optim.msg = optim.msg,
                  exc.prop = exc.prop))
    }
    else{
      #Preferred for minimizations
      return(list(inference = inf.type,
                  #                  mle.or = mle.or,
                  mle = mle$par,
                  max.lik = mle$value,
                  hessian = mle$hessian,
                  Q.extreme = Q.ext,
                  R.level = R.lev,
                  optim.msg = optim.msg,
                  exc.prop = exc.prop))
    }
  }
  
  #Bayesian inference setting
  if(inf.type == "Bayes"){
    #Outcome setting
    optim.msg <- "Something went wrong"
    #Check for the number of posterior samples
    if(is.null(R)){stop("Missing number of replications for MCMC")}
    #Set the number of quantiles to extrapolate
    nQuant <- length(p)
    #Correct specification of the prior
    if(!any(prior == "uniform", prior == "empirical")){
      stop("Specify a prior type within uniform and empirical")
    }
    
    #Empirical Bayes
    #Maximum likelihood estimates for norming constants in the tilde parameterization
    optimfun <- function(par){
      return(gev.lik(par = par, data = data, t = t, p = t.prob, log = TRUE,
                     llik.type = "Gev-Cens", param.type = "tilde"))
    }
    #Compute the maximum consered likelihood of the tilde parameters
    mle.tilde <- optim(par0, optimfun, method = optim.meth, control = control,
                       hessian = hessian)
    #Prior hyperparameters
    hyp <- mle.tilde$par[2:3]
    hyp.ret <- mle.tilde$par[2:3]
    names(hyp.ret) <- c("b(n/k)", "a(n/k)")
    #Optimization
    optimfun <- function(par){
      return(gev.post(par = par, data = data, hyp = hyp, t = t, p = t.prob,
                      llik.type = llik.type, param.type = param.type,
                      prior = prior, log = TRUE))
    }
    #Compute the posterior mode estimates
    #Hessian = TRUE preferred for minimizations
    mode <- optim(par0, optimfun, method = optim.meth, control = control,
                  hessian = hessian)
    optim.msg <- mode$convergence
    names(optim.msg) <- c("Optimization message")
    #Original estimate
    #mode.or <- mode$par
    #Posterior negative inverse hessian
    #sigma <- solve(-mode$hessian)
    #Posterior mode GEV parameters
    #scale <- mode$par[3] * s^(mode$par[1])
    #loc <- mode$par[2] + scale * (1 - s^(-mode$par[1]))/mode$par[1]
    #mode$par[3] <- scale
    #mode$par[2] <- loc
    names(mode$par) <- c("Shape", "Location", "Scale")
    #Posterior chain
    post.mcmc <- adamh(R = R, data = data, hyp = hyp,
                       par0 = par0, k0 = k,
                       t = t, p = p, prior = prior,
                       llik.type = llik.type,
                       param.type = param.type,
                       etastar = etastar)
    #Acceptance rate
    post.acc <- post.mcmc$accepted
    mean.acc <- rep(NA, length(post.mcmc$accepted.prob))
    for(j in c(1:length(post.mcmc$accepted.prob))){
      mean.acc[j] <- mean(post.mcmc$accepted.prob[round(j/2) :j])
    }
    index <- c(1:2000, seq(2001, length(post.mcmc$k.vec), by = 100))
    #Plots of steplength and acceptance probability
    if(val.show){
      #oldpar <- par(mar = c(4.4, 4.8, 0.5, 0.5))
      #on.exit(par(oldpar))
      par(mfrow=c(1, 2))
      plot(cbind(index, post.mcmc$k.vec[index]^2), type = "l", col = 3,
           ylim = c(0, max(post.mcmc$k.vec^2)), ylab = expression(kappa),
           xlab = "Iterations", lwd = 2)
      abline(h = 0, lwd = 2)
      plot(cbind(index, mean.acc[index]), type = "l", col = 2, ylim = c(0,1),
           ylab = "Acceptance Probability", xlab = "Iterations", lwd = 2)
      abline(h = etastar, lwd = 2)
      par(mfrow = c(1, 1))
      #Deduce approximate samples from the parameter's posterior distributions
    }
    #Posterior distribution of the original parameters (after burn-in)
    if(is.null(burn)){
      burn <- round(R/4)
      message(paste0("\n Burn-in set to the 25% of ", R, " by default"))
    } 
    post.sam <- post.mcmc$values[(burn + 1):R,]
    #post.sam.or <- post.sam
    #colnames(post.sam.or) <- c("Shape", "Location", "Scale")
    #Posterior distributions of the GEV parameters
    #scale <- apply(post.sam, 1, function(x) x[3]*s^(x[1]))
    #post.sam[,3] <- scale
    scale <- post.sam[,3]
    #loc <- apply(post.sam, 1, function(x) x[2] + x[3] *(1 - s^(-x[1]))/x[1])
    #post.sam[,2] <- loc
    loc <- post.sam[,2]
    colnames(post.sam) <- c("Shape", "Location", "Scale")
    #Extreme quantile chain
    Q.ext <- matrix(NaN, nrow = nrow(post.sam), ncol= nQuant)
    for(j in 1:nQuant){
      Q.ext[,j] <- loc + scale * ((exc.prop/p[j])^post.sam[,1] - 1)/post.sam[,1]
    }
    #Q.ext <- loc + scale *((s*p)^(-post.sam[,1]))/post.sam[,1]
    #Return level
    R.lev <- loc + scale * ((-log(1 - 1/T.ret))^(-post.sam[,1]) - 1)/post.sam[,1]
    
    #Return object
    if(hessian == FALSE){
      return(list(inference = inf.type,
                  #                  mle.or = mle.or,
                  mle = mle$par,
                  #                  mode.or = mode.or,
                  mode = mode$par,
                  max.post = mode$value,
                  optim.msg = optim.msg,
                  emp.bayes.hyp = hyp.ret,
                  mcmc.acc = post.acc,
                  #                  parameters.or = post.sam.or,
                  parameters = post.sam,
                  Q.extreme = Q.ext,
                  R.level = R.lev,
                  straight.reject = post.mcmc$straight.reject[(burn + 1):R],
                  k.vec = post.mcmc$k.vec,
                  accept.prob = cbind(index, mean.acc[index]),
                  msg = post.mcmc$msg,
                  exc.prop = exc.prop))
    }
    else{
      #Preferred for minimizations
      return(list(inference = inf.type,
                  #                  mle.or = mle.or,
                  mle = mle$par,
                  #                  mode.or = mode.or,
                  mode = mode$par,
                  hessian = mode$hessian,
                  max.post = mode$value,
                  optim.msg = optim.msg,
                  emp.bayes.hyp = hyp.ret,
                  mcmc.acc = post.acc,
                  #                  parameters.or = post.sam.or,
                  parameters = post.sam,
                  Q.extreme = Q.ext,
                  R.level = R.lev,
                  straight.reject = post.mcmc$straight.reject[(burn + 1):R],
                  k.vec = post.mcmc$k.vec,
                  accept.prob = cbind(index, mean.acc[index]),
                  msg = post.mcmc$msg,
                  exc.prop = exc.prop))
    }
  }
}