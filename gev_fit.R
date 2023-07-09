#GEV distribution

library(evd)

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
  cond <- data > t
  #Exceedances
  exceed <- data[cond]
  #Sample size
  n <- length(data)
  #Effective sample size
  k <- sum(cond)
  
  #Likelihood
  if(log == FALSE){
    lik.cens <- (pgev.V(t, par, log = FALSE))^(n - k)
    lik.obs <- prod(dgev.V(exceed, par, log = FALSE))
    lik <- prod(lik.cens, lik.obs)
    if(is.infinite(lik)) return(0)
    return(lik)
  }
  
  #Log-likelihood
  if(log == TRUE){
    llik.cens <- (n - k)*(pgev.V(t, par, log = TRUE))
    llik.obs <- sum(dgev.V(exceed, par, log = TRUE))
    llik <- sum(llik.cens, llik.obs)
    if(is.infinite(llik)) return(Low)
    return(llik)
  }
}

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

#Posterior distribution for (hat, mu, delta)

gev.post <- function(par, data, hyp, t = NULL, p = 0.95, log = FALSE){
  if(log == FALSE) return(prod(gev.lik(par, data, t = t, p = p, log = FALSE),
                               gev.prior(par, hyp, log = FALSE)))
  if(log == TRUE) return(sum(gev.lik(par, data, t = t, p = p, log = TRUE),
                             gev.prior(par, hyp, log = TRUE)))
}

#Adaptive Random Walk Metropolis-Hastings Markov Chain Monte Carlo

#library(mvtnorm)
#library(tictoc)
#library(plyr)
#library(coda)
#library(TeachingDemos)

#MH algorithm
adamh <- function(R, data, hyp, par0, sigma0, k0, t = NULL, p = p,
                  #etastar: desired overall sampler acceptance probability
                  etastar = 0.234,
                  #time = desired interval of generated samples for which we
                  #require the timing
                  time = 500
                  ){
  
  #Required packages
  require(mvtnorm)
  require(tictoc)
  
  #Timing
  if(R %% time == 0){
    timing <- rep(0, R/time)
    time.labs <- as.character(seq(time, R, by = time))
    names(timing) <- time.labs
  } 
  else{
    timing <- rep(0, ceiling(R/time))
    time.labs <- as.character(c(seq(time, R, by = time), R))
    names(timing) <- time.labs
  } 
  tic()
  
  #Quantities related to the adaptive update of k
  zeta0 <- -qnorm(etastar/2)
  #Steplength constant
  a <- ((2*pi)^(1/2))*exp((zeta0^2)/2)/(2*zeta0)
  
  p <- length(par0)
  out <- array(dim = c(R, p))
  #Acceptance vector
  accepted <- 0
  #Initialization
  par <- par0; sigma <- sigma0; k <- k0
  
  cat("---------------------\n 
     CHAINING            \n")
  
  for(i in 1:R){
    #Proposal
    pars <- par + rmvnorm(1, sigma = k*sigma)
    #Acceptance proability
    eta <- min(1, gev.post(pars, data, hyp, t = t, p = p)/gev.post(par, data, hyp, t = t, p = p))
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
    if(i %% time == 0){
      cat(i, "samples: ")
      temp <- toc()
      timing[i / time] <- temp$toc - temp$tic
      if(i != R) tic()
    }
  }
  if(R %% time != 0){
    cat(R, " samples: ")
    temp <- toc()
    timing[length(timing)] <- temp$toc - temp$tic
  }
  cat("---------------------\n")
  list(values = out, accepted = accepted/R, timing = timing)
}

#Optimization: MLE /Bayesian Mode (empirical bayes prior)

fit.gev.inference <- function(data, t = NULL, t.prob = 0.95, llik.type = "Gev-Cens",
                              T.ret = 50, p = 1/length(data),
                              par0 = NULL, optim.meth = "Nelder-Mead", control = NULL,
                              hessian = FALSE,
                              #inf.type = c("Frequent", "Bayes")
                              inf.type = NULL, 
                              #Posterior chain size
                              R = 10^4,
                              #Starting RW steplength k
                              k = 3,
                              #Desired overall sampler acceptance probability
                              etastar = 0.234,
                              #time = desired interval of generated samples for which we
                              #require the timing
                              time = 500, 
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
  #Define the proportion of exceedances
  #Sample size
  n <- length(data)
  #Effective sample size
  k <- sum(data > t)
  exc.prop <- k/n
  names(exc.prop) <- c("Excesses proportion")
  #s = Reciprocal of proportion of exceedances
  s <- 1/exc.prop
  
  #Outcome setting
  optim.msg <- "Something went wrong"
  if(llik.type == "Gev-Cens"){
    optimfun <- function(par){
      return(gev.lik(par = par, data = data, t = t, p = t.prob, log = TRUE))
    }
  }
  #Compute the maximum consered likelihood estimates
  #Hessian = TRUE preferred for minimizations
  mle <- optim(par0, optimfun, method = optim.meth, control = control,
               hessian = hessian)
  #Original estimate
  mle.or <- mle$par
  
  #Define the estimated GEV parameters: GEV fitting for the maximum asymptotic
  #distribution --> I have to refer to this estimates for return values.
  #I need these special estimates in order to compute extreme quantiles.
  scale <- mle$par[3] * s^(mle$par[1])
  loc <- mle$par[2] + scale * (1 - s^(-mle$par[1]))/mle$par[1]
  mle$par[3] <- scale
  mle$par[2] <- loc
  names(mle$par) <- c("Shape", "Location", "Scale")

  #Frequentist inference setting
  if(inf.type == "Frequent"){
    
    optim.msg <- mle$convergence
    names(optim.msg) <- c("Optimization message")
    
    #Define the extreme quantiles
    Q.ext <- loc + scale * ((s*p)^(-mle$par[1]) - 1) / mle$par[1]
    #Taylor expansion from the following...?
    #Q.ext <- loc + scale * ((-log(1 - s*p))^(-mle$par[1]) - 1) / mle$par[1]
    names(Q.ext) <- c("Extreme-Quantile")
    #Define the return-level: G(R.lev) = 1 - 1/T.ex
    #Return level associated with the return period T.ex
    R.lev <- loc + scale * ((-log(1 - 1/T.ret))^(-mle$par[1]) - 1) / mle$par[1]
    names(R.lev) <- c("Return-Level")
    
    #Return object
    if(hessian == FALSE){
      return(list(inference = inf.type,
                  mle.or = mle.or,
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
                  mle.or = mle.or,
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
    #Empirical Bayes
    #Prior hyperparameters
    hyp <- mle$par[2:3]
    hyp.ret <- mle$par
    names(hyp.ret) <- c("gamma", "mu", "delta")
    #Optimization
    if(llik.type == "Gev-Cens"){
      optimfun <- function(par){
        return(gev.post(par = par, data = data, hyp = hyp, t = t, p = t.prob, log = TRUE))
      }
    }
    #Compute the posterior mode estimates
    #Hessian = TRUE preferred for minimizations
    mode <- optim(mle.or, optimfun, method = optim.meth, control = control,
                 hessian = hessian)
    optim.msg <- mode$convergence
    names(optim.msg) <- c("Optimization message")
    #Original estimate
    mode.or <- mode$par
    #Posterior negative inverse hessian
    sigma <- solve(-mode$hessian)
    #Posterior mode GEV parameters
    scale <- mode$par[3] * s^(mode$par[1])
    loc <- mode$par[2] + scale * (1 - s^(-mode$par[1]))/mode$par[1]
    mode$par[3] <- scale
    mode$par[2] <- loc
    names(mode$par) <- c("Shape", "Location", "Scale")
    #Posterior chain
    post.mcmc <- adamh(R = R, data = data, hyp = hyp,
                      par0 = mode.or, sigma0 = sigma, k0 = k,
                      t = t, p = p, etastar = etastar, time = time)
    #Acceptance rate
    post.acc <- post.mcmc$accepted
    #Posterior distribution of the original parameters
    post.sam <- post.mcmc$values
    post.sam.or <- post.sam
    colnames(post.sam.or) <- c("Shape", "Location", "Scale")
    #Posterior distributions of the GEV parameters
    scale <- apply(post.sam, 1, function(x) x[3]*s^(x[1]))
    post.sam[,3] <- scale
    loc <- apply(post.sam, 1, function(x) x[2] + x[3] *(1 - s^(-x[1]))/x[1])
    post.sam[,2] <- loc
    colnames(post.sam) <- c("Shape", "Location", "Scale")
    #Extreme quantile chain
    Q.ext <- loc + scale *((s*p)^(-post.sam[,1]))/post.sam[,1]
    #Return level
    R.lev <- loc + scale * ((-log(1 - 1/T.ret))^(-post.sam[,1]) - 1) / post.sam[,1]
    
    #Return object
    if(hessian == FALSE){
      return(list(inference = inf.type,
                  mle.or = mle.or,
                  mle = mle$par,
                  mode.or = mode.or,
                  mode = mode$par,
                  max.post = mode$value,
                  optim.msg = optim.msg,
                  emp.bayes.hyp = hyp.ret,
                  mcmc.acc = post.acc,
                  parameters.or = post.sam.or,
                  parameters = post.sam,
                  Q.extreme = Q.ext,
                  R.level = R.lev,
                  exc.prop = exc.prop,
                  timing = post.mcmc$timing))
    }
    else{
      #Preferred for minimizations
      return(list(inference = inf.type,
                  mle.or = mle.or,
                  mle = mle$par,
                  mode.or = mode.or,
                  mode = mode$par,
                  hessian = mode$hessian,
                  max.post = mode$value,
                  optim.msg = optim.msg,
                  emp.bayes.hyp = hyp.ret,
                  mcmc.acc = post.acc,
                  parameters.or = post.sam.or,
                  parameters = post.sam,
                  Q.extreme = Q.ext,
                  R.level = R.lev,
                  exc.prop = exc.prop,
                  timing = post.mcmc$timing))
    }
  }
}

#-------------------------------------------------------------------------------

#Simulation setting
#Sample size
n <- 1000
#Exceeding probability used in the estimation
t.prob <- 0.05
#Number of exceedances used in the estimation
k <- n * t.prob
#Block-size for the correspoding block maxima approach
m <- 1/t.prob
#Exceeding probability used for prediction
p <- 1/n
#Return period
T.ret <- 50

#Unit Fréchet example
#Fréchet parameters
shape <- 1 #Inverse of the tail index
loc <- 0   #Location
scale <- 1 #Scale
#Set norming constants
#Centering
bn <- (log(n/k) - log(n/k-1))^(-1/shape)
#Scaling
an <- 1/(shape*(n/k-1)) * (log(n/k)-log(n/k-1))^(-1/shape-1)

#Define the true parameters
true.par <- c(1/shape, bn, an)
#Define true extreme quantile
Q.ext.true <- qfrechet(1 - p, loc, scale, shape)
#Define true return level
R.lev.true <- qgev(1 - 1/T.ret, bn, an, 1/shape)

#Define the estimation setting
start <- c(0.1, 0.1, 1)
set.seed(1)
data <- rfrechet(n, loc, scale, shape)
th <- quantile(data, probs=1-t.prob, type=3)

#Frequentist inference
fit.freq <- fit.gev.inference(data = data, t = th, t.prob = t.prob,
                              llik.type = "Gev-Cens", T.ret = T.ret, p = p,
                              par0 = start, hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq
#Comparison
#Parameters
fit.freq$mle; true.par
#Extreme quantile
fit.freq$Q.extreme; Q.ext.true
#Return level
fit.freq$R.level; R.lev.true

#Bayesian inference
library(TeachingDemos)
library(tictoc)
tic()
fit.bayes <- fit.gev.inference(data = data, t = th, t.prob = t.prob,
                              llik.type = "Gev-Cens", T.ret = T.ret, p = p,
                              par0 = start, hessian = TRUE, inf.type = "Bayes",
                              R = 5000)
toc()
#Results
#fit.bayes
#Comparison
#Parameters
fit.bayes$mle; fit.bayes$mode; true.par
sapply(1:3, function(x) summary(fit.bayes$parameters[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(fit.bayes$parameters[,i]), fit.bayes$parameters[,i], type = "l",
       xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(fit.bayes$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(fit.bayes$parameters[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true
summary(fit.bayes$Q.extreme)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(fit.bayes$Q.extreme), fit.bayes$Q.extreme, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes$Q.extreme, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes$Q.extreme, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes$Q.extreme, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true
summary(fit.bayes$R.level)
#Trace plot
plot(1:length(fit.bayes$R.level), fit.bayes$R.level, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes$R.level, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes$R.level, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes$R.level, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))
par(mfrow = c(1,1))
#Timing for a fixed interval of samples
fit.bayes$timing

#-------------------------------------------------------------------------------

#Timing trials
library(tictoc)
#Overall posterior sample desired
R.tot <- 10000
#Number of posterior subsamples of length R.tot/R.trials
R.trials <- 10
#Total timing
t.tot <- 0
#split timing
t.trials <- rep(0, R.trials)
#Prameters graphical labels
lab <- c(expression(gamma), expression(mu), expression(delta))

#Timing comparison
#Total chain
tic()
set.seed(0)
fit.bayes.tot <- fit.gev.inference(data = data, t = th, t.prob = t.prob,
                               llik.type = "Gev-Cens", T.ret = T.ret, p = p,
                               par0 = start, hessian = TRUE, inf.type = "Bayes",
                               R = R.tot)
cat("\n", R.tot, " samples: ")
time <- toc()
t.tot <- time$toc - time$tic

#Unified chains
tot.par <- NULL
tot.q.ext <- c()
tot.r.lev <- c()
#Starting point
start <- start
for(i in 1:R.trials){
  tic()
  set.seed(i)
  fit.bayes <- fit.gev.inference(data = data, t = th, t.prob = t.prob,
                                 llik.type = "Gev-Cens", T.ret = T.ret, p = p,
                                 par0 = start, hessian = TRUE, inf.type = "Bayes",
                                 R = R.tot/R.trials)
  #Unification
  tot.par <- rbind(tot.par, fit.bayes$parameters)
  tot.q.ext <- c(tot.q.ext, fit.bayes$Q.extreme)
  tot.r.lev <- c(tot.r.lev, fit.bayes$R.level)
  #Starting from the last value generated by the previous iteration...
  start <- as.vector(tot.par[dim(tot.par)[1],])
  cat("\n", i*(R.tot/R.trials), " samples: ")
  time <- toc()
  t.trials[i] <- time$toc - time$tic
}

#Comparisons
par(mfrow = c(2,3))

#Parameters (gamma, mu, delta)
cat("\n Parameters: total chain")
summary(fit.bayes.tot$parameters)
cat("\n Parameters: unified chains")
summary(tot.par)
for(i in 1:3){
  #Total chain
  #Trace plot
  plot(1:length(fit.bayes.tot$parameters[,i]), fit.bayes.tot$parameters[,i], type = "l",
       xlab = lab[i], ylab = "MCMC chain", main = "Trace plot (total chain)")
  abline(h = true.par[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes.tot$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100, total chain)")
  #Histogram
  hist(fit.bayes.tot$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution (total chain)")
  abline(v = true.par[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(fit.bayes.tot$parameters[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  
  #Unified chains
  #Trace plot
  plot(1:length(tot.par[,i]), tot.par[,i], type = "l",
       xlab = lab[i], ylab = "MCMC chain", main = "Trace plot (unified chains)")
  abline(h = true.par[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(tot.par[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100, unified chains)")
  #Histogram
  hist(tot.par[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution (unified chains)")
  abline(v = true.par[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(tot.par[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
}

#Extreme quantile
cat("\n Extreme quantile: total chain")
summary(fit.bayes.tot$Q.extreme)
cat("\n Extreme quantile: unified chains")
summary(tot.q.ext)
#Total chain
#Trace plot
plot(1:length(fit.bayes.tot$Q.extreme), fit.bayes.tot$Q.extreme, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot (total chain)")
abline(h = Q.ext.true, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.tot$Q.extreme, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100, total chain)")
#Histogram
hist(fit.bayes.tot$Q.extreme, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution (total chain)")
abline(v = Q.ext.true, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.tot$Q.extreme, conf = 0.95), col = "red",
       lty = 3, lwd = 2)

#Unified chains
#Trace plot
plot(1:length(tot.q.ext), tot.q.ext, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot (unified chains)")
abline(h = Q.ext.true, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(tot.q.ext, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100, unified chains)")
#Histogram
hist(tot.q.ext, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution (unified chains)")
abline(v = Q.ext.true, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(tot.q.ext, conf = 0.95), col = "red",
       lty = 3, lwd = 2)

#Return level
cat("\n Return level: total chain")
summary(fit.bayes.tot$R.level)
cat("\n Return level: unified chains")
summary(tot.r.lev)
#Total chain
#Trace plot
plot(1:length(fit.bayes.tot$R.level), fit.bayes.tot$R.level, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot (total chain)")
abline(h = R.lev.true, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.tot$R.level, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100, total chain)")
#Histogram
hist(fit.bayes.tot$R.level, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution (total chain)")
abline(v = R.lev.true, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.tot$R.level, conf = 0.95), col = "red",
       lty = 3, lwd = 2)

#Unified chains
#Trace plot
plot(1:length(tot.r.lev), tot.r.lev, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot (unified chains)")
abline(h = Q.ext.true, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(tot.r.lev, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100, unified chains)")
#Histogram
hist(tot.r.lev, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution (unified chains)")
abline(v = Q.ext.true, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(tot.r.lev, conf = 0.95), col = "red",
       lty = 3, lwd = 2)

par(mfrow = c(1,1))

#Timing results
cat("\n", t.tot, "\n", t.trials, "\n", sum(t.trials))