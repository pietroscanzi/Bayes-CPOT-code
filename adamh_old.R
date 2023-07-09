#-------------------------------------------------------------------------------

path <- "C:\\Users\\scanz\\Desktop\\Magistrale\\Extreme Value Theory\\Code\\"
setwd(path)
source(paste0(path ,"functions.R"))

#-------------------------------------------------------------------------------

#Adaptive Random Walk Metropolis-Hastings Markov Chain Monte Carlo

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

#-------------------------------------------------------------------------------

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