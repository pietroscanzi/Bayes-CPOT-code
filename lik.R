#Comparison between:
#1) Classical GEV censored likelihood;
#2) GEV censored likelihood of the maximum.
#Comparison between:
#1) empirical posterior density obtained with the classical GEV likelihood;
#2) empirical posterior density obtained with the GEV likelihood of the maximum.

#-------------------------------------------------------------------------------

path <- path <- "C:\\Users\\scanz\\OneDrive - Università Commerciale Luigi Bocconi\\Desktop/Magistrale\\Extreme Value Theory\\Code\\"
setwd(path)
source(paste0(path ,"functions.R"))

#-------------------------------------------------------------------------------

#Simulation setting
library(LaplacesDemon)
library(TeachingDemos)
library(evd)
#library(extremefit)
library(tictoc)

tic()

#Sample size
n <- 1000
#Exceeding probability used in the estimation
t.prob <- 0.05
#Number of exceedances used in the estimation
k <- n * t.prob
#Block-size for the correspoding block maxima approach
m <- 1/t.prob
#Exceeding probability used for prediction --> extreme quantile
p <- 1/n
#Return period --> return level
T.ret <- 50
#Posterior sample size
R <- 50000
#Burn-in period
burn <- round(R/4)
#Confidence levels
conf.levels <- c(0, 0.5, 0.75, 0.9, 0.95, 0.99)
#Optimization method
method = "Nelder-Mead"
#Optimization control
control <- list(fnscale = -1, maxit = 8e+5)

#-------------------------------------------------------------------------------

#Profile likelihood
#(gamma, mu, \hat{delta}_{(gamma, mu)})
gev.lik.gammu <- function(par, data, t = NULL, p = 0.95, log = FALSE,
                          llik.type = "Max-Gev-Cens", param.type = "tilde"){
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
                          llik.type = "Max-Gev-Cens", param.type = "tilde"){
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
                          llik.type = "Max-Gev-Cens", param.type = "tilde"){
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

#Unit Fréchet (tail index = 1)

#Fréchet parameters
shape.fre <- 1 
loc.fre <- 0   
scale.fre <- 1 
#Set norming constants
#Centering
bn.fre <- (log(n/k) - log(n/k-1))^(-1/shape.fre)
#Scaling
an.fre <- 1/(shape.fre*(n/k-1)) * (log(n/k)-log(n/k-1))^(-1/shape.fre-1)
an2.fre <- (1/shape.fre)*bn.fre
#Define the true parameters
true.par.fre <- c(1/shape.fre, bn.fre, an.fre)

#Define the estimation setting
start.fre <- c(0.1, 0.1, 1)
set.seed(1)
data.fre <- rfrechet(n, loc.fre, scale.fre, shape.fre)
th.fre <- quantile(data.fre, probs = 1 - t.prob, type = 3)

#Likelihood

#Maximum likelihood
mle1.fre <- optim(start.fre, function(x) gev.lik(x, data.fre, t = th.fre,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle1.fre
mle2.fre <- optim(start.fre, function(x) gev.lik(x, data.fre, t = th.fre,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Max-Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle2.fre
mle.bar.fre <- optim(start.fre, function(x) gev.lik(x, data.fre, t = th.fre,
                                                   p = t.prob, log = TRUE,
                                                   llik.type = "Gev-Cens",
                                                   param.type = "bar"),
                     method = method, control = control, hessian = T); mle.bar.fre
#Same estimates
mle1.fre$par; mle2.fre$par; true.par.fre
#Bar parameterization
mle.bar.fre$par
#Standard errors
se1.fre <- sqrt(diag(solve(-mle1.fre$hessian))); se1.fre
se2.fre <- sqrt(diag(solve(-mle2.fre$hessian))); se2.fre
se.bar.fre <- sqrt(diag(solve(-mle.bar.fre$hessian))); se.bar.fre

#Contour plots
par(mfrow = c(3,1))

#gamma, mu
gam.fre <- seq(mle1.fre$par[1] - 3*se1.fre[1], mle1.fre$par[1] + 5*se1.fre[1],
               length = 50)
mu.fre <- seq(mle1.fre$par[2] - 4.5*se1.fre[2], mle1.fre$par[2] + 0.4*se1.fre[2],
              length = 50)
gam.bar.fre <- seq(mle.bar.fre$par[1] - 3*se.bar.fre[1], mle.bar.fre$par[1] + 5*se.bar.fre[1],
               length = 50)
mu.bar.fre <- seq(mle.bar.fre$par[2] - 8*se.bar.fre[2], mle.bar.fre$par[2] + 2*se.bar.fre[2],
              length = 50)
gammu.fre <- expand.grid(gam.fre, mu.fre)
gammu.bar.fre <- expand.grid(gam.bar.fre, mu.bar.fre)
gammu.lik.fre1 <- matrix(apply(gammu.fre,
                               1, function(x) gev.lik.gammu(par = x, data = data.fre,
                                                            t = th.fre,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                        nrow = 50, ncol = 50)
gammu.lik.fre2 <- matrix(apply(gammu.fre,
                               1, function(x) gev.lik.gammu(par = x, data = data.fre,
                                                            t = th.fre,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.bar.fre <- matrix(apply(gammu.bar.fre,
                              1, function(x) gev.lik.gammu(par = x, data = data.fre,
                                                           t = th.fre,
                                                           p = t.prob, log = TRUE,
                                                           llik.type = "Gev-Cens",
                                                           param.type = "bar")),
                        nrow = 50, ncol = 50)
contour(gam.fre, mu.fre, gammu.lik.fre1 - max(gammu.lik.fre1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.fre$par[1], mle1.fre$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.fre, mu.fre, gammu.lik.fre2 - max(gammu.lik.fre2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.fre$par[1], mle1.fre$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.fre, mu.bar.fre, gammu.lik.bar.fre - max(gammu.lik.bar.fre),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(mu)), 
        main = "Classical GEV censored relative log-likelihood (bar parameterization)")
points(mle.bar.fre$par[1], mle.bar.fre$par[2], col = "blue", cex = 1.5, lwd = 2)

#gamma, delta
gam.fre <- seq(mle1.fre$par[1] - 3*se1.fre[1], mle1.fre$par[1] + 5*se1.fre[1],
               length = 50)
del.fre <- seq(mle1.fre$par[3] - 2.5*se1.fre[3], mle1.fre$par[3] + 5*se1.fre[3],
               length = 50)
gam.bar.fre <- seq(mle.bar.fre$par[1] - 3*se.bar.fre[1], mle.bar.fre$par[1] + 5*se.bar.fre[1],
               length = 50)
del.bar.fre <- seq(mle.bar.fre$par[3] - 2.5*se.bar.fre[3], mle.bar.fre$par[3] + 8*se.bar.fre[3],
               length = 50)
gamdel.fre <- expand.grid(gam.fre, del.fre)
gamdel.bar.fre <- expand.grid(gam.bar.fre, del.bar.fre)
gamdel.lik.fre1 <- matrix(apply(gamdel.fre,
                               1, function(x) gev.lik.gamdel(par = x, data = data.fre,
                                                            t = th.fre,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gamdel.lik.fre2 <- matrix(apply(gamdel.fre,
                               1, function(x) gev.lik.gamdel(par = x, data = data.fre,
                                                            t = th.fre,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gamdel.lik.bar.fre <- matrix(apply(gamdel.bar.fre,
                               1, function(x) gev.lik.gamdel(par = x, data = data.fre,
                                                             t = th.fre,
                                                             p = t.prob, log = TRUE,
                                                             llik.type = "Gev-Cens",
                                                             param.type = "bar")),
                         nrow = 50, ncol = 50)
contour(gam.fre, del.fre, gamdel.lik.fre1 - max(gamdel.lik.fre1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.fre$par[1], mle1.fre$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.fre, del.fre, gamdel.lik.fre2 - max(gamdel.lik.fre2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.fre$par[1], mle1.fre$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.fre, del.bar.fre, gamdel.lik.bar.fre - max(gamdel.lik.bar.fre),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.fre$par[1], mle.bar.fre$par[3], col = "blue", cex = 1.5, lwd = 2)

#mu, delta
mu.fre <- seq(mle1.fre$par[2] - se1.fre[2], mle1.fre$par[2] + se1.fre[2],
              length = 50)
del.fre <- seq(mle1.fre$par[3] - 2*se1.fre[3], mle1.fre$par[3] + 5*se1.fre[3],
               length = 50)
mu.bar.fre <- seq(mle.bar.fre$par[2] - 6*se.bar.fre[2], mle.bar.fre$par[2] + 2*se.bar.fre[2],
              length = 50)
del.bar.fre <- seq(mle.bar.fre$par[3] - 1.5*se.bar.fre[3], mle.bar.fre$par[3] + 8*se.bar.fre[3],
               length = 50)
mudel.fre <- expand.grid(mu.fre, del.fre)
mudel.bar.fre <- expand.grid(mu.bar.fre, del.bar.fre)
mudel.lik.fre1 <- matrix(apply(mudel.fre,
                                1, function(x) gev.lik.mudel(par = x, data = data.fre,
                                                             t = th.fre,
                                                             p = t.prob, log = TRUE,
                                                             llik.type = "Gev-Cens",
                                                             param.type = "tilde")),
                          nrow = 50, ncol = 50)
mudel.lik.fre2 <- matrix(apply(mudel.fre,
                                1, function(x) gev.lik.mudel(par = x, data = data.fre,
                                                             t = th.fre,
                                                             p = t.prob, log = TRUE,
                                                             llik.type = "Max-Gev-Cens",
                                                             param.type = "tilde")),
                          nrow = 50, ncol = 50)
mudel.lik.bar.fre <- matrix(apply(mudel.bar.fre,
                               1, function(x) gev.lik.mudel(par = x, data = data.fre,
                                                            t = th.fre,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "bar")),
                         nrow = 50, ncol = 50)
contour(mu.fre, del.fre, mudel.lik.fre1 - max(mudel.lik.fre1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.fre$par[2], mle1.fre$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.fre, del.fre, mudel.lik.fre2 - max(mudel.lik.fre2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.fre$par[2], mle1.fre$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.bar.fre, del.bar.fre, mudel.lik.bar.fre - max(mudel.lik.bar.fre),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(mu)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.fre$par[2], mle.bar.fre$par[3], col = "blue", cex = 1.5, lwd = 2)

#POSTERIOR CHAINS
set.seed(1)
fit.bayes.fre1 <- fit.gev.inference(data = data.fre, t = th.fre, t.prob = t.prob,
                                   llik.type = "Gev-Cens", param.type = "tilde",
                                   T.ret = T.ret, p = p, par0 = start.fre,
                                   hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
pars1 <- fit.bayes.fre1$parameters[-burn,]
fit.bayes.fre2 <- fit.gev.inference(data = data.fre, t = th.fre, t.prob = t.prob,
                                    llik.type = "Max-Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.fre,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars2 <- fit.bayes.fre2$parameters[-burn,]
fit.bayes.bar.fre <- fit.gev.inference(data = data.fre, t = th.fre, t.prob = t.prob,
                                    llik.type = "Gev-Cens", param.type = "bar",
                                    T.ret = T.ret, p = p, par0 = start.fre,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars.bar <- fit.bayes.bar.fre$parameters[-burn,]
#smoothScatter plots
par(mfrow = c(3,1))

#gamma, mu
smoothScatter(pars1[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,1:2])
smoothScatter(pars2[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,1:2])
smoothScatter(pars.bar[,1:2], xlab = expression(bar(gamma)),
              ylab = expression(bar(mu)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,1:2])

#gamma, delta
smoothScatter(pars1[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(1,3)])
smoothScatter(pars2[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(1,3)])
smoothScatter(pars.bar[,c(1,3)], xlab = expression(bar(gamma)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(1,3)])

#mu, delta
smoothScatter(pars1[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(2,3)])
smoothScatter(pars2[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(2,3)])
smoothScatter(pars.bar[,c(2,3)], xlab = expression(bar(mu)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(2,3)])

par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#standard Pareto (tail index = 1)

#standard pareto parameters
shape.par <- 1 
loc.par <- 0   
scale.par <- 1 
#Set norming constants
#Centering
bn.par <- (n/k)^(1/shape.par)
#Scaling
an.par <- (1/shape.par)*((n/k)^(1/shape.par))
an2.par <- (1/shape.par)*((n/k)^(1/shape.par))

#Define the true parameters
true.par.par <- c(1/shape.par, bn.par, an.par)
#Define true extreme quantile
Q.ext.true.par <- extremefit::qpareto(1 - p, shape.par, loc.par, scale.par)
#Define true return level
R.lev.true.par <- qgev(1 - 1/T.ret, bn.par, an.par, 1/shape.par)

#Define the estimation setting
start.par <- c(0.1, 0.1, 1)
set.seed(2)
data.par <- extremefit::rpareto(n, shape.par, loc.par, scale.par)
th.par <- quantile(data.par, probs = 1 - t.prob, type = 3)

#Likelihood

#Maximum likelihood
mle1.par <- optim(start.par, function(x) gev.lik(x, data.par, t = th.par,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle1.par
mle2.par <- optim(start.par, function(x) gev.lik(x, data.par, t = th.par,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Max-Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle2.par
mle.bar.par <- optim(start.par, function(x) gev.lik(x, data.par, t = th.par,
                                                   p = t.prob, log = TRUE,
                                                   llik.type = "Gev-Cens",
                                                   param.type = "bar"),
                     method = method, control = control, hessian = T); mle.bar.par
#Same estimates
mle1.par$par; mle2.par$par; true.par.par
#Bar parameterization
mle.bar.par$par
#Standard errors
se1.par <- sqrt(diag(solve(-mle1.par$hessian))); se1.par
se2.par <- sqrt(diag(solve(-mle2.par$hessian))); se2.par
se.bar.par <- sqrt(diag(solve(-mle.bar.par$hessian))); se.bar.par

#Contour plots
par(mfrow = c(3,1))

#gamma, mu
gam.par <- seq(mle1.par$par[1] - 2*se1.par[1], mle1.par$par[1] + 4*se1.par[1],
               length = 50)
mu.par <- seq(mle1.par$par[2] - 4*se1.par[2], mle1.par$par[2] + 0.5*se1.par[2],
              length = 50)
gam.bar.par <- seq(mle.bar.par$par[1] - 3*se.bar.par[1], mle.bar.par$par[1] + 5*se.bar.par[1],
                   length = 50)
mu.bar.par <- seq(mle.bar.par$par[2] - 6*se.bar.par[2], mle.bar.par$par[2] + 2*se.bar.par[2],
                  length = 50)
gammu.par <- expand.grid(gam.par, mu.par)
gammu.bar.par <- expand.grid(gam.bar.par, mu.bar.par)
gammu.lik.par1 <- matrix(apply(gammu.par,
                               1, function(x) gev.lik.gammu(par = x, data = data.par,
                                                            t = th.par,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.par2 <- matrix(apply(gammu.par,
                               1, function(x) gev.lik.gammu(par = x, data = data.par,
                                                            t = th.par,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.bar.par <- matrix(apply(gammu.bar.par,
                                  1, function(x) gev.lik.gammu(par = x, data = data.par,
                                                               t = th.par,
                                                               p = t.prob, log = TRUE,
                                                               llik.type = "Gev-Cens",
                                                               param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(gam.par, mu.par, gammu.lik.par1 - max(gammu.lik.par1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.par$par[1], mle1.par$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.par, mu.par, gammu.lik.par2 - max(gammu.lik.par2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.par$par[1], mle1.par$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.par, mu.bar.par, gammu.lik.bar.par - max(gammu.lik.bar.par),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(mu)), 
        main = "Classical GEV censored relative log-likelihood (bar parameterization)")
points(mle.bar.par$par[1], mle.bar.par$par[2], col = "blue", cex = 1.5, lwd = 2)

#gamma, delta
gam.par <- seq(mle1.par$par[1] - 2.5*se1.par[1], mle1.par$par[1] + 5*se1.par[1],
               length = 50)
del.par <- seq(mle1.par$par[3] - 2.5*se1.par[3], mle1.par$par[3] + 5*se1.par[3],
               length = 50)
gam.bar.par <- seq(mle.bar.par$par[1] - 2.5*se.bar.par[1], mle.bar.par$par[1] + 4*se.bar.par[1],
                   length = 50)
del.bar.par <- seq(mle.bar.par$par[3] - 2.5*se.bar.par[3], mle.bar.par$par[3] + 8*se.bar.par[3],
                   length = 50)
gamdel.par <- expand.grid(gam.par, del.par)
gamdel.bar.par <- expand.grid(gam.bar.par, del.bar.par)
gamdel.lik.par1 <- matrix(apply(gamdel.par,
                                1, function(x) gev.lik.gamdel(par = x, data = data.par,
                                                              t = th.par,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.par2 <- matrix(apply(gamdel.par,
                                1, function(x) gev.lik.gamdel(par = x, data = data.par,
                                                              t = th.par,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Max-Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.bar.par <- matrix(apply(gamdel.bar.par,
                                   1, function(x) gev.lik.gamdel(par = x, data = data.par,
                                                                 t = th.par,
                                                                 p = t.prob, log = TRUE,
                                                                 llik.type = "Gev-Cens",
                                                                 param.type = "bar")),
                             nrow = 50, ncol = 50)
contour(gam.par, del.par, gamdel.lik.par1 - max(gamdel.lik.par1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.par$par[1], mle1.par$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.par, del.par, gamdel.lik.par2 - max(gamdel.lik.par2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.par$par[1], mle1.par$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.par, del.bar.par, gamdel.lik.bar.par - max(gamdel.lik.bar.par),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.par$par[1], mle.bar.par$par[3], col = "blue", cex = 1.5, lwd = 2)

#mu, delta
mu.par <- seq(mle1.par$par[2] - se1.par[2], mle1.par$par[2] + se1.par[2],
              length = 50)
del.par <- seq(mle1.par$par[3] - 3*se1.par[3], mle1.par$par[3] + 5*se1.par[3],
               length = 50)
mu.bar.par <- seq(mle.bar.par$par[2] - 6*se.bar.par[2], mle.bar.par$par[2] + 2*se.bar.par[2],
                  length = 50)
del.bar.par <- seq(mle.bar.par$par[3] - 1.5*se.bar.par[3], mle.bar.par$par[3] + 8*se.bar.par[3],
                   length = 50)
mudel.par <- expand.grid(mu.par, del.par)
mudel.bar.par <- expand.grid(mu.bar.par, del.bar.par)
mudel.lik.par1 <- matrix(apply(mudel.par,
                               1, function(x) gev.lik.mudel(par = x, data = data.par,
                                                            t = th.par,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
mudel.lik.par2 <- matrix(apply(mudel.par,
                               1, function(x) gev.lik.mudel(par = x, data = data.par,
                                                            t = th.par,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
mudel.lik.bar.par <- matrix(apply(mudel.bar.par,
                                  1, function(x) gev.lik.mudel(par = x, data = data.par,
                                                               t = th.par,
                                                               p = t.prob, log = TRUE,
                                                               llik.type = "Gev-Cens",
                                                               param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(mu.par, del.par, mudel.lik.par1 - max(mudel.lik.par1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.par$par[2], mle1.par$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.par, del.par, mudel.lik.par2 - max(mudel.lik.par2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.par$par[2], mle1.par$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.bar.par, del.bar.par, mudel.lik.bar.par - max(mudel.lik.bar.par),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(mu)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.par$par[2], mle.bar.par$par[3], col = "blue", cex = 1.5, lwd = 2)

#POSTERIOR CHAINS
set.seed(2)
fit.bayes.par1 <- fit.gev.inference(data = data.par, t = th.par, t.prob = t.prob,
                                    llik.type = "Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.par,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars1 <- fit.bayes.par1$parameters[-burn,]
fit.bayes.par2 <- fit.gev.inference(data = data.par, t = th.par, t.prob = t.prob,
                                    llik.type = "Max-Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.par,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars2 <- fit.bayes.par2$parameters[-burn,]
fit.bayes.bar.par <- fit.gev.inference(data = data.par, t = th.par, t.prob = t.prob,
                                       llik.type = "Gev-Cens", param.type = "bar",
                                       T.ret = T.ret, p = p, par0 = start.par,
                                       hessian = TRUE, inf.type = "Bayes",
                                       k = 1, R = R, burn  = burn, prior = "empirical",
                                       val.show = TRUE)
pars.bar <- fit.bayes.bar.par$parameters[-burn,]
#smoothScatter plots
par(mfrow = c(3,1))

#gamma, mu
smoothScatter(pars1[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,1:2])
smoothScatter(pars2[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,1:2])
smoothScatter(pars.bar[,1:2], xlab = expression(bar(gamma)),
              ylab = expression(bar(mu)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,1:2])

#gamma, delta
smoothScatter(pars1[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(1,3)])
smoothScatter(pars2[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(1,3)])
smoothScatter(pars.bar[,c(1,3)], xlab = expression(bar(gamma)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(1,3)])

#mu, delta
smoothScatter(pars1[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(2,3)])
smoothScatter(pars2[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(2,3)])
smoothScatter(pars.bar[,c(2,3)], xlab = expression(bar(mu)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(2,3)])

par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#half Cauchy (tail index = 1)

#half Cauchy parameters
scale.hc <- 1 
#Set norming constants
#Centering
bn.hc <- tan((pi/2)*(1 - (k/n)))
bn.hc2 <- (2*n)/(pi*k)
#Scaling
an.hc <- (pi/2)*(k/n)*(1 + (tan((pi/2)*(1 - (k/n))))^2)
an2.hc <- bn.hc
an3.hc <- (2*n)/(pi*k)

#Define the true parameters
true.par.hc <- c(1/scale.hc, bn.hc, an.hc)
#Define true extreme quantile
Q.ext.true.hc <- qhalfcauchy(1 - p, scale.hc)
#Define true return level
R.lev.true.hc <- qgev(1 - 1/T.ret, bn.hc, an.hc, 1/scale.hc)

#Define the estimation setting
start.hc <- c(0.5, 2, 5)
set.seed(2)
data.hc <- rhalfcauchy(n, scale.hc)
th.hc <- quantile(data.hc, probs = 1 - t.prob, type = 3)

#Likelihood

#Maximum likelihood
mle1.hc <- optim(start.hc, function(x) gev.lik(x, data.hc, t = th.hc,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle1.hc
mle2.hc <- optim(start.hc, function(x) gev.lik(x, data.hc, t = th.hc,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Max-Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle2.hc
mle.bar.hc <- optim(start.hc, function(x) gev.lik(x, data.hc, t = th.hc,
                                                   p = t.prob, log = TRUE,
                                                   llik.type = "Gev-Cens",
                                                   param.type = "bar"),
                     method = method, control = control, hessian = T); mle.bar.hc
#Same estimates
mle1.hc$par; mle2.hc$par; true.par.hc
#Bar parameterization
mle.bar.hc$par
#Standard errors
se1.hc <- sqrt(diag(solve(-mle1.hc$hessian))); se1.hc
se2.hc <- sqrt(diag(solve(-mle2.hc$hessian))); se2.hc
se.bar.hc <- sqrt(diag(solve(-mle.bar.hc$hessian))); se.bar.hc

#Contour plots
par(mfrow = c(3,1))

#gamma, mu
gam.hc <- seq(mle1.hc$par[1] - 3*se1.hc[1], mle1.hc$par[1] + 5*se1.hc[1],
               length = 50)
mu.hc <- seq(mle1.hc$par[2] - 4*se1.hc[2], mle1.hc$par[2] + 0.5*se1.hc[2],
              length = 50)
gam.bar.hc <- seq(mle.bar.hc$par[1] - 2*se.bar.hc[1], mle.bar.hc$par[1] + 4*se.bar.hc[1],
                   length = 50)
mu.bar.hc <- seq(mle.bar.hc$par[2] - 6.5*se.bar.hc[2], mle.bar.hc$par[2] + 2.5*se.bar.hc[2],
                  length = 50)
gammu.hc <- expand.grid(gam.hc, mu.hc)
gammu.bar.hc <- expand.grid(gam.bar.hc, mu.bar.hc)
gammu.lik.hc1 <- matrix(apply(gammu.hc,
                               1, function(x) gev.lik.gammu(par = x, data = data.hc,
                                                            t = th.hc,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.hc2 <- matrix(apply(gammu.hc,
                               1, function(x) gev.lik.gammu(par = x, data = data.hc,
                                                            t = th.hc,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.bar.hc <- matrix(apply(gammu.bar.hc,
                                  1, function(x) gev.lik.gammu(par = x, data = data.hc,
                                                               t = th.hc,
                                                               p = t.prob, log = TRUE,
                                                               llik.type = "Gev-Cens",
                                                               param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(gam.hc, mu.hc, gammu.lik.hc1 - max(gammu.lik.hc1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.hc$par[1], mle1.hc$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.hc, mu.hc, gammu.lik.hc2 - max(gammu.lik.hc2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.hc$par[1], mle1.hc$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.hc, mu.bar.hc, gammu.lik.bar.hc - max(gammu.lik.bar.hc),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(mu)), 
        main = "Classical GEV censored relative log-likelihood (bar parameterization)")
points(mle.bar.hc$par[1], mle.bar.hc$par[2], col = "blue", cex = 1.5, lwd = 2)

#gamma, delta
gam.hc <- seq(mle1.hc$par[1] - 3*se1.hc[1], mle1.hc$par[1] + 5*se1.hc[1],
               length = 50)
del.hc <- seq(mle1.hc$par[3] - 2.5*se1.hc[3], mle1.hc$par[3] + 5*se1.hc[3],
               length = 50)
gam.bar.hc <- seq(mle.bar.hc$par[1] - 2.5*se.bar.hc[1], mle.bar.hc$par[1] + 6*se.bar.hc[1],
                   length = 50)
del.bar.hc <- seq(mle.bar.hc$par[3] - 2.5*se.bar.hc[3], mle.bar.hc$par[3] + 7*se.bar.hc[3],
                   length = 50)
gamdel.hc <- expand.grid(gam.hc, del.hc)
gamdel.bar.hc <- expand.grid(gam.bar.hc, del.bar.hc)
gamdel.lik.hc1 <- matrix(apply(gamdel.hc,
                                1, function(x) gev.lik.gamdel(par = x, data = data.hc,
                                                              t = th.hc,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.hc2 <- matrix(apply(gamdel.hc,
                                1, function(x) gev.lik.gamdel(par = x, data = data.hc,
                                                              t = th.hc,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Max-Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.bar.hc <- matrix(apply(gamdel.bar.hc,
                                   1, function(x) gev.lik.gamdel(par = x, data = data.hc,
                                                                 t = th.hc,
                                                                 p = t.prob, log = TRUE,
                                                                 llik.type = "Gev-Cens",
                                                                 param.type = "bar")),
                             nrow = 50, ncol = 50)
contour(gam.hc, del.hc, gamdel.lik.hc1 - max(gamdel.lik.hc1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.hc$par[1], mle1.hc$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.hc, del.hc, gamdel.lik.hc2 - max(gamdel.lik.hc2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.hc$par[1], mle1.hc$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.hc, del.bar.hc, gamdel.lik.bar.hc - max(gamdel.lik.bar.hc),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.hc$par[1], mle.bar.hc$par[3], col = "blue", cex = 1.5, lwd = 2)

#mu, delta
mu.hc <- seq(mle1.hc$par[2] - se1.hc[2], mle1.hc$par[2] + se1.hc[2],
             length = 50)
del.hc <- seq(mle1.hc$par[3] - 3*se1.hc[3], mle1.hc$par[3] + 5*se1.hc[3],
              length = 50)
mu.bar.hc <- seq(mle.bar.hc$par[2] - 6*se.bar.hc[2], mle.bar.hc$par[2] + 2*se.bar.hc[2],
                 length = 50)
del.bar.hc <- seq(mle.bar.hc$par[3] - 1.5*se.bar.hc[3], mle.bar.hc$par[3] + 8*se.bar.hc[3],
                  length = 50)
mudel.hc <- expand.grid(mu.hc, del.hc)
mudel.bar.hc <- expand.grid(mu.bar.hc, del.bar.hc)
mudel.lik.hc1 <- matrix(apply(mudel.hc,
                               1, function(x) gev.lik.mudel(par = x, data = data.hc,
                                                            t = th.hc,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
mudel.lik.hc2 <- matrix(apply(mudel.hc,
                               1, function(x) gev.lik.mudel(par = x, data = data.hc,
                                                            t = th.hc,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
mudel.lik.bar.hc <- matrix(apply(mudel.bar.hc,
                                  1, function(x) gev.lik.mudel(par = x, data = data.hc,
                                                               t = th.hc,
                                                               p = t.prob, log = TRUE,
                                                               llik.type = "Gev-Cens",
                                                               param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(mu.hc, del.hc, mudel.lik.hc1 - max(mudel.lik.hc1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.hc$par[2], mle1.hc$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.hc, del.hc, mudel.lik.hc2 - max(mudel.lik.hc2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.hc$par[2], mle1.hc$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.bar.hc, del.bar.hc, mudel.lik.bar.hc - max(mudel.lik.bar.hc),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(mu)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.hc$par[2], mle.bar.hc$par[3], col = "blue", cex = 1.5, lwd = 2)

#POSTERIOR CHAINS
set.seed(2)
fit.bayes.hc1 <- fit.gev.inference(data = data.hc, t = th.hc, t.prob = t.prob,
                                    llik.type = "Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.hc,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars1 <- fit.bayes.hc1$parameters[-burn,]
fit.bayes.hc2 <- fit.gev.inference(data = data.hc, t = th.hc, t.prob = t.prob,
                                    llik.type = "Max-Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.hc,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars2 <- fit.bayes.hc2$parameters[-burn,]
fit.bayes.bar.hc <- fit.gev.inference(data = data.hc, t = th.hc, t.prob = t.prob,
                                       llik.type = "Gev-Cens", param.type = "bar",
                                       T.ret = T.ret, p = p, par0 = start.hc,
                                       hessian = TRUE, inf.type = "Bayes",
                                       k = 1, R = R, burn  = burn, prior = "empirical",
                                       val.show = TRUE)
pars.bar <- fit.bayes.bar.hc$parameters[-burn,]
#smoothScatter plots
par(mfrow = c(3,1))

#gamma, mu
smoothScatter(pars1[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,1:2])
smoothScatter(pars2[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,1:2])
smoothScatter(pars.bar[,1:2], xlab = expression(bar(gamma)),
              ylab = expression(bar(mu)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,1:2])

#gamma, delta
smoothScatter(pars1[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(1,3)])
smoothScatter(pars2[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(1,3)])
smoothScatter(pars.bar[,c(1,3)], xlab = expression(bar(gamma)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(1,3)])

#mu, delta
smoothScatter(pars1[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(2,3)])
smoothScatter(pars2[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(2,3)])
smoothScatter(pars.bar[,c(2,3)], xlab = expression(bar(mu)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(2,3)])

par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#standard Gumbel (tail index = 0)

#standard Gumbel parameters
loc.gum <- 0   
scale.gum <- 1 
#Set norming constants
#Centering
bn.gum <- -log(log(n) - log(n - k))
#Scaling
an.gum <- k/((n - k)*(log(n) - log(n - k)))
an2.gum <- bn.gum - (k/n)*integrate(function(x) -log(-log(1 - (1/x))), lower = 1,
                                    upper = (n/k))$value

#Define the true parameters
true.par.gum <- c(0, bn.gum, an.gum)
#Define true extreme quantile
Q.ext.true.gum <- qgumbel(1 - p, loc.gum, scale.gum)
#Define true return level
R.lev.true.gum <- qgev(1 - 1/T.ret, bn.gum, an.gum, 0)

#Define the estimation setting
start.gum <- c(0.1, 0.1, 1)
set.seed(5)
data.gum <- rgumbel(n, loc.gum, scale.gum)
th.gum <- quantile(data.gum, probs = 1 - t.prob, type = 3)

#Likelihood

#Maximum likelihood
mle1.gum <- optim(start.gum, function(x) gev.lik(x, data.gum, t = th.gum,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle1.gum
mle2.gum <- optim(start.gum, function(x) gev.lik(x, data.gum, t = th.gum,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Max-Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle2.gum
mle.bar.gum <- optim(start.gum, function(x) gev.lik(x, data.gum, t = th.gum,
                                                   p = t.prob, log = TRUE,
                                                   llik.type = "Gev-Cens",
                                                   param.type = "bar"),
                     method = method, control = control, hessian = T); mle.bar.gum
#Same estimates
mle1.gum$par; mle2.gum$par; true.par.gum
#Bar parameterization
mle.bar.gum$par
#Standard errors
se1.gum <- sqrt(diag(solve(-mle1.gum$hessian))); se1.gum
se2.gum <- sqrt(diag(solve(-mle2.gum$hessian))); se2.gum
se.bar.gum <- sqrt(diag(solve(-mle.bar.gum$hessian))); se.bar.gum

#Contour plots
par(mfrow = c(3,1))

#gamma, mu
gam.gum <- seq(mle1.gum$par[1] - 3*se1.gum[1], mle1.gum$par[1] + 5*se1.gum[1],
               length = 50)
mu.gum <- seq(mle1.gum$par[2] - 5*se1.gum[2], mle1.gum$par[2] + 5*se1.gum[2],
              length = 50)
gam.bar.gum <- seq(mle.bar.gum$par[1] - 2*se.bar.gum[1], mle.bar.gum$par[1] + 4*se.bar.gum[1],
                   length = 50)
mu.bar.gum <- seq(mle.bar.gum$par[2] - 4*se.bar.gum[2], mle.bar.gum$par[2] + 2*se.bar.gum[2],
                  length = 50)
gammu.gum <- expand.grid(gam.gum, mu.gum)
gammu.bar.gum <- expand.grid(gam.bar.gum, mu.bar.gum)
gammu.lik.gum1 <- matrix(apply(gammu.gum,
                               1, function(x) gev.lik.gammu(par = x, data = data.gum,
                                                            t = th.gum,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.gum2 <- matrix(apply(gammu.gum,
                               1, function(x) gev.lik.gammu(par = x, data = data.gum,
                                                            t = th.gum,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.bar.gum <- matrix(apply(gammu.bar.gum,
                                  1, function(x) gev.lik.gammu(par = x, data = data.gum,
                                                               t = th.gum,
                                                               p = t.prob, log = TRUE,
                                                               llik.type = "Gev-Cens",
                                                               param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(gam.gum, mu.gum, gammu.lik.gum1 - max(gammu.lik.gum1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.gum$par[1], mle1.gum$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.gum, mu.gum, gammu.lik.gum2 - max(gammu.lik.gum2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.gum$par[1], mle1.gum$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.gum, mu.bar.gum, gammu.lik.bar.gum - max(gammu.lik.bar.gum),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(mu)), 
        main = "Classical GEV censored relative log-likelihood (bar parameterization)")
points(mle.bar.gum$par[1], mle.bar.gum$par[2], col = "blue", cex = 1.5, lwd = 2)

#gamma, delta
gam.gum <- seq(mle1.gum$par[1] - 3*se1.gum[1], mle1.gum$par[1] + 5*se1.gum[1],
               length = 50)
del.gum <- seq(mle1.gum$par[3] - 2.5*se1.gum[3], mle1.gum$par[3] + 5*se1.gum[3],
               length = 50)
gam.bar.gum <- seq(mle.bar.gum$par[1] - 3*se.bar.gum[1], mle.bar.gum$par[1] + 5*se.bar.gum[1],
                   length = 50)
del.bar.gum <- seq(mle.bar.gum$par[3] - 2.5*se.bar.gum[3], mle.bar.gum$par[3] + 8*se.bar.gum[3],
                   length = 50)
gamdel.gum <- expand.grid(gam.gum, del.gum)
gamdel.bar.gum <- expand.grid(gam.bar.gum, del.bar.gum)
gamdel.lik.gum1 <- matrix(apply(gamdel.gum,
                                1, function(x) gev.lik.gamdel(par = x, data = data.gum,
                                                              t = th.gum,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.gum2 <- matrix(apply(gamdel.gum,
                                1, function(x) gev.lik.gamdel(par = x, data = data.gum,
                                                              t = th.gum,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Max-Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.bar.gum <- matrix(apply(gamdel.bar.gum,
                                   1, function(x) gev.lik.gamdel(par = x, data = data.gum,
                                                                 t = th.gum,
                                                                 p = t.prob, log = TRUE,
                                                                 llik.type = "Gev-Cens",
                                                                 param.type = "bar")),
                             nrow = 50, ncol = 50)
contour(gam.gum, del.gum, gamdel.lik.gum1 - max(gamdel.lik.gum1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.gum$par[1], mle1.gum$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.gum, del.gum, gamdel.lik.gum2 - max(gamdel.lik.gum2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.gum$par[1], mle1.gum$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.gum, del.bar.gum, gamdel.lik.bar.gum - max(gamdel.lik.bar.gum),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.gum$par[1], mle.bar.gum$par[3], col = "blue", cex = 1.5, lwd = 2)

#mu, delta
mu.gum <- seq(mle1.gum$par[2] - 2*se1.gum[2], mle1.gum$par[2] + 2*se1.gum[2],
              length = 50)
del.gum <- seq(mle1.gum$par[3] - 3*se1.gum[3], mle1.gum$par[3] + 5*se1.gum[3],
               length = 50)
mu.bar.gum <- seq(mle.bar.gum$par[2] - se.bar.gum[2], mle.bar.gum$par[2] + 2*se.bar.gum[2],
                  length = 50)
del.bar.gum <- seq(mle.bar.gum$par[3] - 2.5*se.bar.gum[3], mle.bar.gum$par[3] + se.bar.gum[3],
                   length = 50)
mudel.gum <- expand.grid(mu.gum, del.gum)
mudel.bar.gum <- expand.grid(mu.bar.gum, del.bar.gum)
mudel.lik.gum1 <- matrix(apply(mudel.gum,
                               1, function(x) gev.lik.mudel(par = x, data = data.gum,
                                                            t = th.gum,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
mudel.lik.gum2 <- matrix(apply(mudel.gum,
                               1, function(x) gev.lik.mudel(par = x, data = data.gum,
                                                            t = th.gum,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
mudel.lik.bar.gum <- matrix(apply(mudel.bar.gum,
                                  1, function(x) gev.lik.mudel(par = x, data = data.gum,
                                                               t = th.gum,
                                                               p = t.prob, log = TRUE,
                                                               llik.type = "Gev-Cens",
                                                               param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(mu.gum, del.gum, mudel.lik.gum1 - max(mudel.lik.gum1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.gum$par[2], mle1.gum$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.gum, del.gum, mudel.lik.gum2 - max(mudel.lik.gum2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.gum$par[2], mle1.gum$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.bar.gum, del.bar.gum, mudel.lik.bar.gum - max(mudel.lik.bar.gum),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(mu)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.gum$par[2], mle.bar.gum$par[3], col = "blue", cex = 1.5, lwd = 2)

#POSTERIOR CHAINS
set.seed(5)
fit.bayes.gum1 <- fit.gev.inference(data = data.gum, t = th.gum, t.prob = t.prob,
                                    llik.type = "Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.gum,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars1 <- fit.bayes.gum1$parameters[-burn,]
fit.bayes.gum2 <- fit.gev.inference(data = data.gum, t = th.gum, t.prob = t.prob,
                                    llik.type = "Max-Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.gum,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars2 <- fit.bayes.gum2$parameters[-burn,]
fit.bayes.bar.gum <- fit.gev.inference(data = data.gum, t = th.gum, t.prob = t.prob,
                                       llik.type = "Gev-Cens", param.type = "bar",
                                       T.ret = T.ret, p = p, par0 = start.gum,
                                       hessian = TRUE, inf.type = "Bayes",
                                       k = 1, R = R, burn  = burn, prior = "empirical",
                                       val.show = TRUE)
pars.bar <- fit.bayes.bar.gum$parameters[-burn,]
#smoothScatter plots
par(mfrow = c(3,1))

#gamma, mu
smoothScatter(pars1[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,1:2])
smoothScatter(pars2[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,1:2])
smoothScatter(pars.bar[,1:2], xlab = expression(bar(gamma)),
              ylab = expression(bar(mu)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,1:2])

#gamma, delta
smoothScatter(pars1[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(1,3)])
smoothScatter(pars2[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(1,3)])
smoothScatter(pars.bar[,c(1,3)], xlab = expression(bar(gamma)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(1,3)])

#mu, delta
smoothScatter(pars1[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(2,3)])
smoothScatter(pars2[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(2,3)])
smoothScatter(pars.bar[,c(2,3)], xlab = expression(bar(mu)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(2,3)])

par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#Exponential(1) (tail index = 0)

#Exponential parameters
rate.exp <- 1
#Set norming constants
#Centering
bn.exp <- log(n) - log(k)
#Scaling
an.exp <- 1
an2.exp <- 1

#Define the true parameters
true.par.exp <- c(0, bn.exp, an.exp)
#Define true extreme quantile
Q.ext.true.exp <- qexp(1 - p, rate.exp)
#Define true return level
R.lev.true.exp <- qgev(1 - 1/T.ret, bn.exp, an.exp, 0)

#Define the estimation setting
start.exp <- c(0.1, 1, 1)
set.seed(12)
data.exp <- rexp(n, rate.exp)
th.exp <- quantile(data.exp, probs = 1 - t.prob, type = 3)

#Likelihood

#Maximum likelihood
mle1.exp <- optim(start.exp, function(x) gev.lik(x, data.exp, t = th.exp,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle1.exp
mle2.exp <- optim(start.exp, function(x) gev.lik(x, data.exp, t = th.exp,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Max-Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle2.exp
mle.bar.exp <- optim(start.exp, function(x) gev.lik(x, data.exp, t = th.exp,
                                                   p = t.prob, log = TRUE,
                                                   llik.type = "Gev-Cens",
                                                   param.type = "bar"),
                     method = method, control = control, hessian = T); mle.bar.exp
#Same estimates
mle1.exp$par; mle2.exp$par; true.par.exp
#Bar parameterization
mle.bar.exp$par
#Standard errors
se1.exp <- sqrt(diag(solve(-mle1.exp$hessian))); se1.exp
se2.exp <- sqrt(diag(solve(-mle2.exp$hessian))); se2.exp
se.bar.exp <- sqrt(diag(solve(-mle.bar.exp$hessian))); se.bar.exp

#Contour plots
par(mfrow = c(3,1))

#gamma, mu
gam.exp <- seq(mle1.exp$par[1] - 3*se1.exp[1], mle1.exp$par[1] + 4*se1.exp[1],
               length = 50)
mu.exp <- seq(mle1.exp$par[2] - 5*se1.exp[2], mle1.exp$par[2] + 4*se1.exp[2],
              length = 50)
gam.bar.exp <- seq(mle.bar.exp$par[1] - 3*se.bar.exp[1], mle.bar.exp$par[1] + 5*se.bar.exp[1],
                   length = 50)
mu.bar.exp <- seq(mle.bar.exp$par[2] - 5*se.bar.exp[2], mle.bar.exp$par[2] + 2*se.bar.exp[2],
                  length = 50)
gammu.exp <- expand.grid(gam.exp, mu.exp)
gammu.bar.exp <- expand.grid(gam.bar.exp, mu.bar.exp)
gammu.lik.exp1 <- matrix(apply(gammu.exp,
                               1, function(x) gev.lik.gammu(par = x, data = data.exp,
                                                            t = th.exp,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.exp2 <- matrix(apply(gammu.exp,
                               1, function(x) gev.lik.gammu(par = x, data = data.exp,
                                                            t = th.exp,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.bar.exp <- matrix(apply(gammu.bar.exp,
                                  1, function(x) gev.lik.gammu(par = x, data = data.exp,
                                                               t = th.exp,
                                                               p = t.prob, log = TRUE,
                                                               llik.type = "Gev-Cens",
                                                               param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(gam.exp, mu.exp, gammu.lik.exp1 - max(gammu.lik.exp1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.exp$par[1], mle1.exp$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.exp, mu.exp, gammu.lik.exp2 - max(gammu.lik.exp2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.exp$par[1], mle1.exp$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.exp, mu.bar.exp, gammu.lik.bar.exp - max(gammu.lik.bar.exp),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(mu)), 
        main = "Classical GEV censored relative log-likelihood (bar parameterization)")
points(mle.bar.exp$par[1], mle.bar.exp$par[2], col = "blue", cex = 1.5, lwd = 2)

#gamma, delta
gam.exp <- seq(mle1.exp$par[1] - 3*se1.exp[1], mle1.exp$par[1] + 5*se1.exp[1],
               length = 50)
del.exp <- seq(mle1.exp$par[3] - 2.5*se1.exp[3], mle1.exp$par[3] + 5*se1.exp[3],
               length = 50)
gam.bar.exp <- seq(mle.bar.exp$par[1] - 3*se.bar.exp[1], mle.bar.exp$par[1] + 5*se.bar.exp[1],
                   length = 50)
del.bar.exp <- seq(mle.bar.exp$par[3] - 2.5*se.bar.exp[3], mle.bar.exp$par[3] + 8*se.bar.exp[3],
                   length = 50)
gamdel.exp <- expand.grid(gam.exp, del.exp)
gamdel.bar.exp <- expand.grid(gam.bar.exp, del.bar.exp)
gamdel.lik.exp1 <- matrix(apply(gamdel.exp,
                                1, function(x) gev.lik.gamdel(par = x, data = data.exp,
                                                              t = th.exp,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.exp2 <- matrix(apply(gamdel.exp,
                                1, function(x) gev.lik.gamdel(par = x, data = data.exp,
                                                              t = th.exp,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Max-Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.bar.exp <- matrix(apply(gamdel.bar.exp,
                                   1, function(x) gev.lik.gamdel(par = x, data = data.exp,
                                                                 t = th.exp,
                                                                 p = t.prob, log = TRUE,
                                                                 llik.type = "Gev-Cens",
                                                                 param.type = "bar")),
                             nrow = 50, ncol = 50)
contour(gam.exp, del.exp, gamdel.lik.exp1 - max(gamdel.lik.exp1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.exp$par[1], mle1.exp$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.exp, del.exp, gamdel.lik.exp2 - max(gamdel.lik.exp2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.exp$par[1], mle1.exp$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.exp, del.bar.exp, gamdel.lik.bar.exp - max(gamdel.lik.bar.exp),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.exp$par[1], mle.bar.exp$par[3], col = "blue", cex = 1.5, lwd = 2)

#mu, delta
mu.exp <- seq(mle1.exp$par[2] - 2*se1.exp[2], mle1.exp$par[2] + 2*se1.exp[2],
              length = 50)
del.exp <- seq(mle1.exp$par[3] - 3*se1.exp[3], mle1.exp$par[3] + 5*se1.exp[3],
               length = 50)
mu.bar.exp <- seq(mle.bar.exp$par[2] - se.bar.exp[2], mle.bar.exp$par[2] + 2*se.bar.exp[2],
                  length = 50)
del.bar.exp <- seq(mle.bar.exp$par[3] - 2.5*se.bar.exp[3], mle.bar.exp$par[3] + se.bar.exp[3],
                   length = 50)
mudel.exp <- expand.grid(mu.exp, del.exp)
mudel.bar.exp <- expand.grid(mu.bar.exp, del.bar.exp)
mudel.lik.exp1 <- matrix(apply(mudel.exp,
                               1, function(x) gev.lik.mudel(par = x, data = data.exp,
                                                            t = th.exp,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
mudel.lik.exp2 <- matrix(apply(mudel.exp,
                               1, function(x) gev.lik.mudel(par = x, data = data.exp,
                                                            t = th.exp,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
mudel.lik.bar.exp <- matrix(apply(mudel.bar.exp,
                                  1, function(x) gev.lik.mudel(par = x, data = data.exp,
                                                               t = th.exp,
                                                               p = t.prob, log = TRUE,
                                                               llik.type = "Gev-Cens",
                                                               param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(mu.exp, del.exp, mudel.lik.exp1 - max(mudel.lik.exp1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.exp$par[2], mle1.exp$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.exp, del.exp, mudel.lik.exp2 - max(mudel.lik.exp2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.exp$par[2], mle1.exp$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.bar.exp, del.bar.exp, mudel.lik.bar.exp - max(mudel.lik.bar.exp),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(mu)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.exp$par[2], mle.bar.exp$par[3], col = "blue", cex = 1.5, lwd = 2)

#POSTERIOR CHAINS
set.seed(12)
fit.bayes.exp1 <- fit.gev.inference(data = data.exp, t = th.exp, t.prob = t.prob,
                                    llik.type = "Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.exp,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars1 <- fit.bayes.exp1$parameters[-burn,]
fit.bayes.exp2 <- fit.gev.inference(data = data.exp, t = th.exp, t.prob = t.prob,
                                    llik.type = "Max-Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.exp,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars2 <- fit.bayes.exp2$parameters[-burn,]
fit.bayes.bar.exp <- fit.gev.inference(data = data.exp, t = th.exp, t.prob = t.prob,
                                       llik.type = "Gev-Cens", param.type = "bar",
                                       T.ret = T.ret, p = p, par0 = start.exp,
                                       hessian = TRUE, inf.type = "Bayes",
                                       k = 1, R = R, burn  = burn, prior = "empirical",
                                       val.show = TRUE)
pars.bar <- fit.bayes.bar.exp$parameters[-burn,]
#smoothScatter plots
par(mfrow = c(3,1))

#gamma, mu
smoothScatter(pars1[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,1:2])
smoothScatter(pars2[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,1:2])
smoothScatter(pars.bar[,1:2], xlab = expression(bar(gamma)),
              ylab = expression(bar(mu)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,1:2])

#gamma, delta
smoothScatter(pars1[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(1,3)])
smoothScatter(pars2[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(1,3)])
smoothScatter(pars.bar[,c(1,3)], xlab = expression(bar(gamma)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(1,3)])

#mu, delta
smoothScatter(pars1[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(2,3)])
smoothScatter(pars2[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(2,3)])
smoothScatter(pars.bar[,c(2,3)], xlab = expression(bar(mu)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(2,3)])

par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#Gamma(2,2) (tail index = 0)

#Gamma parameters
shape.gam <- 2 
rate.gam <- 2
#Set norming constants
#Centering
bn.gam <- qgamma(1 - (k/n), shape = shape.gam, rate = rate.gam)
#Scaling
#h <- 1e-10
an.gam <- k/(n*dgamma(bn.gam, shape = shape.gam, rate = rate.gam))
an2.gam <- bn.gam - (k/n)*integrate(function(x) qgamma(1 - (1/x), shape = shape.gam,
                                                       rate = rate.gam), lower = 1,
                                    upper = (n/k))$value

#Define the true parameters
true.par.gam <- c(0, bn.gam, an.gam)
#Define true extreme quantile
Q.ext.true.gam <- qgamma(1 - p, shape.gam, rate.gam)
#Define true return level
R.lev.true.gam <- qgev(1 - 1/T.ret, bn.gam, an.gam, 0)

#Define the estimation setting
start.gam <- c(0.1, 0.1, 1)
set.seed(16)
data.gam <- rgamma(n, shape.gam, rate.gam)
th.gam <- quantile(data.gam, probs = 1 - t.prob, type = 3)

#Likelihood

#Maximum likelihood
mle1.gam <- optim(start.gam, function(x) gev.lik(x, data.gam, t = th.gam,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle1.gam
mle2.gam <- optim(start.gam, function(x) gev.lik(x, data.gam, t = th.gam,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Max-Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle2.gam
mle.bar.gam <- optim(start.gam, function(x) gev.lik(x, data.gam, t = th.gam,
                                                   p = t.prob, log = TRUE,
                                                   llik.type = "Gev-Cens",
                                                   param.type = "bar"),
                     method = method, control = control, hessian = T); mle.bar.gam
#Same estimates
mle1.gam$par; mle2.gam$par; true.par.gam
#Bar parameterization
mle.bar.gam$par
#Standard errors
se1.gam <- sqrt(diag(solve(-mle1.gam$hessian))); se1.gam
se2.gam <- sqrt(diag(solve(-mle2.gam$hessian))); se2.gam
se.bar.gam <- sqrt(diag(solve(-mle.bar.gam$hessian))); se.bar.gam

#Contour plots
par(mfrow = c(3,1))

#gamma, mu
gam.gam <- seq(mle1.gam$par[1] - 3*se1.gam[1], mle1.gam$par[1] + 4*se1.gam[1],
               length = 50)
mu.gam <- seq(mle1.gam$par[2] - 5*se1.gam[2], mle1.gam$par[2] + 4*se1.gam[2],
              length = 50)
gam.bar.gam <- seq(mle.bar.gam$par[1] - 3*se.bar.gam[1], mle.bar.gam$par[1] + 5*se.bar.gam[1],
                   length = 50)
mu.bar.gam <- seq(mle.bar.gam$par[2] - 5*se.bar.gam[2], mle.bar.gam$par[2] + 2*se.bar.gam[2],
                  length = 50)
gammu.gam <- expand.grid(gam.gam, mu.gam)
gammu.bar.gam <- expand.grid(gam.bar.gam, mu.bar.gam)
gammu.lik.gam1 <- matrix(apply(gammu.gam,
                               1, function(x) gev.lik.gammu(par = x, data = data.gam,
                                                            t = th.gam,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.gam2 <- matrix(apply(gammu.gam,
                               1, function(x) gev.lik.gammu(par = x, data = data.gam,
                                                            t = th.gam,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.bar.gam <- matrix(apply(gammu.bar.gam,
                                  1, function(x) gev.lik.gammu(par = x, data = data.gam,
                                                               t = th.gam,
                                                               p = t.prob, log = TRUE,
                                                               llik.type = "Gev-Cens",
                                                               param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(gam.gam, mu.gam, gammu.lik.gam1 - max(gammu.lik.gam1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.gam$par[1], mle1.gam$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.gam, mu.gam, gammu.lik.gam2 - max(gammu.lik.gam2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.gam$par[1], mle1.gam$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.gam, mu.bar.gam, gammu.lik.bar.gam - max(gammu.lik.bar.gam),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(mu)), 
        main = "Classical GEV censored relative log-likelihood (bar parameterization)")
points(mle.bar.gam$par[1], mle.bar.gam$par[2], col = "blue", cex = 1.5, lwd = 2)

#gamma, delta
gam.gam <- seq(mle1.gam$par[1] - 3*se1.gam[1], mle1.gam$par[1] + 5*se1.gam[1],
               length = 50)
del.gam <- seq(mle1.gam$par[3] - 2.5*se1.gam[3], mle1.gam$par[3] + 5*se1.gam[3],
               length = 50)
gam.bar.gam <- seq(mle.bar.gam$par[1] - 3*se.bar.gam[1], mle.bar.gam$par[1] + 5*se.bar.gam[1],
                   length = 50)
del.bar.gam <- seq(mle.bar.gam$par[3] - 2.5*se.bar.gam[3], mle.bar.gam$par[3] + 8*se.bar.gam[3],
                   length = 50)
gamdel.gam <- expand.grid(gam.gam, del.gam)
gamdel.bar.gam <- expand.grid(gam.bar.gam, del.bar.gam)
gamdel.lik.gam1 <- matrix(apply(gamdel.gam,
                                1, function(x) gev.lik.gamdel(par = x, data = data.gam,
                                                              t = th.gam,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.gam2 <- matrix(apply(gamdel.gam,
                                1, function(x) gev.lik.gamdel(par = x, data = data.gam,
                                                              t = th.gam,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Max-Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.bar.gam <- matrix(apply(gamdel.bar.gam,
                                   1, function(x) gev.lik.gamdel(par = x, data = data.gam,
                                                                 t = th.gam,
                                                                 p = t.prob, log = TRUE,
                                                                 llik.type = "Gev-Cens",
                                                                 param.type = "bar")),
                             nrow = 50, ncol = 50)
contour(gam.gam, del.gam, gamdel.lik.gam1 - max(gamdel.lik.gam1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.gam$par[1], mle1.gam$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.gam, del.gam, gamdel.lik.gam2 - max(gamdel.lik.gam2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.gam$par[1], mle1.gam$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.gam, del.bar.gam, gamdel.lik.bar.gam - max(gamdel.lik.bar.gam),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.gam$par[1], mle.bar.gam$par[3], col = "blue", cex = 1.5, lwd = 2)

#mu, delta
mu.gam <- seq(mle1.gam$par[2] - 2*se1.gam[2], mle1.gam$par[2] + 2*se1.gam[2],
              length = 50)
del.gam <- seq(mle1.gam$par[3] - 3*se1.gam[3], mle1.gam$par[3] + 5*se1.gam[3],
               length = 50)
mu.bar.gam <- seq(mle.bar.gam$par[2] - se.bar.gam[2], mle.bar.gam$par[2] + 2*se.bar.gam[2],
                  length = 50)
del.bar.gam <- seq(mle.bar.gam$par[3] - 2.5*se.bar.gam[3], mle.bar.gam$par[3] + se.bar.gam[3],
                   length = 50)
mudel.gam <- expand.grid(mu.gam, del.gam)
mudel.bar.gam <- expand.grid(mu.bar.gam, del.bar.gam)
mudel.lik.gam1 <- matrix(apply(mudel.gam,
                               1, function(x) gev.lik.mudel(par = x, data = data.gam,
                                                            t = th.gam,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
mudel.lik.gam2 <- matrix(apply(mudel.gam,
                               1, function(x) gev.lik.mudel(par = x, data = data.gam,
                                                            t = th.gam,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
mudel.lik.bar.gam <- matrix(apply(mudel.bar.gam,
                                  1, function(x) gev.lik.mudel(par = x, data = data.gam,
                                                               t = th.gam,
                                                               p = t.prob, log = TRUE,
                                                               llik.type = "Gev-Cens",
                                                               param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(mu.gam, del.gam, mudel.lik.gam1 - max(mudel.lik.gam1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.gam$par[2], mle1.gam$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.gam, del.gam, mudel.lik.gam2 - max(mudel.lik.gam2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.gam$par[2], mle1.gam$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.bar.gam, del.bar.gam, mudel.lik.bar.gam - max(mudel.lik.bar.gam),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(mu)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.gam$par[2], mle.bar.gam$par[3], col = "blue", cex = 1.5, lwd = 2)

#POSTERIOR CHAINS
set.seed(11)
fit.bayes.gam1 <- fit.gev.inference(data = data.gam, t = th.gam, t.prob = t.prob,
                                    llik.type = "Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.gam,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars1 <- fit.bayes.gam1$parameters[-burn,]
fit.bayes.gam2 <- fit.gev.inference(data = data.gam, t = th.gam, t.prob = t.prob,
                                    llik.type = "Max-Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.gam,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars2 <- fit.bayes.gam2$parameters[-burn,]
fit.bayes.bar.gam <- fit.gev.inference(data = data.gam, t = th.gam, t.prob = t.prob,
                                       llik.type = "Gev-Cens", param.type = "bar",
                                       T.ret = T.ret, p = p, par0 = start.gam,
                                       hessian = TRUE, inf.type = "Bayes",
                                       k = 1, R = R, burn  = burn, prior = "empirical",
                                       val.show = TRUE)
pars.bar <- fit.bayes.bar.gam$parameters[-burn,]
#smoothScatter plots
par(mfrow = c(3,1))

#gamma, mu
smoothScatter(pars1[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,1:2])
smoothScatter(pars2[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,1:2])
smoothScatter(pars.bar[,1:2], xlab = expression(bar(gamma)),
              ylab = expression(bar(mu)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,1:2])

#gamma, delta
smoothScatter(pars1[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(1,3)])
smoothScatter(pars2[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(1,3)])
smoothScatter(pars.bar[,c(1,3)], xlab = expression(bar(gamma)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(1,3)])

#mu, delta
smoothScatter(pars1[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(2,3)])
smoothScatter(pars2[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(2,3)])
smoothScatter(pars.bar[,c(2,3)], xlab = expression(bar(mu)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(2,3)])

par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#Power Law distribution (tail index = -1/3)
#F(x) = 1 - K(x* - x)^(alpha);
#x* = 5: end-point of the distribution;
#alpha = 3: shape parameter;
#K = 1/9 positive constant.
K <- 1/9

#Power Law cdf
pl.cdf <- function(x, xstar, alpha, K){
  #K <- 1/9
  cdf <- 1 - K*(xstar - x)^alpha
  return(cdf)
}

#Power law pdf
pl.pdf <- function(x, xstar, alpha, K){
  #K <- 1/9
  pdf <- K*alpha*(xstar - x)^(alpha - 1)
  return(pdf)
}

#Power law quantile function
pl.qdf <- function(p, xstar, alpha, K){
  #K <- 1/9
  qdf <- xstar - ((1 - p)/K)^(1/alpha)
  return(qdf)
}

#Power Law parameters
xstar.pl <- 5 
alpha.pl <- 3  
#Set norming constants
#Centering
bn.pl <- xstar.pl - (k/(n*K))^(1/alpha.pl)
#Scaling
an.pl <- (1/alpha.pl)*((k/(n*K))^(1/alpha.pl))

#Define the true parameters
true.par.pl <- c(-1/alpha.pl, bn.pl, an.pl)
#Define true extreme quantile
Q.ext.true.pl <- pl.qdf(1 - p, xstar.pl, alpha.pl, K = K)
#Define true return level
R.lev.true.pl <- qgev(1 - 1/T.ret, bn.pl, an.pl, -1/alpha.pl)

#Define the estimation setting
start.pl <- c(0.1, 2, 1)
set.seed(15)
data.pl <- pl.qdf(runif(n), xstar.pl, alpha.pl, K = K)
th.pl <- quantile(data.pl, probs = 1 - t.prob, type = 3)

#Likelihood

#Maximum likelihood
mle1.pl <- optim(start.pl, function(x) gev.lik(x, data.pl, t = th.pl,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle1.pl
mle2.pl <- optim(start.pl, function(x) gev.lik(x, data.pl, t = th.pl,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Max-Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle2.pl
mle.bar.pl <- optim(start.pl, function(x) gev.lik(x, data.pl, t = th.pl,
                                                   p = t.prob, log = TRUE,
                                                   llik.type = "Gev-Cens",
                                                   param.type = "bar"),
                     method = method, control = control, hessian = T); mle.bar.pl
#Same estimates
mle1.pl$par; mle2.pl$par; true.par.pl
#Bar parameterization
mle.bar.pl$par
#Standard errors
#se1.pl <- sqrt(diag(solve(-mle1.pl$hessian))); se1.pl
#se2.pl <- sqrt(diag(solve(-mle2.pl$hessian))); se2.pl
se.bar.pl <- sqrt(diag(solve(-mle.bar.pl$hessian))); se.bar.pl
se1.pl <- se2.pl <- se.bar.pl

#Contour plots
par(mfrow = c(3,1))

#gamma, mu
gam.pl <- seq(mle1.pl$par[1] - 0.05*se1.pl[1], mle1.pl$par[1] + 0.05*se1.pl[1],
               length = 50)
mu.pl <- seq(mle1.pl$par[2] - 0.05*se1.pl[2], mle1.pl$par[2] + 0.05*se1.pl[2],
              length = 50)
gam.bar.pl <- seq(mle.bar.pl$par[1] - 3*se.bar.pl[1], mle.bar.pl$par[1] + 4.5*se.bar.pl[1],
                   length = 50)
mu.bar.pl <- seq(mle.bar.pl$par[2] - 4*se.bar.pl[2], mle.bar.pl$par[2] + 2*se.bar.pl[2],
                  length = 50)
gammu.pl <- expand.grid(gam.pl, mu.pl)
gammu.bar.pl <- expand.grid(gam.bar.pl, mu.bar.pl)
gammu.lik.pl1 <- matrix(apply(gammu.pl,
                               1, function(x) gev.lik.gammu(par = x, data = data.pl,
                                                            t = th.pl,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.pl2 <- matrix(apply(gammu.pl,
                               1, function(x) gev.lik.gammu(par = x, data = data.pl,
                                                            t = th.pl,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.bar.pl <- matrix(apply(gammu.bar.pl,
                                  1, function(x) gev.lik.gammu(par = x, data = data.pl,
                                                               t = th.pl,
                                                               p = t.prob, log = TRUE,
                                                               llik.type = "Gev-Cens",
                                                               param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(gam.pl, mu.pl, gammu.lik.pl1 - max(gammu.lik.pl1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.pl$par[1], mle1.pl$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.pl, mu.pl, gammu.lik.pl2 - max(gammu.lik.pl2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.pl$par[1], mle1.pl$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.pl, mu.bar.pl, gammu.lik.bar.pl - max(gammu.lik.bar.pl),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(mu)), 
        main = "Classical GEV censored relative log-likelihood (bar parameterization)")
points(mle.bar.pl$par[1], mle.bar.pl$par[2], col = "blue", cex = 1.5, lwd = 2)

#gamma, delta
gam.pl <- seq(mle1.pl$par[1] - 0.05*se1.pl[1], mle1.pl$par[1] + 0.05*se1.pl[1],
               length = 50)
del.pl <- seq(mle1.pl$par[3] - 0.01*se1.pl[3], mle1.pl$par[3] + 0.01*se1.pl[3],
               length = 50)
gam.bar.pl <- seq(mle.bar.pl$par[1] - 3*se.bar.pl[1], mle.bar.pl$par[1] + 6*se.bar.pl[1],
                   length = 50)
del.bar.pl <- seq(mle.bar.pl$par[3] - 4*se.bar.pl[3], mle.bar.pl$par[3] + 8*se.bar.pl[3],
                   length = 50)
gamdel.pl <- expand.grid(gam.pl, del.pl)
gamdel.bar.pl <- expand.grid(gam.bar.pl, del.bar.pl)
gamdel.lik.pl1 <- matrix(apply(gamdel.pl,
                                1, function(x) gev.lik.gamdel(par = x, data = data.pl,
                                                              t = th.pl,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.pl2 <- matrix(apply(gamdel.pl,
                                1, function(x) gev.lik.gamdel(par = x, data = data.pl,
                                                              t = th.pl,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Max-Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.bar.pl <- matrix(apply(gamdel.bar.pl,
                                   1, function(x) gev.lik.gamdel(par = x, data = data.pl,
                                                                 t = th.pl,
                                                                 p = t.prob, log = TRUE,
                                                                 llik.type = "Gev-Cens",
                                                                 param.type = "bar")),
                             nrow = 50, ncol = 50)
contour(gam.pl, del.pl, gamdel.lik.pl1 - max(gamdel.lik.pl1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.pl$par[1], mle1.pl$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.pl, del.pl, gamdel.lik.pl2 - max(gamdel.lik.pl2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.pl$par[1], mle1.pl$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.pl, del.bar.pl, gamdel.lik.bar.pl - max(gamdel.lik.bar.pl),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.pl$par[1], mle.bar.pl$par[3], col = "blue", cex = 1.5, lwd = 2)

#mu, delta
# mu.pl <- seq(mle1.pl$par[2] - 1e-15*se1.pl[2], mle1.pl$par[2] + 1e-15*se1.pl[2],
#               length = 50)
# del.pl <- seq(mle1.pl$par[3] - 1e-15*se1.pl[3], mle1.pl$par[3] + 1e-15*se1.pl[3],
#                length = 50)
# mu.bar.pl <- seq(mle.bar.pl$par[2] - 4*se.bar.pl[2], mle.bar.pl$par[2] + 2*se.bar.pl[2],
#                   length = 50)
# del.bar.pl <- seq(mle.bar.pl$par[3] - 4*se.bar.pl[3], mle.bar.pl$par[3] + 8*se.bar.pl[3],
#                    length = 50)
# mudel.pl <- expand.grid(mu.pl, del.pl)
# mudel.bar.pl <- expand.grid(mu.bar.pl, del.bar.pl)
# mudel.lik.pl1 <- matrix(apply(mudel.pl,
#                                1, function(x) gev.lik.mudel(par = x, data = data.pl,
#                                                             t = th.pl,
#                                                             p = t.prob, log = TRUE,
#                                                             llik.type = "Gev-Cens",
#                                                             param.type = "tilde")),
#                          nrow = 50, ncol = 50)
# mudel.lik.pl2 <- matrix(apply(mudel.pl,
#                                1, function(x) gev.lik.mudel(par = x, data = data.pl,
#                                                             t = th.pl,
#                                                             p = t.prob, log = TRUE,
#                                                             llik.type = "Max-Gev-Cens",
#                                                             param.type = "tilde")),
#                          nrow = 50, ncol = 50)
# mudel.lik.bar.pl <- matrix(apply(mudel.bar.pl,
#                                   1, function(x) gev.lik.mudel(par = x, data = data.pl,
#                                                                t = th.pl,
#                                                                p = t.prob, log = TRUE,
#                                                                llik.type = "Gev-Cens",
#                                                                param.type = "bar")),
#                             nrow = 50, ncol = 50)
# contour(mu.pl, del.pl, mudel.lik.pl1 - max(mudel.lik.pl1),
#         levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
#         xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
#         main = "Classical GEV censored relative log-likelihood")
# points(mle1.pl$par[2], mle1.pl$par[3], col = "blue", cex = 1.5, lwd = 2)
# contour(mu.pl, del.pl, mudel.lik.pl2 - max(mudel.lik.pl2),
#         levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
#         xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
#         main = "GEV censored relative log-likelihood of the maximum")
# points(mle1.pl$par[2], mle1.pl$par[3], col = "blue", cex = 1.5, lwd = 2)
# contour(mu.bar.pl, del.bar.pl, mudel.lik.bar.pl - max(mudel.lik.bar.pl),
#         levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
#         xlab = expression(bar(mu)), ylab = expression(bar(delta)), 
#         main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
# points(mle.bar.pl$par[2], mle.bar.pl$par[3], col = "blue", cex = 1.5, lwd = 2)

#POSTERIOR CHAINS
set.seed(15)
fit.bayes.pl1 <- fit.gev.inference(data = data.pl, t = th.pl, t.prob = t.prob,
                                    llik.type = "Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.pl,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars1 <- fit.bayes.pl1$parameters[-burn,]
fit.bayes.pl2 <- fit.gev.inference(data = data.pl, t = th.pl, t.prob = t.prob,
                                    llik.type = "Max-Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.pl,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars2 <- fit.bayes.pl2$parameters[-burn,]
fit.bayes.bar.pl <- fit.gev.inference(data = data.pl, t = th.pl, t.prob = t.prob,
                                       llik.type = "Gev-Cens", param.type = "bar",
                                       T.ret = T.ret, p = p, par0 = start.pl,
                                       hessian = TRUE, inf.type = "Bayes",
                                       k = 1, R = R, burn  = burn, prior = "empirical",
                                       val.show = TRUE)
pars.bar <- fit.bayes.bar.pl$parameters[-burn,]
#smoothScatter plots
par(mfrow = c(3,1))

#gamma, mu
smoothScatter(pars1[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,1:2])
smoothScatter(pars2[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,1:2])
smoothScatter(pars.bar[,1:2], xlab = expression(bar(gamma)),
              ylab = expression(bar(mu)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,1:2])

#gamma, delta
smoothScatter(pars1[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(1,3)])
smoothScatter(pars2[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(1,3)])
smoothScatter(pars.bar[,c(1,3)], xlab = expression(bar(gamma)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(1,3)])

#mu, delta
smoothScatter(pars1[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(2,3)])
smoothScatter(pars2[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(2,3)])
smoothScatter(pars.bar[,c(2,3)], xlab = expression(bar(mu)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(2,3)])

par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#Reverse Weibull (shape = 3) --> tail index = -1/3

#Reverse Weibull parameters
shape.rw <- 3 
#Set norming constants
#Centering
bn.rw <- -(-log(1 - k/n))^(1/shape.rw)
#Scaling
an.rw <- (k/(shape.rw*(n - k)))*(-log(1 - k/n))^(1/shape.rw - 1)
#an.rw2 <- 

#Define the true parameters
true.par.rw <- c(-1/shape.rw, bn.rw, an.rw)
#Define true extreme quantile
Q.ext.true.rw <- qrweibull(1 - p, shape = shape.rw)
#Define true return level
R.lev.true.rw <- qgev(1 - 1/T.ret, bn.rw, an.rw, -1/shape.rw)

#Define the estimation setting
start.rw <- c(0.1, 0.1, 1)
set.seed(13)
data.rw <- rrweibull(n, shape = shape.rw)
th.rw <- quantile(data.rw, probs = 1 - t.prob, type = 3)

#Likelihood

#Maximum likelihood
mle1.rw <- optim(start.rw, function(x) gev.lik(x, data.rw, t = th.rw,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle1.rw
mle2.rw <- optim(start.rw, function(x) gev.lik(x, data.rw, t = th.rw,
                                                p = t.prob, log = TRUE,
                                                llik.type = "Max-Gev-Cens",
                                                param.type = "tilde"),
                  method = method, control = control, hessian = T); mle2.rw
mle.bar.rw <- optim(start.rw, function(x) gev.lik(x, data.rw, t = th.rw,
                                                   p = t.prob, log = TRUE,
                                                   llik.type = "Gev-Cens",
                                                   param.type = "bar"),
                     method = method, control = control, hessian = T); mle.bar.rw
#Same estimates
mle1.rw$par; mle2.rw$par; true.par.rw
#Bar parameterization
mle.bar.rw$par
#Standard errors
se1.rw <- sqrt(diag(solve(-mle1.rw$hessian))); se1.rw
se2.rw <- sqrt(diag(solve(-mle2.rw$hessian))); se2.rw
se.bar.rw <- sqrt(diag(solve(-mle.bar.rw$hessian))); se.bar.rw

#Contour plots
par(mfrow = c(3,1))

#gamma, mu
gam.rw <- seq(mle1.rw$par[1] - 5*se1.rw[1], mle1.rw$par[1] + 4.5*se1.rw[1],
               length = 50)
mu.rw <- seq(mle1.rw$par[2] - 5*se1.rw[2], mle1.rw$par[2] + 5*se1.rw[2],
              length = 50)
gam.bar.rw <- seq(mle.bar.rw$par[1] - 3*se.bar.rw[1], mle.bar.rw$par[1] + 4*se.bar.rw[1],
                   length = 50)
mu.bar.rw <- seq(mle.bar.rw$par[2] - 4*se.bar.rw[2], mle.bar.rw$par[2] + 2*se.bar.rw[2],
                  length = 50)
gammu.rw <- expand.grid(gam.rw, mu.rw)
gammu.bar.rw <- expand.grid(gam.bar.rw, mu.bar.rw)
gammu.lik.rw1 <- matrix(apply(gammu.rw,
                               1, function(x) gev.lik.gammu(par = x, data = data.rw,
                                                            t = th.rw,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.rw2 <- matrix(apply(gammu.rw,
                               1, function(x) gev.lik.gammu(par = x, data = data.rw,
                                                            t = th.rw,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
gammu.lik.bar.rw <- matrix(apply(gammu.bar.rw,
                                  1, function(x) gev.lik.gammu(par = x, data = data.rw,
                                                               t = th.rw,
                                                               p = t.prob, log = TRUE,
                                                               llik.type = "Gev-Cens",
                                                               param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(gam.rw, mu.rw, gammu.lik.rw1 - max(gammu.lik.rw1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.rw$par[1], mle1.rw$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.rw, mu.rw, gammu.lik.rw2 - max(gammu.lik.rw2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.rw$par[1], mle1.rw$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.rw, mu.bar.rw, gammu.lik.bar.rw - max(gammu.lik.bar.rw),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(mu)), 
        main = "Classical GEV censored relative log-likelihood (bar parameterization)")
points(mle.bar.rw$par[1], mle.bar.rw$par[2], col = "blue", cex = 1.5, lwd = 2)

#gamma, delta
gam.rw <- seq(mle1.rw$par[1] - 2*se1.rw[1], mle1.rw$par[1] + se1.rw[1],
               length = 50)
del.rw <- seq(mle1.rw$par[3] - 2*se1.rw[3], mle1.rw$par[3] + se1.rw[3],
               length = 50)
gam.bar.rw <- seq(mle.bar.rw$par[1] - 4*se.bar.rw[1], mle.bar.rw$par[1] + 3.5*se.bar.rw[1],
                   length = 50)
del.bar.rw <- seq(mle.bar.rw$par[3] - 1.5*se.bar.rw[3], mle.bar.rw$par[3] + 8*se.bar.rw[3],
                   length = 50)
gamdel.rw <- expand.grid(gam.rw, del.rw)
gamdel.bar.rw <- expand.grid(gam.bar.rw, del.bar.rw)
gamdel.lik.rw1 <- matrix(apply(gamdel.rw,
                                1, function(x) gev.lik.gamdel(par = x, data = data.rw,
                                                              t = th.rw,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.rw2 <- matrix(apply(gamdel.rw,
                                1, function(x) gev.lik.gamdel(par = x, data = data.rw,
                                                              t = th.rw,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Max-Gev-Cens",
                                                              param.type = "tilde")),
                          nrow = 50, ncol = 50)
gamdel.lik.bar.rw <- matrix(apply(gamdel.bar.rw,
                                   1, function(x) gev.lik.gamdel(par = x, data = data.rw,
                                                                 t = th.rw,
                                                                 p = t.prob, log = TRUE,
                                                                 llik.type = "Gev-Cens",
                                                                 param.type = "bar")),
                             nrow = 50, ncol = 50)
contour(gam.rw, del.rw, gamdel.lik.rw1 - max(gamdel.lik.rw1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.rw$par[1], mle1.rw$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.rw, del.rw, gamdel.lik.rw2 - max(gamdel.lik.rw2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.rw$par[1], mle1.rw$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.rw, del.bar.rw, gamdel.lik.bar.rw - max(gamdel.lik.bar.rw),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.rw$par[1], mle.bar.rw$par[3], col = "blue", cex = 1.5, lwd = 2)

#mu, delta
mu.rw <- seq(mle1.rw$par[2] - 4*se1.rw[2], mle1.rw$par[2] + 4*se1.rw[2],
              length = 50)
del.rw <- seq(mle1.rw$par[3] - se1.rw[3], mle1.rw$par[3] + 0.5*se1.rw[3],
               length = 50)
mu.bar.rw <- seq(mle.bar.rw$par[2] - 0.05*se.bar.rw[2], mle.bar.rw$par[2] + se.bar.rw[2],
                  length = 50)
del.bar.rw <- seq(mle.bar.rw$par[3] - 0.05*se.bar.rw[3], mle.bar.rw$par[3] + se.bar.rw[3],
                   length = 50)
mudel.rw <- expand.grid(mu.rw, del.rw)
mudel.bar.rw <- expand.grid(mu.bar.rw, del.bar.rw)
mudel.lik.rw1 <- matrix(apply(mudel.rw,
                               1, function(x) gev.lik.mudel(par = x, data = data.rw,
                                                            t = th.rw,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
mudel.lik.rw2 <- matrix(apply(mudel.rw,
                               1, function(x) gev.lik.mudel(par = x, data = data.rw,
                                                            t = th.rw,
                                                            p = t.prob, log = TRUE,
                                                            llik.type = "Max-Gev-Cens",
                                                            param.type = "tilde")),
                         nrow = 50, ncol = 50)
mudel.lik.bar.rw <- matrix(apply(mudel.bar.rw,
                                  1, function(x) gev.lik.mudel(par = x, data = data.rw,
                                                               t = th.rw,
                                                               p = t.prob, log = TRUE,
                                                               llik.type = "Gev-Cens",
                                                               param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(mu.rw, del.rw, mudel.lik.rw1 - max(mudel.lik.rw1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.rw$par[2], mle1.rw$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.rw, del.rw, mudel.lik.rw2 - max(mudel.lik.rw2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.rw$par[2], mle1.rw$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.bar.rw, del.bar.rw, mudel.lik.bar.rw - max(mudel.lik.bar.rw),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(mu)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.rw$par[2], mle.bar.rw$par[3], col = "blue", cex = 1.5, lwd = 2)

#POSTERIOR CHAINS
set.seed(13)
fit.bayes.rw1 <- fit.gev.inference(data = data.rw, t = th.rw, t.prob = t.prob,
                                    llik.type = "Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.rw,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars1 <- fit.bayes.rw1$parameters[-burn,]
fit.bayes.rw2 <- fit.gev.inference(data = data.rw, t = th.rw, t.prob = t.prob,
                                    llik.type = "Max-Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p, par0 = start.rw,
                                    hessian = TRUE, inf.type = "Bayes",
                                    k = 1, R = R, burn  = burn, prior = "empirical",
                                    val.show = TRUE)
pars2 <- fit.bayes.rw2$parameters[-burn,]
fit.bayes.bar.rw <- fit.gev.inference(data = data.rw, t = th.rw, t.prob = t.prob,
                                       llik.type = "Gev-Cens", param.type = "bar",
                                       T.ret = T.ret, p = p, par0 = start.rw,
                                       hessian = TRUE, inf.type = "Bayes",
                                       k = 1, R = R, burn  = burn, prior = "empirical",
                                       val.show = TRUE)
pars.bar <- fit.bayes.bar.rw$parameters[-burn,]
#smoothScatter plots
par(mfrow = c(3,1))

#gamma, mu
smoothScatter(pars1[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,1:2])
smoothScatter(pars2[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,1:2])
smoothScatter(pars.bar[,1:2], xlab = expression(bar(gamma)),
              ylab = expression(bar(mu)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,1:2])

#gamma, delta
smoothScatter(pars1[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(1,3)])
smoothScatter(pars2[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(1,3)])
smoothScatter(pars.bar[,c(1,3)], xlab = expression(bar(gamma)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(1,3)])

#mu, delta
smoothScatter(pars1[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(2,3)])
smoothScatter(pars2[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(2,3)])
smoothScatter(pars.bar[,c(2,3)], xlab = expression(bar(mu)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(2,3)])

par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#Beta(1,3) --> tail index = -1/3

#Beta parameters
shape1.be <- 1
shape2.be <- 3
#Set norming constants
#Centering
bn.be <- qbeta(1 -(k/n), shape1 = shape1.be, shape2 = shape2.be)
#Scaling
an.be <- k/(n*dbeta(bn.be, shape1 = shape1.be, shape2 = shape2.be))

#Define the true parameters
true.par.be <- c(-1/shape2.be, bn.be, an.be)
#Define true extreme quantile
Q.ext.true.be <- qbeta(1 - p, shape1 = shape1.be, shape2 = shape2.be)
#Define true return level
R.lev.true.be <- qgev(1 - 1/T.ret, bn.be, an.be, -1/shape2.be)

#Define the estimation setting
start.be <- c(0.1, 0.1, 1)
set.seed(15)
data.be <- rbeta(n, shape1 = shape1.be, shape2 = shape2.be)
th.be <- quantile(data.be, probs = 1 - t.prob, type = 3)

#Likelihood

#Maximum likelihood
mle1.be <- optim(start.be, function(x) gev.lik(x, data.be, t = th.be,
                                               p = t.prob, log = TRUE,
                                               llik.type = "Gev-Cens",
                                               param.type = "tilde"),
                 method = method, control = control, hessian = T); mle1.be
mle2.be <- optim(start.be, function(x) gev.lik(x, data.be, t = th.be,
                                               p = t.prob, log = TRUE,
                                               llik.type = "Max-Gev-Cens",
                                               param.type = "tilde"),
                 method = method, control = control, hessian = T); mle2.be
mle.bar.be <- optim(start.be, function(x) gev.lik(x, data.be, t = th.be,
                                                  p = t.prob, log = TRUE,
                                                  llik.type = "Gev-Cens",
                                                  param.type = "bar"),
                    method = method, control = control, hessian = T); mle.bar.be
#Same estimates
mle1.be$par; mle2.be$par; true.par.be
#Bar parameterization
mle.bar.be$par
#Standard errors
se1.be <- sqrt(diag(solve(-mle1.be$hessian))); se1.be
se2.be <- sqrt(diag(solve(-mle2.be$hessian))); se2.be
se.bar.be <- sqrt(diag(solve(-mle.bar.be$hessian))); se.bar.be

#Contour plots
par(mfrow = c(3,1))

#gamma, mu
gam.be <- seq(mle1.be$par[1] - 5*se1.be[1], mle1.be$par[1] + 4.5*se1.be[1],
              length = 50)
mu.be <- seq(mle1.be$par[2] - 5*se1.be[2], mle1.be$par[2] + 5*se1.be[2],
             length = 50)
gam.bar.be <- seq(mle.bar.be$par[1] - 3*se.bar.be[1], mle.bar.be$par[1] + 4*se.bar.be[1],
                  length = 50)
mu.bar.be <- seq(mle.bar.be$par[2] - 4*se.bar.be[2], mle.bar.be$par[2] + 2*se.bar.be[2],
                 length = 50)
gammu.be <- expand.grid(gam.be, mu.be)
gammu.bar.be <- expand.grid(gam.bar.be, mu.bar.be)
gammu.lik.be1 <- matrix(apply(gammu.be,
                              1, function(x) gev.lik.gammu(par = x, data = data.be,
                                                           t = th.be,
                                                           p = t.prob, log = TRUE,
                                                           llik.type = "Gev-Cens",
                                                           param.type = "tilde")),
                        nrow = 50, ncol = 50)
gammu.lik.be2 <- matrix(apply(gammu.be,
                              1, function(x) gev.lik.gammu(par = x, data = data.be,
                                                           t = th.be,
                                                           p = t.prob, log = TRUE,
                                                           llik.type = "Max-Gev-Cens",
                                                           param.type = "tilde")),
                        nrow = 50, ncol = 50)
gammu.lik.bar.be <- matrix(apply(gammu.bar.be,
                                 1, function(x) gev.lik.gammu(par = x, data = data.be,
                                                              t = th.be,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Gev-Cens",
                                                              param.type = "bar")),
                           nrow = 50, ncol = 50)
contour(gam.be, mu.be, gammu.lik.be1 - max(gammu.lik.be1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.be$par[1], mle1.be$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.be, mu.be, gammu.lik.be2 - max(gammu.lik.be2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.be$par[1], mle1.be$par[2], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.be, mu.bar.be, gammu.lik.bar.be - max(gammu.lik.bar.be),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(mu)), 
        main = "Classical GEV censored relative log-likelihood (bar parameterization)")
points(mle.bar.be$par[1], mle.bar.be$par[2], col = "blue", cex = 1.5, lwd = 2)

#gamma, delta
gam.be <- seq(mle1.be$par[1] - 2*se1.be[1], mle1.be$par[1] + se1.be[1],
              length = 50)
del.be <- seq(mle1.be$par[3] - 2*se1.be[3], mle1.be$par[3] + se1.be[3],
              length = 50)
gam.bar.be <- seq(mle.bar.be$par[1] - 4*se.bar.be[1], mle.bar.be$par[1] + 3.5*se.bar.be[1],
                  length = 50)
del.bar.be <- seq(mle.bar.be$par[3] - 1.5*se.bar.be[3], mle.bar.be$par[3] + 8*se.bar.be[3],
                  length = 50)
gamdel.be <- expand.grid(gam.be, del.be)
gamdel.bar.be <- expand.grid(gam.bar.be, del.bar.be)
gamdel.lik.be1 <- matrix(apply(gamdel.be,
                               1, function(x) gev.lik.gamdel(par = x, data = data.be,
                                                             t = th.be,
                                                             p = t.prob, log = TRUE,
                                                             llik.type = "Gev-Cens",
                                                             param.type = "tilde")),
                         nrow = 50, ncol = 50)
gamdel.lik.be2 <- matrix(apply(gamdel.be,
                               1, function(x) gev.lik.gamdel(par = x, data = data.be,
                                                             t = th.be,
                                                             p = t.prob, log = TRUE,
                                                             llik.type = "Max-Gev-Cens",
                                                             param.type = "tilde")),
                         nrow = 50, ncol = 50)
gamdel.lik.bar.be <- matrix(apply(gamdel.bar.be,
                                  1, function(x) gev.lik.gamdel(par = x, data = data.be,
                                                                t = th.be,
                                                                p = t.prob, log = TRUE,
                                                                llik.type = "Gev-Cens",
                                                                param.type = "bar")),
                            nrow = 50, ncol = 50)
contour(gam.be, del.be, gamdel.lik.be1 - max(gamdel.lik.be1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.be$par[1], mle1.be$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.be, del.be, gamdel.lik.be2 - max(gamdel.lik.be2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.be$par[1], mle1.be$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(gam.bar.be, del.bar.be, gamdel.lik.bar.be - max(gamdel.lik.bar.be),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.be$par[1], mle.bar.be$par[3], col = "blue", cex = 1.5, lwd = 2)

#mu, delta
mu.be <- seq(mle1.be$par[2] - 4*se1.be[2], mle1.be$par[2] + 4*se1.be[2],
             length = 50)
del.be <- seq(mle1.be$par[3] - se1.be[3], mle1.be$par[3] + 0.5*se1.be[3],
              length = 50)
mu.bar.be <- seq(mle.bar.be$par[2] - 0.05*se.bar.be[2], mle.bar.be$par[2] + se.bar.be[2],
                 length = 50)
del.bar.be <- seq(mle.bar.be$par[3] - 0.05*se.bar.be[3], mle.bar.be$par[3] + se.bar.be[3],
                  length = 50)
mudel.be <- expand.grid(mu.be, del.be)
mudel.bar.be <- expand.grid(mu.bar.be, del.bar.be)
mudel.lik.be1 <- matrix(apply(mudel.be,
                              1, function(x) gev.lik.mudel(par = x, data = data.be,
                                                           t = th.be,
                                                           p = t.prob, log = TRUE,
                                                           llik.type = "Gev-Cens",
                                                           param.type = "tilde")),
                        nrow = 50, ncol = 50)
mudel.lik.be2 <- matrix(apply(mudel.be,
                              1, function(x) gev.lik.mudel(par = x, data = data.be,
                                                           t = th.be,
                                                           p = t.prob, log = TRUE,
                                                           llik.type = "Max-Gev-Cens",
                                                           param.type = "tilde")),
                        nrow = 50, ncol = 50)
mudel.lik.bar.be <- matrix(apply(mudel.bar.be,
                                 1, function(x) gev.lik.mudel(par = x, data = data.be,
                                                              t = th.be,
                                                              p = t.prob, log = TRUE,
                                                              llik.type = "Gev-Cens",
                                                              param.type = "bar")),
                           nrow = 50, ncol = 50)
contour(mu.be, del.be, mudel.lik.be1 - max(mudel.lik.be1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "Classical GEV censored relative log-likelihood")
points(mle1.be$par[2], mle1.be$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.be, del.be, mudel.lik.be2 - max(mudel.lik.be2),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "GEV censored relative log-likelihood of the maximum")
points(mle1.be$par[2], mle1.be$par[3], col = "blue", cex = 1.5, lwd = 2)
contour(mu.bar.be, del.bar.be, mudel.lik.bar.be - max(mudel.lik.bar.be),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(mu)), ylab = expression(bar(delta)), 
        main = "Classical GEV censored relative log-likelihood (bar reparameterization)")
points(mle.bar.be$par[2], mle.bar.be$par[3], col = "blue", cex = 1.5, lwd = 2)

#POSTERIOR CHAINS
set.seed(15)
fit.bayes.be1 <- fit.gev.inference(data = data.be, t = th.be, t.prob = t.prob,
                                   llik.type = "Gev-Cens", param.type = "tilde",
                                   T.ret = T.ret, p = p, par0 = start.be,
                                   hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
pars1 <- fit.bayes.be1$parameters[-burn,]
fit.bayes.be2 <- fit.gev.inference(data = data.be, t = th.be, t.prob = t.prob,
                                   llik.type = "Max-Gev-Cens", param.type = "tilde",
                                   T.ret = T.ret, p = p, par0 = start.be,
                                   hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
pars2 <- fit.bayes.be2$parameters[-burn,]
fit.bayes.bar.be <- fit.gev.inference(data = data.be, t = th.be, t.prob = t.prob,
                                      llik.type = "Gev-Cens", param.type = "bar",
                                      T.ret = T.ret, p = p, par0 = start.be,
                                      hessian = TRUE, inf.type = "Bayes",
                                      k = 1, R = R, burn  = burn, prior = "empirical",
                                      val.show = TRUE)
pars.bar <- fit.bayes.bar.be$parameters[-burn,]
#smoothScatter plots
par(mfrow = c(3,1))

#gamma, mu
smoothScatter(pars1[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,1:2])
smoothScatter(pars2[,1:2], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(mu)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,1:2])
smoothScatter(pars.bar[,1:2], xlab = expression(bar(gamma)),
              ylab = expression(bar(mu)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,1:2])

#gamma, delta
smoothScatter(pars1[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(1,3)])
smoothScatter(pars2[,c(1,3)], xlab = expression(tilde(gamma)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(1,3)])
smoothScatter(pars.bar[,c(1,3)], xlab = expression(bar(gamma)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(1,3)])

#mu, delta
smoothScatter(pars1[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "Classical GEV censored posterior density")
cor(pars1[,c(2,3)])
smoothScatter(pars2[,c(2,3)], xlab = expression(tilde(mu)),
              ylab = expression(tilde(delta)),  
              main = "GEV censored posterior density of the maximum")
cor(pars2[,c(2,3)])
smoothScatter(pars.bar[,c(2,3)], xlab = expression(bar(mu)),
              ylab = expression(bar(delta)),  
              main = "Classical GEV censored posterior density (bar reparameterization)")
cor(pars.bar[,c(2,3)])

par(mfrow = c(1,1))

toc()

#-------------------------------------------------------------------------------

save.image(file = "lik.RData")

#-------------------------------------------------------------------------------