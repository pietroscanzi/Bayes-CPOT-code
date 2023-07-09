#Bayesian censored Peaks over threshold complete analysis

#Power-Law model

#-------------------------------------------------------------------------------

path <- "C:\\Users\\scanz\\OneDrive - Università Commerciale Luigi Bocconi\\Desktop/Magistrale\\Extreme Value Theory\\Code\\"
setwd(path)
source(paste0(path ,"functions.R"))

#-------------------------------------------------------------------------------

#Simulation setting
library(LaplacesDemon)
library(TeachingDemos)
library(evd)
library(tictoc)
library(mvtnorm)
library(plyr)
library(coda)

#Sample size
n <- 10000
#Exceeding probability used in the estimation
t.prob <- 0.05
#Number of exceedances used in the estimation
k <- n * t.prob
#Block-size for the correspoding block maxima approach
m <- 1/t.prob
#Inverse of the proportion of exceedances: s = n/k
s <- n/k
#Exceeding probability used for prediction --> extreme quantile
p <- 1/n
#Return period --> return level
T.ret <- 50
#Posterior sample size
R <- 60000
#Burn-in period
burn <- 10000
#Number of MCMC chains for diagnostics
nchains <- 5

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

#Define the true parameters: tilde parameterization
true.par.pl.tilde <- c(-1/alpha.pl, bn.pl, an.pl); true.par.pl.tilde
#Define the true parameters: bar parameterization
delta.bar <- s^(-(-1/alpha.pl))*an.pl
mu.bar <- bn.pl - an.pl*(1 - s^(-(-1/alpha.pl)))/(-1/alpha.pl)
true.par.pl.bar <- c(-1/alpha.pl, mu.bar, delta.bar); true.par.pl.bar
#Define true extreme quantile
Q.ext.true.pl <- pl.qdf(1 - p, xstar.pl, alpha.pl, K = K)
#Define true return level
R.lev.true.pl <- qgev(1 - 1/T.ret, bn.pl, an.pl, -1/alpha.pl)

#Define the estimation setting
start.pl <- c(0.1, 2, 1)
set.seed(3)
data.pl <- pl.qdf(runif(n), xstar.pl, alpha.pl, K = K)
th.pl <- quantile(data.pl, probs = 1 - t.prob, type = 3)
#Define threshold exceedances
excess.pl <- data.pl[data.pl > th.pl]

#Inspect the observed sample
#Probability density function
hist(data.pl, probability = T, nclass = 1000, 
     xlab = "Sample", ylab = "Histogram", main = "",
     xlim = c(0, xstar.pl + 2))
curve(pl.pdf(x, xstar = xstar.pl, alpha = alpha.pl, K = K),
      from = 0, to = xstar.pl,
      add = TRUE, col = "red", lwd = 4, lty = 3)
abline(v = th.pl, col = "blue", lty = "dashed", lwd = 3)
curve(dgev.tilde.V(x, par = true.par.pl.tilde, s = s), add = TRUE,
      col = "orange", lwd = 4, lty = 5)
curve(dgev.V(x, par = true.par.pl.bar), add = TRUE, col = "green", lwd = 3,
      lty = 6)
legend("topright", legend = c("True density", "Censoring threshold", 
                              "Gev tilde density", "Gev bar density"),
       lty = c(3, 2, 5, 6), lwd = c(4, 3, 4, 3),
       col = c("red", "blue", "orange", "green"))

#Cumulative density function
plot(ecdf(data.pl), xlim = c(0, xstar.pl + 1), xlab = "x", ylab = "c.d.f.",
     main = "", lwd = 3)
curve(pl.cdf(x, xstar = xstar.pl, alpha = alpha.pl, K = K),
      from = 0, to = xstar.pl,
      add = TRUE, col = "red", lwd = 3, lty = 3)
abline(v = th.pl, col = "blue", lty = "dashed", lwd = 2)
curve(pgev.tilde.V(x, par = true.par.pl.tilde, s = s), add = TRUE,
      col = "orange", lwd = 3, lty = 5)
curve(pgev.V(x, par = true.par.pl.bar), add = TRUE, col = "green", lwd = 2,
      lty = 6)
legend("topleft", legend = c("Empirical c.d.f.", "True c.d.f.", "Censoring threshold", 
                                 "Gev tilde c.d.f.", "Gev bar c.d.f."),
       lty = c(1, 3, 2, 5, 6), lwd = c(3, 3, 2, 3, 2),
       col = c("black", "red", "blue", "orange", "green"))

#Focus on threshold exceedances
#Probability density function
hist(excess.pl, probability = T, nclass = 1000, 
     xlab = "t.exc", ylab = "Histogram", main = "Power-Law model",
     xlim = c(0, xstar.pl + 2))
abline(v = th.pl, col = "blue", lty = "dashed", lwd = 2)
curve(pl.pdf(x, xstar = xstar.pl, alpha = alpha.pl, K = K),
      from = 0, to = xstar.pl,
      add = TRUE, col = "red", lwd = 2, lty = 3)
curve(dgev.tilde.V(x, par = true.par.pl.tilde, s = s), add = TRUE,
      col = "orange", lwd = 4, lty = 5)
curve(dgev.V(x, par = true.par.pl.bar), add = TRUE, col = "green", lwd = 2,
      lty = 6)
legend("topleft", legend = c("True density", "Censoring threshold", 
                              "Gev tilde density", "Gev bar density"),
       lty = c(3, 2, 5, 6), lwd = c(2, 2, 4, 2),
       col = c("red", "blue", "orange", "green"))

#Cumulative density function
plot(ecdf(excess.pl), xlim = c(0, xstar.pl + 1), xlab = "t.exc", ylab = "cdf",
     main = "Power-Law model", lwd = 2)
curve(pl.cdf(x, xstar = xstar.pl, alpha = alpha.pl, K = K),
      from = 0, to = xstar.pl,
      add = TRUE, col = "red", lwd = 3, lty = 3)
abline(v = th.pl, col = "blue", lty = "dashed", lwd = 2)
curve(pgev.tilde.V(x, par = true.par.pl.tilde, s = s), add = TRUE,
      col = "orange", lwd = 3, lty = 5)
curve(pgev.V(x, par = true.par.pl.bar), add = TRUE, col = "green", lwd = 2,
      lty = 6)
legend("topleft", legend = c("Empirical cdf", "True cdf", "Censoring threshold", 
                                 "Gev tilde cdf", "Gev bar cdf"),
       lty = c(1, 3, 2, 5, 6), lwd = c(2, 3, 2, 3, 2),
       col = c("black", "red", "blue", "orange", "green"))

#Kolmogorov-Smirnov tests
ks.test(data.pl, pl.pdf, xstar = xstar.pl, alpha = alpha.pl, K = K)
ks.test(data.pl, pgev.tilde.V, par = true.par.pl.tilde, s = s)
ks.test(data.pl, pgev.tilde.V, par = true.par.pl.tilde, s = s,
        alternative = "greater")
ks.test(data.pl, pgev.V, par = true.par.pl.bar)
ks.test(data.pl, pgev.V, par = true.par.pl.bar, alternative = "greater")

ks.test(excess.pl, pgev.tilde.V, par = true.par.pl.tilde, s = s)
ks.test(excess.pl, pgev.tilde.V, par = true.par.pl.tilde, s = s,
        alternative = "greater")
ks.test(excess.pl, pgev.V, par = true.par.pl.bar)
ks.test(excess.pl, pgev.V, par = true.par.pl.bar, alternative = "greater")

#-------------------------------------------------------------------------------

#Frequentist inference
fit.freq.pl <- fit.gev.inference(data = data.pl, t = th.pl, t.prob = t.prob,
                                 llik.type = "Max-Gev-Cens", param.type = "tilde",
                                 T.ret = T.ret, p = p, par0 = start.pl,
                                 hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.pl
#Parameters
fit.freq.pl$mle; true.par.pl.tilde
#Standard errors
sqrt(diag(solve(-fit.freq.pl$hessian)))
#Extreme quantile
fit.freq.pl$Q.extreme; Q.ext.true.pl
#Return level
fit.freq.pl$R.level; R.lev.true.pl

#-------------------------------------------------------------------------------

#Bayesian inference
tic()
fit.bayes.pl <- fit.gev.inference(data = data.pl, t = th.pl, t.prob = t.prob,
                                  llik.type = "Max-Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p, par0 = start.pl,
                                  hessian = TRUE, inf.type = "Bayes",
                                  k = 1, R = R, burn  = burn, prior = "empirical",
                                  val.show = TRUE)
toc()
#411.56 sec elapsed
#Results
fit.bayes.pl$mcmc.acc
#Parameters
fit.bayes.pl$mle; fit.bayes.pl$mode; true.par.pl.tilde
sapply(1:3, function(x) summary(fit.bayes.pl$parameters[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(fit.bayes.pl$parameters[,i]), fit.bayes.pl$parameters[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.pl.tilde[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes.pl$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(fit.bayes.pl$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.pl.tilde[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(fit.bayes.pl$parameters[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.pl
summary(fit.bayes.pl$Q.extreme)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(fit.bayes.pl$Q.extreme), fit.bayes.pl$Q.extreme, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.pl, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.pl$Q.extreme, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.pl$Q.extreme, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.pl, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.pl$Q.extreme, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.pl
summary(fit.bayes.pl$R.level)
#Trace plot
plot(1:length(fit.bayes.pl$R.level), fit.bayes.pl$R.level, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.pl, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.pl$R.level, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.pl$R.level, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.pl, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.pl$R.level, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))
par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#Convergence diagnostics
set.seed(3)
starts.pl <- rmvnorm(n = nchains, mean = fit.bayes.pl$mode, 
                      sigma = solve(-fit.bayes.pl$hessian)); starts.pl
res.list.pl <- list()
tic()
for(i in 1:nchains){
  inference <- fit.gev.inference(data = data.pl, t = th.pl, t.prob = t.prob,
                                 llik.type = "Max-Gev-Cens", param.type = "tilde",
                                 T.ret = T.ret, p = p, par0 = starts.pl[i,],
                                 hessian = TRUE, inf.type = "Bayes",
                                 k = 1, R = R, burn = 0, prior = "empirical",
                                 val.show = FALSE)
  res.list.pl[[i]] <- cbind(inference$parameters, inference$Q.extreme,
                             inference$R.level)
  colnames(res.list.pl[[i]]) <- c("Shape", "Location", "Scale", "Q.ext", "R.lev")
  cat(i, "")
}
toc()
#5219.72 sec elapsed
save.image(file = "powerlaw_complete.RData")

res.mcmc.pl <- llply(res.list.pl, function(x) mcmc(window(x, start = burn + 1),
                                                     start = burn + 1))
res.mcmc.pl <- mcmc.list(res.mcmc.pl)
summary(res.mcmc.pl)
plot(res.mcmc.pl)
acfplot(res.mcmc.pl)
#Potential scale reduction factor
gelman.diag(res.mcmc.pl)
gelman.plot(res.mcmc.pl)
#Effective sample size
effectiveSize(res.mcmc.pl)

#-------------------------------------------------------------------------------

#Final sample
resf.pl <- NULL
for(i in 1:nchains){
  resf.pl <- rbind(resf.pl, res.mcmc.pl[[i]])
}
dim(resf.pl)

#Parameters
colMeans(resf.pl)[1:3]; true.par.pl.tilde
sapply(1:3, function(x) summary(resf.pl[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(resf.pl[,i]), resf.pl[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.pl.tilde[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(resf.pl[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(resf.pl[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.pl.tilde[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(resf.pl[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.pl
summary(resf.pl[,4])
par(mfrow = c(3,1))
#Trace plot
plot(1:length(resf.pl[,4]), resf.pl[,4], type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.pl, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(resf.pl[,4], lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(resf.pl[,4], nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.pl, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(resf.pl[,4], conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.pl
summary(resf.pl[,5])
#Trace plot
plot(1:length(resf.pl[,5]), resf.pl[,5], type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.pl, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(resf.pl[,5], lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(resf.pl[,5], nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.pl, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(resf.pl[,5], conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))
par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

save.image(file = "powerlaw_complete.RData")

#-------------------------------------------------------------------------------