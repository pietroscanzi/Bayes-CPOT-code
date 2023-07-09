#Bayesian censored Peaks over threshold complete analysis

#Exponential model

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

#Exponential(1) (tail index = 0)

#Exponential parameters
rate.exp <- 1
#Set norming constants
#Centering
bn.exp <- log(n) - log(k)
#Scaling
an.exp <- 1
an2.exp <- 1

#Define the true parameters: tilde parameterization
true.par.exp.tilde <- c(0, bn.exp, an.exp); true.par.exp.tilde
#Define the true parameters: bar parameterization
delta.bar <- an.exp
mu.bar <- bn.exp - (an.exp*log(s))
true.par.exp.bar <- c(0, mu.bar, delta.bar); true.par.exp.bar
#Define true extreme quantile
Q.ext.true.exp <- qexp(1 - p, rate.exp)
#Define true return level
R.lev.true.exp <- qgev(1 - 1/T.ret, bn.exp, an.exp, 0)

#Define the estimation setting
start.exp <- c(0.1, 1, 1)
set.seed(2)
data.exp <- rexp(n, rate.exp)
th.exp <- quantile(data.exp, probs = 1 - t.prob, type = 3)
#Define threshold exceedances
excess.exp <- data.exp[data.exp > th.exp]

#Inspect the observed sample
#Probability density function
hist(data.exp, probability = T, nclass = 1000, 
     xlab = "Sample", ylab = "Histogram", main = "",
     xlim = c(0, 10))
curve(dexp(x, rate = rate.exp), add = TRUE, col = "red",
      lwd = 4, lty = 3)
abline(v = th.exp, col = "blue", lty = "dashed", lwd = 3)
curve(dgev.tilde.V(x, par = true.par.exp.tilde, s = s), add = TRUE,
      col = "orange", lwd = 4, lty = 5)
curve(dgev.V(x, par = true.par.exp.bar), add = TRUE, col = "green", lwd = 3,
      lty = 6)
legend("topright", legend = c("True density", "Censoring threshold", 
                              "Gev tilde density", "Gev bar density"),
       lty = c(3, 2, 5, 6), lwd = c(4, 3, 4, 3),
       col = c("red", "blue", "orange", "green"))

#Cumulative density function
plot(ecdf(data.exp), xlim = c(0, 10), xlab = "x", ylab = "c.d.f.",
     main = "", lwd = 3)
curve(pexp(x, rate = rate.exp), add = TRUE, col = "red",
      lwd = 3, lty = 3)
abline(v = th.exp, col = "blue", lty = "dashed", lwd = 2)
curve(pgev.tilde.V(x, par = true.par.exp.tilde, s = s), add = TRUE,
      col = "orange", lwd = 3, lty = 5)
curve(pgev.V(x, par = true.par.exp.bar), add = TRUE, col = "green", lwd = 2,
      lty = 6)
legend("bottomright", legend = c("Empirical c.d.f.", "True c.d.f.", "Censoring threshold", 
                                 "Gev tilde c.d.f.", "Gev bar c.d.f."),
       lty = c(1, 3, 2, 5, 6), lwd = c(3, 3, 2, 3, 2),
       col = c("black", "red", "blue", "orange", "green"))

#Focus on threshold exceedances
#Probability density function
hist(excess.exp, probability = T, nclass = 1000, 
     xlab = "t.exc", ylab = "Histogram", main = "Exponential model",
     xlim = c(th.exp - 1, 10))
abline(v = th.exp, col = "blue", lty = "dashed", lwd = 2)
curve(dexp(x, rate = rate.exp), add = TRUE, col = "red",
      lwd = 2, lty = 3)
curve(dgev.tilde.V(x, par = true.par.exp.tilde, s = s), add = TRUE,
      col = "orange", lwd = 4, lty = 5)
curve(dgev.V(x, par = true.par.exp.bar), add = TRUE, col = "green", lwd = 2,
      lty = 6)
legend("topright", legend = c("True density", "Censoring threshold", 
                              "Gev tilde density", "Gev bar density"),
       lty = c(3, 2, 5, 6), lwd = c(2, 2, 4, 2),
       col = c("red", "blue", "orange", "green"))

#Cumulative density function
plot(ecdf(excess.exp), xlim = c(th.exp -1, 10), xlab = "t.exc", ylab = "cdf",
     main = "Exponential model", lwd = 2)
curve(pexp(x, rate = rate.exp), add = TRUE, col = "red",
      lwd = 3, lty = 3)
abline(v = th.exp, col = "blue", lty = "dashed", lwd = 2)
curve(pgev.tilde.V(x, par = true.par.exp.tilde, s = s), add = TRUE,
      col = "orange", lwd = 3, lty = 5)
curve(pgev.V(x, par = true.par.exp.bar), add = TRUE, col = "green", lwd = 2,
      lty = 6)
legend("bottomright", legend = c("Empirical cdf", "True cdf", "Censoring threshold", 
                                 "Gev tilde cdf", "Gev bar cdf"),
       lty = c(1, 3, 2, 5, 6), lwd = c(2, 3, 2, 3, 2),
       col = c("black", "red", "blue", "orange", "green"))

#Kolmogorov-Smirnov tests
ks.test(data.exp, pexp, rate = rate.exp)
ks.test(data.exp, pgev.tilde.V, par = true.par.exp.tilde, s = s)
ks.test(data.exp, pgev.tilde.V, par = true.par.exp.tilde, s = s,
        alternative = "greater")
ks.test(data.exp, pgev.V, par = true.par.exp.bar)
ks.test(data.exp, pgev.V, par = true.par.exp.bar, alternative = "greater")

ks.test(excess.exp, pgev.tilde.V, par = true.par.exp.tilde, s = s)
ks.test(excess.exp, pgev.tilde.V, par = true.par.exp.tilde, s = s,
        alternative = "greater")
ks.test(excess.exp, pgev.V, par = true.par.exp.bar)
ks.test(excess.exp, pgev.V, par = true.par.exp.bar, alternative = "greater")

#-------------------------------------------------------------------------------

#Frequentist inference
fit.freq.exp <- fit.gev.inference(data = data.exp, t = th.exp, t.prob = t.prob,
                                  llik.type = "Max-Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p, par0 = start.exp,
                                  hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.exp
#Parameters
fit.freq.exp$mle; true.par.exp.tilde
#Standard errors
sqrt(diag(solve(-fit.freq.exp$hessian)))
#Extreme quantile
fit.freq.exp$Q.extreme; Q.ext.true.exp
#Return level
fit.freq.exp$R.level; R.lev.true.exp

#-------------------------------------------------------------------------------

#Bayesian inference
tic()
fit.bayes.exp <- fit.gev.inference(data = data.exp, t = th.exp, t.prob = t.prob,
                                   llik.type = "Max-Gev-Cens", param.type = "tilde",
                                   T.ret = T.ret, p = p, par0 = start.exp,
                                   hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
toc()
#438.39 sec elapsed
#Results
fit.bayes.exp$mcmc.acc
#Parameters
fit.bayes.exp$mle; fit.bayes.exp$mode; true.par.exp.tilde
sapply(1:3, function(x) summary(fit.bayes.exp$parameters[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(fit.bayes.exp$parameters[,i]), fit.bayes.exp$parameters[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.exp.tilde[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes.exp$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(fit.bayes.exp$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.exp.tilde[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(fit.bayes.exp$parameters[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.exp
summary(fit.bayes.exp$Q.extreme)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(fit.bayes.exp$Q.extreme), fit.bayes.exp$Q.extreme, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.exp, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.exp$Q.extreme, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.exp$Q.extreme, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.exp, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.exp$Q.extreme, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.exp
summary(fit.bayes.exp$R.level)
#Trace plot
plot(1:length(fit.bayes.exp$R.level), fit.bayes.exp$R.level, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.exp, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.exp$R.level, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.exp$R.level, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.exp, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.exp$R.level, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))
par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#Convergence diagnostics
set.seed(2)
starts.exp <- rmvnorm(n = nchains, mean = fit.bayes.exp$mode, 
                      sigma = solve(-fit.bayes.exp$hessian)); starts.exp
res.list.exp <- list()
tic()
for(i in 1:nchains){
  inference <- fit.gev.inference(data = data.exp, t = th.exp, t.prob = t.prob,
                                 llik.type = "Max-Gev-Cens", param.type = "tilde",
                                 T.ret = T.ret, p = p, par0 = starts.exp[i,],
                                 hessian = TRUE, inf.type = "Bayes",
                                 k = 1, R = R, burn = 0, prior = "empirical",
                                 val.show = FALSE)
  res.list.exp[[i]] <- cbind(inference$parameters, inference$Q.extreme,
                             inference$R.level)
  colnames(res.list.exp[[i]]) <- c("Shape", "Location", "Scale", "Q.ext", "R.lev")
  cat(i, "")
}
toc()
#6331.18 sec elapsed
save.image(file = "exponential_complete.RData")

res.mcmc.exp <- llply(res.list.exp, function(x) mcmc(window(x, start = burn + 1),
                                                     start = burn + 1))
res.mcmc.exp <- mcmc.list(res.mcmc.exp)
summary(res.mcmc.exp)
plot(res.mcmc.exp)
acfplot(res.mcmc.exp)
#Potential scale reduction factor
gelman.diag(res.mcmc.exp)
gelman.plot(res.mcmc.exp)
#Effective sample size
effectiveSize(res.mcmc.exp)

#-------------------------------------------------------------------------------

#Final sample
resf.exp <- NULL
for(i in 1:nchains){
  resf.exp <- rbind(resf.exp, res.mcmc.exp[[i]])
}
dim(resf.exp)

#Parameters
colMeans(resf.exp)[1:3]; true.par.exp.tilde
sapply(1:3, function(x) summary(resf.exp[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(resf.exp[,i]), resf.exp[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.exp.tilde[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(resf.exp[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(resf.exp[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.exp.tilde[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(resf.exp[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.exp
summary(resf.exp[,4])
par(mfrow = c(3,1))
#Trace plot
plot(1:length(resf.exp[,4]), resf.exp[,4], type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.exp, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(resf.exp[,4], lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(resf.exp[,4], nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.exp, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(resf.exp[,4], conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.exp
summary(resf.exp[,5])
#Trace plot
plot(1:length(resf.exp[,5]), resf.exp[,5], type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.exp, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(resf.exp[,5], lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(resf.exp[,5], nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.exp, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(resf.exp[,5], conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))
par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

save.image(file = "exponential_complete.RData")

#-------------------------------------------------------------------------------