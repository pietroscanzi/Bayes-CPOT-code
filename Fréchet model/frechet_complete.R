#Bayesian censored Peaks over threshold complete analysis

#Fréchet model

#-------------------------------------------------------------------------------

path <- "C:\\Users\\scanz\\OneDrive - Università Commerciale Luigi Bocconi\\Desktop\\Magistrale\\Extreme Value Theory\\Code\\"
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

#Define the true parameters: tilde parameterization
true.par.fre.tilde <- c(1/shape.fre, bn.fre, an.fre); true.par.fre.tilde
#Define the true parameters: bar parameterization
delta.bar <- s^(-(1/shape.fre))*an.fre
mu.bar <- bn.fre - an.fre*(1 - s^(-(1/shape.fre)))/(1/shape.fre)
true.par.fre.bar <- c(1/shape.fre, mu.bar, delta.bar); true.par.fre.bar
#Define true extreme quantile
Q.ext.true.fre <- qfrechet(1 - p, loc.fre, scale.fre, shape.fre)
#Define true return level
R.lev.true.fre <- qgev(1 - 1/T.ret, bn.fre, an.fre, 1/shape.fre)

#Define the estimation setting
start.fre <- c(0.1, 0.1, 1)
set.seed(6)
data.fre <- rfrechet(n, loc.fre, scale.fre, shape.fre)
th.fre <- quantile(data.fre, probs = 1 - t.prob, type = 3)
#Define threshold exceedances
excess.fre <- data.fre[data.fre > th.fre]

#Inspect the observed sample
#Probability density function
hist(data.fre, probability = T, nclass = 1000, 
     xlab = "Sample", ylab = "Histogram",
     main = "",
     xlim = c(0, 80), ylim = c(0, 0.45))
curve(dfrechet(x, loc.fre, scale.fre, shape.fre), add = TRUE, col = "red",
      lwd = 4, lty = 3)
abline(v = th.fre, col = "blue", lty = "dashed", lwd = 3)
curve(dgev.tilde.V(x, par = true.par.fre.tilde, s = s), add = TRUE,
      col = "orange", lwd = 4, lty = 5)
curve(dgev.V(x, par = true.par.fre.bar), add = TRUE, col = "green", lwd = 3,
      lty = 6)
legend("topright", legend = c("True density", "Censoring threshold", 
                              "Gev tilde density", "Gev bar density"),
       lty = c(3, 2, 5, 6), lwd = c(4, 3, 4, 3),
       col = c("red", "blue", "orange", "green"))

#Cumulative density function
plot(ecdf(data.fre), xlim = c(0, 80), xlab = "x", ylab = "c.d.f.",
     main = "",
     lwd = 3)
curve(pfrechet(x, loc.fre, scale.fre, shape.fre), add = TRUE, col = "red",
      lwd = 3, lty = 4)
abline(v = th.fre, col = "blue", lty = "dashed", lwd = 4)
curve(pgev.tilde.V(x, par = true.par.fre.tilde, s = s), add = TRUE,
      col = "orange", lwd = 3, lty = 5)
curve(pgev.V(x, par = true.par.fre.bar), add = TRUE, col = "green", lwd = 3,
      lty = 6)
legend("bottomright", legend = c("Empirical c.d.f.", "True c.d.f.", "Censoring threshold", 
                              "Gev tilde c.d.f.", "Gev bar c.d.f."),
       lty = c(1, 3, 2, 5, 6), lwd = c(3, 4, 4, 3, 2),
       col = c("black", "red", "blue", "orange", "green"))

#Focus on threshold exceedances
#Probability density function
hist(excess.fre, probability = T, nclass = 1000, 
     xlab = "t.exc", ylab = "Histogram", main = "Fréchet model",
     xlim = c(th.fre - 5, 300))
abline(v = th.fre, col = "blue", lty = "dashed", lwd = 2)
curve(dfrechet(x, loc.fre, scale.fre, shape.fre), add = TRUE, col = "red",
      lwd = 2, lty = 3)
curve(dgev.tilde.V(x, par = true.par.fre.tilde, s = s), add = TRUE,
      col = "orange", lwd = 4, lty = 5)
curve(dgev.V(x, par = true.par.fre.bar), add = TRUE, col = "green", lwd = 2,
      lty = 6)
legend("topright", legend = c("True density", "Censoring threshold", 
                              "Gev tilde density", "Gev bar density"),
       lty = c(3, 2, 5, 6), lwd = c(2, 2, 4, 2),
       col = c("red", "blue", "orange", "green"))

#Cumulative density function
plot(ecdf(excess.fre), xlim = c(th.fre -5, 300), xlab = "t.exc", ylab = "cdf",
     main = "Fréchet model", lwd = 2)
curve(pfrechet(x, loc.fre, scale.fre, shape.fre), add = TRUE, col = "red",
      lwd = 3, lty = 3)
abline(v = th.fre, col = "blue", lty = "dashed", lwd = 2)
curve(pgev.tilde.V(x, par = true.par.fre.tilde, s = s), add = TRUE,
      col = "orange", lwd = 3, lty = 5)
curve(pgev.V(x, par = true.par.fre.bar), add = TRUE, col = "green", lwd = 2,
      lty = 6)
legend("bottomright", legend = c("Empirical cdf", "True cdf", "Censoring threshold", 
                                 "Gev tilde cdf", "Gev bar cdf"),
       lty = c(1, 3, 2, 5, 6), lwd = c(2, 3, 2, 3, 2),
       col = c("black", "red", "blue", "orange", "green"))

#Kolmogorov-Smirnov tests
ks.test(data.fre, pfrechet, loc = loc.fre, scale = scale.fre, shape = shape.fre)
ks.test(data.fre, pgev.tilde.V, par = true.par.fre.tilde, s = s)
ks.test(data.fre, pgev.tilde.V, par = true.par.fre.tilde, s = s,
        alternative = "greater")
ks.test(data.fre, pgev.V, par = true.par.fre.bar)
ks.test(data.fre, pgev.V, par = true.par.fre.bar, alternative = "greater")

ks.test(excess.fre, pgev.tilde.V, par = true.par.fre.tilde, s = s)
ks.test(excess.fre, pgev.tilde.V, par = true.par.fre.tilde, s = s,
        alternative = "greater")
ks.test(excess.fre, pgev.V, par = true.par.fre.bar)
ks.test(excess.fre, pgev.V, par = true.par.fre.bar, alternative = "greater")

#-------------------------------------------------------------------------------

#Frequentist inference
fit.freq.fre <- fit.gev.inference(data = data.fre, t = th.fre, t.prob = t.prob,
                                  llik.type = "Max-Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p, par0 = start.fre,
                                  hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.fre
#Parameters
fit.freq.fre$mle; true.par.fre.tilde
#Standard errors
sqrt(diag(solve(-fit.freq.fre$hessian)))
#Extreme quantile
fit.freq.fre$Q.extreme; Q.ext.true.fre
#Return level
fit.freq.fre$R.level; R.lev.true.fre

#-------------------------------------------------------------------------------

#Bayesian inference
tic()
fit.bayes.fre <- fit.gev.inference(data = data.fre, t = th.fre, t.prob = t.prob,
                                   llik.type = "Max-Gev-Cens", param.type = "tilde",
                                   T.ret = T.ret, p = p, par0 = start.fre,
                                   hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
toc()
#342.22 sec elapsed
#Results
fit.bayes.fre$mcmc.acc
#Parameters
fit.bayes.fre$mle; fit.bayes.fre$mode; true.par.fre.tilde
sapply(1:3, function(x) summary(fit.bayes.fre$parameters[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(fit.bayes.fre$parameters[,i]), fit.bayes.fre$parameters[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.fre.tilde[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes.fre$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(fit.bayes.fre$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.fre.tilde[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(fit.bayes.fre$parameters[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.fre
summary(fit.bayes.fre$Q.extreme)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(fit.bayes.fre$Q.extreme), fit.bayes.fre$Q.extreme, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.fre, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.fre$Q.extreme, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.fre$Q.extreme, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.fre, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.fre$Q.extreme, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.fre
summary(fit.bayes.fre$R.level)
#Trace plot
plot(1:length(fit.bayes.fre$R.level), fit.bayes.fre$R.level, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.fre, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.fre$R.level, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.fre$R.level, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.fre, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.fre$R.level, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))
par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#Convergence diagnostics
set.seed(7)
starts.fre <- rmvnorm(n = nchains, mean = fit.bayes.fre$mode, 
                      sigma = solve(-fit.bayes.fre$hessian)); starts.fre
res.list.fre <- list()
tic()
for(i in 1:nchains){
  inference <- fit.gev.inference(data = data.fre, t = th.fre, t.prob = t.prob,
                                         llik.type = "Max-Gev-Cens", param.type = "tilde",
                                         T.ret = T.ret, p = p, par0 = starts.fre[i,],
                                         hessian = TRUE, inf.type = "Bayes",
                                         k = 1, R = R, burn = 0, prior = "empirical",
                                         val.show = FALSE)
  res.list.fre[[i]] <- cbind(inference$parameters, inference$Q.extreme,
                             inference$R.level)
  colnames(res.list.fre[[i]]) <- c("Shape", "Location", "Scale", "Q.ext", "R.lev")
  cat(i, "")
}
toc()
#946.49 sec elapsed
save.image(file = "frechet_complete.RData")

res.mcmc.fre <- llply(res.list.fre, function(x) mcmc(window(x, start = burn + 1),
                                                     start = burn + 1))
res.mcmc.fre <- mcmc.list(res.mcmc.fre)
summary(res.mcmc.fre)
plot(res.mcmc.fre)
acfplot(res.mcmc.fre)
#Potential Scale Reduction Factor
gelman.diag(res.mcmc.fre)
gelman.plot(res.mcmc.fre)
#Effective sample size
effectiveSize(res.mcmc.fre)

#-------------------------------------------------------------------------------

#Final sample
resf.fre <- NULL
for(i in 1:nchains){
  resf.fre <- rbind(resf.fre, res.mcmc.fre[[i]])
}
dim(resf.fre)

#Parameters
colMeans(resf.fre)[1:3]; true.par.fre.tilde
sapply(1:3, function(x) summary(resf.fre[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(resf.fre[,i]), resf.fre[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.fre.tilde[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(resf.fre[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(resf.fre[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.fre.tilde[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(resf.fre[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.fre
summary(resf.fre[,4])
par(mfrow = c(3,1))
#Trace plot
plot(1:length(resf.fre[,4]), resf.fre[,4], type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.fre, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(resf.fre[,4], lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(resf.fre[,4], nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.fre, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(resf.fre[,4], conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.fre
summary(resf.fre[,5])
#Trace plot
plot(1:length(resf.fre[,5]), resf.fre[,5], type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.fre, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(resf.fre[,5], lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(resf.fre[,5], nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.fre, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(resf.fre[,5], conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))
par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

save.image(file = "frechet_complete.RData")

#-------------------------------------------------------------------------------