#-------------------------------------------------------------------------------

#Trials of censored-likelihood threshold Max-Gev inference methods for different
#classes of models belonging to the Gev family:
#Frechét domain: unit-Fréchet, standard Pareto, half Cauchy (tail index = 1)
#Gumbel domain: standard Gumbel, Exponential(1), Gamma(2,2) (tail index = 0)
#Reverse Weibull domain: Power-law, Reverse Weibull, Beta(1,3) (tail index = -1/3)

#Power-law distribution: F(x) = 1 - K(x* - x)^(alpha);
#x* = 5: end-point of the distribution;
#alpha = 3: shape parameter;
#K = 1/9 positive constant;

#-------------------------------------------------------------------------------

path <- "C:\\Users\\scanz\\OneDrive - Università Commerciale Luigi Bocconi\\Desktop\\Magistrale\\Extreme Value Theory\\Code\\"
setwd(path)
source(paste0(path ,"functions.R"))

#-------------------------------------------------------------------------------

#Simulation setting
library(LaplacesDemon)
library(TeachingDemos)
library(evd)
#library(extremefit)
library(tictoc)

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
#Define true extreme quantile
Q.ext.true.fre <- qfrechet(1 - p, loc.fre, scale.fre, shape.fre)
#Define true return level
R.lev.true.fre <- qgev(1 - 1/T.ret, bn.fre, an.fre, 1/shape.fre)

#Define the estimation setting
start.fre <- c(0.1, 0.1, 1)
set.seed(1)
data.fre <- rfrechet(n, loc.fre, scale.fre, shape.fre)
th.fre <- quantile(data.fre, probs = 1 - t.prob, type = 3)

#Frequentist inference
fit.freq.fre <- fit.gev.inference(data = data.fre, t = th.fre, t.prob = t.prob,
                                  llik.type = "Max-Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p, par0 = start.fre,
                                  hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.fre
#Comparison
#Parameters
fit.freq.fre$mle; true.par.fre
#Extreme quantile
fit.freq.fre$Q.extreme; Q.ext.true.fre
#Return level
fit.freq.fre$R.level; R.lev.true.fre

#Bayesian inference
tic()
fit.bayes.fre <- fit.gev.inference(data = data.fre, t = th.fre, t.prob = t.prob,
                                   llik.type = "Max-Gev-Cens", param.type = "tilde",
                                   T.ret = T.ret, p = p, par0 = start.fre,
                                   hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
toc()
#58.8 sec elapsed
#Results
fit.bayes.fre$mcmc.acc
#Comparison
#Parameters
fit.bayes.fre$mle; fit.bayes.fre$mode; true.par.fre
sapply(1:3, function(x) summary(fit.bayes.fre$parameters[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(fit.bayes.fre$parameters[,i]), fit.bayes.fre$parameters[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.fre[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes.fre$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(fit.bayes.fre$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.fre[i], col = "blue", lty = 2, lwd = 2)
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

#Frequentist inference
fit.freq.par <- fit.gev.inference(data = data.par, t = th.par, t.prob = t.prob,
                                  llik.type = "Max-Gev-Cens",  param.type = "tilde",
                                  T.ret = T.ret, p = p, par0 = start.par,
                                  hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.par
#Comparison
#Parameters
fit.freq.par$mle; true.par.par
#Extreme quantile
fit.freq.par$Q.extreme; Q.ext.true.par
#Return level
fit.freq.par$R.level; R.lev.true.par

#Bayesian inference
tic()
fit.bayes.par <- fit.gev.inference(data = data.par, t = th.par, t.prob = t.prob,
                                   llik.type = "Max-Gev-Cens", param.type = "tilde",
                                   T.ret = T.ret, p = p, par0 = start.par,
                                   hessian = TRUE, inf.type = "Bayes", k = 1,
                                   R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
toc()
#57.35 sec elapsed
#Results
fit.bayes.par$mcmc.acc
#Comparison
#Parameters
fit.bayes.par$mle; fit.bayes.par$mode; true.par.par
sapply(1:3, function(x) summary(fit.bayes.par$parameters[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(fit.bayes.par$parameters[,i]), fit.bayes.par$parameters[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.par[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes.par$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(fit.bayes.par$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.par[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(fit.bayes.par$parameters[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.par
summary(fit.bayes.par$Q.extreme)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(fit.bayes.par$Q.extreme), fit.bayes.par$Q.extreme, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.par, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.par$Q.extreme, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.par$Q.extreme, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.par, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.par$Q.extreme, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.par
summary(fit.bayes.par$R.level)
#Trace plot
plot(1:length(fit.bayes.par$R.level), fit.bayes.par$R.level, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.par, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.par$R.level, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.par$R.level, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.par, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.par$R.level, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))
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

#Frequentist inference
fit.freq.hc <- fit.gev.inference(data = data.hc, t = th.hc, t.prob = t.prob,
                                 llik.type = "Max-Gev-Cens", param.type = "tilde",
                                 T.ret = T.ret, p = p, par0 = start.hc,
                                 hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.hc
#Comparison
#Parameters
fit.freq.hc$mle; true.par.hc
#Extreme quantile
fit.freq.hc$Q.extreme; Q.ext.true.hc
#Return level
fit.freq.hc$R.level; R.lev.true.hc

#Bayesian inference
tic()
fit.bayes.hc <- fit.gev.inference(data = data.hc, t = th.hc, t.prob = t.prob,
                                  llik.type = "Max-Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p, par0 = start.hc,
                                  hessian = TRUE, inf.type = "Bayes", k = 1, R = R,
                                  burn  = burn, prior = "empirical",
                                  val.show = TRUE)
toc()
#55.14 sec elapsed
#Results
fit.bayes.hc$mcmc.acc
#Comparison
#Parameters
fit.bayes.hc$mle; fit.bayes.hc$mode; true.par.hc
sapply(1:3, function(x) summary(fit.bayes.hc$parameters[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(fit.bayes.hc$parameters[,i]), fit.bayes.hc$parameters[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.hc[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes.hc$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(fit.bayes.hc$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.hc[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(fit.bayes.hc$parameters[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.hc
summary(fit.bayes.hc$Q.extreme)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(fit.bayes.hc$Q.extreme), fit.bayes.hc$Q.extreme, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.hc, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.hc$Q.extreme, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.hc$Q.extreme, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.hc, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.hc$Q.extreme, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.hc
summary(fit.bayes.hc$R.level)
#Trace plot
plot(1:length(fit.bayes.hc$R.level), fit.bayes.hc$R.level, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.hc, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.hc$R.level, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.hc$R.level, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.hc, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.hc$R.level, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))
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

#Frequentist inference
fit.freq.gum <- fit.gev.inference(data = data.gum, t = th.gum, t.prob = t.prob,
                                  llik.type = "Max-Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p, par0 = start.gum,
                                  hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.gum
#Comparison
#Parameters
fit.freq.gum$mle; true.par.gum
#Extreme quantile
fit.freq.gum$Q.extreme; Q.ext.true.gum
#Return level
fit.freq.gum$R.level; R.lev.true.gum

#Bayesian inference
tic()
fit.bayes.gum <- fit.gev.inference(data = data.gum, t = th.gum, t.prob = t.prob,
                                   llik.type = "Max-Gev-Cens", param.type = "tilde",
                                   T.ret = T.ret, p = p, par0 = start.gum,
                                   hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
toc()
#55.36 sec elapsed
#Results
fit.bayes.gum$mcmc.acc
#Comparison
#Parameters
fit.bayes.gum$mle; fit.bayes.gum$mode; true.par.gum
sapply(1:3, function(x) summary(fit.bayes.gum$parameters[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(fit.bayes.gum$parameters[,i]), fit.bayes.gum$parameters[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.gum[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes.gum$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(fit.bayes.gum$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.gum[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(fit.bayes.gum$parameters[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.gum
summary(fit.bayes.gum$Q.extreme)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(fit.bayes.gum$Q.extreme), fit.bayes.gum$Q.extreme, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.gum, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.gum$Q.extreme, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.gum$Q.extreme, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.gum, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.gum$Q.extreme, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.gum
summary(fit.bayes.gum$R.level)
#Trace plot
plot(1:length(fit.bayes.gum$R.level), fit.bayes.gum$R.level, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.gum, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.gum$R.level, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.gum$R.level, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.gum, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.gum$R.level, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))
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
set.seed(11)
data.exp <- rexp(n, rate.exp)
th.exp <- quantile(data.exp, probs = 1 - t.prob, type = 3)

#Frequentist inference
fit.freq.exp <- fit.gev.inference(data = data.exp, t = th.exp, t.prob = t.prob,
                                  llik.type = "Max-Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p, par0 = start.exp,
                                  hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.exp
#Comparison
#Parameters
fit.freq.exp$mle; true.par.exp
#Extreme quantile
fit.freq.exp$Q.extreme; Q.ext.true.exp
#Return level
fit.freq.exp$R.level; R.lev.true.exp

#Bayesian inference
tic()
fit.bayes.exp <- fit.gev.inference(data = data.exp, t = th.exp, t.prob = t.prob,
                                   llik.type = "Max-Gev-Cens", param.type = "tilde",
                                   T.ret = T.ret, p = p, par0 = start.exp,
                                   hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
toc()
#51.92 sec elapsed
#Results
fit.bayes.exp$mcmc.acc
#Comparison
#Parameters
fit.bayes.exp$mle; fit.bayes.exp$mode; true.par.exp
sapply(1:3, function(x) summary(fit.bayes.exp$parameters[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(fit.bayes.exp$parameters[,i]), fit.bayes.exp$parameters[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.exp[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes.exp$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(fit.bayes.exp$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.exp[i], col = "blue", lty = 2, lwd = 2)
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

#Frequentist inference
fit.freq.gam <- fit.gev.inference(data = data.gam, t = th.gam, t.prob = t.prob,
                                  llik.type = "Max-Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p, par0 = start.gam,
                                  hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.gam
#Comparison
#Parameters
fit.freq.gam$mle; true.par.gam
#Extreme quantile
fit.freq.gam$Q.extreme; Q.ext.true.gam
#Return level
fit.freq.gam$R.level; R.lev.true.gam

#Bayesian inference
tic()
fit.bayes.gam <- fit.gev.inference(data = data.gam, t = th.gam, t.prob = t.prob,
                                   llik.type = "Max-Gev-Cens", param.type = "tilde",
                                   T.ret = T.ret, p = p, par0 = start.gam,
                                   hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
toc()
#59.13 sec elapsed
#Results
fit.bayes.gam$mcmc.acc
#Comparison
#Parameters
fit.bayes.gam$mle; fit.bayes.gam$mode; true.par.gam
sapply(1:3, function(x) summary(fit.bayes.gam$parameters[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(fit.bayes.gam$parameters[,i]), fit.bayes.gam$parameters[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.gam[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes.gam$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(fit.bayes.gam$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.gam[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(fit.bayes.gam$parameters[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.gam
summary(fit.bayes.gam$Q.extreme)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(fit.bayes.gam$Q.extreme), fit.bayes.gam$Q.extreme, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.gam, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.gam$Q.extreme, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.gam$Q.extreme, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.gam, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.gam$Q.extreme, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.gam
summary(fit.bayes.gam$R.level)
#Trace plot
plot(1:length(fit.bayes.gam$R.level), fit.bayes.gam$R.level, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.gam, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.gam$R.level, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.gam$R.level, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.gam, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.gam$R.level, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))
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
set.seed(12)
data.pl <- pl.qdf(runif(n), xstar.pl, alpha.pl, K = K)
th.pl <- quantile(data.pl, probs = 1 - t.prob, type = 3)

#Frequentist inference
fit.freq.pl <- fit.gev.inference(data = data.pl, t = th.pl, t.prob = t.prob,
                                 llik.type = "Max-Gev-Cens", param.type = "tilde",
                                 T.ret = T.ret, p = p, par0 = start.pl,
                                 hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.pl
#Comparison
#Parameters
fit.freq.pl$mle; true.par.pl
#Extreme quantile
fit.freq.pl$Q.extreme; Q.ext.true.pl
#Return level
fit.freq.pl$R.level; R.lev.true.pl

#Bayesian inference
tic()
fit.bayes.pl <- fit.gev.inference(data = data.pl, t = th.pl, t.prob = t.prob,
                                  llik.type = "Max-Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p, par0 = start.pl,
                                  hessian = TRUE, inf.type = "Bayes",
                                  k = 1, R = R, burn  = burn, prior = "empirical",
                                  val.show = TRUE)
toc()
#56.58 sec elapsed
#Results
fit.bayes.pl$mcmc.acc
#Comparison
#Parameters
fit.bayes.pl$mle; fit.bayes.pl$mode; true.par.pl
sapply(1:3, function(x) summary(fit.bayes.pl$parameters[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(fit.bayes.pl$parameters[,i]), fit.bayes.pl$parameters[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.pl[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes.pl$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(fit.bayes.pl$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.pl[i], col = "blue", lty = 2, lwd = 2)
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

#Frequentist inference
fit.freq.rw <- fit.gev.inference(data = data.rw, t = th.rw, t.prob = t.prob,
                                 llik.type = "Max-Gev-Cens", param.type = "tilde",
                                 T.ret = T.ret, p = p, par0 = start.rw,
                                 hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.rw
#Comparison
#Parameters
fit.freq.rw$mle; true.par.rw
#Extreme quantile
fit.freq.rw$Q.extreme; Q.ext.true.rw
#Return level
fit.freq.rw$R.level; R.lev.true.rw

#Bayesian inference
tic()
fit.bayes.rw <- fit.gev.inference(data = data.rw, t = th.rw, t.prob = t.prob,
                                  llik.type = "Max-Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p, par0 = start.rw,
                                  hessian = TRUE, inf.type = "Bayes",
                                  k = 1, R = R, burn  = burn, prior = "empirical",
                                  val.show = TRUE)
toc()
#64.71 sec elapsed
#Results
fit.bayes.rw$mcmc.acc
#Comparison
#Parameters
fit.bayes.rw$mle; fit.bayes.rw$mode; true.par.rw
sapply(1:3, function(x) summary(fit.bayes.rw$parameters[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(fit.bayes.rw$parameters[,i]), fit.bayes.rw$parameters[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.rw[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes.rw$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(fit.bayes.rw$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.rw[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(fit.bayes.rw$parameters[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.rw
summary(fit.bayes.rw$Q.extreme)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(fit.bayes.rw$Q.extreme), fit.bayes.rw$Q.extreme, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.rw, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.rw$Q.extreme, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.rw$Q.extreme, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.rw, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.rw$Q.extreme, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.rw
summary(fit.bayes.rw$R.level)
#Trace plot
plot(1:length(fit.bayes.rw$R.level), fit.bayes.rw$R.level, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.rw, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.rw$R.level, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.rw$R.level, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.rw, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.rw$R.level, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))
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
set.seed(11)
data.be <- rbeta(n, shape1 = shape1.be, shape2 = shape2.be)
th.be <- quantile(data.be, probs = 1 - t.prob, type = 3)

#Frequentist inference
fit.freq.be <- fit.gev.inference(data = data.be, t = th.be, t.prob = t.prob,
                                 llik.type = "Max-Gev-Cens", param.type = "tilde",
                                 T.ret = T.ret, p = p, par0 = start.be,
                                 hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.be
#Comparison
#Parameters
fit.freq.be$mle; true.par.be
#Extreme quantile
fit.freq.be$Q.extreme; Q.ext.true.be
#Return level
fit.freq.be$R.level; R.lev.true.be

#Bayesian inference
tic()
fit.bayes.be <- fit.gev.inference(data = data.be, t = th.be, t.prob = t.prob,
                                  llik.type = "Max-Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p, par0 = start.be,
                                  hessian = TRUE, inf.type = "Bayes",
                                  k = 1, R = R, burn  = burn, prior = "empirical",
                                  val.show = TRUE)
toc()
#45.57 sec elapsed
#Results
fit.bayes.be$mcmc.acc
#Comparison
#Parameters
fit.bayes.be$mle; fit.bayes.be$mode; true.par.be
sapply(1:3, function(x) summary(fit.bayes.be$parameters[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(fit.bayes.be$parameters[,i]), fit.bayes.be$parameters[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.be[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(fit.bayes.be$parameters[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(fit.bayes.be$parameters[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.be[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(fit.bayes.be$parameters[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.be
summary(fit.bayes.be$Q.extreme)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(fit.bayes.be$Q.extreme), fit.bayes.be$Q.extreme, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.be, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.be$Q.extreme, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.be$Q.extreme, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.be, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.be$Q.extreme, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.be
summary(fit.bayes.be$R.level)
#Trace plot
plot(1:length(fit.bayes.be$R.level), fit.bayes.be$R.level, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.be, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(fit.bayes.be$R.level, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(fit.bayes.be$R.level, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.be, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(fit.bayes.be$R.level, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))
par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#Memory usage

sort(sapply(ls(),function(x){object.size(get(x))}))
mem <- sum(sort(sapply(ls(),function(x){object.size(get(x))})))
if(mem < 10^9) cat("\n Weight of created and used objects:", mem/(10^6), "Mbytes")
if(mem >= 10^9) cat("\n Weight of created and used objects:", mem/(10^9), "Gbytes")

#-------------------------------------------------------------------------------

save.image("models_max.RData")

#-------------------------------------------------------------------------------