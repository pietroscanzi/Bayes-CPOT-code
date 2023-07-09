#-------------------------------------------------------------------------------

#Trials of censored-likelihood threshold Gev inference methods for different
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
#Proportion of exceedances
s <- n/k
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
                                  llik.type = "Gev-Cens", param.type = "bar",
                                  T.ret = T.ret, p = p, par0 = start.fre,
                                  hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.fre
#Tilde parameterization
scale.fre <- fit.freq.fre$mle[3]*s^(fit.freq.fre$mle[1])
loc.fre <- fit.freq.fre$mle[2] + 
            scale.fre*(1 - s^(-fit.freq.fre$mle[1]))/fit.freq.fre$mle[1]
mle.fre <- c(fit.freq.fre$mle[1], loc.fre, scale.fre)
names(mle.fre) <- c("Shape", "Location", "Scale")
Q.ext.fre <- mle.fre[2] + mle.fre[3]*((s*p)^(-mle.fre[1]) - 1)/mle.fre[1]
names(Q.ext.fre) <- "Extreme-Quantile"
R.lev.fre <- mle.fre[2] + 
              mle.fre[3]*((-log(1 - 1/T.ret))^(-mle.fre[1]) - 1)/mle.fre[1]
names(R.lev.fre) <- "Return-Level"
#Comparison
#Parameters
mle.fre; true.par.fre
#Extreme quantile
Q.ext.fre; Q.ext.true.fre
#Return level
R.lev.fre; R.lev.true.fre

#Bayesian inference
tic()
fit.bayes.fre <- fit.gev.inference(data = data.fre, t = th.fre, t.prob = t.prob,
                                   llik.type = "Gev-Cens", param.type = "bar",
                                   T.ret = T.ret, p = p, par0 = start.fre,
                                   hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
toc()
#82.51 sec elapsed
#Results
fit.bayes.fre$mcmc.acc
#Tilde parameterization
mode.tilde.fre <- fit.bayes.fre$mode
mode.tilde.fre[3] <- fit.bayes.fre$mode[3]*s^(fit.bayes.fre$mode[1])
mode.tilde.fre[2] <- mode.tilde.fre[2] + 
              mode.tilde.fre[3]*(1 - s^(-mode.tilde.fre[1]))/mode.tilde.fre[1]
params.tilde.fre <- fit.bayes.fre$parameters
scale.post.fre <- apply(fit.bayes.fre$parameters, 1, function(x) x[3]*s^(x[1]))
params.tilde.fre[,3] <- scale.post.fre
loc.post.fre <- apply(params.tilde.fre, 1,
                      function(x) x[2] + x[3]*(1 - s^(-x[1]))/x[1])
params.tilde.fre[,2] <- loc.post.fre
Q.ext.post.fre <- params.tilde.fre[,2] + 
  params.tilde.fre[,3]*((s*p)^(-params.tilde.fre[,1]) - 1)/params.tilde.fre[,1]
R.lev.post.fre <- params.tilde.fre[,2] +
  params.tilde.fre[,3]*((-log(1 - 1/T.ret))^(-params.tilde.fre[,1]) - 1)/params.tilde.fre[,1]
#Comparison
#Parameters
mle.fre; mode.tilde.fre; true.par.fre
sapply(1:3, function(x) summary(params.tilde.fre[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(params.tilde.fre[,i]), params.tilde.fre[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.fre[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(params.tilde.fre[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(params.tilde.fre[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.fre[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(params.tilde.fre[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.fre
summary(Q.ext.post.fre)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(Q.ext.post.fre), Q.ext.post.fre, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.fre, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(Q.ext.post.fre, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(Q.ext.post.fre, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.fre, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(Q.ext.post.fre, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.fre
summary(R.lev.post.fre)
#Trace plot
plot(1:length(R.lev.post.fre), R.lev.post.fre, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.fre, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(R.lev.post.fre, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(R.lev.post.fre, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.fre, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(R.lev.post.fre, conf = 0.95), col = "red",
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
                                  llik.type = "Gev-Cens",  param.type = "bar",
                                  T.ret = T.ret, p = p, par0 = start.par,
                                  hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.par
#Tilde parameterization
scale.par <- fit.freq.par$mle[3]*s^(fit.freq.par$mle[1])
loc.par <- fit.freq.par$mle[2] + 
  scale.par*(1 - s^(-fit.freq.par$mle[1]))/fit.freq.par$mle[1]
mle.par <- c(fit.freq.par$mle[1], loc.par, scale.par)
names(mle.par) <- c("Shape", "Location", "Scale")
Q.ext.par <- mle.par[2] + mle.par[3]*((s*p)^(-mle.par[1]) - 1)/mle.par[1]
names(Q.ext.par) <- "Extreme-Quantile"
R.lev.par <- mle.par[2] + 
  mle.par[3]*((-log(1 - 1/T.ret))^(-mle.par[1]) - 1)/mle.par[1]
names(R.lev.par) <- "Return-Level"
#Comparison
#Parameters
mle.par; true.par.par
#Extreme quantile
Q.ext.par; Q.ext.true.par
#Return level
R.lev.par; R.lev.true.par

#Bayesian inference
tic()
fit.bayes.par <- fit.gev.inference(data = data.par, t = th.par, t.prob = t.prob,
                                   llik.type = "Gev-Cens", param.type = "bar",
                                   T.ret = T.ret, p = p, par0 = start.par,
                                   hessian = TRUE, inf.type = "Bayes", k = 1,
                                   R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
toc()
#82.84 sec elapsed
#Results
fit.bayes.par$mcmc.acc
#Tilde parameterization
mode.tilde.par <- fit.bayes.par$mode
mode.tilde.par[3] <- fit.bayes.par$mode[3]*s^(fit.bayes.par$mode[1])
mode.tilde.par[2] <- mode.tilde.par[2] + 
  mode.tilde.par[3]*(1 - s^(-mode.tilde.par[1]))/mode.tilde.par[1]
params.tilde.par <- fit.bayes.par$parameters
scale.post.par <- apply(fit.bayes.par$parameters, 1, function(x) x[3]*s^(x[1]))
params.tilde.par[,3] <- scale.post.par
loc.post.par <- apply(params.tilde.par, 1,
                      function(x) x[2] + x[3]*(1 - s^(-x[1]))/x[1])
params.tilde.par[,2] <- loc.post.par
Q.ext.post.par <- params.tilde.par[,2] + 
  params.tilde.par[,3]*((s*p)^(-params.tilde.par[,1]) - 1)/params.tilde.par[,1]
R.lev.post.par <- params.tilde.par[,2] +
  params.tilde.par[,3]*((-log(1 - 1/T.ret))^(-params.tilde.par[,1]) - 1)/params.tilde.par[,1]
#Comparison
#Parameters
mle.par; mode.tilde.par; true.par.par
sapply(1:3, function(x) summary(params.tilde.par[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(params.tilde.par[,i]), params.tilde.par[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.par[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(params.tilde.par[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(params.tilde.par[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.par[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(params.tilde.par[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.par
summary(Q.ext.post.par)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(Q.ext.post.par), Q.ext.post.par, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.par, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(Q.ext.post.par, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(Q.ext.post.par, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.par, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(Q.ext.post.par, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.par
summary(R.lev.post.par)
#Trace plot
plot(1:length(R.lev.post.par), R.lev.post.par, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.par, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(R.lev.post.par, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(R.lev.post.par, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.par, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(R.lev.post.par, conf = 0.95), col = "red",
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
                                 llik.type = "Gev-Cens", param.type = "bar",
                                 T.ret = T.ret, p = p, par0 = start.hc,
                                 hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.hc
#Tilde parameterization
scale.hc <- fit.freq.hc$mle[3]*s^(fit.freq.hc$mle[1])
loc.hc <- fit.freq.hc$mle[2] + 
  scale.hc*(1 - s^(-fit.freq.hc$mle[1]))/fit.freq.hc$mle[1]
mle.hc <- c(fit.freq.hc$mle[1], loc.hc, scale.hc)
names(mle.hc) <- c("Shape", "Location", "Scale")
Q.ext.hc <- mle.hc[2] + mle.hc[3]*((s*p)^(-mle.hc[1]) - 1)/mle.hc[1]
names(Q.ext.hc) <- "Extreme-Quantile"
R.lev.hc <- mle.hc[2] + 
  mle.hc[3]*((-log(1 - 1/T.ret))^(-mle.hc[1]) - 1)/mle.hc[1]
names(R.lev.hc) <- "Return-Level"
#Comparison
#Parameters
mle.hc; true.par.hc
#Extreme quantile
Q.ext.hc; Q.ext.true.hc
#Return level
R.lev.hc; R.lev.true.hc

#Bayesian inference
tic()
fit.bayes.hc <- fit.gev.inference(data = data.hc, t = th.hc, t.prob = t.prob,
                                  llik.type = "Gev-Cens", param.type = "bar",
                                  T.ret = T.ret, p = p, par0 = start.hc,
                                  hessian = TRUE, inf.type = "Bayes", k = 1, R = R,
                                  burn  = burn, prior = "empirical",
                                  val.show = TRUE)
toc()
#82.37 sec elapsed
#Results
fit.bayes.hc$mcmc.acc
#Tilde parameterization
mode.tilde.hc <- fit.bayes.hc$mode
mode.tilde.hc[3] <- fit.bayes.hc$mode[3]*s^(fit.bayes.hc$mode[1])
mode.tilde.hc[2] <- mode.tilde.hc[2] + 
  mode.tilde.hc[3]*(1 - s^(-mode.tilde.hc[1]))/mode.tilde.hc[1]
params.tilde.hc <- fit.bayes.hc$parameters
scale.post.hc <- apply(fit.bayes.hc$parameters, 1, function(x) x[3]*s^(x[1]))
params.tilde.hc[,3] <- scale.post.hc
loc.post.hc <- apply(params.tilde.hc, 1,
                      function(x) x[2] + x[3]*(1 - s^(-x[1]))/x[1])
params.tilde.hc[,2] <- loc.post.hc
Q.ext.post.hc <- params.tilde.hc[,2] + 
  params.tilde.hc[,3]*((s*p)^(-params.tilde.hc[,1]) - 1)/params.tilde.hc[,1]
R.lev.post.hc <- params.tilde.hc[,2] +
  params.tilde.hc[,3]*((-log(1 - 1/T.ret))^(-params.tilde.hc[,1]) - 1)/params.tilde.hc[,1]
#Comparison
#Parameters
mle.hc; mode.tilde.hc; true.par.hc
sapply(1:3, function(x) summary(params.tilde.hc[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(params.tilde.hc[,i]), params.tilde.hc[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.hc[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(params.tilde.hc[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(params.tilde.hc[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.hc[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(params.tilde.hc[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.hc
summary(Q.ext.post.hc)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(Q.ext.post.hc), Q.ext.post.hc, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.hc, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(Q.ext.post.hc, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(Q.ext.post.hc, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.hc, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(Q.ext.post.hc, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.hc
summary(R.lev.post.hc)
#Trace plot
plot(1:length(R.lev.post.hc), R.lev.post.hc, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.hc, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(R.lev.post.hc, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(R.lev.post.hc, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.hc, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(R.lev.post.hc, conf = 0.95), col = "red",
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
                                  llik.type = "Gev-Cens", param.type = "bar",
                                  T.ret = T.ret, p = p, par0 = start.gum,
                                  hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.gum
#Tilde parameterization
scale.gum <- fit.freq.gum$mle[3]*s^(fit.freq.gum$mle[1])
loc.gum <- fit.freq.gum$mle[2] + 
  scale.gum*(1 - s^(-fit.freq.gum$mle[1]))/fit.freq.gum$mle[1]
mle.gum <- c(fit.freq.gum$mle[1], loc.gum, scale.gum)
names(mle.gum) <- c("Shape", "Location", "Scale")
Q.ext.gum <- mle.gum[2] + mle.gum[3]*((s*p)^(-mle.gum[1]) - 1)/mle.gum[1]
names(Q.ext.gum) <- "Extreme-Quantile"
R.lev.gum <- mle.gum[2] + 
  mle.gum[3]*((-log(1 - 1/T.ret))^(-mle.gum[1]) - 1)/mle.gum[1]
names(R.lev.gum) <- "Return-Level"
#Comparison
#Parameters
mle.gum; true.par.gum
#Extreme quantile
Q.ext.gum; Q.ext.true.gum
#Return level
R.lev.gum; R.lev.true.gum

#Bayesian inference
tic()
fit.bayes.gum <- fit.gev.inference(data = data.gum, t = th.gum, t.prob = t.prob,
                                   llik.type = "Gev-Cens", param.type = "bar",
                                   T.ret = T.ret, p = p, par0 = start.gum,
                                   hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
toc()
#72.67 sec elapsed
#Results
fit.bayes.gum$mcmc.acc
#Tilde parameterization
mode.tilde.gum <- fit.bayes.gum$mode
mode.tilde.gum[3] <- fit.bayes.gum$mode[3]*s^(fit.bayes.gum$mode[1])
mode.tilde.gum[2] <- mode.tilde.gum[2] + 
  mode.tilde.gum[3]*(1 - s^(-mode.tilde.gum[1]))/mode.tilde.gum[1]
params.tilde.gum <- fit.bayes.gum$parameters
scale.post.gum <- apply(fit.bayes.gum$parameters, 1, function(x) x[3]*s^(x[1]))
params.tilde.gum[,3] <- scale.post.gum
loc.post.gum <- apply(params.tilde.gum, 1,
                     function(x) x[2] + x[3]*(1 - s^(-x[1]))/x[1])
params.tilde.gum[,2] <- loc.post.gum
Q.ext.post.gum <- params.tilde.gum[,2] + 
  params.tilde.gum[,3]*((s*p)^(-params.tilde.gum[,1]) - 1)/params.tilde.gum[,1]
R.lev.post.gum <- params.tilde.gum[,2] +
  params.tilde.gum[,3]*((-log(1 - 1/T.ret))^(-params.tilde.gum[,1]) - 1)/params.tilde.gum[,1]
#Comparison
#Parameters
mle.gum; mode.tilde.gum; true.par.gum
sapply(1:3, function(x) summary(params.tilde.gum[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(params.tilde.gum[,i]), params.tilde.gum[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.gum[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(params.tilde.gum[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(params.tilde.gum[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.gum[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(params.tilde.gum[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.gum
summary(Q.ext.post.gum)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(Q.ext.post.gum), Q.ext.post.gum, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.gum, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(Q.ext.post.gum, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(Q.ext.post.gum, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.gum, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(Q.ext.post.gum, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.gum
summary(R.lev.post.gum)
#Trace plot
plot(1:length(R.lev.post.gum), R.lev.post.gum, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.gum, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(R.lev.post.gum, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(R.lev.post.gum, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.gum, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(R.lev.post.gum, conf = 0.95), col = "red",
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
                                  llik.type = "Gev-Cens", param.type = "bar",
                                  T.ret = T.ret, p = p, par0 = start.exp,
                                  hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.exp
#Tilde parameterization
scale.exp <- fit.freq.exp$mle[3]*s^(fit.freq.exp$mle[1])
loc.exp <- fit.freq.exp$mle[2] + 
  scale.exp*(1 - s^(-fit.freq.exp$mle[1]))/fit.freq.exp$mle[1]
mle.exp <- c(fit.freq.exp$mle[1], loc.exp, scale.exp)
names(mle.exp) <- c("Shape", "Location", "Scale")
Q.ext.exp <- mle.exp[2] + mle.exp[3]*((s*p)^(-mle.exp[1]) - 1)/mle.exp[1]
names(Q.ext.exp) <- "Extreme-Quantile"
R.lev.exp <- mle.exp[2] + 
  mle.exp[3]*((-log(1 - 1/T.ret))^(-mle.exp[1]) - 1)/mle.exp[1]
names(R.lev.exp) <- "Return-Level"
#Comparison
#Parameters
mle.exp; true.par.exp
#Extreme quantile
Q.ext.exp; Q.ext.true.exp
#Return level
R.lev.exp; R.lev.true.exp

#Bayesian inference
tic()
fit.bayes.exp <- fit.gev.inference(data = data.exp, t = th.exp, t.prob = t.prob,
                                   llik.type = "Gev-Cens", param.type = "bar",
                                   T.ret = T.ret, p = p, par0 = start.exp,
                                   hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
toc()
#73.91 sec elapsed
#Results
fit.bayes.exp$mcmc.acc
#Tilde parameterization
mode.tilde.exp <- fit.bayes.exp$mode
mode.tilde.exp[3] <- fit.bayes.exp$mode[3]*s^(fit.bayes.exp$mode[1])
mode.tilde.exp[2] <- mode.tilde.exp[2] + 
  mode.tilde.exp[3]*(1 - s^(-mode.tilde.exp[1]))/mode.tilde.exp[1]
params.tilde.exp <- fit.bayes.exp$parameters
scale.post.exp <- apply(fit.bayes.exp$parameters, 1, function(x) x[3]*s^(x[1]))
params.tilde.exp[,3] <- scale.post.exp
loc.post.exp <- apply(params.tilde.exp, 1,
                     function(x) x[2] + x[3]*(1 - s^(-x[1]))/x[1])
params.tilde.exp[,2] <- loc.post.exp
Q.ext.post.exp <- params.tilde.exp[,2] + 
  params.tilde.exp[,3]*((s*p)^(-params.tilde.exp[,1]) - 1)/params.tilde.exp[,1]
R.lev.post.exp <- params.tilde.exp[,2] +
  params.tilde.exp[,3]*((-log(1 - 1/T.ret))^(-params.tilde.exp[,1]) - 1)/params.tilde.exp[,1]
#Comparison
#Parameters
mle.exp; mode.tilde.exp; true.par.exp
sapply(1:3, function(x) summary(params.tilde.exp[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(params.tilde.exp[,i]), params.tilde.exp[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.exp[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(params.tilde.exp[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(params.tilde.exp[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.exp[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(params.tilde.exp[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.exp
summary(Q.ext.post.exp)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(Q.ext.post.exp), Q.ext.post.exp, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.exp, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(Q.ext.post.exp, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(Q.ext.post.exp, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.exp, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(Q.ext.post.exp, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.exp
summary(R.lev.post.exp)
#Trace plot
plot(1:length(R.lev.post.exp), R.lev.post.exp, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.exp, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(R.lev.post.exp, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(R.lev.post.exp, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.exp, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(R.lev.post.exp, conf = 0.95), col = "red",
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
set.seed(2)
data.gam <- rgamma(n, shape.gam, rate.gam)
th.gam <- quantile(data.gam, probs = 1 - t.prob, type = 3)

#Frequentist inference
fit.freq.gam <- fit.gev.inference(data = data.gam, t = th.gam, t.prob = t.prob,
                                  llik.type = "Gev-Cens", param.type = "bar",
                                  T.ret = T.ret, p = p, par0 = start.gam,
                                  hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.gam
#Tilde parameterization
scale.gam <- fit.freq.gam$mle[3]*s^(fit.freq.gam$mle[1])
loc.gam <- fit.freq.gam$mle[2] + 
  scale.gam*(1 - s^(-fit.freq.gam$mle[1]))/fit.freq.gam$mle[1]
mle.gam <- c(fit.freq.gam$mle[1], loc.gam, scale.gam)
names(mle.gam) <- c("Shape", "Location", "Scale")
Q.ext.gam <- mle.gam[2] + mle.gam[3]*((s*p)^(-mle.gam[1]) - 1)/mle.gam[1]
names(Q.ext.gam) <- "Extreme-Quantile"
R.lev.gam <- mle.gam[2] + 
  mle.gam[3]*((-log(1 - 1/T.ret))^(-mle.gam[1]) - 1)/mle.gam[1]
names(R.lev.gam) <- "Return-Level"
#Comparison
#Parameters
mle.gam; true.par.gam
#Extreme quantile
Q.ext.gam; Q.ext.true.gam
#Return level
R.lev.gam; R.lev.true.gam

#Bayesian inference
tic()
fit.bayes.gam <- fit.gev.inference(data = data.gam, t = th.gam, t.prob = t.prob,
                                   llik.type = "Gev-Cens", param.type = "bar",
                                   T.ret = T.ret, p = p, par0 = start.gam,
                                   hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)
toc()
#69.2 sec elapsed
#Results
fit.bayes.gam$mcmc.acc
#Tilde parameterization
mode.tilde.gam <- fit.bayes.gam$mode
mode.tilde.gam[3] <- fit.bayes.gam$mode[3]*s^(fit.bayes.gam$mode[1])
mode.tilde.gam[2] <- mode.tilde.gam[2] + 
  mode.tilde.gam[3]*(1 - s^(-mode.tilde.gam[1]))/mode.tilde.gam[1]
params.tilde.gam <- fit.bayes.gam$parameters
scale.post.gam <- apply(fit.bayes.gam$parameters, 1, function(x) x[3]*s^(x[1]))
params.tilde.gam[,3] <- scale.post.gam
loc.post.gam <- apply(params.tilde.gam, 1,
                     function(x) x[2] + x[3]*(1 - s^(-x[1]))/x[1])
params.tilde.gam[,2] <- loc.post.gam
Q.ext.post.gam <- params.tilde.gam[,2] + 
  params.tilde.gam[,3]*((s*p)^(-params.tilde.gam[,1]) - 1)/params.tilde.gam[,1]
R.lev.post.gam <- params.tilde.gam[,2] +
  params.tilde.gam[,3]*((-log(1 - 1/T.ret))^(-params.tilde.gam[,1]) - 1)/params.tilde.gam[,1]
#Comparison
#Parameters
mle.gam; mode.tilde.gam; true.par.gam
sapply(1:3, function(x) summary(params.tilde.gam[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(params.tilde.gam[,i]), params.tilde.gam[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.gam[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(params.tilde.gam[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(params.tilde.gam[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.gam[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(params.tilde.gam[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.gam
summary(Q.ext.post.gam)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(Q.ext.post.gam), Q.ext.post.gam, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.gam, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(Q.ext.post.gam, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(Q.ext.post.gam, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.gam, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(Q.ext.post.gam, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.gam
summary(R.lev.post.gam)
#Trace plot
plot(1:length(R.lev.post.gam), R.lev.post.gam, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.gam, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(R.lev.post.gam, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(R.lev.post.gam, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.gam, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(R.lev.post.gam, conf = 0.95), col = "red",
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
                                 llik.type = "Gev-Cens", param.type = "bar",
                                 T.ret = T.ret, p = p, par0 = start.pl,
                                 hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.pl
#Tilde parameterization
scale.pl <- fit.freq.pl$mle[3]*s^(fit.freq.pl$mle[1])
loc.pl <- fit.freq.pl$mle[2] + 
  scale.pl*(1 - s^(-fit.freq.pl$mle[1]))/fit.freq.pl$mle[1]
mle.pl <- c(fit.freq.pl$mle[1], loc.pl, scale.pl)
names(mle.pl) <- c("Shape", "Location", "Scale")
Q.ext.pl <- mle.pl[2] + mle.pl[3]*((s*p)^(-mle.pl[1]) - 1)/mle.pl[1]
names(Q.ext.pl) <- "Extreme-Quantile"
R.lev.pl <- mle.pl[2] + 
  mle.pl[3]*((-log(1 - 1/T.ret))^(-mle.pl[1]) - 1)/mle.pl[1]
names(R.lev.pl) <- "Return-Level"
#Comparison
#Parameters
mle.pl; true.par.pl
#Extreme quantile
Q.ext.pl; Q.ext.true.pl
#Return level
R.lev.pl; R.lev.true.pl

#Bayesian inference
tic()
fit.bayes.pl <- fit.gev.inference(data = data.pl, t = th.pl, t.prob = t.prob,
                                  llik.type = "Gev-Cens", param.type = "bar",
                                  T.ret = T.ret, p = p, par0 = start.pl,
                                  hessian = TRUE, inf.type = "Bayes",
                                  k = 1, R = R, burn  = burn, prior = "empirical",
                                  val.show = TRUE)
toc()
#73.54 sec elapsed
#Results
fit.bayes.pl$mcmc.acc
#Tilde parameterization
mode.tilde.pl <- fit.bayes.pl$mode
mode.tilde.pl[3] <- fit.bayes.pl$mode[3]*s^(fit.bayes.pl$mode[1])
mode.tilde.pl[2] <- mode.tilde.pl[2] + 
  mode.tilde.pl[3]*(1 - s^(-mode.tilde.pl[1]))/mode.tilde.pl[1]
params.tilde.pl <- fit.bayes.pl$parameters
scale.post.pl <- apply(fit.bayes.pl$parameters, 1, function(x) x[3]*s^(x[1]))
params.tilde.pl[,3] <- scale.post.pl
loc.post.pl <- apply(params.tilde.pl, 1,
                     function(x) x[2] + x[3]*(1 - s^(-x[1]))/x[1])
params.tilde.pl[,2] <- loc.post.pl
Q.ext.post.pl <- params.tilde.pl[,2] + 
  params.tilde.pl[,3]*((s*p)^(-params.tilde.pl[,1]) - 1)/params.tilde.pl[,1]
R.lev.post.pl <- params.tilde.pl[,2] +
  params.tilde.pl[,3]*((-log(1 - 1/T.ret))^(-params.tilde.pl[,1]) - 1)/params.tilde.pl[,1]
#Comparison
#Parameters
mle.pl; mode.tilde.pl; true.par.pl
sapply(1:3, function(x) summary(params.tilde.pl[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(params.tilde.pl[,i]), params.tilde.pl[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.pl[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(params.tilde.pl[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(params.tilde.pl[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.pl[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(params.tilde.pl[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.pl
summary(Q.ext.post.pl)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(Q.ext.post.pl), Q.ext.post.pl, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.pl, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(Q.ext.post.pl, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(Q.ext.post.pl, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.pl, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(Q.ext.post.pl, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.pl
summary(R.lev.post.pl)
#Trace plot
plot(1:length(R.lev.post.pl), R.lev.post.pl, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.pl, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(R.lev.post.pl, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(R.lev.post.pl, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.pl, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(R.lev.post.pl, conf = 0.95), col = "red",
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
start.rw <- c(-0.46, -0.37, 0.13)
set.seed(50)
data.rw <- rrweibull(n, shape = shape.rw)
th.rw <- quantile(data.rw, probs = 1 - t.prob, type = 3)

#Frequentist inference
fit.freq.rw <- fit.gev.inference(data = data.rw, t = th.rw, t.prob = t.prob,
                                 llik.type = "Gev-Cens", param.type = "bar",
                                 T.ret = T.ret, p = p, par0 = start.rw,
                                 hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.rw
#Tilde parameterization
scale.rw <- fit.freq.rw$mle[3]*s^(fit.freq.rw$mle[1])
loc.rw <- fit.freq.rw$mle[2] + 
  scale.rw*(1 - s^(-fit.freq.rw$mle[1]))/fit.freq.rw$mle[1]
mle.rw <- c(fit.freq.rw$mle[1], loc.rw, scale.rw)
names(mle.rw) <- c("Shape", "Location", "Scale")
Q.ext.rw <- mle.rw[2] + mle.rw[3]*((s*p)^(-mle.rw[1]) - 1)/mle.rw[1]
names(Q.ext.rw) <- "Extreme-Quantile"
R.lev.rw <- mle.rw[2] + 
  mle.rw[3]*((-log(1 - 1/T.ret))^(-mle.rw[1]) - 1)/mle.rw[1]
names(R.lev.rw) <- "Return-Level"
#Comparison
#Parameters
mle.rw; true.par.rw
#Extreme quantile
Q.ext.rw; Q.ext.true.rw
#Return level
R.lev.rw; R.lev.true.rw

#Bayesian inference
tic()
fit.bayes.rw <- fit.gev.inference(data = data.rw, t = th.rw, t.prob = t.prob,
                                  llik.type = "Gev-Cens", param.type = "bar",
                                  T.ret = T.ret, p = p, par0 = start.rw,
                                  hessian = TRUE, inf.type = "Bayes",
                                  k = 1, R = R, burn  = burn, prior = "empirical",
                                  val.show = TRUE)
toc()
#71.11 sec elapsed
#Results
fit.bayes.rw$mcmc.acc
#Tilde parameterization
mode.tilde.rw <- fit.bayes.rw$mode
mode.tilde.rw[3] <- fit.bayes.rw$mode[3]*s^(fit.bayes.rw$mode[1])
mode.tilde.rw[2] <- mode.tilde.rw[2] + 
  mode.tilde.rw[3]*(1 - s^(-mode.tilde.rw[1]))/mode.tilde.rw[1]
params.tilde.rw <- fit.bayes.rw$parameters
scale.post.rw <- apply(fit.bayes.rw$parameters, 1, function(x) x[3]*s^(x[1]))
params.tilde.rw[,3] <- scale.post.rw
loc.post.rw <- apply(params.tilde.rw, 1,
                     function(x) x[2] + x[3]*(1 - s^(-x[1]))/x[1])
params.tilde.rw[,2] <- loc.post.rw
Q.ext.post.rw <- params.tilde.rw[,2] + 
  params.tilde.rw[,3]*((s*p)^(-params.tilde.rw[,1]) - 1)/params.tilde.rw[,1]
R.lev.post.rw <- params.tilde.rw[,2] +
  params.tilde.rw[,3]*((-log(1 - 1/T.ret))^(-params.tilde.rw[,1]) - 1)/params.tilde.rw[,1]
#Comparison
#Parameters
mle.rw; mode.tilde.rw; true.par.rw
sapply(1:3, function(x) summary(params.tilde.rw[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(params.tilde.rw[,i]), params.tilde.rw[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.rw[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(params.tilde.rw[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(params.tilde.rw[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.rw[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(params.tilde.rw[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.rw
summary(Q.ext.post.rw)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(Q.ext.post.rw), Q.ext.post.rw, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.rw, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(Q.ext.post.rw, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(Q.ext.post.rw, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.rw, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(Q.ext.post.rw, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.rw
summary(R.lev.post.rw)
#Trace plot
plot(1:length(R.lev.post.rw), R.lev.post.rw, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.rw, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(R.lev.post.rw, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(R.lev.post.rw, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.rw, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(R.lev.post.rw, conf = 0.95), col = "red",
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
set.seed(18)
data.be <- rbeta(n, shape1 = shape1.be, shape2 = shape2.be)
th.be <- quantile(data.be, probs = 1 - t.prob, type = 3)

#Frequentist inference
fit.freq.be <- fit.gev.inference(data = data.be, t = th.be, t.prob = t.prob,
                                 llik.type = "Gev-Cens", param.type = "bar",
                                 T.ret = T.ret, p = p, par0 = start.be,
                                 hessian = TRUE, inf.type = "Frequent")
#Results
fit.freq.be
#Tilde parameterization
scale.be <- fit.freq.be$mle[3]*s^(fit.freq.be$mle[1])
loc.be <- fit.freq.be$mle[2] + 
  scale.be*(1 - s^(-fit.freq.be$mle[1]))/fit.freq.be$mle[1]
mle.be <- c(fit.freq.be$mle[1], loc.be, scale.be)
names(mle.be) <- c("Shape", "Location", "Scale")
Q.ext.be <- mle.be[2] + mle.be[3]*((s*p)^(-mle.be[1]) - 1)/mle.be[1]
names(Q.ext.be) <- "Extreme-Quantile"
R.lev.be <- mle.be[2] + 
  mle.be[3]*((-log(1 - 1/T.ret))^(-mle.be[1]) - 1)/mle.be[1]
names(R.lev.be) <- "Return-Level"
#Comparison
#Parameters
mle.be; true.par.be
#Extreme quantile
Q.ext.be; Q.ext.true.be
#Return level
R.lev.be; R.lev.true.be

#Bayesian inference
tic()
fit.bayes.be <- fit.gev.inference(data = data.be, t = th.be, t.prob = t.prob,
                                  llik.type = "Gev-Cens", param.type = "bar",
                                  T.ret = T.ret, p = p, par0 = start.be,
                                  hessian = TRUE, inf.type = "Bayes",
                                  k = 1, R = R, burn  = burn, prior = "empirical",
                                  val.show = TRUE)
toc()
#131.47 sec elapsed
#Results
fit.bayes.be$mcmc.acc
#Tilde parameterization
mode.tilde.be <- fit.bayes.be$mode
mode.tilde.be[3] <- fit.bayes.be$mode[3]*s^(fit.bayes.be$mode[1])
mode.tilde.be[2] <- mode.tilde.be[2] + 
  mode.tilde.be[3]*(1 - s^(-mode.tilde.be[1]))/mode.tilde.be[1]
params.tilde.be <- fit.bayes.be$parameters
scale.post.be <- apply(fit.bayes.be$parameters, 1, function(x) x[3]*s^(x[1]))
params.tilde.be[,3] <- scale.post.be
loc.post.be <- apply(params.tilde.be, 1,
                     function(x) x[2] + x[3]*(1 - s^(-x[1]))/x[1])
params.tilde.be[,2] <- loc.post.be
Q.ext.post.be <- params.tilde.be[,2] + 
  params.tilde.be[,3]*((s*p)^(-params.tilde.be[,1]) - 1)/params.tilde.be[,1]
R.lev.post.be <- params.tilde.be[,2] +
  params.tilde.be[,3]*((-log(1 - 1/T.ret))^(-params.tilde.be[,1]) - 1)/params.tilde.be[,1]
#Comparison
#Parameters
mle.be; mode.tilde.be; true.par.be
sapply(1:3, function(x) summary(params.tilde.be[,x]))
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(params.tilde.be[,i]), params.tilde.be[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.be[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(params.tilde.be[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(params.tilde.be[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.be[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(params.tilde.be[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))

#Extreme quantile
Q.ext.true.be
summary(Q.ext.post.be)
par(mfrow = c(3,1))
#Trace plot
plot(1:length(Q.ext.post.be), Q.ext.post.be, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.be, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(Q.ext.post.be, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(Q.ext.post.be, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.be, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(Q.ext.post.be, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#legend("topright", legend = c("True value", "95% hpd interval"),
#       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2))

#Return level
R.lev.true.be
summary(R.lev.post.be)
#Trace plot
plot(1:length(R.lev.post.be), R.lev.post.be, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.be, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(R.lev.post.be, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(R.lev.post.be, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.be, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(R.lev.post.be, conf = 0.95), col = "red",
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

save.image("models_bar.RData")

#-------------------------------------------------------------------------------