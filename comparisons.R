#Comparisons between functions

#-------------------------------------------------------------------------------

path <- "C:\\Users\\scanz\\Desktop\\Magistrale\\Extreme Value Theory\\Code\\"
setwd(path)

#Pietro Scanzi
source(paste0(path ,"functions.R"))
#Simone Padoan
source(paste0(path, "cens_loglik_SP.R"))

library(LaplacesDemon)
library(TeachingDemos)
library(evd)
#library(extremefit)
library(tictoc)
library(data.table)

#-------------------------------------------------------------------------------

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
#Posterior sample size
R <- 50000
#Burn-in period
burn <- 10000

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
set.seed(8)
data.fre <- rfrechet(n, loc.fre, scale.fre, shape.fre)
th.fre <- quantile(data.fre, probs = 1 - t.prob, type = 3)

#-------------------------------------------------------------------------------

#Frequentist inference

#PS
fit.freq.fre.ps <- fit.gev.inference(data = data.fre, t = th.fre, t.prob = t.prob,
                                  llik.type = "Gev-Cens", T.ret = T.ret, p = p,
                                  par0 = start.fre, hessian = TRUE,
                                  inf.type = "Frequent")
#SP
fit.freq.fre.sp <- FitCensLLik(data = data.fre, th = th.fre, param0 = start.fre, 
                               llik.type = "Gev-Cens", p = p, T = T.ret, 
                               t.prob = t.prob, method = "frequentist")

#Comparison
fit.freq.fre.ps$mle; fit.freq.fre.sp$estimate; true.par.fre
fit.freq.fre.ps$Q.extreme; fit.freq.fre.sp$Q.extreme; Q.ext.true.fre
fit.freq.fre.ps$R.level; fit.freq.fre.sp$R.level; R.lev.true.fre
#
par.ps <- fit.freq.fre.ps$mle; par.sp <- fit.freq.fre.sp$estimate
Q.ext.ps <- par.ps[2] + par.ps[3] *(50^par.ps[1] - 1)/par.ps[1]; Q.ext.ps
Q.ext.sp <- par.sp[1] + par.sp[2] *(50^par.sp[3] - 1)/par.sp[3]; Q.ext.sp
R.lev.ps <- par.ps[2] + par.ps[3] *((-log(1 - 1/T.ret))^(-par.ps[1]) - 1)/par.ps[1]; R.lev.ps
R.lev.sp <- par.sp[1] + par.sp[2] *((-log(1 - 1/T.ret))^(-par.sp[3]) - 1)/par.sp[3]; R.lev.sp
#

#-------------------------------------------------------------------------------

#Bayesian inference
set.seed(3)
#PS
fit.bayes.fre.ps <- fit.gev.inference(data = data.fre, t = th.fre, t.prob = t.prob,
                                   llik.type = "Gev-Cens", T.ret = T.ret, p = p,
                                   par0 = start.fre, hessian = TRUE, inf.type = "Bayes",
                                   k = 1, R = R, burn  = burn, prior = "empirical",
                                   val.show = TRUE)

#SP
fit.bayes.fre.sp <- FitCensLLik(data = data.fre, th = th.fre, param0 = start.fre, 
                                llik.type = "Gev-Cens", p = p, T = T.ret, 
                                t.prob = t.prob, method = "bayesian", sig0 = 1,
                                nsim = R, burn = burn, prior = "empirical",
                                val.show = TRUE)

#Comparison
parb.ps <- fit.bayes.fre.ps$parameters; parb.sp <- fit.bayes.fre.sp$post_sample
parb.sp <- cbind(parb.sp[,3], parb.sp[,1:2])
qextb.ps <- fit.bayes.fre.ps$Q.extreme; qextb.sp <- fit.bayes.fre.sp$Q.extreme
rlevb.ps <- fit.bayes.fre.ps$R.level; rlevb.sp <- fit.bayes.fre.sp$R.level
#
sapply(1:3, function(x) summary(parb.ps[,x]))
sapply(1:3, function(x) summary(parb.sp[,x]))
colMeans(parb.ps); colMeans(parb.sp); true.par.fre
#
summary(qextb.ps); summary(qextb.sp)
mean(qextb.ps); mean(qextb.sp); Q.ext.true.fre
#
summary(rlevb.ps); summary(rlevb.sp)
mean(rlevb.ps); mean(rlevb.sp); Q.ext.true.fre
#
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(parb.ps[,i]), parb.ps[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.fre[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(parb.ps[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(parb.ps[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.fre[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(parb.ps[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))
#
lab <- c(expression(gamma), expression(mu), expression(delta))
par(mfrow = c(3,3))
for(i in 1:3){
  #Trace plot
  plot(1:length(parb.sp[,i]), parb.sp[,i],
       type = "l", xlab = lab[i], ylab = "MCMC chain", main = "Trace plot")
  abline(h = true.par.fre[i], col = "green", lty = "dashed", lwd = 3)
  #ACF plot
  acf(parb.sp[,i], lag.max = 100,
      xlab = lab[i], ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
  #Histogram
  hist(parb.sp[,i], nclass = 50, probability = TRUE,
       xlab = lab[i], main = "Posterior distribution")
  abline(v = true.par.fre[i], col = "blue", lty = 2, lwd = 2)
  abline(v = emp.hpd(parb.sp[,i], conf = 0.95), col = "red",
         lty = 3, lwd = 2)
  #legend("topright", legend = c("True value", "95% hpd interval"),
  #       col = c("blue", "red"), lty = c(2,3), lwd = rep(2,2), cex = 0.75)
}
par(mfrow = c(1,1))
#
par(mfrow = c(3,1))
#Trace plot
plot(1:length(qextb.ps), qextb.ps, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.fre, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(qextb.ps, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(qextb.ps, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.fre, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(qextb.ps, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#
plot(1:length(qextb.sp), qextb.sp, type = "l",
     xlab = "Extreme quantile", ylab = "MCMC chain", main = "Trace plot")
abline(h = Q.ext.true.fre, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(qextb.sp, lag.max = 100, xlab = "Extreme quantile",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(qextb.sp, nclass = 50, probability = TRUE,
     xlab = "Extreme quantile", main = "Posterior distribution")
abline(v = Q.ext.true.fre, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(qextb.sp, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#
#Trace plot
plot(1:length(rlevb.ps), rlevb.ps, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.fre, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(rlevb.ps, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(rlevb.ps, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.fre, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(rlevb.ps, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
#
#Trace plot
plot(1:length(rlevb.sp), rlevb.sp, type = "l",
     xlab = "Return level", ylab = "MCMC chain", main = "Trace plot")
abline(h = R.lev.true.fre, col = "green", lty = "dashed", lwd = 3)
#ACF plot
acf(rlevb.sp, lag.max = 100, xlab = "Return level",
    ylab = "ACF", main = "Autocorrelation plot (lag.max = 100)")
#Histogram
hist(rlevb.sp, nclass = 50, probability = TRUE,
     xlab = "Return level", main = "Posterior distribution")
abline(v = R.lev.true.fre, col = "blue", lty = 2, lwd = 2)
abline(v = emp.hpd(rlevb.sp, conf = 0.95), col = "red",
       lty = 3, lwd = 2)
par(mfrow = c(1,1))

#-------------------------------------------------------------------------------

#Simulation study

#Sample size
n <- c(800, 1800, 5450, 23400)
#Number of exceedances used in the estimation
k <- c(20, 30, 50, 100)
#Number of possible simulation configurations
nconf <- length(n)
#Exceeding probability used in the estimation
t.prob <- k/n
#Exceeding probability used for prediction
p <- 1/n
#Return period
T.ret <- 50
#Posterior sample size
R <- 10000
#Burn-in period
burn <- 2500
#Number of simulations
Nsim <- 100

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
#Define the estimation setting
start.fre <- c(0.1, 0.1, 1)

#Return objects
lab.par <- c("gamma", "mu", "delta")
lab.ci <- c("gamma.low", "gamma.up", "mu.low", "mu.up", "delta.low", "delta.up",
            "Q.ext.low", "Q.ext.up", "R.lev.low", "R.lev.up")
lab.cov <- c("gamma", "mu", "delta", "Q.ext", "R.lev")
true.par.fre <- matrix(0, nrow = nconf, ncol = 3)
colnames(true.par.fre) <- lab.par
Q.ext.true.fre <- rep(0, nconf)
R.lev.true.fre <- rep(0, nconf)
#
quant.ci.fre.ps <- array(0, dim = c(Nsim, 10, nconf))
colnames(quant.ci.fre.ps) <- lab.ci
quant.cov.fre.ps <- array(0, dim = c(Nsim, 5, nconf))
colnames(quant.cov.fre.ps) <- lab.cov
wald.ci.fre.ps <- array(0, dim = c(Nsim, 10, nconf))
colnames(wald.ci.fre.ps) <- lab.ci
wald.cov.fre.ps <- array(0, dim = c(Nsim, 5, nconf))
colnames(wald.cov.fre.ps) <- lab.cov
#
quant.ci.fre.sp <- array(0, dim = c(Nsim, 10, nconf))
colnames(quant.ci.fre.sp) <- lab.ci
quant.cov.fre.sp <- array(0, dim = c(Nsim, 5, nconf))
colnames(quant.cov.fre.sp) <- lab.cov
wald.ci.fre.sp <- array(0, dim = c(Nsim, 10, nconf))
colnames(wald.ci.fre.sp) <- lab.ci
wald.cov.fre.sp <- array(0, dim = c(Nsim, 5, nconf))
colnames(wald.cov.fre.sp) <- lab.cov

#Simulation loop
set.seed(1)
tic()
for(s in 1:nconf){
  
  cat("\n", "Starting with n =", n[s], "and k =", k[s], "\n")
  
  #Define the true parameters
  true.par.fre[s,] <- c(1/shape.fre, bn.fre[s], an.fre[s])
  #Define true extreme quantile
  Q.ext.true.fre[s] <- qfrechet(1 - p[s], loc.fre, scale.fre, shape.fre)
  #Define true return level
  R.lev.true.fre[s] <- qgev(1 - 1/T.ret, bn.fre[s], an.fre[s], 1/shape.fre)
  
  #Single configuration loop
  for(i in 1:Nsim){
    
    #Sample
    data.fre <- rfrechet(n[s], loc.fre, scale.fre, shape.fre)
    #Threshold
    th.fre <- quantile(data.fre, probs = 1 - t.prob[s], type = 3)
    
    fit.bayes.fre.ps <- fit.gev.inference(data = data.fre, t = th.fre,
                                          t.prob = t.prob[s], llik.type = "Gev-Cens",
                                          T.ret = T.ret, p = p[s],  par0 = start.fre,
                                          hessian = TRUE, inf.type = "Bayes",
                                          k = 1, R = R, burn  = burn, prior = "empirical")
    fit.bayes.fre.sp <- FitCensLLik(data = data.fre, th = th.fre, param0 = start.fre, 
                                    llik.type = "Gev-Cens", p = p[s], T = T.ret, 
                                    t.prob = t.prob[s], method = "bayesian", sig0 = 1,
                                    nsim = R, burn = burn, prior = "empirical",
                                    val.show = FALSE)
    #
    params.ps <- fit.bayes.fre.ps$parameters
    params.sp <- fit.bayes.fre.sp$post_sample
    params.sp <- cbind(params.sp[,3], params.sp[,1:2])
    Q.ext.ps <- fit.bayes.fre.ps$Q.extreme
    Q.ext.sp <- fit.bayes.fre.sp$Q.extreme
    R.lev.ps <- fit.bayes.fre.ps$R.level
    R.lev.sp <- fit.bayes.fre.sp$R.level
    #
    quant.ci.fre.ps[i,1:2,s] <- as.numeric(quantile(params.ps[,1], c(0.025, 0.975)))
    quant.cov.fre.ps[i,1,s] <- as.numeric(true.par.fre[s,1] %between% quant.ci.fre.ps[i,1:2,s])
    wald.ci.fre.ps[i,1:2,s] <- mean(params.ps[,1]) + c(-1,1)*qnorm(0.975)*sd(params.ps[,1])
    wald.cov.fre.ps[i,1,s] <- as.numeric(true.par.fre[s,1] %between% wald.ci.fre.ps[i,1:2,s])
    quant.ci.fre.ps[i,3:4,s] <- as.numeric(quantile(params.ps[,2], c(0.025, 0.975)))
    quant.cov.fre.ps[i,2,s] <- as.numeric(true.par.fre[s,2] %between% quant.ci.fre.ps[i,3:4,s])
    wald.ci.fre.ps[i,3:4,s] <- mean(params.ps[,2]) + c(-1,1)*qnorm(0.975)*sd(params.ps[,2])
    wald.cov.fre.ps[i,2,s] <- as.numeric(true.par.fre[s,2] %between% wald.ci.fre.ps[i,3:4,s])
    quant.ci.fre.ps[i,5:6,s] <- as.numeric(quantile(params.ps[,3], c(0.025, 0.975)))
    quant.cov.fre.ps[i,3,s] <- as.numeric(true.par.fre[s,3] %between% quant.ci.fre.ps[i,5:6,s])
    wald.ci.fre.ps[i,5:6,s] <- mean(params.ps[,3]) + c(-1,1)*qnorm(0.975)*sd(params.ps[,3])
    wald.cov.fre.ps[i,3,s] <- as.numeric(true.par.fre[s,3] %between% wald.ci.fre.ps[i,5:6,s])
    #
    quant.ci.fre.ps[i,7:8,s] <- as.numeric(quantile(Q.ext.ps, c(0.025, 0.975)))
    quant.cov.fre.ps[i,4,s] <- as.numeric(Q.ext.true.fre[s] %between% quant.ci.fre.ps[i,7:8,s])
    wald.ci.fre.ps[i,7:8,s] <- mean(Q.ext.ps) + c(-1,1)*qnorm(0.975)*sd(Q.ext.ps)
    wald.cov.fre.ps[i,4,s] <- as.numeric(Q.ext.true.fre[s] %between% wald.ci.fre.ps[i,7:8,s])
    #
    quant.ci.fre.ps[i,9:10,s] <- as.numeric(quantile(R.lev.ps, c(0.025, 0.975)))
    quant.cov.fre.ps[i,5,s] <- as.numeric(R.lev.true.fre[s] %between% quant.ci.fre.ps[i,9:10,s])
    wald.ci.fre.ps[i,9:10,s] <- mean(R.lev.ps) + c(-1,1)*qnorm(0.975)*sd(R.lev.ps)
    wald.cov.fre.ps[i,5,s] <- as.numeric(R.lev.true.fre[s] %between% wald.ci.fre.ps[i,9:10,s])
    
    #
    quant.ci.fre.sp[i,1:2,s] <- as.numeric(quantile(params.sp[,1], c(0.025, 0.975)))
    quant.cov.fre.sp[i,1,s] <- as.numeric(true.par.fre[s,1] %between% quant.ci.fre.sp[i,1:2,s])
    wald.ci.fre.sp[i,1:2,s] <- mean(params.sp[,1]) + c(-1,1)*qnorm(0.975)*sd(params.sp[,1])
    wald.cov.fre.sp[i,1,s] <- as.numeric(true.par.fre[s,1] %between% wald.ci.fre.sp[i,1:2,s])
    quant.ci.fre.sp[i,3:4,s] <- as.numeric(quantile(params.sp[,2], c(0.025, 0.975)))
    quant.cov.fre.sp[i,2,s] <- as.numeric(true.par.fre[s,2] %between% quant.ci.fre.sp[i,3:4,s])
    wald.ci.fre.sp[i,3:4,s] <- mean(params.sp[,2]) + c(-1,1)*qnorm(0.975)*sd(params.sp[,2])
    wald.cov.fre.sp[i,2,s] <- as.numeric(true.par.fre[s,2] %between% wald.ci.fre.sp[i,3:4,s])
    quant.ci.fre.sp[i,5:6,s] <- as.numeric(quantile(params.sp[,3], c(0.025, 0.975)))
    quant.cov.fre.sp[i,3,s] <- as.numeric(true.par.fre[s,3] %between% quant.ci.fre.sp[i,5:6,s])
    wald.ci.fre.sp[i,5:6,s] <- mean(params.sp[,3]) + c(-1,1)*qnorm(0.975)*sd(params.sp[,3])
    wald.cov.fre.sp[i,3,s] <- as.numeric(true.par.fre[s,3] %between% wald.ci.fre.sp[i,5:6,s])
    #
    quant.ci.fre.sp[i,7:8,s] <- as.numeric(quantile(Q.ext.sp, c(0.025, 0.975)))
    quant.cov.fre.sp[i,4,s] <- as.numeric(Q.ext.true.fre[s] %between% quant.ci.fre.sp[i,7:8,s])
    wald.ci.fre.sp[i,7:8,s] <- mean(Q.ext.sp) + c(-1,1)*qnorm(0.975)*sd(Q.ext.sp)
    wald.cov.fre.sp[i,4,s] <- as.numeric(Q.ext.true.fre[s] %between% wald.ci.fre.sp[i,7:8,s])
    #
    quant.ci.fre.sp[i,9:10,s] <- as.numeric(quantile(R.lev.sp, c(0.025, 0.975)))
    quant.cov.fre.sp[i,5,s] <- as.numeric(R.lev.true.fre[s] %between% quant.ci.fre.sp[i,9:10,s])
    wald.ci.fre.sp[i,9:10,s] <- mean(R.lev.sp) + c(-1,1)*qnorm(0.975)*sd(R.lev.sp)
    wald.cov.fre.sp[i,5,s] <- as.numeric(R.lev.true.fre[s] %between% wald.ci.fre.sp[i,9:10,s])
    
    #Loop indicator 
    cat(i, " ")
  }
}
toc()

#Empirical coverage
lab.sim <- c("n = 800", "n = 1800", "n = 5450", "n = 23400")
quant.emp.cov.fre.ps <- colMeans(quant.cov.fre.ps)
colnames(quant.emp.cov.fre.ps) <- lab.sim
wald.emp.cov.fre.ps <- colMeans(wald.cov.fre.ps)
colnames(wald.emp.cov.fre.ps) <- lab.sim
#
quant.emp.cov.fre.sp <- colMeans(quant.cov.fre.sp)
colnames(quant.emp.cov.fre.sp) <- lab.sim
wald.emp.cov.fre.sp <- colMeans(wald.cov.fre.sp)
colnames(wald.emp.cov.fre.sp) <- lab.sim

#Results
quant.emp.cov.fre.ps
quant.emp.cov.fre.sp
wald.emp.cov.fre.ps
wald.emp.cov.fre.sp
#-------------------------------------------------------------------------------