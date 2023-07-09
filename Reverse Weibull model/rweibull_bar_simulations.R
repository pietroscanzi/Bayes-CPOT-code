#Reverse Weibull posterior simulations (bar parameterization)

#-------------------------------------------------------------------------------

path <- "C:\\Users\\scanz\\Desktop\\Magistrale\\Extreme Value Theory\\Code\\"
setwd(path)
source(paste0(path ,"functions_simulations.R"))

#-------------------------------------------------------------------------------

#Simulation setting
library(LaplacesDemon)
library(TeachingDemos)
library(evd)
#library(extremefit)
library(tictoc)
library(data.table)

#Sample size
n <- c(800, 1800, 5450, 23400)
#Number of exceedances used in the estimation
k <- c(20, 30, 50, 100)
#Proportion of exceedances
s.p <- n/k
#Number of possible simulation configurations
nconf <- length(n)
#Exceeding probability used in the estimation
t.prob <- k/n
#Exceeding probability used for prediction
p <- 1/n
#Return period
T.ret <- 50
#Posterior sample size
R <- 50000
#Burn-in period
burn <- 10000
#Number of simulations
Nsim <- 1000

#-------------------------------------------------------------------------------

#Reverse Weibull (shape = 3) --> tail index = -1/3

#Reverse Weibull parameters
shape.rw <- 3 
#Set norming constants
#Centering
bn.rw <- -(-log(1 - k/n))^(1/shape.rw)
#Scaling
#h <- 1e-10
an.rw <- (k/(shape.rw*(n - k)))*(-log(1 - k/n))^(1/shape.rw - 1)
#an.rw2 <- 
#Define the estimation setting
start.rw <- c(-0.40, -1.07, 0.41)

#Return objects
lab.par <- c("gamma", "mu", "delta")
lab.ci <- c("gamma.low", "gamma.up", "mu.low", "mu.up", "delta.low", "delta.up",
            "Q.ext.low", "Q.ext.up", "R.lev.low", "R.lev.up")
lab.cov <- c("gamma", "mu", "delta", "Q.ext", "R.lev")
true.par.rw <- matrix(0, nrow = nconf, ncol = 3)
colnames(true.par.rw) <- lab.par
Q.ext.true.rw <- rep(0, nconf)
R.lev.true.rw <- rep(0, nconf)
quant.ci.rw <- array(0, dim = c(Nsim, 10, nconf))
colnames(quant.ci.rw) <- lab.ci
quant.cov.rw <- array(0, dim = c(Nsim, 5, nconf))
colnames(quant.cov.rw) <- lab.cov
wald.ci.rw <- array(0, dim = c(Nsim, 10, nconf))
colnames(wald.ci.rw) <- lab.ci
wald.cov.rw <- array(0, dim = c(Nsim, 5, nconf))
colnames(wald.cov.rw) <- lab.cov
hpd.ci.rw <- array(0, dim = c(Nsim, 10, nconf))
colnames(hpd.ci.rw) <- lab.ci
hpd.cov.rw <- array(0, dim = c(Nsim, 5, nconf))
colnames(hpd.cov.rw) <- lab.cov

#Simulation loop
set.seed(50)
for(s in 1:nconf){
  
  cat("\n", "Starting with n =", n[s], "and k =", k[s], "\n")
  
  #Define the true parameters
  true.par.rw[s,] <- c(-1/shape.rw, bn.rw[s], an.rw[s])
  #Define true extreme quantile
  Q.ext.true.rw[s] <- qrweibull(1 - p[s], shape = shape.rw)
  #Define true return level
  R.lev.true.rw[s] <- qgev(1 - 1/T.ret, bn.rw[s], an.rw[s], -1/shape.rw)
  
  #Single configuration loop
  for(i in 1:Nsim){
    
    #Sample
    data.rw <- rrweibull(n[s], shape = shape.rw)
    #Threshold
    th.rw <- quantile(data.rw, probs = 1 - t.prob[s], type = 3)
    
    fit.bayes.rw <- fit.gev.bayes(data = data.rw, t = th.rw, t.prob = t.prob[s],
                                  llik.type = "Gev-Cens", param.type = "bar",
                                  T.ret = T.ret, p = p[s],
                                  par0 = start.rw, k = 1, R = R, burn  = burn,
                                  prior = "empirical", etastar = 0.234)
    
    params.bar <- fit.bayes.rw$parameters
    params <- params.bar
    params[,3] <- apply(params.bar, 1, function(x) x[3]*s.p[s]^(x[1]))
    params[,2] <- apply(params, 1, function(x) x[2] + x[3]*(1 - s.p[s]^(-x[1]))/x[1])
    quant.ci.rw[i,1:2,s] <- as.numeric(quantile(params[,1], c(0.025, 0.975), na.rm = T))
    quant.cov.rw[i,1,s] <- as.numeric(true.par.rw[s,1] %between% quant.ci.rw[i,1:2,s])
    wald.ci.rw[i,1:2,s] <- mean(params[,1], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,1], na.rm = T)
    wald.cov.rw[i,1,s] <- as.numeric(true.par.rw[s,1] %between% wald.ci.rw[i,1:2,s])
    hpd.ci.rw[i,1:2,s] <- emp.hpd(as.numeric(na.omit(params[,1])), conf = 0.95)
    hpd.cov.rw[i,1,s] <- as.numeric(true.par.rw[s,1] %between% hpd.ci.rw[i,1:2,s])
    quant.ci.rw[i,3:4,s] <- as.numeric(quantile(params[,2], c(0.025, 0.975), na.rm = T))
    quant.cov.rw[i,2,s] <- as.numeric(true.par.rw[s,2] %between% quant.ci.rw[i,3:4,s])
    wald.ci.rw[i,3:4,s] <- mean(params[,2], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,2], na.rm = T)
    wald.cov.rw[i,2,s] <- as.numeric(true.par.rw[s,2] %between% wald.ci.rw[i,3:4,s])
    hpd.ci.rw[i,3:4,s] <- emp.hpd(as.numeric(na.omit(params[,2])), conf = 0.95)
    hpd.cov.rw[i,2,s] <- as.numeric(true.par.rw[s,2] %between% hpd.ci.rw[i,3:4,s])
    quant.ci.rw[i,5:6,s] <- as.numeric(quantile(params[,3], c(0.025, 0.975), na.rm = T))
    quant.cov.rw[i,3,s] <- as.numeric(true.par.rw[s,3] %between% quant.ci.rw[i,5:6,s])
    wald.ci.rw[i,5:6,s] <- mean(params[,3], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,3], na.rm = T)
    wald.cov.rw[i,3,s] <- as.numeric(true.par.rw[s,3] %between% wald.ci.rw[i,5:6,s])
    hpd.ci.rw[i,5:6,s] <- emp.hpd(as.numeric(na.omit(params[,3])), conf = 0.95)
    hpd.cov.rw[i,3,s] <- as.numeric(true.par.rw[s,3] %between% hpd.ci.rw[i,5:6,s])
    
    Q.ext <- params[,2] + params[,3]*((s.p[s]*p[s])^(-params[,1]) - 1)/params[,1]
    quant.ci.rw[i,7:8,s] <- as.numeric(quantile(Q.ext, c(0.025, 0.975), na.rm = T))
    quant.cov.rw[i,4,s] <- as.numeric(Q.ext.true.rw[s] %between% quant.ci.rw[i,7:8,s])
    wald.ci.rw[i,7:8,s] <- mean(Q.ext, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(Q.ext, na.rm = T)
    wald.cov.rw[i,4,s] <- as.numeric(Q.ext.true.rw[s] %between% wald.ci.rw[i,7:8,s])
    hpd.ci.rw[i,7:8,s] <- emp.hpd(as.numeric(na.omit(Q.ext)), conf = 0.95)
    hpd.cov.rw[i,4,s] <- as.numeric(Q.ext.true.rw[s] %between% hpd.ci.rw[i,7:8,s])
    
    R.lev <- params[,2] + params[,3]*((-log(1 - 1/T.ret))^(-params[,1]) - 1)/params[,1]
    quant.ci.rw[i,9:10,s] <- as.numeric(quantile(R.lev, c(0.025, 0.975), na.rm = T))
    quant.cov.rw[i,5,s] <- as.numeric(R.lev.true.rw[s] %between% quant.ci.rw[i,9:10,s])
    wald.ci.rw[i,9:10,s] <- mean(R.lev, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(R.lev, na.rm = T)
    wald.cov.rw[i,5,s] <- as.numeric(R.lev.true.rw[s] %between% wald.ci.rw[i,9:10,s])
    hpd.ci.rw[i,9:10,s] <- emp.hpd(as.numeric(na.omit(R.lev)), conf = 0.95)
    hpd.cov.rw[i,5,s] <- as.numeric(R.lev.true.rw[s] %between% hpd.ci.rw[i,9:10,s])
    
    #Loop indicator 
    cat(i, " ")
  }
}

#Empirical coverage
lab.sim <- c("n = 800", "n = 1800", "n = 5450", "n = 23400")
quant.emp.cov.rw <- colMeans(quant.cov.rw)
colnames(quant.emp.cov.rw) <- lab.sim
wald.emp.cov.rw <- colMeans(wald.cov.rw)
colnames(wald.emp.cov.rw) <- lab.sim
hpd.emp.cov.rw <- colMeans(hpd.cov.rw)
colnames(hpd.emp.cov.rw) <- lab.sim

#Save
save(true.par.rw, Q.ext.true.rw, R.lev.true.rw, quant.ci.rw, quant.cov.rw,
     wald.ci.rw, wald.cov.rw, hpd.ci.rw, hpd.cov.rw, quant.emp.cov.rw,
     wald.emp.cov.rw, hpd.emp.cov.rw,
     file = "rweibull_bar_simulations.RData")

#-------------------------------------------------------------------------------