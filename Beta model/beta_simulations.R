#Beta posterior simulations

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

#Beta(1,3) --> tail index = -1/3

#Beta parameters
shape1.be <- 1
shape2.be <- 3
#Set norming constants
#Centering
bn.be <- qbeta(1 -(k/n), shape1 = shape1.be, shape2 = shape2.be)
#Scaling
an.be <- k/(n*dbeta(bn.be, shape1 = shape1.be, shape2 = shape2.be))
#Define the estimation setting
start.be <- c(-0.36, 0.63, 0.13)

#Return objects
lab.par <- c("gamma", "mu", "delta")
lab.ci <- c("gamma.low", "gamma.up", "mu.low", "mu.up", "delta.low", "delta.up",
            "Q.ext.low", "Q.ext.up", "R.lev.low", "R.lev.up")
lab.cov <- c("gamma", "mu", "delta", "Q.ext", "R.lev")
true.par.be <- matrix(0, nrow = nconf, ncol = 3)
colnames(true.par.be) <- lab.par
Q.ext.true.be <- rep(0, nconf)
R.lev.true.be <- rep(0, nconf)
quant.ci.be <- array(0, dim = c(Nsim, 10, nconf))
colnames(quant.ci.be) <- lab.ci
quant.cov.be <- array(0, dim = c(Nsim, 5, nconf))
colnames(quant.cov.be) <- lab.cov
wald.ci.be <- array(0, dim = c(Nsim, 10, nconf))
colnames(wald.ci.be) <- lab.ci
wald.cov.be <- array(0, dim = c(Nsim, 5, nconf))
colnames(wald.cov.be) <- lab.cov
hpd.ci.be <- array(0, dim = c(Nsim, 10, nconf))
colnames(hpd.ci.be) <- lab.ci
hpd.cov.be <- array(0, dim = c(Nsim, 5, nconf))
colnames(hpd.cov.be) <- lab.cov

#Simulation loop
set.seed(15)
for(s in 1:nconf){
  
  cat("\n", "Starting with n =", n[s], "and k =", k[s], "\n")
  
  #Define the true parameters
  true.par.be[s,] <- c(-1/shape2.be, bn.be[s], an.be[s])
  #Define true extreme quantile
  Q.ext.true.be[s] <- qbeta(1 - p[s], shape1 = shape1.be, shape2 = shape2.be)
  #Define true return level
  R.lev.true.be[s] <- qgev(1 - 1/T.ret, bn.be[s], an.be[s], -1/shape2.be)
  
  #Single configuration loop
  for(i in 1:Nsim){
    
    #Sample
    data.be <- rbeta(n[s], shape1 = shape1.be, shape2 = shape2.be)
    #Threshold
    th.be <- quantile(data.be, probs = 1 - t.prob[s], type = 3)
    
    fit.bayes.be <- fit.gev.bayes(data = data.be, t = th.be, t.prob = t.prob[s],
                                  llik.type = "Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p[s],
                                  par0 = start.be, k = 1, R = R, burn  = burn,
                                  prior = "empirical", etastar = 0.234)
    
    params <- fit.bayes.be$parameters
    quant.ci.be[i,1:2,s] <- as.numeric(quantile(params[,1], c(0.025, 0.975), na.rm = T))
    quant.cov.be[i,1,s] <- as.numeric(true.par.be[s,1] %between% quant.ci.be[i,1:2,s])
    wald.ci.be[i,1:2,s] <- mean(params[,1], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,1], na.rm = T)
    wald.cov.be[i,1,s] <- as.numeric(true.par.be[s,1] %between% wald.ci.be[i,1:2,s])
    hpd.ci.be[i,1:2,s] <- emp.hpd(as.numeric(na.omit(params[,1])), conf = 0.95)
    hpd.cov.be[i,1,s] <- as.numeric(true.par.be[s,1] %between% hpd.ci.be[i,1:2,s])
    quant.ci.be[i,3:4,s] <- as.numeric(quantile(params[,2], c(0.025, 0.975), na.rm = T))
    quant.cov.be[i,2,s] <- as.numeric(true.par.be[s,2] %between% quant.ci.be[i,3:4,s])
    wald.ci.be[i,3:4,s] <- mean(params[,2], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,2], na.rm = T)
    wald.cov.be[i,2,s] <- as.numeric(true.par.be[s,2] %between% wald.ci.be[i,3:4,s])
    hpd.ci.be[i,3:4,s] <- emp.hpd(as.numeric(na.omit(params[,2])), conf = 0.95)
    hpd.cov.be[i,2,s] <- as.numeric(true.par.be[s,2] %between% hpd.ci.be[i,3:4,s])
    quant.ci.be[i,5:6,s] <- as.numeric(quantile(params[,3], c(0.025, 0.975), na.rm = T))
    quant.cov.be[i,3,s] <- as.numeric(true.par.be[s,3] %between% quant.ci.be[i,5:6,s])
    wald.ci.be[i,5:6,s] <- mean(params[,3], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,3], na.rm = T)
    wald.cov.be[i,3,s] <- as.numeric(true.par.be[s,3] %between% wald.ci.be[i,5:6,s])
    hpd.ci.be[i,5:6,s] <- emp.hpd(as.numeric(na.omit(params[,3])), conf = 0.95)
    hpd.cov.be[i,3,s] <- as.numeric(true.par.be[s,3] %between% hpd.ci.be[i,5:6,s])
    Q.ext <- fit.bayes.be$Q.extreme
    quant.ci.be[i,7:8,s] <- as.numeric(quantile(Q.ext, c(0.025, 0.975), na.rm = T))
    quant.cov.be[i,4,s] <- as.numeric(Q.ext.true.be[s] %between% quant.ci.be[i,7:8,s])
    wald.ci.be[i,7:8,s] <- mean(Q.ext, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(Q.ext, na.rm = T)
    wald.cov.be[i,4,s] <- as.numeric(Q.ext.true.be[s] %between% wald.ci.be[i,7:8,s])
    hpd.ci.be[i,7:8,s] <- emp.hpd(as.numeric(na.omit(Q.ext)), conf = 0.95)
    hpd.cov.be[i,4,s] <- as.numeric(Q.ext.true.be[s] %between% hpd.ci.be[i,7:8,s])
    R.lev <- fit.bayes.be$R.level
    quant.ci.be[i,9:10,s] <- as.numeric(quantile(R.lev, c(0.025, 0.975), na.rm = T))
    quant.cov.be[i,5,s] <- as.numeric(R.lev.true.be[s] %between% quant.ci.be[i,9:10,s])
    wald.ci.be[i,9:10,s] <- mean(R.lev, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(R.lev, na.rm = T)
    wald.cov.be[i,5,s] <- as.numeric(R.lev.true.be[s] %between% wald.ci.be[i,9:10,s])
    hpd.ci.be[i,9:10,s] <- emp.hpd(as.numeric(na.omit(R.lev)), conf = 0.95)
    hpd.cov.be[i,5,s] <- as.numeric(R.lev.true.be[s] %between% hpd.ci.be[i,9:10,s])
    
    #Loop indicator 
    cat(i, " ")
  }
}

#Empirical coverage
lab.sim <- c("n = 800", "n = 1800", "n = 5450", "n = 23400")
quant.emp.cov.be <- colMeans(quant.cov.be)
colnames(quant.emp.cov.be) <- lab.sim
wald.emp.cov.be <- colMeans(wald.cov.be)
colnames(wald.emp.cov.be) <- lab.sim
hpd.emp.cov.be <- colMeans(hpd.cov.be)
colnames(hpd.emp.cov.be) <- lab.sim

#Save
save(true.par.be, Q.ext.true.be, R.lev.true.be, quant.ci.be, quant.cov.be,
     wald.ci.be, wald.cov.be, hpd.ci.be, hpd.cov.be, quant.emp.cov.be,
     wald.emp.cov.be, hpd.emp.cov.be,
     file = "beta_simulations.RData")

#-------------------------------------------------------------------------------