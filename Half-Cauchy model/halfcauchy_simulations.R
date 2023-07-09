#Half-Cauchy posterior simulations

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
#Define the estimation setting
start.hc <- c(0.81, 13.85, 13.43)

#Return objects
lab.par <- c("gamma", "mu", "delta")
lab.ci <- c("gamma.low", "gamma.up", "mu.low", "mu.up", "delta.low", "delta.up",
            "Q.ext.low", "Q.ext.up", "R.lev.low", "R.lev.up")
lab.cov <- c("gamma", "mu", "delta", "Q.ext", "R.lev")
true.par.hc <- matrix(0, nrow = nconf, ncol = 3)
colnames(true.par.hc) <- lab.par
Q.ext.true.hc <- rep(0, nconf)
R.lev.true.hc <- rep(0, nconf)
quant.ci.hc <- array(0, dim = c(Nsim, 10, nconf))
colnames(quant.ci.hc) <- lab.ci
quant.cov.hc <- array(0, dim = c(Nsim, 5, nconf))
colnames(quant.cov.hc) <- lab.cov
wald.ci.hc <- array(0, dim = c(Nsim, 10, nconf))
colnames(wald.ci.hc) <- lab.ci
wald.cov.hc <- array(0, dim = c(Nsim, 5, nconf))
colnames(wald.cov.hc) <- lab.cov
hpd.ci.hc <- array(0, dim = c(Nsim, 10, nconf))
colnames(hpd.ci.hc) <- lab.ci
hpd.cov.hc <- array(0, dim = c(Nsim, 5, nconf))
colnames(hpd.cov.hc) <- lab.cov

#Simulation loop
set.seed(3)
for(s in 1:nconf){
  
  cat("\n", "Starting with n =", n[s], "and k =", k[s], "\n")
  
  #Define the true parameters
  true.par.hc[s,] <- c(1/scale.hc, bn.hc[s], an.hc[s])
  #Define true extreme quantile
  Q.ext.true.hc[s] <- qhalfcauchy(1 - p[s], scale.hc)
  #Define true return level
  R.lev.true.hc[s] <- qgev(1 - 1/T.ret, bn.hc[s], an.hc[s], 1/scale.hc)
  
  #Single configuration loop
  for(i in 1:Nsim){
    
    #Sample
    data.hc <- rhalfcauchy(n[s], scale.hc)
    #Threshold
    th.hc <- quantile(data.hc, probs = 1 - t.prob[s], type = 3)
    
    fit.bayes.hc <- fit.gev.bayes(data = data.hc, t = th.hc, t.prob = t.prob[s],
                                   llik.type = "Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p[s],
                                   par0 = start.hc, k = 1, R = R, burn  = burn,
                                   prior = "empirical", etastar = 0.234)
    
    params <- fit.bayes.hc$parameters
    quant.ci.hc[i,1:2,s] <- as.numeric(quantile(params[,1], c(0.025, 0.975), na.rm = T))
    quant.cov.hc[i,1,s] <- as.numeric(true.par.hc[s,1] %between% quant.ci.hc[i,1:2,s])
    wald.ci.hc[i,1:2,s] <- mean(params[,1], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,1], na.rm = T)
    wald.cov.hc[i,1,s] <- as.numeric(true.par.hc[s,1] %between% wald.ci.hc[i,1:2,s])
    hpd.ci.hc[i,1:2,s] <- emp.hpd(as.numeric(na.omit(params[,1])), conf = 0.95)
    hpd.cov.hc[i,1,s] <- as.numeric(true.par.hc[s,1] %between% hpd.ci.hc[i,1:2,s])
    quant.ci.hc[i,3:4,s] <- as.numeric(quantile(params[,2], c(0.025, 0.975), na.rm = T))
    quant.cov.hc[i,2,s] <- as.numeric(true.par.hc[s,2] %between% quant.ci.hc[i,3:4,s])
    wald.ci.hc[i,3:4,s] <- mean(params[,2], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,2], na.rm = T)
    wald.cov.hc[i,2,s] <- as.numeric(true.par.hc[s,2] %between% wald.ci.hc[i,3:4,s])
    hpd.ci.hc[i,3:4,s] <- emp.hpd(as.numeric(na.omit(params[,2])), conf = 0.95)
    hpd.cov.hc[i,2,s] <- as.numeric(true.par.hc[s,2] %between% hpd.ci.hc[i,3:4,s])
    quant.ci.hc[i,5:6,s] <- as.numeric(quantile(params[,3], c(0.025, 0.975), na.rm = T))
    quant.cov.hc[i,3,s] <- as.numeric(true.par.hc[s,3] %between% quant.ci.hc[i,5:6,s])
    wald.ci.hc[i,5:6,s] <- mean(params[,3], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,3], na.rm = T)
    wald.cov.hc[i,3,s] <- as.numeric(true.par.hc[s,3] %between% wald.ci.hc[i,5:6,s])
    hpd.ci.hc[i,5:6,s] <- emp.hpd(as.numeric(na.omit(params[,3])), conf = 0.95)
    hpd.cov.hc[i,3,s] <- as.numeric(true.par.hc[s,3] %between% hpd.ci.hc[i,5:6,s])
    Q.ext <- fit.bayes.hc$Q.extreme
    quant.ci.hc[i,7:8,s] <- as.numeric(quantile(Q.ext, c(0.025, 0.975), na.rm = T))
    quant.cov.hc[i,4,s] <- as.numeric(Q.ext.true.hc[s] %between% quant.ci.hc[i,7:8,s])
    wald.ci.hc[i,7:8,s] <- mean(Q.ext, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(Q.ext, na.rm = T)
    wald.cov.hc[i,4,s] <- as.numeric(Q.ext.true.hc[s] %between% wald.ci.hc[i,7:8,s])
    hpd.ci.hc[i,7:8,s] <- emp.hpd(as.numeric(na.omit(Q.ext)), conf = 0.95)
    hpd.cov.hc[i,4,s] <- as.numeric(Q.ext.true.hc[s] %between% hpd.ci.hc[i,7:8,s])
    R.lev <- fit.bayes.hc$R.level
    quant.ci.hc[i,9:10,s] <- as.numeric(quantile(R.lev, c(0.025, 0.975), na.rm = T))
    quant.cov.hc[i,5,s] <- as.numeric(R.lev.true.hc[s] %between% quant.ci.hc[i,9:10,s])
    wald.ci.hc[i,9:10,s] <- mean(R.lev, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(R.lev, na.rm = T)
    wald.cov.hc[i,5,s] <- as.numeric(R.lev.true.hc[s] %between% wald.ci.hc[i,9:10,s])
    hpd.ci.hc[i,9:10,s] <- emp.hpd(as.numeric(na.omit(R.lev)), conf = 0.95)
    hpd.cov.hc[i,5,s] <- as.numeric(R.lev.true.hc[s] %between% hpd.ci.hc[i,9:10,s])
    
    #Loop indicator 
    cat(i, " ")
  }
}

#Empirical coverage
lab.sim <- c("n = 800", "n = 1800", "n = 5450", "n = 23400")
quant.emp.cov.hc <- colMeans(quant.cov.hc)
colnames(quant.emp.cov.hc) <- lab.sim
wald.emp.cov.hc <- colMeans(wald.cov.hc)
colnames(wald.emp.cov.hc) <- lab.sim
hpd.emp.cov.hc <- colMeans(hpd.cov.hc)
colnames(hpd.emp.cov.hc) <- lab.sim

#Save
save(true.par.hc, Q.ext.true.hc, R.lev.true.hc, quant.ci.hc, quant.cov.hc,
     wald.ci.hc, wald.cov.hc, hpd.ci.hc, hpd.cov.hc, quant.emp.cov.hc,
     wald.emp.cov.hc, hpd.emp.cov.hc,
     file = "halfcauchy_simulations.RData")

#-------------------------------------------------------------------------------