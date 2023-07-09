#Gamma posterior simulations (bar parameterization)

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

#Unit Fréchet (tail index = 1)

#Gamma parameters
shape.gam <- 2 
rate.gam <- 2
#Set norming constants
#Centering
bn.gam <- qgamma(1 - (k/n), shape = shape.gam, rate = rate.gam)
#Scaling
#h <- 1e-10
an.gam <- k/(n*dgamma(bn.gam, shape = shape.gam, rate = rate.gam))
#an2.gam <- bn.gam - (k/n)*integrate(function(x) qgamma(1 - (1/x), shape = shape.gam,
#                                                       rate = rate.gam), lower = 1,
#                                    upper = (n/k))$value
#Define the estimation setting
start.gam <- c(0.14, 1.21, 0.34)

#Return objects
lab.par <- c("gamma", "mu", "delta")
lab.ci <- c("gamma.low", "gamma.up", "mu.low", "mu.up", "delta.low", "delta.up",
            "Q.ext.low", "Q.ext.up", "R.lev.low", "R.lev.up")
lab.cov <- c("gamma", "mu", "delta", "Q.ext", "R.lev")
true.par.gam <- matrix(0, nrow = nconf, ncol = 3)
colnames(true.par.gam) <- lab.par
Q.ext.true.gam <- rep(0, nconf)
R.lev.true.gam <- rep(0, nconf)
quant.ci.gam <- array(0, dim = c(Nsim, 10, nconf))
colnames(quant.ci.gam) <- lab.ci
quant.cov.gam <- array(0, dim = c(Nsim, 5, nconf))
colnames(quant.cov.gam) <- lab.cov
wald.ci.gam <- array(0, dim = c(Nsim, 10, nconf))
colnames(wald.ci.gam) <- lab.ci
wald.cov.gam <- array(0, dim = c(Nsim, 5, nconf))
colnames(wald.cov.gam) <- lab.cov
hpd.ci.gam <- array(0, dim = c(Nsim, 10, nconf))
colnames(hpd.ci.gam) <- lab.ci
hpd.cov.gam <- array(0, dim = c(Nsim, 5, nconf))
colnames(hpd.cov.gam) <- lab.cov

#Simulation loop
set.seed(9)
for(s in 1:nconf){
  
  cat("\n", "Starting with n =", n[s], "and k =", k[s], "\n")
  
  #Define the true parameters
  true.par.gam[s,] <- c(0, bn.gam[s], an.gam[s])
  #Define true extreme quantile
  Q.ext.true.gam[s] <- qgamma(1 - p[s], shape.gam, rate.gam)
  #Define true return level
  R.lev.true.gam[s] <- qgev(1 - 1/T.ret, bn.gam[s], an.gam[s], 0)
  
  #Single configuration loop
  for(i in 1:Nsim){
    
    #Sample
    data.gam <- rgamma(n[s], shape.gam, rate.gam)
    #Threshold
    th.gam <- quantile(data.gam, probs = 1 - t.prob[s], type = 3)
    
    fit.bayes.gam <- fit.gev.bayes(data = data.gam, t = th.gam, t.prob = t.prob[s],
                                   llik.type = "Gev-Cens", param.type = "bar",
                                   T.ret = T.ret, p = p[s],
                                   par0 = start.gam, k = 1, R = R, burn  = burn,
                                   prior = "empirical", etastar = 0.234)
    
    params.bar <- fit.bayes.gam$parameters
    params <- params.bar
    params[,3] <- apply(params.bar, 1, function(x) x[3]*s.p[s]^(x[1]))
    params[,2] <- apply(params, 1, function(x) x[2] + x[3]*(1 - s.p[s]^(-x[1]))/x[1])
    quant.ci.gam[i,1:2,s] <- as.numeric(quantile(params[,1], c(0.025, 0.975), na.rm = T))
    quant.cov.gam[i,1,s] <- as.numeric(true.par.gam[s,1] %between% quant.ci.gam[i,1:2,s])
    wald.ci.gam[i,1:2,s] <- mean(params[,1], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,1], na.rm = T)
    wald.cov.gam[i,1,s] <- as.numeric(true.par.gam[s,1] %between% wald.ci.gam[i,1:2,s])
    hpd.ci.gam[i,1:2,s] <- emp.hpd(as.numeric(na.omit(params[,1])), conf = 0.95)
    hpd.cov.gam[i,1,s] <- as.numeric(true.par.gam[s,1] %between% hpd.ci.gam[i,1:2,s])
    quant.ci.gam[i,3:4,s] <- as.numeric(quantile(params[,2], c(0.025, 0.975), na.rm = T))
    quant.cov.gam[i,2,s] <- as.numeric(true.par.gam[s,2] %between% quant.ci.gam[i,3:4,s])
    wald.ci.gam[i,3:4,s] <- mean(params[,2], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,2], na.rm = T)
    wald.cov.gam[i,2,s] <- as.numeric(true.par.gam[s,2] %between% wald.ci.gam[i,3:4,s])
    hpd.ci.gam[i,3:4,s] <- emp.hpd(as.numeric(na.omit(params[,2])), conf = 0.95)
    hpd.cov.gam[i,2,s] <- as.numeric(true.par.gam[s,2] %between% hpd.ci.gam[i,3:4,s])
    quant.ci.gam[i,5:6,s] <- as.numeric(quantile(params[,3], c(0.025, 0.975), na.rm = T))
    quant.cov.gam[i,3,s] <- as.numeric(true.par.gam[s,3] %between% quant.ci.gam[i,5:6,s])
    wald.ci.gam[i,5:6,s] <- mean(params[,3], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,3], na.rm = T)
    wald.cov.gam[i,3,s] <- as.numeric(true.par.gam[s,3] %between% wald.ci.gam[i,5:6,s])
    hpd.ci.gam[i,5:6,s] <- emp.hpd(as.numeric(na.omit(params[,3])), conf = 0.95)
    hpd.cov.gam[i,3,s] <- as.numeric(true.par.gam[s,3] %between% hpd.ci.gam[i,5:6,s])
    
    Q.ext <- params[,2] + params[,3]*((s.p[s]*p[s])^(-params[,1]) - 1)/params[,1]
    quant.ci.gam[i,7:8,s] <- as.numeric(quantile(Q.ext, c(0.025, 0.975), na.rm = T))
    quant.cov.gam[i,4,s] <- as.numeric(Q.ext.true.gam[s] %between% quant.ci.gam[i,7:8,s])
    wald.ci.gam[i,7:8,s] <- mean(Q.ext, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(Q.ext, na.rm = T)
    wald.cov.gam[i,4,s] <- as.numeric(Q.ext.true.gam[s] %between% wald.ci.gam[i,7:8,s])
    hpd.ci.gam[i,7:8,s] <- emp.hpd(as.numeric(na.omit(Q.ext)), conf = 0.95)
    hpd.cov.gam[i,4,s] <- as.numeric(Q.ext.true.gam[s] %between% hpd.ci.gam[i,7:8,s])
    
    R.lev <- params[,2] + params[,3]*((-log(1 - 1/T.ret))^(-params[,1]) - 1)/params[,1]
    quant.ci.gam[i,9:10,s] <- as.numeric(quantile(R.lev, c(0.025, 0.975), na.rm = T))
    quant.cov.gam[i,5,s] <- as.numeric(R.lev.true.gam[s] %between% quant.ci.gam[i,9:10,s])
    wald.ci.gam[i,9:10,s] <- mean(R.lev, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(R.lev, na.rm = T)
    wald.cov.gam[i,5,s] <- as.numeric(R.lev.true.gam[s] %between% wald.ci.gam[i,9:10,s])
    hpd.ci.gam[i,9:10,s] <- emp.hpd(as.numeric(na.omit(R.lev)), conf = 0.95)
    hpd.cov.gam[i,5,s] <- as.numeric(R.lev.true.gam[s] %between% hpd.ci.gam[i,9:10,s])
    
    #Loop indicator 
    cat(i, " ")
  }
}

#Empirical coverage
lab.sim <- c("n = 800", "n = 1800", "n = 5450", "n = 23400")
quant.emp.cov.gam <- colMeans(quant.cov.gam)
colnames(quant.emp.cov.gam) <- lab.sim
wald.emp.cov.gam <- colMeans(wald.cov.gam)
colnames(wald.emp.cov.gam) <- lab.sim
hpd.emp.cov.gam <- colMeans(hpd.cov.gam)
colnames(hpd.emp.cov.gam) <- lab.sim

#Save
save(true.par.gam, Q.ext.true.gam, R.lev.true.gam, quant.ci.gam, quant.cov.gam,
     wald.ci.gam, wald.cov.gam, hpd.ci.gam, hpd.cov.gam, quant.emp.cov.gam,
     wald.emp.cov.gam, hpd.emp.cov.gam,
     file = "gamma_bar_simulations.RData")

#-------------------------------------------------------------------------------