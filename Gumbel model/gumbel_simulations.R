#Gumbel posterior simulations

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

#standard Gumbel (tail index = 0)

#standard Gumbel parameters
loc.gum <- 0   
scale.gum <- 1 
#Set norming constants
#Centering
bn.gum <- -log(log(n) - log(n - k))
#Scaling
an.gum <- k/((n - k)*(log(n) - log(n - k)))
#an2.gum <- bn.gum - (k/n)*integrate(function(x) -log(-log(1 - (1/x))), lower = 1,
#                                    upper = (n/k))$value
#Define the estimation setting
start.gum <- c(-0.037, 2.91, 1.17)

#Return objects
lab.par <- c("gamma", "mu", "delta")
lab.ci <- c("gamma.low", "gamma.up", "mu.low", "mu.up", "delta.low", "delta.up",
            "Q.ext.low", "Q.ext.up", "R.lev.low", "R.lev.up")
lab.cov <- c("gamma", "mu", "delta", "Q.ext", "R.lev")
true.par.gum <- matrix(0, nrow = nconf, ncol = 3)
colnames(true.par.gum) <- lab.par
Q.ext.true.gum <- rep(0, nconf)
R.lev.true.gum <- rep(0, nconf)
quant.ci.gum <- array(0, dim = c(Nsim, 10, nconf))
colnames(quant.ci.gum) <- lab.ci
quant.cov.gum <- array(0, dim = c(Nsim, 5, nconf))
colnames(quant.cov.gum) <- lab.cov
wald.ci.gum <- array(0, dim = c(Nsim, 10, nconf))
colnames(wald.ci.gum) <- lab.ci
wald.cov.gum <- array(0, dim = c(Nsim, 5, nconf))
colnames(wald.cov.gum) <- lab.cov
hpd.ci.gum <- array(0, dim = c(Nsim, 10, nconf))
colnames(hpd.ci.gum) <- lab.ci
hpd.cov.gum <- array(0, dim = c(Nsim, 5, nconf))
colnames(hpd.cov.gum) <- lab.cov

#Simulation loop
set.seed(2)
for(s in 1:nconf){
  
  cat("\n", "Starting with n =", n[s], "and k =", k[s], "\n")
  
  #Define the true parameters
  true.par.gum[s,] <- c(0, bn.gum[s], an.gum[s])
  #Define true extreme quantile
  Q.ext.true.gum[s] <- qgumbel(1 - p[s], loc.gum, scale.gum)
  #Define true return level
  R.lev.true.gum[s] <- qgev(1 - 1/T.ret, bn.gum[s], an.gum[s], 0)
  
  #Single configuration loop
  for(i in 1:Nsim){
    
    #Sample
    data.gum <- rgumbel(n[s], loc.gum, scale.gum)
    #Threshold
    th.gum <- quantile(data.gum, probs = 1 - t.prob[s], type = 3)
    
    fit.bayes.gum <- fit.gev.bayes(data = data.gum, t = th.gum, t.prob = t.prob[s],
                                   llik.type = "Gev-Cens", param.type = "tilde",
                                   T.ret = T.ret, p = p[s],
                                   par0 = start.gum, k = 1, R = R, burn  = burn,
                                   prior = "empirical", etastar = 0.234)
    
    params <- fit.bayes.gum$parameters
    quant.ci.gum[i,1:2,s] <- as.numeric(quantile(params[,1], c(0.025, 0.975), na.rm = T))
    quant.cov.gum[i,1,s] <- as.numeric(true.par.gum[s,1] %between% quant.ci.gum[i,1:2,s])
    wald.ci.gum[i,1:2,s] <- mean(params[,1], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,1], na.rm = T)
    wald.cov.gum[i,1,s] <- as.numeric(true.par.gum[s,1] %between% wald.ci.gum[i,1:2,s])
    hpd.ci.gum[i,1:2,s] <- emp.hpd(as.numeric(na.omit(params[,1])), conf = 0.95)
    hpd.cov.gum[i,1,s] <- as.numeric(true.par.gum[s,1] %between% hpd.ci.gum[i,1:2,s])
    quant.ci.gum[i,3:4,s] <- as.numeric(quantile(params[,2], c(0.025, 0.975), na.rm = T))
    quant.cov.gum[i,2,s] <- as.numeric(true.par.gum[s,2] %between% quant.ci.gum[i,3:4,s])
    wald.ci.gum[i,3:4,s] <- mean(params[,2], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,2], na.rm = T)
    wald.cov.gum[i,2,s] <- as.numeric(true.par.gum[s,2] %between% wald.ci.gum[i,3:4,s])
    hpd.ci.gum[i,3:4,s] <- emp.hpd(as.numeric(na.omit(params[,2])), conf = 0.95)
    hpd.cov.gum[i,2,s] <- as.numeric(true.par.gum[s,2] %between% hpd.ci.gum[i,3:4,s])
    quant.ci.gum[i,5:6,s] <- as.numeric(quantile(params[,3], c(0.025, 0.975), na.rm = T))
    quant.cov.gum[i,3,s] <- as.numeric(true.par.gum[s,3] %between% quant.ci.gum[i,5:6,s])
    wald.ci.gum[i,5:6,s] <- mean(params[,3], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,3], na.rm = T)
    wald.cov.gum[i,3,s] <- as.numeric(true.par.gum[s,3] %between% wald.ci.gum[i,5:6,s])
    hpd.ci.gum[i,5:6,s] <- emp.hpd(as.numeric(na.omit(params[,3])), conf = 0.95)
    hpd.cov.gum[i,3,s] <- as.numeric(true.par.gum[s,3] %between% hpd.ci.gum[i,5:6,s])
    Q.ext <- fit.bayes.gum$Q.extreme
    quant.ci.gum[i,7:8,s] <- as.numeric(quantile(Q.ext, c(0.025, 0.975), na.rm = T))
    quant.cov.gum[i,4,s] <- as.numeric(Q.ext.true.gum[s] %between% quant.ci.gum[i,7:8,s])
    wald.ci.gum[i,7:8,s] <- mean(Q.ext, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(Q.ext, na.rm = T)
    wald.cov.gum[i,4,s] <- as.numeric(Q.ext.true.gum[s] %between% wald.ci.gum[i,7:8,s])
    hpd.ci.gum[i,7:8,s] <- emp.hpd(as.numeric(na.omit(Q.ext)), conf = 0.95)
    hpd.cov.gum[i,4,s] <- as.numeric(Q.ext.true.gum[s] %between% hpd.ci.gum[i,7:8,s])
    R.lev <- fit.bayes.gum$R.level
    quant.ci.gum[i,9:10,s] <- as.numeric(quantile(R.lev, c(0.025, 0.975), na.rm = T))
    quant.cov.gum[i,5,s] <- as.numeric(R.lev.true.gum[s] %between% quant.ci.gum[i,9:10,s])
    wald.ci.gum[i,9:10,s] <- mean(R.lev, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(R.lev, na.rm = T)
    wald.cov.gum[i,5,s] <- as.numeric(R.lev.true.gum[s] %between% wald.ci.gum[i,9:10,s])
    hpd.ci.gum[i,9:10,s] <- emp.hpd(as.numeric(na.omit(R.lev)), conf = 0.95)
    hpd.cov.gum[i,5,s] <- as.numeric(R.lev.true.gum[s] %between% hpd.ci.gum[i,9:10,s])
    
    #Loop indicator 
    cat(i, " ")
  }
}

#Empirical coverage
lab.sim <- c("n = 800", "n = 1800", "n = 5450", "n = 23400")
quant.emp.cov.gum <- colMeans(quant.cov.gum)
colnames(quant.emp.cov.gum) <- lab.sim
wald.emp.cov.gum <- colMeans(wald.cov.gum)
colnames(wald.emp.cov.gum) <- lab.sim
hpd.emp.cov.gum <- colMeans(hpd.cov.gum)
colnames(hpd.emp.cov.gum) <- lab.sim

#Save
save(true.par.gum, Q.ext.true.gum, R.lev.true.gum, quant.ci.gum, quant.cov.gum,
     wald.ci.gum, wald.cov.gum, hpd.ci.gum, hpd.cov.gum, quant.emp.cov.gum,
     wald.emp.cov.gum, hpd.emp.cov.gum,
     file = "gumbel_simulations.RData")

#-------------------------------------------------------------------------------