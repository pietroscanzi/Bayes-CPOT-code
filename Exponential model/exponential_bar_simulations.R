#Exponential posterior simulations (bar parameterization)

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

#Exponential(1) (tail index = 0)

#Exponential parameters
rate.exp <- 1
#Set norming constants
#Centering
bn.exp <- log(n) - log(k)
#Scaling
an.exp <- rep(1, nconf)
an2.exp <- rep(1, nconf)
#Define the estimation setting
start.exp <- c(0.0027, 0.35, 0.89)

#Return objects
lab.par <- c("gamma", "mu", "delta")
lab.ci <- c("gamma.low", "gamma.up", "mu.low", "mu.up", "delta.low", "delta.up",
            "Q.ext.low", "Q.ext.up", "R.lev.low", "R.lev.up")
lab.cov <- c("gamma", "mu", "delta", "Q.ext", "R.lev")
true.par.exp <- matrix(0, nrow = nconf, ncol = 3)
colnames(true.par.exp) <- lab.par
Q.ext.true.exp <- rep(0, nconf)
R.lev.true.exp <- rep(0, nconf)
quant.ci.exp <- array(0, dim = c(Nsim, 10, nconf))
colnames(quant.ci.exp) <- lab.ci
quant.cov.exp <- array(0, dim = c(Nsim, 5, nconf))
colnames(quant.cov.exp) <- lab.cov
wald.ci.exp <- array(0, dim = c(Nsim, 10, nconf))
colnames(wald.ci.exp) <- lab.ci
wald.cov.exp <- array(0, dim = c(Nsim, 5, nconf))
colnames(wald.cov.exp) <- lab.cov
hpd.ci.exp <- array(0, dim = c(Nsim, 10, nconf))
colnames(hpd.ci.exp) <- lab.ci
hpd.cov.exp <- array(0, dim = c(Nsim, 5, nconf))
colnames(hpd.cov.exp) <- lab.cov

#Simulation loop
set.seed(5)
for(s in 1:nconf){
  
  cat("\n", "Starting with n =", n[s], "and k =", k[s], "\n")
  
  #Define the true parameters
  true.par.exp[s,] <- c(0, bn.exp[s], an.exp[s])
  #Define true extreme quantile
  Q.ext.true.exp[s] <- qexp(1 - p[s], rate.exp)
  #Define true return level
  R.lev.true.exp[s] <- qgev(1 - 1/T.ret, bn.exp[s], an.exp[s], 0)
  
  #Single configuration loop
  for(i in 1:Nsim){
    
    #Sample
    data.exp <- rexp(n[s], rate.exp)
    #Threshold
    th.exp <- quantile(data.exp, probs = 1 - t.prob[s], type = 3)
    
    fit.bayes.exp <- fit.gev.bayes(data = data.exp, t = th.exp, t.prob = t.prob[s],
                                   llik.type = "Gev-Cens", param.type = "bar",
                                   T.ret = T.ret, p = p[s],
                                   par0 = start.exp, k = 1, R = R, burn  = burn,
                                   prior = "empirical", etastar = 0.234)
    
    params.bar <- fit.bayes.exp$parameters
    params <- params.bar
    params[,3] <- apply(params.bar, 1, function(x) x[3]*s.p[s]^(x[1]))
    params[,2] <- apply(params, 1, function(x) x[2] + x[3]*(1 - s.p[s]^(-x[1]))/x[1])
    quant.ci.exp[i,1:2,s] <- as.numeric(quantile(params[,1], c(0.025, 0.975), na.rm = T))
    quant.cov.exp[i,1,s] <- as.numeric(true.par.exp[s,1] %between% quant.ci.exp[i,1:2,s])
    wald.ci.exp[i,1:2,s] <- mean(params[,1], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,1], na.rm = T)
    wald.cov.exp[i,1,s] <- as.numeric(true.par.exp[s,1] %between% wald.ci.exp[i,1:2,s])
    hpd.ci.exp[i,1:2,s] <- emp.hpd(as.numeric(na.omit(params[,1])), conf = 0.95)
    hpd.cov.exp[i,1,s] <- as.numeric(true.par.exp[s,1] %between% hpd.ci.exp[i,1:2,s])
    quant.ci.exp[i,3:4,s] <- as.numeric(quantile(params[,2], c(0.025, 0.975), na.rm = T))
    quant.cov.exp[i,2,s] <- as.numeric(true.par.exp[s,2] %between% quant.ci.exp[i,3:4,s])
    wald.ci.exp[i,3:4,s] <- mean(params[,2], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,2], na.rm = T)
    wald.cov.exp[i,2,s] <- as.numeric(true.par.exp[s,2] %between% wald.ci.exp[i,3:4,s])
    hpd.ci.exp[i,3:4,s] <- emp.hpd(as.numeric(na.omit(params[,2])), conf = 0.95)
    hpd.cov.exp[i,2,s] <- as.numeric(true.par.exp[s,2] %between% hpd.ci.exp[i,3:4,s])
    quant.ci.exp[i,5:6,s] <- as.numeric(quantile(params[,3], c(0.025, 0.975), na.rm = T))
    quant.cov.exp[i,3,s] <- as.numeric(true.par.exp[s,3] %between% quant.ci.exp[i,5:6,s])
    wald.ci.exp[i,5:6,s] <- mean(params[,3], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,3], na.rm = T)
    wald.cov.exp[i,3,s] <- as.numeric(true.par.exp[s,3] %between% wald.ci.exp[i,5:6,s])
    hpd.ci.exp[i,5:6,s] <- emp.hpd(as.numeric(na.omit(params[,3])), conf = 0.95)
    hpd.cov.exp[i,3,s] <- as.numeric(true.par.exp[s,3] %between% hpd.ci.exp[i,5:6,s])
    
    Q.ext <- params[,2] + params[,3]*((s.p[s]*p[s])^(-params[,1]) - 1)/params[,1]
    quant.ci.exp[i,7:8,s] <- as.numeric(quantile(Q.ext, c(0.025, 0.975), na.rm = T))
    quant.cov.exp[i,4,s] <- as.numeric(Q.ext.true.exp[s] %between% quant.ci.exp[i,7:8,s])
    wald.ci.exp[i,7:8,s] <- mean(Q.ext, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(Q.ext, na.rm = T)
    wald.cov.exp[i,4,s] <- as.numeric(Q.ext.true.exp[s] %between% wald.ci.exp[i,7:8,s])
    hpd.ci.exp[i,7:8,s] <- emp.hpd(as.numeric(na.omit(Q.ext)), conf = 0.95)
    hpd.cov.exp[i,4,s] <- as.numeric(Q.ext.true.exp[s] %between% hpd.ci.exp[i,7:8,s])
    
    R.lev <- params[,2] + params[,3]*((-log(1 - 1/T.ret))^(-params[,1]) - 1)/params[,1]
    quant.ci.exp[i,9:10,s] <- as.numeric(quantile(R.lev, c(0.025, 0.975), na.rm = T))
    quant.cov.exp[i,5,s] <- as.numeric(R.lev.true.exp[s] %between% quant.ci.exp[i,9:10,s])
    wald.ci.exp[i,9:10,s] <- mean(R.lev, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(R.lev, na.rm = T)
    wald.cov.exp[i,5,s] <- as.numeric(R.lev.true.exp[s] %between% wald.ci.exp[i,9:10,s])
    hpd.ci.exp[i,9:10,s] <- emp.hpd(as.numeric(na.omit(R.lev)), conf = 0.95)
    hpd.cov.exp[i,5,s] <- as.numeric(R.lev.true.exp[s] %between% hpd.ci.exp[i,9:10,s])
    
    #Loop indicator 
    cat(i, " ")
  }
}

#Empirical coverage
lab.sim <- c("n = 800", "n = 1800", "n = 5450", "n = 23400")
quant.emp.cov.exp <- colMeans(quant.cov.exp)
colnames(quant.emp.cov.exp) <- lab.sim
wald.emp.cov.exp <- colMeans(wald.cov.exp)
colnames(wald.emp.cov.exp) <- lab.sim
hpd.emp.cov.exp <- colMeans(hpd.cov.exp)
colnames(hpd.emp.cov.exp) <- lab.sim

#Save
save(true.par.exp, Q.ext.true.exp, R.lev.true.exp, quant.ci.exp, quant.cov.exp,
     wald.ci.exp, wald.cov.exp, hpd.ci.exp, hpd.cov.exp, quant.emp.cov.exp,
     wald.emp.cov.exp, hpd.emp.cov.exp,
     file = "exponential_bar_simulations.RData")

#-------------------------------------------------------------------------------