#Pareto posterior simulations

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
#Define the estimation setting
start.par <- c(0.81, 21.79, 21.09)

#Return objects
lab.par <- c("gamma", "mu", "delta")
lab.ci <- c("gamma.low", "gamma.up", "mu.low", "mu.up", "delta.low", "delta.up",
            "Q.ext.low", "Q.ext.up", "R.lev.low", "R.lev.up")
lab.cov <- c("gamma", "mu", "delta", "Q.ext", "R.lev")
true.par.par <- matrix(0, nrow = nconf, ncol = 3)
colnames(true.par.par) <- lab.par
Q.ext.true.par <- rep(0, nconf)
R.lev.true.par <- rep(0, nconf)
quant.ci.par <- array(0, dim = c(Nsim, 10, nconf))
colnames(quant.ci.par) <- lab.ci
quant.cov.par <- array(0, dim = c(Nsim, 5, nconf))
colnames(quant.cov.par) <- lab.cov
wald.ci.par <- array(0, dim = c(Nsim, 10, nconf))
colnames(wald.ci.par) <- lab.ci
wald.cov.par <- array(0, dim = c(Nsim, 5, nconf))
colnames(wald.cov.par) <- lab.cov
hpd.ci.par <- array(0, dim = c(Nsim, 10, nconf))
colnames(hpd.ci.par) <- lab.ci
hpd.cov.par <- array(0, dim = c(Nsim, 5, nconf))
colnames(hpd.cov.par) <- lab.cov

#Simulation loop
set.seed(2)
for(s in 1:nconf){
  
  cat("\n", "Starting with n =", n[s], "and k =", k[s], "\n")
  
  #Define the true parameters
  true.par.par[s,] <- c(1/shape.par, bn.par[s], an.par[s])
  #Define true extreme quantile
  Q.ext.true.par[s] <- extremefit::qpareto(1 - p[s], shape.par, loc.par, scale.par)
  #Define true return level
  R.lev.true.par[s] <- qgev(1 - 1/T.ret, bn.par[s], an.par[s], 1/shape.par)
  
  #Single configuration loop
  for(i in 1:Nsim){
    
    #Sample
    data.par <- extremefit::rpareto(n[s], shape.par, loc.par, scale.par)
    #Threshold
    th.par <- quantile(data.par, probs = 1 - t.prob[s], type = 3)
    
    fit.bayes.par <- fit.gev.bayes(data = data.par, t = th.par, t.prob = t.prob[s],
                                   llik.type = "Max-Gev-Cens", param.type = "tilde",
                                   T.ret = T.ret, p = p[s],
                                   par0 = start.par, k = 1, R = R, burn  = burn,
                                   prior = "empirical", etastar = 0.234)
    
    params <- fit.bayes.par$parameters
    quant.ci.par[i,1:2,s] <- as.numeric(quantile(params[,1], c(0.025, 0.975), na.rm = T))
    quant.cov.par[i,1,s] <- as.numeric(true.par.par[s,1] %between% quant.ci.par[i,1:2,s])
    wald.ci.par[i,1:2,s] <- mean(params[,1], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,1], na.rm = T)
    wald.cov.par[i,1,s] <- as.numeric(true.par.par[s,1] %between% wald.ci.par[i,1:2,s])
    hpd.ci.par[i,1:2,s] <- emp.hpd(as.numeric(na.omit(params[,1])), conf = 0.95)
    hpd.cov.par[i,1,s] <- as.numeric(true.par.par[s,1] %between% hpd.ci.par[i,1:2,s])
    quant.ci.par[i,3:4,s] <- as.numeric(quantile(params[,2], c(0.025, 0.975), na.rm = T))
    quant.cov.par[i,2,s] <- as.numeric(true.par.par[s,2] %between% quant.ci.par[i,3:4,s])
    wald.ci.par[i,3:4,s] <- mean(params[,2], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,2], na.rm = T)
    wald.cov.par[i,2,s] <- as.numeric(true.par.par[s,2] %between% wald.ci.par[i,3:4,s])
    hpd.ci.par[i,3:4,s] <- emp.hpd(as.numeric(na.omit(params[,2])), conf = 0.95)
    hpd.cov.par[i,2,s] <- as.numeric(true.par.par[s,2] %between% hpd.ci.par[i,3:4,s])
    quant.ci.par[i,5:6,s] <- as.numeric(quantile(params[,3], c(0.025, 0.975), na.rm = T))
    quant.cov.par[i,3,s] <- as.numeric(true.par.par[s,3] %between% quant.ci.par[i,5:6,s])
    wald.ci.par[i,5:6,s] <- mean(params[,3], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,3], na.rm = T)
    wald.cov.par[i,3,s] <- as.numeric(true.par.par[s,3] %between% wald.ci.par[i,5:6,s])
    hpd.ci.par[i,5:6,s] <- emp.hpd(as.numeric(na.omit(params[,3])), conf = 0.95)
    hpd.cov.par[i,3,s] <- as.numeric(true.par.par[s,3] %between% hpd.ci.par[i,5:6,s])
    Q.ext <- fit.bayes.par$Q.extreme
    quant.ci.par[i,7:8,s] <- as.numeric(quantile(Q.ext, c(0.025, 0.975), na.rm = T))
    quant.cov.par[i,4,s] <- as.numeric(Q.ext.true.par[s] %between% quant.ci.par[i,7:8,s])
    wald.ci.par[i,7:8,s] <- mean(Q.ext, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(Q.ext, na.rm = T)
    wald.cov.par[i,4,s] <- as.numeric(Q.ext.true.par[s] %between% wald.ci.par[i,7:8,s])
    hpd.ci.par[i,7:8,s] <- emp.hpd(as.numeric(na.omit(Q.ext)), conf = 0.95)
    hpd.cov.par[i,4,s] <- as.numeric(Q.ext.true.par[s] %between% hpd.ci.par[i,7:8,s])
    R.lev <- fit.bayes.par$R.level
    quant.ci.par[i,9:10,s] <- as.numeric(quantile(R.lev, c(0.025, 0.975), na.rm = T))
    quant.cov.par[i,5,s] <- as.numeric(R.lev.true.par[s] %between% quant.ci.par[i,9:10,s])
    wald.ci.par[i,9:10,s] <- mean(R.lev, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(R.lev, na.rm = T)
    wald.cov.par[i,5,s] <- as.numeric(R.lev.true.par[s] %between% wald.ci.par[i,9:10,s])
    hpd.ci.par[i,9:10,s] <- emp.hpd(as.numeric(na.omit(R.lev)), conf = 0.95)
    hpd.cov.par[i,5,s] <- as.numeric(R.lev.true.par[s] %between% hpd.ci.par[i,9:10,s])
    
    #Loop indicator 
    cat(i, " ")
  }
}

#Empirical coverage
lab.sim <- c("n = 800", "n = 1800", "n = 5450", "n = 23400")
quant.emp.cov.par <- colMeans(quant.cov.par)
colnames(quant.emp.cov.par) <- lab.sim
wald.emp.cov.par <- colMeans(wald.cov.par)
colnames(wald.emp.cov.par) <- lab.sim
hpd.emp.cov.par <- colMeans(hpd.cov.par)
colnames(hpd.emp.cov.par) <- lab.sim

#Save
save(true.par.par, Q.ext.true.par, R.lev.true.par, quant.ci.par, quant.cov.par,
     wald.ci.par, wald.cov.par, hpd.ci.par, hpd.cov.par, quant.emp.cov.par,
     wald.emp.cov.par, hpd.emp.cov.par,
     file = "pareto_max_simulations.RData")

#-------------------------------------------------------------------------------