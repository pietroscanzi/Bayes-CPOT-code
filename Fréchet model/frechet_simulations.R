#Fréchet posterior simulations

#-------------------------------------------------------------------------------

path <- "C:\\Users\\scanz\\OneDrive - Università Commerciale Luigi Bocconi\\Desktop/Magistrale\\Extreme Value Theory\\Code\\"
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
start.fre <- c(0.90, 19.41, 21.32)

#Return objects
lab.par <- c("gamma", "mu", "delta")
lab.ci <- c("gamma.low", "gamma.up", "mu.low", "mu.up", "delta.low", "delta.up",
            "Q.ext.low", "Q.ext.up", "R.lev.low", "R.lev.up")
lab.cov <- c("gamma", "mu", "delta", "Q.ext", "R.lev")
true.par.fre <- matrix(0, nrow = nconf, ncol = 3)
colnames(true.par.fre) <- lab.par
Q.ext.true.fre <- rep(0, nconf)
R.lev.true.fre <- rep(0, nconf)
quant.ci.fre <- array(0, dim = c(Nsim, 10, nconf))
colnames(quant.ci.fre) <- lab.ci
quant.cov.fre <- array(0, dim = c(Nsim, 5, nconf))
colnames(quant.cov.fre) <- lab.cov
wald.ci.fre <- array(0, dim = c(Nsim, 10, nconf))
colnames(wald.ci.fre) <- lab.ci
wald.cov.fre <- array(0, dim = c(Nsim, 5, nconf))
colnames(wald.cov.fre) <- lab.cov
hpd.ci.fre <- array(0, dim = c(Nsim, 10, nconf))
colnames(hpd.ci.fre) <- lab.ci
hpd.cov.fre <- array(0, dim = c(Nsim, 5, nconf))
colnames(hpd.cov.fre) <- lab.cov

#Simulation loop
set.seed(1)
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
    
    fit.bayes.fre <- fit.gev.bayes(data = data.fre, t = th.fre, t.prob = t.prob[s],
                                    llik.type = "Gev-Cens", param.type = "tilde",
                                    T.ret = T.ret, p = p[s],
                                    par0 = start.fre, k = 1, R = R, burn  = burn,
                                    prior = "empirical", etastar = 0.234)
    
    params <- fit.bayes.fre$parameters
    quant.ci.fre[i,1:2,s] <- as.numeric(quantile(params[,1], c(0.025, 0.975), na.rm = T))
    quant.cov.fre[i,1,s] <- as.numeric(true.par.fre[s,1] %between% quant.ci.fre[i,1:2,s])
    wald.ci.fre[i,1:2,s] <- mean(params[,1], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,1], na.rm = T)
    wald.cov.fre[i,1,s] <- as.numeric(true.par.fre[s,1] %between% wald.ci.fre[i,1:2,s])
    hpd.ci.fre[i,1:2,s] <- emp.hpd(as.numeric(na.omit(params[,1])), conf = 0.95)
    hpd.cov.fre[i,1,s] <- as.numeric(true.par.fre[s,1] %between% hpd.ci.fre[i,1:2,s])
    quant.ci.fre[i,3:4,s] <- as.numeric(quantile(params[,2], c(0.025, 0.975), na.rm = T))
    quant.cov.fre[i,2,s] <- as.numeric(true.par.fre[s,2] %between% quant.ci.fre[i,3:4,s])
    wald.ci.fre[i,3:4,s] <- mean(params[,2], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,2], na.rm = T)
    wald.cov.fre[i,2,s] <- as.numeric(true.par.fre[s,2] %between% wald.ci.fre[i,3:4,s])
    hpd.ci.fre[i,3:4,s] <- emp.hpd(as.numeric(na.omit(params[,2])), conf = 0.95)
    hpd.cov.fre[i,2,s] <- as.numeric(true.par.fre[s,2] %between% hpd.ci.fre[i,3:4,s])
    quant.ci.fre[i,5:6,s] <- as.numeric(quantile(params[,3], c(0.025, 0.975), na.rm = T))
    quant.cov.fre[i,3,s] <- as.numeric(true.par.fre[s,3] %between% quant.ci.fre[i,5:6,s])
    wald.ci.fre[i,5:6,s] <- mean(params[,3], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,3], na.rm = T)
    wald.cov.fre[i,3,s] <- as.numeric(true.par.fre[s,3] %between% wald.ci.fre[i,5:6,s])
    hpd.ci.fre[i,5:6,s] <- emp.hpd(as.numeric(na.omit(params[,3])), conf = 0.95)
    hpd.cov.fre[i,3,s] <- as.numeric(true.par.fre[s,3] %between% hpd.ci.fre[i,5:6,s])
    Q.ext <- fit.bayes.fre$Q.extreme
    quant.ci.fre[i,7:8,s] <- as.numeric(quantile(Q.ext, c(0.025, 0.975), na.rm = T))
    quant.cov.fre[i,4,s] <- as.numeric(Q.ext.true.fre[s] %between% quant.ci.fre[i,7:8,s])
    wald.ci.fre[i,7:8,s] <- mean(Q.ext, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(Q.ext, na.rm = T)
    wald.cov.fre[i,4,s] <- as.numeric(Q.ext.true.fre[s] %between% wald.ci.fre[i,7:8,s])
    hpd.ci.fre[i,7:8,s] <- emp.hpd(as.numeric(na.omit(Q.ext)), conf = 0.95)
    hpd.cov.fre[i,4,s] <- as.numeric(Q.ext.true.fre[s] %between% hpd.ci.fre[i,7:8,s])
    R.lev <- fit.bayes.fre$R.level
    quant.ci.fre[i,9:10,s] <- as.numeric(quantile(R.lev, c(0.025, 0.975), na.rm = T))
    quant.cov.fre[i,5,s] <- as.numeric(R.lev.true.fre[s] %between% quant.ci.fre[i,9:10,s])
    wald.ci.fre[i,9:10,s] <- mean(R.lev, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(R.lev, na.rm = T)
    wald.cov.fre[i,5,s] <- as.numeric(R.lev.true.fre[s] %between% wald.ci.fre[i,9:10,s])
    hpd.ci.fre[i,9:10,s] <- emp.hpd(as.numeric(na.omit(R.lev)), conf = 0.95)
    hpd.cov.fre[i,5,s] <- as.numeric(R.lev.true.fre[s] %between% hpd.ci.fre[i,9:10,s])
    
    #Loop indicator 
    cat(i, " ")
  }
}

#Empirical coverage
lab.sim <- c("n = 800", "n = 1800", "n = 5450", "n = 23400")
quant.emp.cov.fre <- colMeans(quant.cov.fre)
colnames(quant.emp.cov.fre) <- lab.sim
wald.emp.cov.fre <- colMeans(wald.cov.fre)
colnames(wald.emp.cov.fre) <- lab.sim
hpd.emp.cov.fre <- colMeans(hpd.cov.fre)
colnames(hpd.emp.cov.fre) <- lab.sim

#Save
save(true.par.fre, Q.ext.true.fre, R.lev.true.fre, quant.ci.fre, quant.cov.fre,
     wald.ci.fre, wald.cov.fre, hpd.ci.fre, hpd.cov.fre, quant.emp.cov.fre,
     wald.emp.cov.fre, hpd.emp.cov.fre,
     file = "frechet_simulations.RData")

#-------------------------------------------------------------------------------