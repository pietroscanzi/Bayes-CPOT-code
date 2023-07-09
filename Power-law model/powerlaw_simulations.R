#Power Law posterior simulations

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
#Define the estimation setting
start.pl <- c(-0.28, 4.21, 0.21)

#Return objects
lab.par <- c("gamma", "mu", "delta")
lab.ci <- c("gamma.low", "gamma.up", "mu.low", "mu.up", "delta.low", "delta.up",
            "Q.ext.low", "Q.ext.up", "R.lev.low", "R.lev.up")
lab.cov <- c("gamma", "mu", "delta", "Q.ext", "R.lev")
true.par.pl <- matrix(0, nrow = nconf, ncol = 3)
colnames(true.par.pl) <- lab.par
Q.ext.true.pl <- rep(0, nconf)
R.lev.true.pl <- rep(0, nconf)
quant.ci.pl <- array(0, dim = c(Nsim, 10, nconf))
colnames(quant.ci.pl) <- lab.ci
quant.cov.pl <- array(0, dim = c(Nsim, 5, nconf))
colnames(quant.cov.pl) <- lab.cov
wald.ci.pl <- array(0, dim = c(Nsim, 10, nconf))
colnames(wald.ci.pl) <- lab.ci
wald.cov.pl <- array(0, dim = c(Nsim, 5, nconf))
colnames(wald.cov.pl) <- lab.cov
hpd.ci.pl <- array(0, dim = c(Nsim, 10, nconf))
colnames(hpd.ci.pl) <- lab.ci
hpd.cov.pl <- array(0, dim = c(Nsim, 5, nconf))
colnames(hpd.cov.pl) <- lab.cov

#Simulation loop
set.seed(7)
for(s in 1:nconf){
  
  cat("\n", "Starting with n =", n[s], "and k =", k[s], "\n")
  
  #Define the true parameters
  true.par.pl[s,] <- c(-1/alpha.pl, bn.pl[s], an.pl[s])
  #Define true extreme quantile
  Q.ext.true.pl[s] <- pl.qdf(1 - p[s], xstar.pl, alpha.pl, K = K)
  #Define true return level
  R.lev.true.pl[s] <- qgev(1 - 1/T.ret, bn.pl[s], an.pl[s], -1/alpha.pl)
  
  #Single configuration loop
  for(i in 1:Nsim){
    
    #Sample
    data.pl <- pl.qdf(runif(n[s]), xstar.pl, alpha.pl, K = K)
    #Threshold
    th.pl <- quantile(data.pl, probs = 1 - t.prob[s], type = 3)
    
    fit.bayes.pl <- fit.gev.bayes(data = data.pl, t = th.pl, t.prob = t.prob[s],
                                   llik.type = "Gev-Cens", param.type = "tilde",
                                  T.ret = T.ret, p = p[s],
                                   par0 = start.pl, k = 1, R = R, burn  = burn,
                                   prior = "empirical", etastar = 0.234)
    
    params <- fit.bayes.pl$parameters
    quant.ci.pl[i,1:2,s] <- as.numeric(quantile(params[,1], c(0.025, 0.975), na.rm = T))
    quant.cov.pl[i,1,s] <- as.numeric(true.par.pl[s,1] %between% quant.ci.pl[i,1:2,s])
    wald.ci.pl[i,1:2,s] <- mean(params[,1], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,1], na.rm = T)
    wald.cov.pl[i,1,s] <- as.numeric(true.par.pl[s,1] %between% wald.ci.pl[i,1:2,s])
    hpd.ci.pl[i,1:2,s] <- emp.hpd(as.numeric(na.omit(params[,1])), conf = 0.95)
    hpd.cov.pl[i,1,s] <- as.numeric(true.par.pl[s,1] %between% hpd.ci.pl[i,1:2,s])
    quant.ci.pl[i,3:4,s] <- as.numeric(quantile(params[,2], c(0.025, 0.975), na.rm = T))
    quant.cov.pl[i,2,s] <- as.numeric(true.par.pl[s,2] %between% quant.ci.pl[i,3:4,s])
    wald.ci.pl[i,3:4,s] <- mean(params[,2], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,2], na.rm = T)
    wald.cov.pl[i,2,s] <- as.numeric(true.par.pl[s,2] %between% wald.ci.pl[i,3:4,s])
    hpd.ci.pl[i,3:4,s] <- emp.hpd(as.numeric(na.omit(params[,2])), conf = 0.95)
    hpd.cov.pl[i,2,s] <- as.numeric(true.par.pl[s,2] %between% hpd.ci.pl[i,3:4,s])
    quant.ci.pl[i,5:6,s] <- as.numeric(quantile(params[,3], c(0.025, 0.975), na.rm = T))
    quant.cov.pl[i,3,s] <- as.numeric(true.par.pl[s,3] %between% quant.ci.pl[i,5:6,s])
    wald.ci.pl[i,5:6,s] <- mean(params[,3], na.rm = T) + c(-1,1)*qnorm(0.975)*sd(params[,3], na.rm = T)
    wald.cov.pl[i,3,s] <- as.numeric(true.par.pl[s,3] %between% wald.ci.pl[i,5:6,s])
    hpd.ci.pl[i,5:6,s] <- emp.hpd(as.numeric(na.omit(params[,3])), conf = 0.95)
    hpd.cov.pl[i,3,s] <- as.numeric(true.par.pl[s,3] %between% hpd.ci.pl[i,5:6,s])
    Q.ext <- fit.bayes.pl$Q.extreme
    quant.ci.pl[i,7:8,s] <- as.numeric(quantile(Q.ext, c(0.025, 0.975), na.rm = T))
    quant.cov.pl[i,4,s] <- as.numeric(Q.ext.true.pl[s] %between% quant.ci.pl[i,7:8,s])
    wald.ci.pl[i,7:8,s] <- mean(Q.ext, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(Q.ext, na.rm = T)
    wald.cov.pl[i,4,s] <- as.numeric(Q.ext.true.pl[s] %between% wald.ci.pl[i,7:8,s])
    hpd.ci.pl[i,7:8,s] <- emp.hpd(as.numeric(na.omit(Q.ext)), conf = 0.95)
    hpd.cov.pl[i,4,s] <- as.numeric(Q.ext.true.pl[s] %between% hpd.ci.pl[i,7:8,s])
    R.lev <- fit.bayes.pl$R.level
    quant.ci.pl[i,9:10,s] <- as.numeric(quantile(R.lev, c(0.025, 0.975), na.rm = T))
    quant.cov.pl[i,5,s] <- as.numeric(R.lev.true.pl[s] %between% quant.ci.pl[i,9:10,s])
    wald.ci.pl[i,9:10,s] <- mean(R.lev, na.rm = T) + c(-1,1)*qnorm(0.975)*sd(R.lev, na.rm = T)
    wald.cov.pl[i,5,s] <- as.numeric(R.lev.true.pl[s] %between% wald.ci.pl[i,9:10,s])
    hpd.ci.pl[i,9:10,s] <- emp.hpd(as.numeric(na.omit(R.lev)), conf = 0.95)
    hpd.cov.pl[i,5,s] <- as.numeric(R.lev.true.pl[s] %between% hpd.ci.pl[i,9:10,s])
    
    #Loop indicator 
    cat(i, " ")
  }
}

#Empirical coverage
lab.sim <- c("n = 800", "n = 1800", "n = 5450", "n = 23400")
quant.emp.cov.pl <- colMeans(quant.cov.pl)
colnames(quant.emp.cov.pl) <- lab.sim
wald.emp.cov.pl <- colMeans(wald.cov.pl)
colnames(wald.emp.cov.pl) <- lab.sim
hpd.emp.cov.pl <- colMeans(hpd.cov.pl)
colnames(hpd.emp.cov.pl) <- lab.sim

#Save
save(true.par.pl, Q.ext.true.pl, R.lev.true.pl, quant.ci.pl, quant.cov.pl,
     wald.ci.pl, wald.cov.pl, hpd.ci.pl, hpd.cov.pl, quant.emp.cov.pl,
     wald.emp.cov.pl, hpd.emp.cov.pl,
     file = "powerlaw_simulations.RData")

#-------------------------------------------------------------------------------