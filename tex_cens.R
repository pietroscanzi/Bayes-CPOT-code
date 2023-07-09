#Threshold exceedance via censored GEV likelihood

#Function computing the negative log-likelihood
nlogL <- function(par, data){
  #GEV distribution parameters
  gamma <- par[1]; mu <- par[2]; delta <- par[3]
  #Condition on the data for proper definition of the GEV distribution
#  if (sum(1 + gamma*((data - mu)/delta) < 0) > 0) return(Inf)
  #Empirical threshold
  t <- as.numeric(quantile(data, 0.95))
  #Indicator of censoring
  cond <- as.numeric(data > t)
  
  #gamma = 0 case --> Gumbel distribution
  if(gamma == 0){
    #Censored negative log-likelihood
    nlogL.cens <- sum(cond == 0)*(exp(-((t - mu)/delta)))
    #Observed negative log-likelihood
    nlogL.obs <- sum(cond)*log(delta) + sum((data[cond] - mu)/delta) +
      sum(exp(-((data[cond] - mu)/delta)))
    return(nlogL.cens + nlogL.obs)
  }
  
  #gamma != 0 case
  if(gamma != 0){
    #Censored negative log-likelihood
    nlogL.cens <- sum(cond == 0)*((1 + gamma*((t - mu)/delta))^(-1/gamma))
    #Observed negative log-likelihood
    nlogL.obs <- sum(cond)*log(delta) +
      ((gamma + 1)/gamma)*sum(log(1 + gamma*((data[cond] - mu)/delta))) +
      sum((1 + gamma*((data[cond] - mu)/delta))^(-1/gamma))
    return(nlogL.cens + nlogL.obs)
  }
}
#gamma \in (-Inf, Inf), mu \in (-Inf, Inf), delta \in (0, Inf)

#Data
n <- 10000


#Domain of attraction: Gumbel --> gamma = 0
#Exp(1) --> a_n = 1, b_n = log(n) (Coles p.51)
set.seed(1)
sam1 <- rexp(n, rate = 1)
hist(sam1, nclass = 50, probability = T,
     xlab = expression(y), ylab = "density", main = "Exp(1)")
curve(dexp(x, 1), col = "red", add = T)
#I try to sketch the form of the log-likelihood for fixed gamma = 0 
mu.val <- seq(-1, 1, length = 100)
delta.val <- seq(1e-5, 2, length = 100)
point.val <- expand.grid(mu.val, delta.val)
llik <- apply(point.val, 1, function(x) -nlogL(c(0, x[1], x[2]), sam1))
llik.mat <- matrix(llik, nrow = length(mu.val), ncol = length(delta.val))
contour(mu.val, delta.val, llik.mat, nlevels = 5,
        xlab = expression(mu), ylab = expression(delta), main = "log-likelihood")
#Maximum likelihood?
max(llik)
point.val[which.max(llik),]

#Optimization for MLE
mle1 <- nlminb(start = c(1,1,1), function(x) nlogL(x, sam1),
               lower = c(-Inf, -Inf, 1e-5), upper = rep(Inf, 3))
mle1
#Really unstable, convergence (when found) is different for near points.

mle2 <- optim(par = c(1,1,1), function(x) nlogL(x, sam1), method = "L-BFGS-B",
              lower = c(-Inf, -Inf, 1e-5), upper = rep(Inf, 3))
#Can't allow for infinite values of the likelihood


#Domain of attraction Fréchet --> gamma = 1

#Fre(0,1,1) --> a_n = n, b_n = 0 (Coles p.51)
#install.packages("evd")
library(evd)
set.seed(3)
sam2 <- rfrechet(n, shape = 1)
summary(sam2)
hist(sam2, nclass = 50, probability = T,
     xlab = expression(y), ylab = "density", main = "Fre(0,1,1)")
curve(dfrechet(x, shape = 1), col = "red", add = T)
ks.test(sam2, pfrechet, shape = 1)
#I try to sketch the form of the log-likelihood for fixed gamma = 1 
mu.val <- seq(-1, 1, length = 100)
delta.val <- seq(1e-5, 2, length = 100)
point.val <- expand.grid(mu.val, delta.val)
llik <- apply(point.val, 1, function(x) -nlogL(c(1, x[1], x[2]), sam2))
#The likelihood produces some NAs...
llik.mat <- matrix(llik, nrow = length(mu.val), ncol = length(delta.val))
contour(mu.val, delta.val, llik.mat, nlevels = 5,
        xlab = expression(mu), ylab = expression(delta), main = "log-likelihood")
#Maximum likelihood?
summary(llik)

#Optimization for MLE
mle1 <- nlminb(start = c(1,0,1), function(x) nlogL(x, sam2),
               lower = c(-Inf, -Inf, 1e-5), upper = rep(Inf, 3))
mle1