#########################################################
#
# START COMMANDS TO CHECK THE GEV CDF AND PDF
#
#########################################################
library(evd)
source("cens_loglik_SP.R")

x <- seq(-10,10, length=200)
# GUMBEL CASE
mu <- 1
sigma <- 1
gam <- 0

plot(x, pgev(x, mu, sigma, gam), type="l", lwd=2, ylab="", main="Gumbel CDF")
lines(x, gevcdf(x, mu, sigma, gam), col=2, lwd=2, lty=2)

plot(x, dgev(x, mu, sigma, gam), type="l", lwd=2, ylab="", main="Gumbel PDF")
lines(x, gevpdf(x, mu, sigma, gam), col=2, lwd=2, lty=2)

# FRECHET CASE
mu <- 1
sigma <- 1
gam <- 1

plot(x, pgev(x, mu, sigma, gam), type="l", lwd=2, ylab="", main="Frechet CDF")
lines(x, gevcdf(x, mu, sigma, gam), col=2, lwd=2, lty=2)

plot(x, dgev(x, mu, sigma, gam), type="l", lwd=2, ylab="", main="Frechet PDF")
lines(x, gevpdf(x, mu, sigma, gam), col=2, lwd=2, lty=2)


# REVERSE-WEIBULL CASE
mu <- 1
sigma <- 1
gam <- -1

plot(x, pgev(x, mu, sigma, gam), type="l", lwd=2, ylab="", main="Reverse-Weibull CDF")
lines(x, gevcdf(x, mu, sigma, gam), col=2, lwd=2, lty=2)

plot(x, dgev(x, mu, sigma, gam), type="l", lwd=2, ylab="", main="Reverse-Weibull PDF")
lines(x, gevpdf(x, mu, sigma, gam), col=2, lwd=2, lty=2)
#
#########################################################
#
# END COMMANDS TO CHECK THE GEV CDF AND PDF
#
#########################################################

#########################################################
#
# START COMMANDS TO CHECK THE CENSORED GEV LOG-LIKELIHOOD
#
#########################################################
library(evd)
library(crch)
source("cens_loglik_SP.R")

### DEFINE SIMULATION SETTING
# sample size
n <- 1000
# exceeding probability used in the estimation
t.prob <- 0.05
# number of exceedances used in the estimation
k <- n * t.prob
# block-size for the correspoding block maxima approach
m <- 1/t.prob
# exceeding probability used for prediction
p <- 1/n
# return period
T <- 50

### 1) UNIT-FRECHET EXAMPLE
### FRECHET'S PARAMETERS
loc <- 0   # LOCATION
scale <- 1 # SCALE
shape <- 1 # INVERSE OF TAIL INDEX
### SET NORMING CONSTANTS
# CENTERING
bn <- (log(n/k) - log(n/k-1))^(-1/shape)
# SCALING
an <- 1/(shape*(n/k-1)) * (log(n/k)-log(n/k-1))^(-1/shape-1)

# DEFINES THE TRUE PARAMETERS
true.par <- c(bn, an, 1/shape)

# DEFINES THE ESTIMATION SETTING
start <- c(.1, 1, .1)
set.seed(1)
data <- rfrechet(n, loc, scale, shape)
th <- quantile(data, probs=1-t.prob, type=3)

# GEV CENSORED LOG-LIKEIHOOD METHOD
fit <- FitCensLLik(data, th, start, llik.type="Gev-Cens", p=p)
# TRUE PARAMETERS
true.par
# TRUE EXTREME QUANTILE
qfrechet(1-p, loc, scale, shape)
# TRUE RETURN LEVEL
qgev(1-1/T, bn, an, 1/shape)
# ESTIMATES
fit # IT SEEMS IT IS WORKING WELL...

# BAYESIAN ESTIMATION
sig0 <- 1
nsim <- 5e+4
burn <- 1e+4
method <- "bayesian"

fit.bayes <- FitCensLLik(data, th, start, p=p, method=method, sig0=sig0, nsim=nsim,
                         prior="empirical", burn=burn, val.show = TRUE)
fit.bayes$mle
colMeans(fit.bayes$post_sample)
sqrt(diag(var(fit.bayes$post_sample)))
colMeans(fit.bayes$post_sample)-qnorm(.975)*sqrt(diag(var(fit.bayes$post_sample)))
colMeans(fit.bayes$post_sample)+qnorm(.975)*sqrt(diag(var(fit.bayes$post_sample)))

par(mfrow=c(2,3))
plot(fit.bayes$post_sample[,1], type="l")
plot(fit.bayes$post_sample[,2], type="l")
plot(fit.bayes$post_sample[,3], type="l")

hist(fit.bayes$post_sample[,1])
hist(fit.bayes$post_sample[,2])
hist(fit.bayes$post_sample[,3])

#########################################################
#
# END COMMANDS TO CHECK THE CENSORED GEV LOG-LIKELIHOOD
#
#########################################################
