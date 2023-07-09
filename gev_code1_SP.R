
# N.B. we should avoid to use "gam" instead of "gamma"
# for the gamma paratemer since that gamma is also a
# function

gevcdf.root <- function(x, mu, sigma, gam)
{
  # define zero neighbourhood
  eps <- .Machine$double.eps^.3
  # set the three domains
  # GUMBEL
  if(abs(gam)<=eps) return (exp(-exp(-(x-mu)/sigma)))
  else{
    # FRECHET
    if(gam>0){
      if(x<=(mu-sigma/gam)) return(0)
      else return(exp(-(1+gam*(x-mu)/sigma)^(-1/gam)))
    }
    # REVERSE WEIBULL
    if(gam<0){
      if(x<=(mu-sigma/gam)) return(exp(-(1+gam*(x-mu)/sigma)^(-1/gam)))
      else return(1)
    }
  }
}

gevpdf.root <- function(x, mu, sigma, gam)
{
  # define zero neighbourhood
  eps <- .Machine$double.eps^.3
  # set the three domains
  # GUMBEL
  if(abs(gam)<=eps){
    z <- (x-mu)/sigma
    return (exp(-exp(-z)-z)/sigma)
  }
  else{
    # FRECHET
    if(gam>0){
      if(x<=(mu-sigma/gam)) return(0)
      else{
        z <- 1 + gam*(x-mu)/sigma
        return(1/sigma * exp(-z^(-1/gam)) * z^(-1/gam-1))
      }
    }
    # REVERSE WEIBULL
    if(gam<0){
      if(x<=(mu-sigma/gam)){
        z <- 1 + gam*(x-mu)/sigma
        return(1/sigma * exp(-z^(-1/gam)) * z^(-1/gam-1))
      }
      else return(0)
    }
  }
}
gevcdf <- Vectorize(gevcdf.root, "x")
gevpdf <- Vectorize(gevpdf.root, "x")


# CHECK
library(evd)

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
