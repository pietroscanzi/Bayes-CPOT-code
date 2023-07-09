

gev.cens.llik <- function(param, data, u){
  # set the smallest value
  LOW <- -1e300
  # check the validity of the data
  if(!is.vector(data)){stop("Data must be a vector")}
  # check the validity of the parameters
  if(length(param) != 3){stop("Wrong length of parameter vector")}
  # parameters definition
  mu <- param[1]
  sigma <- param[2]
  gam <- param[3]
  # check the validity of the scale
  # N.B. THIS SHOULD BE THE ONLY CHECK
  if(sigma<0){return(LOW)}
  # derive the exceedances
  I <- data>u
  exceeds <- data[I]
  # set sample size and effective sample size
  n <- length(data)
  k <- sum(I)
  # likelihood concerning the exceedances
  up <- sum(log(gevpdf(exceeds, mu, sigma, gam)))
  # likelihood concerning the not exceedances
  low <- (n-k) * log(gevcdf(u, mu, sigma, gam))
  res <- up+low
  if(is.infinite(res))return(LOW)
  return(res)
}



FitCensLLik<- function(data, u=NULL, param0=NULL, optim.meth="Nelder-Mead",
                       control=NULL, p=0.9, ...){

  if(is.null(param0)){stop("Need to specify param0 parameter")} # param0 is the starting value for both methods

  if(is.null(u)){
    u <- quantile(data, probs=p, type=3)
    message("U set to 90% quantile by default \n")
  }

  if(is.null(optim.meth)){stop("Need to specify an optimisation method in frequentist setting")}

  if(!("control" %in% names(sys.call())) ){
    control <- list(fnscale=-1, maxit = 8e+5)
  }else{
    if(!("fnscale" %in% names(control))){
      control[["fnscale"]] <- -1
      control[["maxit"]] <- 8e+5
    }
  }
  # COMPUTE MAXIMUM LIKELIHOOD ESTIMATES
  est <- optim(param0, gev.cens.llik, method=optim.meth, control=control, data=data, u=u)

  return(list(est=est$par,
              msg=est$message,
              exc.p=sum(data>u)/length(data)))
}


# CHECK

n <- 1000

data <- rfrechet(n, 1, 1, 1)

#Shape parameter must be set to a positive value
#Location and scale must be set to reasonable small values

param0 <- c(5, 5, 10)

u <- quantile(data, probs=0.95, type=3)

res <- FitCensLLik(data, u, param0)

res
