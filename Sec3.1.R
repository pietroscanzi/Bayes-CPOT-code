#
#
# program code for section 3.1, univariate example
# adaptive MCMC, this appears in the paper
#
# 	written by
# 	Y Fan
#	7th April 2012
#
#################################################################
## The following code could be adjusted for any univariate example. The only
## necessary change is to change the function that gives the target distribution.
## You may well want to change the number of iterations in the Markov chain (niter)
## and the starting value for theta (the log of standard deviation of the proposal distribution).
## Less frequently you might want to change the target acceptance proability to something other
## than 0.44.
##
## The last lines in the program are optional and draw diagnostic plots that monitor
## the Robins Monro search.
##
## For more complicated problems that include univariate RWMH, you might want to copy and paste
## parts of this code into your program so as to estimate the variance of a proposal distribution,
## rather than having to specify that variance. (Usually, little of your code need be deleted.)
#################################################################################################
#


niter=1000                      # niter is the number of iterations in the Markov chain.
                                # Change it to the value you want.
                               
output<-rep(NA, niter)          
output[1]=1
accept=0
theta=0                         # theta=ln(sigma) is the quantity that the Robbins Monro estimates.
                                # Here it is set to 0 (sigma=1) but it should be set to your best
                                # guess of its appropriate value.
                          

calc.target<-function(x) {
#set the distribution of target distribution here,
# That is, you must change this function to give your target distribution.
# (currently the target distribution is a gamma(5,1) distribution.)
#
return(dgamma(x, 5, 1))
}

#function to update x using RWMH
updateX<- function(x, theta) {
        sigma=exp(theta)
	x.prop<-rnorm(1, mean=x, sd=sigma)
	pi.old<- calc.target(x)
	pi.new<-calc.target(x.prop)
	
	accept=min(1, pi.new/pi.old)
	is.ACC <- (runif(1)<accept)
	x<-ifelse(is.ACC, x.prop, x)
	return(list(x=x, acc=accept, is.ACC=is.ACC))
} 



##########################################
## The remainder of the code is mostly for the Robbins Monro procedure and
## gives the algorithm in the paper. (The most likely variation that
## you might want is to change the target acceptance probability from
## 0.44 to another value between 0 and 1.)
#
####################################
# Begin
# Initialisation relevant to RM
#####################################
#the optimal acceptance probability for the univariate case
pstar=0.44  
#set n0, the number of iterations before the start of search
n0=round(5/(pstar*(1-pstar)))
#iMax, is the max number of iterations before the last restart
iMax=100


Numbig=0
Numsmall=0
theta.start<-theta
theta.vec <-theta
theta.restart<-999
numRS<-2        #counts the iteration number at each restart

####################################
# End
# Initialisation relevant to RM
#####################################

i=2
#MCMC algorithm with RM
while (i <= niter) {
	dummy<-updateX(output[i-1], theta)
	output[i]<- dummy$x
	accept<-c(accept, dummy$acc)
        numRS=numRS+1	

######################################
# Begin
# updating sigma using RM section
######################################
	if (i > n0) {
		theta<- theta + (min(1, dummy$acc) - pstar)/((i+n0)*pstar*(1-pstar))
                theta.vec<-c(theta.vec, theta)				
		if (i <= iMax && (Numbig<5 || Numsmall<5)) {
			Toobig<-(exp(theta) > (3*exp(theta.start)))
			Toosmall<-(exp(theta) < (exp(theta.start)/3))	
	
		if (Toobig || Toosmall){
			#restart the algorithm
				cat("restart the program at", numRS, "th iteration", "\n")
                                theta.restart<-c(theta.restart, numRS)
				Numbig<- Numbig + Toobig
				Numsmall <- Numsmall + Toosmall
				i<-n0
				theta.start<-theta
				}	
		}
	}

######################################
# End 
# updating sigma using RM section
######################################

	i=i+1
}


##############################################################################################
## Diagnostic Plotting - The remaining lines are optional and the program will run without them
###############################################################################################
#
#
#

accept<-accept[-1]
sigma.restart<-sigma.restart[-1]
#calculate running mean acceptance rate
meanacc<-rep(NA, length(accept))
        for (i in c(1:length(accept))) {
        meanacc[i]=     mean(accept[round(i/2) :i])
            }


#begin plot
############
par(mfrow=c(2,2))
plot(output, type="l")
title("trace plot of MCMC")
hist(output, probability=T)
xgrid=seq(0, 20, length=100)
lines(xgrid, dgamma(xgrid, 5, 1))
plot(sigma.vec, type="l", col=3, ylim=c(0, max(sigma.vec)), ylab="sigma")
points(sigma.restart, sigma.vec[sigma.restart], pch="X")
plot(meanacc, type="l", col=2, ylim=c(0,1), ylab="acceptance probability")
points(sigma.restart-1, meanacc[sigma.restart-1], pch="X")
abline(h=0.44)

