#Comparison between quantile, wald and hpd intervals

#Reverse Weibull model

load("C:\\Users\\scanz\\Desktop\\Magistrale\\Extreme Value Theory\\Code\\Reverse Weibull model\\rweibull_simulations.RData")

#Sample size
n <- c(800, 1800, 5450, 23400)
#Number of exceedances used in the estimation
k <- c(20, 30, 50, 100)
#Number of possible simulation configurations
nconf <- length(n)

#True parameters
true.par.rw
#True extreme quantile
Q.ext.true.rw
#True return level
R.lev.true.rw

par(mfrow = c(3,2))

#Tail index gamma
true.par.rw[,1]
for(i in 1:nconf){
  plot(density(quant.ci.rw[,1,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.rw[,1,i]), median(quant.ci.rw[,1,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.rw[,2,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.rw[,2,i]), median(quant.ci.rw[,2,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.rw[,1,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.rw[,1,i]), median(wald.ci.rw[,1,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.rw[,2,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.rw[,2,i]), median(wald.ci.rw[,2,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.rw[,1,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.rw[,1,i]), median(hpd.ci.rw[,1,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.rw[,2,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.rw[,2,i]), median(hpd.ci.rw[,2,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.rw[1,]
wald.emp.cov.rw[1,]
hpd.emp.cov.rw[1,]

#Location parameter mu
true.par.rw[,2]
for(i in 1:nconf){
  plot(density(quant.ci.rw[,3,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.rw[,3,i]), median(quant.ci.rw[,3,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.rw[,4,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.rw[,4,i]), median(quant.ci.rw[,4,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.rw[,3,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.rw[,3,i]), median(wald.ci.rw[,3,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.rw[,4,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.rw[,4,i]), median(wald.ci.rw[,4,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.rw[,3,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.rw[,3,i]), median(hpd.ci.rw[,3,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.rw[,4,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.rw[,4,i]), median(hpd.ci.rw[,4,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.rw[2,]
wald.emp.cov.rw[2,]
hpd.emp.cov.rw[2,]

#Scale parameter delta
true.par.rw[,3]
for(i in 1:nconf){
  plot(density(quant.ci.rw[,5,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.rw[,5,i]), median(quant.ci.rw[,5,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.rw[,6,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.rw[,6,i]), median(quant.ci.rw[,6,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.rw[,5,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.rw[,5,i]), median(wald.ci.rw[,5,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.rw[,6,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.rw[,6,i]), median(wald.ci.rw[,6,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.rw[,5,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.rw[,5,i]), median(hpd.ci.rw[,5,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.rw[,6,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.rw[,6,i]), median(hpd.ci.rw[,6,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.rw[3,]
wald.emp.cov.rw[3,]
hpd.emp.cov.rw[3,]

#Wald intervals are stochastically smaller than quantile intervals and this 
#affect coverage probabilities.

#Extreme quantile
Q.ext.true.rw
for(i in 1:nconf){
  plot(density(quant.ci.rw[,7,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.rw[,7,i]), median(quant.ci.rw[,7,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.rw[,8,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.rw[,8,i]), median(quant.ci.rw[,8,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.rw[,7,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.rw[,7,i]), median(wald.ci.rw[,7,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.rw[,8,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.rw[,8,i]), median(wald.ci.rw[,8,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.rw[,7,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.rw[,7,i]), median(hpd.ci.rw[,7,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.rw[,8,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.rw[,8,i]), median(hpd.ci.rw[,8,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.rw[4,]
wald.emp.cov.rw[4,]
hpd.emp.cov.rw[4,]

#Return level
R.lev.true.rw
for(i in 1:nconf){
  plot(density(quant.ci.rw[,9,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.rw[,9,i]), median(quant.ci.rw[,9,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.rw[,10,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.rw[,10,i]), median(quant.ci.rw[,10,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.rw[,9,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.rw[,9,i]), median(wald.ci.rw[,9,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.rw[,10,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.rw[,10,i]), median(wald.ci.rw[,10,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.rw[,9,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.rw[,9,i]), median(hpd.ci.rw[,9,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.rw[,10,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.rw[,10,i]), median(hpd.ci.rw[,10,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.rw[5,]
wald.emp.cov.rw[5,]
hpd.emp.cov.rw[5,]

par(mfrow = c(1,1))