#Comparison between quantile, wald and hpd intervals

#Power-law model

load("C:\\Users\\scanz\\Desktop\\Magistrale\\Extreme Value Theory\\Code\\Power-law model\\powerlaw_simulations.RData")

#Sample size
n <- c(800, 1800, 5450, 23400)
#Number of exceedances used in the estimation
k <- c(20, 30, 50, 100)
#Number of possible simulation configurations
nconf <- length(n)

#True parameters
true.par.pl
#True extreme quantile
Q.ext.true.pl
#True return level
R.lev.true.pl

par(mfrow = c(3,2))

#Tail index ùma
true.par.pl[,1]
for(i in 1:nconf){
  plot(density(quant.ci.pl[,1,i], na.rm = T), xlab = "Lower bound (quantile)",
       ylab = "Density", main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.pl[,1,i], na.rm = T), median(quant.ci.pl[,1,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.pl[,2,i], na.rm = T), xlab = "Upper bound (quantile)",
       ylab = "Density", main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.pl[,2,i], na.rm = T), median(quant.ci.pl[,2,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.pl[,1,i], na.rm = T), xlab = "Lower bound (wald)",
       ylab = "Density", main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.pl[,1,i], na.rm = T), median(wald.ci.pl[,1,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.pl[,2,i], na.rm = T), xlab = "Upper bound (wald)",
       ylab = "Density", main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.pl[,2,i], na.rm = T), median(wald.ci.pl[,2,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.pl[,1,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(plma ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.pl[,1,i]), median(hpd.ci.pl[,1,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.pl[,2,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(plma ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.pl[,2,i]), median(hpd.ci.pl[,2,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.pl[1,]
wald.emp.cov.pl[1,]
hpd.emp.cov.pl[1,]

#Location parameter mu
true.par.pl[,2]
for(i in 1:nconf){
  plot(density(quant.ci.pl[,3,i], na.rm = T), xlab = "Lower bound (quantile)",
       ylab = "Density", main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.pl[,3,i], na.rm = T), median(quant.ci.pl[,3,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.pl[,4,i], na.rm = T), xlab = "Upper bound (quantile)",
       ylab = "Density", main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.pl[,4,i], na.rm = T), median(quant.ci.pl[,4,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.pl[,3,i], na.rm = T), xlab = "Lower bound (wald)",
       ylab = "Density", main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.pl[,3,i], na.rm = T), median(wald.ci.pl[,3,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.pl[,4,i], na.rm = T), xlab = "Upper bound (wald)",
       ylab = "Density", main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.pl[,4,i], na.rm = T), median(wald.ci.pl[,4,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.pl[,3,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.pl[,3,i]), median(hpd.ci.pl[,3,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.pl[,4,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.pl[,4,i]), median(hpd.ci.pl[,4,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.pl[2,]
wald.emp.cov.pl[2,]
hpd.emp.cov.pl[2,]

#Scale parameter delta
true.par.pl[,3]
for(i in 1:nconf){
  plot(density(quant.ci.pl[,5,i], na.rm = T), xlab = "Lower bound (quantile)",
       ylab = "Density", main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.pl[,5,i], na.rm = T), median(quant.ci.pl[,5,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.pl[,6,i], na.rm = T), xlab = "Upper bound (quantile)",
       ylab = "Density", main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.pl[,6,i], na.rm = T), median(quant.ci.pl[,6,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.pl[,5,i], na.rm = T), xlab = "Lower bound (wald)",
       ylab = "Density", main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.pl[,5,i], na.rm = T), median(wald.ci.pl[,5,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.pl[,6,i], na.rm = T), xlab = "Upper bound (wald)",
       ylab = "Density", main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.pl[,6,i], na.rm = T), median(wald.ci.pl[,6,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.pl[,5,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.pl[,5,i]), median(hpd.ci.pl[,5,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.pl[,6,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.pl[,6,i]), median(hpd.ci.pl[,6,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.pl[3,]
wald.emp.cov.pl[3,]
hpd.emp.cov.pl[3,]

#Wald intervals are stochastically smaller than quantile intervals and this 
#affect coverage probabilities.

#Extreme quantile
Q.ext.true.pl
for(i in 1:nconf){
  plot(density(quant.ci.pl[,7,i], na.rm = T), xlab = "Lower bound (quantile)",
       ylab = "Density", main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.pl[,7,i], na.rm = T), median(quant.ci.pl[,7,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.pl[,8,i], na.rm = T), xlab = "Upper bound (quantile)",
       ylab = "Density", main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.pl[,8,i], na.rm = T), median(quant.ci.pl[,8,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.pl[,7,i], na.rm = T), xlab = "Lower bound (wald)",
       ylab = "Density", main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.pl[,7,i], na.rm = T), median(wald.ci.pl[,7,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.pl[,8,i], na.rm = T), xlab = "Upper bound (wald)",
       ylab = "Density", main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.pl[,8,i], na.rm = T), median(wald.ci.pl[,8,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.pl[,7,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.pl[,7,i]), median(hpd.ci.pl[,7,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.pl[,8,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.pl[,8,i]), median(hpd.ci.pl[,8,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.pl[4,]
wald.emp.cov.pl[4,]
hpd.emp.cov.pl[4,]

#Return level
R.lev.true.pl
for(i in 1:nconf){
  plot(density(quant.ci.pl[,9,i], na.rm = T), xlab = "Lower bound (quantile)",
       ylab = "Density", main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.pl[,9,i], na.rm = T), median(quant.ci.pl[,9,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.pl[,10,i], na.rm = T), xlab = "Upper bound (quantile)",
       ylab = "Density", main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.pl[,10,i], na.rm = T), median(quant.ci.pl[,10,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.pl[,9,i], na.rm = T), xlab = "Lower bound (wald)",
       ylab = "Density", main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.pl[,9,i], na.rm = T), median(wald.ci.pl[,9,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.pl[,10,i], na.rm = T), xlab = "Upper bound (wald)",
       ylab = "Density", main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.pl[,10,i], na.rm = T), median(wald.ci.pl[,10,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.pl[,9,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.pl[,9,i]), median(hpd.ci.pl[,9,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.pl[,10,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.pl[,10,i]), median(hpd.ci.pl[,10,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.pl[5,]
wald.emp.cov.pl[5,]
hpd.emp.cov.pl[5,]

par(mfrow = c(1,1))