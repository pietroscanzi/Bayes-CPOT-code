#Comparison between quantile and wald intervals

#Exponential model

load("C:\\Users\\scanz\\Desktop\\Magistrale\\Extreme Value Theory\\Code\\Exponential model\\exponential_simulations.RData")

#Sample size
n <- c(800, 1800, 5450, 23400)
#Number of exceedances used in the estimation
k <- c(20, 30, 50, 100)
#Number of possible simulation configurations
nconf <- length(n)

#True parameters
true.par.exp
#True extreme quantile
Q.ext.true.exp
#True return level
R.lev.true.exp

par(mfrow = c(3,2))

#Tail index gamma
true.par.exp[,1]
for(i in 1:nconf){
  plot(density(quant.ci.exp[,1,i], na.rm = T), xlab = "Lower bound (quantile)",
       ylab = "Density", main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.exp[,1,i], na.rm = T), median(quant.ci.exp[,1,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.exp[,2,i], na.rm = T), xlab = "Upper bound (quantile)",
       ylab = "Density", main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.exp[,2,i], na.rm = T), median(quant.ci.exp[,2,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.exp[,1,i], na.rm = T), xlab = "Lower bound (wald)",
       ylab = "Density", main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.exp[,1,i], na.rm = T), median(wald.ci.exp[,1,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.exp[,2,i], na.rm = T), xlab = "Upper bound (wald)",
       ylab = "Density", main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.exp[,2,i], na.rm = T), median(wald.ci.exp[,2,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.exp[,1,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.exp[,1,i]), median(hpd.ci.exp[,1,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.exp[,2,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.exp[,2,i]), median(hpd.ci.exp[,2,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.exp[1,]
wald.emp.cov.exp[1,]
hpd.emp.cov.exp[1,]

#Location parameter mu
true.par.exp[,2]
for(i in 1:nconf){
  plot(density(quant.ci.exp[,3,i], na.rm = T), xlab = "Lower bound (quantile)",
       ylab = "Density", main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.exp[,3,i], na.rm = T), median(quant.ci.exp[,3,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.exp[,4,i], na.rm = T), xlab = "Upper bound (quantile)",
       ylab = "Density", main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.exp[,4,i], na.rm = T), median(quant.ci.exp[,4,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.exp[,3,i], na.rm = T), xlab = "Lower bound (wald)",
       ylab = "Density", main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.exp[,3,i], na.rm = T), median(wald.ci.exp[,3,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.exp[,4,i], na.rm = T), xlab = "Upper bound (wald)",
       ylab = "Density", main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.exp[,4,i], na.rm = T), median(wald.ci.exp[,4,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.exp[,3,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.exp[,3,i]), median(hpd.ci.exp[,3,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.exp[,4,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.exp[,4,i]), median(hpd.ci.exp[,4,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.exp[2,]
wald.emp.cov.exp[2,]
hpd.emp.cov.exp[2,]

#Scale parameter delta
true.par.exp[,3]
for(i in 1:nconf){
  plot(density(quant.ci.exp[,5,i], na.rm = T), xlab = "Lower bound (quantile)",
       ylab = "Density", main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.exp[,5,i], na.rm = T), median(quant.ci.exp[,5,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.exp[,6,i], na.rm = T), xlab = "Upper bound (quantile)",
       ylab = "Density", main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.exp[,6,i], na.rm = T), median(quant.ci.exp[,6,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.exp[,5,i], na.rm = T), xlab = "Lower bound (wald)",
       ylab = "Density", main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.exp[,5,i], na.rm = T), median(wald.ci.exp[,5,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.exp[,6,i], na.rm = T), xlab = "Upper bound (wald)",
       ylab = "Density", main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.exp[,6,i], na.rm = T), median(wald.ci.exp[,6,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.exp[,5,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.exp[,5,i]), median(hpd.ci.exp[,5,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.exp[,6,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.exp[,6,i]), median(hpd.ci.exp[,6,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.exp[3,]
wald.emp.cov.exp[3,]
hpd.emp.cov.exp[3,]

#Wald intervals are stochastically smaller than quantile intervals and this 
#affect coverage probabilities.

#Extreme quantile
Q.ext.true.exp
for(i in 1:nconf){
  plot(density(quant.ci.exp[,7,i], na.rm = T), xlab = "Lower bound (quantile)",
       ylab = "Density", main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.exp[,7,i], na.rm = T), median(quant.ci.exp[,7,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.exp[,8,i], na.rm = T), xlab = "Upper bound (quantile)",
       ylab = "Density", main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.exp[,8,i], na.rm = T), median(quant.ci.exp[,8,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.exp[,7,i], na.rm = T), xlab = "Lower bound (wald)",
       ylab = "Density", main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.exp[,7,i], na.rm = T), median(wald.ci.exp[,7,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.exp[,8,i], na.rm = T), xlab = "Upper bound (wald)",
       ylab = "Density", main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.exp[,8,i], na.rm = T), median(wald.ci.exp[,8,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.exp[,7,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.exp[,7,i]), median(hpd.ci.exp[,7,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.exp[,8,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.exp[,8,i]), median(hpd.ci.exp[,8,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.exp[4,]
wald.emp.cov.exp[4,]
hpd.emp.cov.exp[4,]

#Return level
R.lev.true.exp
for(i in 1:nconf){
  plot(density(quant.ci.exp[,9,i], na.rm = T), xlab = "Lower bound (quantile)",
       ylab = "Density", main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.exp[,9,i], na.rm = T), median(quant.ci.exp[,9,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.exp[,10,i], na.rm = T), xlab = "Upper bound (quantile)",
       ylab = "Density", main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.exp[,10,i], na.rm = T), median(quant.ci.exp[,10,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.exp[,9,i], na.rm = T), xlab = "Lower bound (wald)",
       ylab = "Density", main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.exp[,9,i], na.rm = T), median(wald.ci.exp[,9,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.exp[,10,i], na.rm = T), xlab = "Upper bound (wald)",
       ylab = "Density", main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.exp[,10,i], na.rm = T), median(wald.ci.exp[,10,i], na.rm = T)),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.exp[,9,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.exp[,9,i]), median(hpd.ci.exp[,9,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.exp[,10,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.exp[,10,i]), median(hpd.ci.exp[,10,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.exp[5,]
wald.emp.cov.exp[5,]
hpd.emp.cov.exp[5,]

par(mfrow = c(1,1))