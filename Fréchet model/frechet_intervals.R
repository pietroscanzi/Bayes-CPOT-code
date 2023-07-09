#Comparison between quantile, wald and hpd intervals

#Frechet model

load("C:\\Users\\scanz\\Desktop\\Magistrale\\Extreme Value Theory\\Code\\Fréchet model\\frechet_simulations.RData")

#Sample size
n <- c(800, 1800, 5450, 23400)
#Number of exceedances used in the estimation
k <- c(20, 30, 50, 100)
#Number of possible simulation configurations
nconf <- length(n)

#True parameters
true.par.fre
#True extreme quantile
Q.ext.true.fre
#True return level
R.lev.true.fre

par(mfrow = c(3,2))

#Tail index gamma
true.par.fre[,1]
for(i in 1:nconf){
  plot(density(quant.ci.fre[,1,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.fre[,1,i]), median(quant.ci.fre[,1,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
                               lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.fre[,2,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.fre[,2,i]), median(quant.ci.fre[,2,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.fre[,1,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.fre[,1,i]), median(wald.ci.fre[,1,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.fre[,2,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.fre[,2,i]), median(wald.ci.fre[,2,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.fre[,1,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.fre[,1,i]), median(hpd.ci.fre[,1,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.fre[,2,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.fre[,2,i]), median(hpd.ci.fre[,2,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.fre[1,]
wald.emp.cov.fre[1,]
hpd.emp.cov.fre[1,]

#Location parameter mu
true.par.fre[,2]
for(i in 1:nconf){
  plot(density(quant.ci.fre[,3,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.fre[,3,i]), median(quant.ci.fre[,3,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.fre[,4,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.fre[,4,i]), median(quant.ci.fre[,4,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.fre[,3,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.fre[,3,i]), median(wald.ci.fre[,3,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.fre[,4,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.fre[,4,i]), median(wald.ci.fre[,4,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.fre[,3,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.fre[,3,i]), median(hpd.ci.fre[,3,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.fre[,4,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.fre[,4,i]), median(hpd.ci.fre[,4,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.fre[2,]
wald.emp.cov.fre[2,]
hpd.emp.cov.fre[2,]

#Scale parameter delta
true.par.fre[,3]
for(i in 1:nconf){
  plot(density(quant.ci.fre[,5,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.fre[,5,i]), median(quant.ci.fre[,5,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.fre[,6,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.fre[,6,i]), median(quant.ci.fre[,6,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.fre[,5,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.fre[,5,i]), median(wald.ci.fre[,5,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.fre[,6,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.fre[,6,i]), median(wald.ci.fre[,6,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.fre[,5,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.fre[,5,i]), median(hpd.ci.fre[,5,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.fre[,6,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.fre[,6,i]), median(hpd.ci.fre[,6,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.fre[3,]
wald.emp.cov.fre[3,]
hpd.emp.cov.fre[3,]

#Wald intervals are stochastically smaller than quantile intervals and this 
#affect coverage probabilities.

#Extreme quantile
Q.ext.true.fre
for(i in 1:nconf){
  plot(density(quant.ci.fre[,7,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.fre[,7,i]), median(quant.ci.fre[,7,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.fre[,8,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.fre[,8,i]), median(quant.ci.fre[,8,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.fre[,7,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.fre[,7,i]), median(wald.ci.fre[,7,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.fre[,8,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.fre[,8,i]), median(wald.ci.fre[,8,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.fre[,7,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.fre[,7,i]), median(hpd.ci.fre[,7,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.fre[,8,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.fre[,8,i]), median(hpd.ci.fre[,8,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.fre[4,]
wald.emp.cov.fre[4,]
hpd.emp.cov.fre[4,]

#Return level
R.lev.true.fre
for(i in 1:nconf){
  plot(density(quant.ci.fre[,9,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.fre[,9,i]), median(quant.ci.fre[,9,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.fre[,10,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.fre[,10,i]), median(quant.ci.fre[,10,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.fre[,9,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.fre[,9,i]), median(wald.ci.fre[,9,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.fre[,10,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.fre[,10,i]), median(wald.ci.fre[,10,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.fre[,9,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.fre[,9,i]), median(hpd.ci.fre[,9,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.fre[,10,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.fre[,10,i]), median(hpd.ci.fre[,10,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.fre[5,]
wald.emp.cov.fre[5,]
hpd.emp.cov.fre[5,]

par(mfrow = c(1,1))