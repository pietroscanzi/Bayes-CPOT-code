#Comparison between quantile, wald and hpd intervals

#Gumbel model

load("C:\\Users\\scanz\\Desktop\\Magistrale\\Extreme Value Theory\\Code\\Gumbel model\\gumbel_simulations.RData")

#Sample size
n <- c(800, 1800, 5450, 23400)
#Number of exceedances used in the estimation
k <- c(20, 30, 50, 100)
#Number of possible simulation configurations
nconf <- length(n)

#True parameters
true.par.gum
#True extreme quantile
Q.ext.true.gum
#True return level
R.lev.true.gum

par(mfrow = c(3,2))

#Tail index gamma
true.par.gum[,1]
for(i in 1:nconf){
  plot(density(quant.ci.gum[,1,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.gum[,1,i]), median(quant.ci.gum[,1,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.gum[,2,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.gum[,2,i]), median(quant.ci.gum[,2,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.gum[,1,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.gum[,1,i]), median(wald.ci.gum[,1,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.gum[,2,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.gum[,2,i]), median(wald.ci.gum[,2,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.gum[,1,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.gum[,1,i]), median(hpd.ci.gum[,1,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.gum[,2,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(gamma ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.gum[,2,i]), median(hpd.ci.gum[,2,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.gum[1,]
wald.emp.cov.gum[1,]
hpd.emp.cov.gum[1,]

#Location parameter mu
true.par.gum[,2]
for(i in 1:nconf){
  plot(density(quant.ci.gum[,3,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.gum[,3,i]), median(quant.ci.gum[,3,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.gum[,4,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.gum[,4,i]), median(quant.ci.gum[,4,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.gum[,3,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.gum[,3,i]), median(wald.ci.gum[,3,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.gum[,4,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.gum[,4,i]), median(wald.ci.gum[,4,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.gum[,3,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.gum[,3,i]), median(hpd.ci.gum[,3,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.gum[,4,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(mu ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.gum[,4,i]), median(hpd.ci.gum[,4,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.gum[2,]
wald.emp.cov.gum[2,]
hpd.emp.cov.gum[2,]

#Scale parameter delta
true.par.gum[,3]
for(i in 1:nconf){
  plot(density(quant.ci.gum[,5,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.gum[,5,i]), median(quant.ci.gum[,5,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.gum[,6,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.gum[,6,i]), median(quant.ci.gum[,6,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.gum[,5,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.gum[,5,i]), median(wald.ci.gum[,5,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.gum[,6,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.gum[,6,i]), median(wald.ci.gum[,6,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.gum[,5,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.gum[,5,i]), median(hpd.ci.gum[,5,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.gum[,6,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(delta ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.gum[,6,i]), median(hpd.ci.gum[,6,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.gum[3,]
wald.emp.cov.gum[3,]
hpd.emp.cov.gum[3,]

#Wald intervals are stochastically smaller than quantile intervals and this 
#affect coverage probabilities.

#Extreme quantile
Q.ext.true.gum
for(i in 1:nconf){
  plot(density(quant.ci.gum[,7,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.gum[,7,i]), median(quant.ci.gum[,7,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.gum[,8,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.gum[,8,i]), median(quant.ci.gum[,8,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.gum[,7,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.gum[,7,i]), median(wald.ci.gum[,7,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.gum[,8,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.gum[,8,i]), median(wald.ci.gum[,8,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.gum[,7,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.gum[,7,i]), median(hpd.ci.gum[,7,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.gum[,8,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(Q.ext ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.gum[,8,i]), median(hpd.ci.gum[,8,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.gum[4,]
wald.emp.cov.gum[4,]
hpd.emp.cov.gum[4,]

#Return level
R.lev.true.gum
for(i in 1:nconf){
  plot(density(quant.ci.gum[,9,i]), xlab = "Lower bound (quantile)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.gum[,9,i]), median(quant.ci.gum[,9,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(quant.ci.gum[,10,i]), xlab = "Upper bound (quantile)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(quant.ci.gum[,10,i]), median(quant.ci.gum[,10,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.gum[,9,i]), xlab = "Lower bound (wald)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.gum[,9,i]), median(wald.ci.gum[,9,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(wald.ci.gum[,10,i]), xlab = "Upper bound (wald)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(wald.ci.gum[,10,i]), median(wald.ci.gum[,10,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.gum[,9,i]), xlab = "Lower bound (hpd)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.gum[,9,i]), median(hpd.ci.gum[,9,i])),
         col = c("blue", "green"), lty = 2)
  legend("topleft", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
  plot(density(hpd.ci.gum[,10,i]), xlab = "Upper bound (hpd)", ylab = "Density",
       main = bquote(R.lev ~ (n == .(n[i]))))
  abline(v = c(mean(hpd.ci.gum[,10,i]), median(hpd.ci.gum[,10,i])),
         col = c("blue", "green"), lty = 2)
  legend("topright", legend = c("mean", "median"), col = c("blue", "green"),
         lty = c(2,2), cex = 0.6)
}
#Coverage probabilities
quant.emp.cov.gum[5,]
wald.emp.cov.gum[5,]
hpd.emp.cov.gum[5,]

par(mfrow = c(1,1))