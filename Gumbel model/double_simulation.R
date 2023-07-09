#Double simulation code

lab.cov <- c("gamma", "mu", "delta", "Q.ext", "R.lev")
lab.sim <- c("n = 800", "n = 1800", "n = 5450", "n = 23400")

quant.ci.gum.tot <- array(0, dim = c(1000,10,4))
quant.ci.gum.tot[1:500,,] <- quant.ci.gum1
quant.ci.gum.tot[501:1000,,] <- quant.ci.gum
rm(quant.ci.gum1, quant.ci.gum)
quant.ci.gum <- quant.ci.gum.tot
rm(quant.ci.gum.tot)

wald.ci.gum.tot <- array(0, dim = c(1000,10,4))
wald.ci.gum.tot[1:500,,] <- wald.ci.gum1
wald.ci.gum.tot[501:1000,,] <- wald.ci.gum
rm(wald.ci.gum1, wald.ci.gum)
wald.ci.gum <- wald.ci.gum.tot
rm(wald.ci.gum.tot)

quant.cov.gum.tot <- array(0, dim = c(1000,5,4))
quant.cov.gum.tot[1:500,,] <- quant.cov.gum1
quant.cov.gum.tot[501:1000,,] <- quant.cov.gum
rm(quant.cov.gum1, quant.cov.gum)
quant.cov.gum <- quant.cov.gum.tot
rm(quant.cov.gum.tot)

wald.cov.gum.tot <- array(0, dim = c(1000,5,4))
wald.cov.gum.tot[1:500,,] <- wald.cov.gum1
wald.cov.gum.tot[501:1000,,] <- wald.cov.gum
rm(wald.cov.gum1, wald.cov.gum)
wald.cov.gum <- wald.cov.gum.tot
rm(wald.cov.gum.tot)

quant.emp.cov.gum <- colMeans(quant.cov.gum)
row.names(quant.emp.cov.gum) <- lab.cov
colnames(quant.emp.cov.gum) <- lab.sim
rm(quant.emp.cov.gum1)
wald.emp.cov.gum <- colMeans(wald.cov.gum)
row.names(wald.emp.cov.gum) <- lab.cov
colnames(wald.emp.cov.gum) <- lab.sim
rm(wald.emp.cov.gum1)

save.image("C:/Users/scanz/Desktop/Magistrale/Extreme Value Theory/Code/Gumbel model/gumbel_simulations.RData")

xtable(quant.emp.cov.gum)
xtable(wald.emp.cov.gum)
