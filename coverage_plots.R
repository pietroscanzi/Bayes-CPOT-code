par(mfrow = c(1,3))

quant.emp.cov.fre
min(quant.emp.cov.fre); max(quant.emp.cov.fre)
plot(c(40,60,109,234), quant.emp.cov.fre[1,], type = "l", col = 2, lwd = 2, 
     xlab = "s", ylab = "Coperture empiriche al 95%", main = "",
     ylim = c(0.8, 1)) 
abline(h = 0.95, lwd = 2, lty = 2, col = "blue")
lines(c(40,60,109,234), quant.emp.cov.fre[2,], type = "l", col = 3, lwd = 2,
     add = TRUE)
lines(c(40,60,109,234), quant.emp.cov.fre[3,], type = "l", col = "cyan", lwd = 2,
      add = TRUE)
lines(c(40,60,109,234), quant.emp.cov.fre[4,], type = "l", col = "orange", lwd = 2,
      add = TRUE)
lines(c(40,60,109,234), quant.emp.cov.fre[5,], type = "l", col = "grey", lwd = 2,
      add = TRUE)

wald.emp.cov.exp
min(wald.emp.cov.exp); max(wald.emp.cov.exp)
plot(c(40,60,109,234), wald.emp.cov.exp[1,], type = "l", col = 2, lwd = 2, 
     xlab = "s", ylab = "Coperture empiriche al 95%", main = "",
     ylim = c(0.8, 1)) 
abline(h = 0.95, lwd = 2, lty = 2, col = "blue")
lines(c(40,60,109,234), wald.emp.cov.exp[2,], type = "l", col = 3, lwd = 2,
      add = TRUE)
lines(c(40,60,109,234), wald.emp.cov.exp[3,], type = "l", col = "cyan", lwd = 2,
      add = TRUE)
lines(c(40,60,109,234), wald.emp.cov.exp[4,], type = "l", col = "orange", lwd = 2,
      add = TRUE)
lines(c(40,60,109,234), wald.emp.cov.exp[5,], type = "l", col = "grey", lwd = 2,
      add = TRUE)

hpd.emp.cov.pl
min(hpd.emp.cov.pl); max(hpd.emp.cov.pl)
plot(c(40,60,109,234), hpd.emp.cov.pl[1,], type = "l", col = 2, lwd = 2, 
     xlab = "s", ylab = "Coperture empiriche al 95%", main = "",
     ylim = c(0.8, 1)) 
abline(h = 0.95, lwd = 2, lty = 2, col = "blue")
lines(c(40,60,109,234), hpd.emp.cov.pl[2,], type = "l", col = 3, lwd = 2,
      add = TRUE)
lines(c(40,60,109,234), hpd.emp.cov.pl[3,], type = "l", col = "cyan", lwd = 2,
      add = TRUE)
lines(c(40,60,109,234), hpd.emp.cov.pl[4,], type = "l", col = "orange", lwd = 2,
      add = TRUE)
lines(c(40,60,109,234), hpd.emp.cov.pl[5,], type = "l", col = "grey", lwd = 2,
      add = TRUE)

par(mfrow = c(1,1))
