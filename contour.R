par(mfrow = c(3,1))

#Frechet model

#tilde
contour(gam.fre, mu.fre, gammu.lik.fre1 - max(gammu.lik.fre1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "", lwd = 2)
#points(mle1.fre$par[1], mle1.fre$par[2], col = "blue", cex = 2, lwd = 2)
contour(gam.fre, del.fre, gamdel.lik.fre1 - max(gamdel.lik.fre1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "", lwd = 2)
#points(mle1.fre$par[1], mle1.fre$par[3], col = "blue", cex = 2, lwd = 2)
contour(mu.fre, del.fre, mudel.lik.fre1 - max(mudel.lik.fre1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "", lwd = 2)
#points(mle1.fre$par[2], mle1.fre$par[3], col = "blue", cex = 2, lwd = 2)

#bar
contour(gam.bar.fre, mu.bar.fre, gammu.lik.bar.fre - max(gammu.lik.bar.fre),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(mu)), 
        main = "", lwd = 2)
#points(mle.bar.fre$par[1], mle.bar.fre$par[2], col = "blue", cex = 2, lwd = 2)
contour(gam.bar.fre, del.bar.fre, gamdel.lik.bar.fre - max(gamdel.lik.bar.fre),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(delta)), 
        main = "", lwd = 2)
#points(mle.bar.fre$par[1], mle.bar.fre$par[3], col = "blue", cex = 2, lwd = 2)
contour(mu.bar.fre, del.bar.fre, mudel.lik.bar.fre - max(mudel.lik.bar.fre),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(mu)), ylab = expression(bar(delta)), 
        main = "", lwd = 2)
#points(mle.bar.fre$par[2], mle.bar.fre$par[3], col = "blue", cex = 2, lwd = 2)


#Exponential model

#tilde
contour(gam.exp, mu.exp, gammu.lik.exp1 - max(gammu.lik.exp1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(mu)), 
        main = "")
contour(gam.exp, del.exp, gamdel.lik.exp1 - max(gamdel.lik.exp1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(gamma)), ylab = expression(tilde(delta)), 
        main = "")
contour(mu.exp, del.exp, mudel.lik.exp1 - max(mudel.lik.exp1),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(tilde(mu)), ylab = expression(tilde(delta)), 
        main = "")

#bar
contour(gam.bar.exp, mu.bar.exp, gammu.lik.bar.exp - max(gammu.lik.bar.exp),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(mu)), 
        main = "")
contour(gam.bar.exp, del.bar.exp, gamdel.lik.bar.exp - max(gamdel.lik.bar.exp),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(gamma)), ylab = expression(bar(delta)), 
        main = "")
contour(mu.bar.exp, del.bar.exp, mudel.lik.bar.exp - max(mudel.lik.bar.exp),
        levels = -qchisq(conf.levels, 2)/2, labels = as.character(conf.levels),
        xlab = expression(bar(mu)), ylab = expression(bar(delta)), 
        main = "")

par(mfrow = c(1,1))