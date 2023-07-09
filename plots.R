#-------------------------------------------------------------------------------

#Plots of GEV and GP distributions

#-------------------------------------------------------------------------------

path <- "C:\\Users\\scanz\\OneDrive - Università Commerciale Luigi Bocconi\\Desktop/Magistrale\\Extreme Value Theory\\Code\\"
setwd(path)
source(paste0(path ,"functions.R"))

#-------------------------------------------------------------------------------

#Generalized Extreme Value Distribution plots

#GEV cdf
plot(function(x) pgev.V(x, par = c(1, 0, 1), log = FALSE),
     from = -10, to = 10, col = 2, lty = 2, lwd = 4.5, ylim = c(0, 1.1),
     main = "", xlab = "", ylab = "")
plot(function(x) pgev.V(x, par = c(0, 0, 1), log = FALSE),
     from = -10, to = 10, col = 4, lty = 3, lwd = 4.5, add = TRUE)
plot(function(x) pgev.V(x, par = c(-1/3, 0, 1), log = FALSE),
     from = -10, to = 10, col = 7, lty = 4, lwd = 4.5, add = TRUE)
grid()
legend("topleft", legend = c(expression(gamma == 1),
                             expression(gamma == 0),
                             expression(gamma == -1/3)),
                             lty = c(2,3,4), lwd = rep(4.5,3), col = c(2, 4, 7),
       bty = "n", cex = 1.5)

#GEV pdf
plot(function(x) dgev.V(x, par = c(1, 0, 1), log = FALSE),
     from = -10, to = 10, col = 2, lty = 2, lwd = 4.5, ylim = c(0, 0.6),
     main = "", xlab = "", ylab = "")
plot(function(x) dgev.V(x, par = c(0, 0, 1), log = FALSE),
     from = -10, to = 10, col = 4, lty = 3, lwd = 4.5, add = TRUE)
plot(function(x) dgev.V(x, par = c(-1/3, 0, 1), log = FALSE),
     from = -10, to = 10, col = 7, lty = 4, lwd = 4.5, add = TRUE)
grid()
legend("topleft", legend = c(expression(gamma == 1),
                             expression(gamma == 0),
                             expression(gamma == -1/3)),
       lty = c(2,3,4), lwd = rep(4.5,3), col = c(2, 4, 7),
       bty = "n", cex = 1.5)

#-------------------------------------------------------------------------------

#Generalized Pareto distribution plots

#GEV cdf
plot(function(x) pgp.V(x, par = c(1, 1), log = FALSE),
     from = -2, to = 15, col = 2, lty = 2, lwd = 4.5, ylim = c(0, 1.1),
     main = "", xlab = "", ylab = "")
plot(function(x) pgp.V(x, par = c(0, 1), log = FALSE),
     from = -2, to = 15, col = 4, lty = 3, lwd = 4.5, add = TRUE)
plot(function(x) pgp.V(x, par = c(-1/3, 1), log = FALSE),
     from = -2, to = 15, col = 7, lty = 4, lwd = 4.5, add = TRUE)
grid()
legend(10, 0.3, legend = c(expression(gamma == 1),
                             expression(gamma == 0),
                             expression(gamma == -1/3)),
       lty = c(2,3,4), lwd = rep(4.5,3), col = c(2, 4, 7),
       bty = "n", cex = 1.5)

#GEV pdf
plot(function(x) dgp.V(x, par = c(1, 1), log = FALSE),
     from = -2, to = 15, col = 2, lty = 2, lwd = 4.5, ylim = c(0, 1),
     main = "", xlab = "", ylab = "")
plot(function(x) dgp.V(x, par = c(0, 1), log = FALSE),
     from = -2, to = 15, col = 4, lty = 3, lwd = 4.5, add = TRUE)
plot(function(x) dgp.V(x, par = c(-1/3, 1), log = FALSE),
     from = -2, to = 15, col = 7, lty = 4, lwd = 4.5, add = TRUE)
grid()
legend(10, 1, legend = c(expression(gamma == 1),
                             expression(gamma == 0),
                             expression(gamma == -1/3)),
       lty = c(2,3,4), lwd = rep(4.5,3), col = c(2, 4, 7),
       bty = "n", cex = 1.5)

#-------------------------------------------------------------------------------