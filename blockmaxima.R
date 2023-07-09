#HURDAT2
#BAYESIAN BLOCK MAXIMA

setwd("C:\\Users\\scanz\\Desktop\\Magistrale\\Extreme Value Theory")
load("AL.rda")
hur <- read.csv("hurdat2-1851-2020-020922.txt", header = F)

#Wind
hur.w <- as.integer(na.omit(hur$V7))
hur.w <- hur.w[which(hur.w != -99)]
length(hur.w)
summary(hur.w)
wind <- as.integer(na.omit(AL$Wind))
length(wind)
summary(wind)
#Sembrano essere le stesse osservazioni del dataset AL con massimo valore
#per wind pari a 165 --> in contrasto con il paper.

#Year 
year <- substr(as.Date(AL$DateTime), 1, 4)
str(year)
idx <- unique(year)
wind.max <- rep(NA, length(idx))
for(i in 1:length(idx)){
  wind.max[i] <- max(wind[year == idx[i]], na.rm = TRUE)
}

#Block maxima: year maximum wind speed
wind.max
summary(wind.max)
n <- length(wind.max)
plot(idx, wind.max, type = "l", xlab = "Year", ylab = "Wind speed")
hist(wind.max, nclass = 15, probability = T)
lines(density(wind.max), col = "red", lwd = 2)