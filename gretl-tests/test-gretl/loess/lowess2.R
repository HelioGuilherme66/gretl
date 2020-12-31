sindat <- read.table("sin_noise.csv", header=TRUE)
attach(sindat)
ymod <- loess(y ~ x, span=0.75)
yhat <- predict(ymod, x)
yhat

