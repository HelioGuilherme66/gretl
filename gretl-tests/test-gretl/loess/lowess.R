period <- 120
x <- 1:120
y <- sin(2*pi*x/period) + runif(length(x),-1,1)
plot(x,y, main="Sine Curve + 'Uniform' Noise")
mtext("showing loess smoothing (local regression smoothing)")
y.loess <- loess(y ~ x, span=0.75, data.frame(x=x, y=y))
y.predict <- predict(y.loess, data.frame(x=x))
lines(x,y.predict)
