library(DiceKriging)



# -------------------------------
# A 1D example with penalized MLE
# -------------------------------
# from Fang K.-T., Li R. and Sudjianto A. (2006), "Design and Modeling for 
# Computer Experiments", Chapman & Hall, pages 145-152

n <- 8; d <- 1
x <- sample(0:20,6)
y <- sin(x)
t <- seq(0,20, length=1000)


# add a small nugget effect, to avoid numerical problems
epsilon <- 1e-3
constant_kernel <- function(x,y){
  0.1
}
model <- km(coef.trend = 5, design=data.frame(x=x), response=data.frame(y=y), 
            kernel = constant_kernel, optim="gen", nugget=epsilon)

p <- predict(model, data.frame(x=t), "UK")

plot(t, p$mean, type="l", xlab="x", ylab="y", 
     main="GP with constant kernel", ylim = c(-1,1), xlim = c(0,20))
points(x, y, col="red", pch=19)
lines(t, sin(t), lty=2, col="blue")
legend(0, -0.5, legend=c("Sine Curve", "Sample", "Fitted Curve"), 
       pch=c(-1,19,-1), lty=c(2,-1,1), col=c("blue","red","black"))


# sample from mvr normal distribution

sample = MASS::mvrnorm(n=100, mu = rep(0,5), Sigma = matrix(rep(1,25), nrow = 5) )

plot(sample)
