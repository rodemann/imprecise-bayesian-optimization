library(DiceKriging)
library(TestFunctions)


library(sjPlot)  # for plotting
library(ggplot2) # to access ggplot-themes

# load sample data set

set_theme(
  geom.outline.color = "antiquewhite4", 
  geom.outline.size = 1, 
  geom.label.size = 2,
  geom.label.color = "grey50",
  title.color = "red", 
  title.size = 1.5, 
  axis.angle.x = 45, 
  axis.textcolor = "blue", 
  base = theme_bw()
)

# -------------------------------
# A 1D example
# -------------------------------
# from Fang K.-T., Li R. and Sudjianto A. (2006), "Design and Modeling for 
# Computer Experiments", Chapman & Hall, pages 145-152

intervall = -30:30
n <- 10; d <- 1
x <- runif(n, min = -30, max= -10)
x <- append(x, runif(n, min = 10, max = 30))

x = as.vector(x)
y <- sapply(x, TF_ackley)

t <- seq(min(intervall), max(intervall), length=1000)
ground_truth = sapply(t, TF_ackley)

# add a small nugget effect, to avoid numerical problems
epsilon <- 1e-5
constant_kernel <- function(x,y){
  0.1
}



mean_fun_constant = 10000
{
model <- km(formula = ~1 , coef.trend = mean_fun_constant, 
            design=data.frame(x=x), response=data.frame(y=y), 
            covtype = "powexp", optim="gen", nugget=epsilon)

p <- predict(model, data.frame(x=t), "UK")

plot(t, p$mean, type="l", xlab="x", ylab="y", xlim = c(min(intervall)+1, max(intervall)-1), ylim = c(0,25), 
     main="GP with constant mean function")
points(x, y, col="red", pch=19)
lines(t, ground_truth, lty=2, col="blue")
legend(10, 6.5, legend=c("Sine Curve", "Sample", "Fitted Curve"), 
       pch=c(-1,16,-1), lty=c(2,-1,1), col=c("blue","red","black"))

}


