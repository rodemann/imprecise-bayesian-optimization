library(mlrMBO)
library(ggplot2)
library(xtable)

source("_Explore_Exploit_Measures/xplxpl-jr.R")


# visualize MBO on Branin Function to explain difference of var-based and distance-based xplxpl measures

obj.fun = makeBraninFunction()
# visualize the function
autoplot(obj.fun, render.levels = TRUE, show.optimum = TRUE)




iters = 20

ctrl = makeMBOControl(y.name = "target", propose.points = 1L, store.model.at = 1:(iters+1))
ctrl = setMBOControlTermination(ctrl, iters = iters)
ctrl = setMBOControlInfill(ctrl, opt = "focussearchSavepts", crit = makeMBOInfillCritEI())


res = mbo(obj.fun, control = ctrl, show.info = FALSE)
xplxpl.res = xplxpl(res)
xtable(xplxpl.res)

# visualize
autoplot(obj.fun, render.levels = TRUE, show.optimum = TRUE) + geom_text(data = as.data.frame(res$opt.path), mapping = aes(label = dob), color = "white")
