

library(mlrMBO)
library(ggplot2)
library(xtable)

source("_Explore_Exploit_Measures/xplxpl-jr.R")


# visualize MBO on Branin Function to explain difference of var-based and distance-based xplxpl measures

obj.fun.list = list( makeBraninFunction(),
makeBBOBFunction(dimension = 4, fid = 2, iid = 3),
makeBBOBFunction(dimension = 2, fid = 3, iid = 6),
makeBBOBFunction(dimension = 6, fid = 20, iid = 2),
makeBBOBFunction(dimension = 7, fid = 17, iid = 4),
makeBBOBFunction(dimension = 3, fid = 5, iid = 16),
makeAckleyFunction(3),
makeAckleyFunction(5),
makeBirdFunction(),
makeEngvallFunction(),
makeGiuntaFunction(),
makeZettlFunction(),
makePriceN4Function(),
makeCarromTableFunction(),
makeBraninFunction(),
makeSchwefelFunction(dimensions = 3),
makeInvertedVincentFunction(dimensions = 4),
makeFreudensteinRothFunction(),
makeShekelFunction(m = 5),
makeKeaneFunction()

)


iters = 40
ctrl = makeMBOControl(y.name = "target", propose.points = 1L, store.model.at = 1:(iters+1))
ctrl = setMBOControlTermination(ctrl, iters = iters)
ctrl = setMBOControlInfill(ctrl, opt = "focussearchSavepts", crit = makeMBOInfillCritEI())

runs = length(obj.fun.list)
cor.matrix.list = matrix(NA, nrow = runs, ncol = 6*6)

for (run in 1:runs) {
obj.fun = obj.fun.list[[run]]
res = mbo(obj.fun, control = ctrl, show.info = FALSE)
xpl.xpl.res = xplxpl(res)
xpl.xpl.res = xpl.xpl.res[,-7]

cor.matrix.list[run,] = as.vector(cor(xpl.xpl.res))

}

mean.cor = apply(cor.matrix.list, 2, mean, na.rm = TRUE)
mean.cor = matrix(mean.cor, nrow = 6)
colnames(mean.cor) = c("Standard Error Ratio", "Standard Error Distribution Value", "Pairwise Distance", "Minimal Distance",
                       "Maximal Distance", "Mean Distance")
rownames(mean.cor) = c("Standard Error Ratio", "Standard Error Distribution Value", "Pairwise Distance", "Minimal Distance",
                       "Maximal Distance", "Mean Distance")
mean.cor = mean.cor * lower.tri(mean.cor)
xtable(mean.cor)

