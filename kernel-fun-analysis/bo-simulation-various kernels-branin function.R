# This is a sensitivity analysis for Bayesian Optimization with Gaussian Processes 
# (Kriging) as surrogate model with regard to different kernel functions

library(mlrMBO)
library(smoof)
library(doParallel)

cov_types = c("gauss", "matern5_2", "matern3_2", "exp", "powexp")
mbo_runs = list()


###########
# branin 2D
###########

set.seed(1)
configureMlr(show.learner.output = FALSE)
obj.fun = makeBraninFunction()
ctrl = makeMBOControl(propose.points = 1L)
ctrl = setMBOControlTermination(ctrl, iters = 20L)
ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI(),
                           opt = "focussearch", opt.focussearch.points = 2000L)

for(i in seq_along(cov_types)){
  
lrn = makeLearner("regr.km", predict.type = "se", covtype = cov_types[i])

mbo_runs[[i]] = list()

for (j in 1:10) {
  
design = generateDesign(10L, getParamSet(obj.fun), fun = lhs::maximinLHS)
mbo_runs[[i]][[j]] = exampleRun(obj.fun, design = design, learner = lrn, control = ctrl,
                 points.per.dim = 50L, show.info = TRUE)
}
}

mbo_results = lapply(mbo_runs[[3]], `[[`,"mbo.res")
proposed_y = lapply(mbo_results, `[[`,"y")
data.frame(Proposals = unlist(proposed_y), Mean = mean(unlist(proposed_y)))


mbo_runs[[3]][[1]][["mbo.res"]]$y

