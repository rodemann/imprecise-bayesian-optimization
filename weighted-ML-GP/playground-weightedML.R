library(mlrMBO)
library(DiceKriging)

source("weighted-ML-GP/BO_weighted_ML.R")
source("_Explore_Exploit_Measures/xplxplMBO-jr.R")

fun = smoof::makeAckleyFunction(1)
parameter_set <- getParamSet(fun)
design <- generateDesign(n = 3L, par.set = parameter_set, fun = lhs::randomLHS)

# set Control Argument of BO 
# iters = Budget for classic bo
ctrl <- makeMBOControl()
ctrl <- setMBOControlTermination(ctrl, iters = 15L)
ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI(), opt = "focussearch")

lrn <- makeLearner("regr.km", predict.type = "se", optim.method = "gen", 
                           control = list(trace = FALSE), config = list(on.par.without.desc = "warn"))
lrn = setHyperPars(learner = lrn, nugget.estim = FALSE)


res_weightedML = BO_weighted_ML(fun, design, ctrl, learner = NULL, show.info = FALSE)

res_unweightedML = mbo(fun, design, ctrl, learner = NULL, show.info = FALSE)
# 
# 
# res_weightedML$res_mbo$models$`2`$learner.model
# res_unweightedML$models$`2`$learner.model
# 
# res_weightedML$xplxpl_prop_total
# 
# 
res_weightedML$res_mbo$y
res_unweightedML$y

#source("weighted-ML-GP/km.R")

#mbo(fun, design, ctrl, learner = NULL) # learner defaults to GP in case of numeric function

#debugonce(logLikFun)
#res_xplxpl <- xplxplMBO(fun, design, ctrl, learner = NULL)


## km
# a 16-points factorial design, and the corresponding response
# d <- 2; n <- 16
# design.fact <- expand.grid(x1=seq(0,1,length=4))
# y <- apply(design.fact, 1, branin)
# 
# # kriging model 1 : matern5_2 covariance structure, no trend, no nugget effect
# #debugonce(DiceKriging:::kmEstimate)
# #test = "hi"
# source("weighted-ML-GP/alter-km.R")
# debugonce(kmEstimate)
# m1 <- km(design=design.fact, response=y)
# 
