source("imprecise-bayes-opt-parallel/imp_BO_univariate.R")
source("imprecise-bayes-opt-parallel/imp_BO_univariate_batch.R")

library(smoof)

obj.fun.list = list(makeAckleyFunction(1),
                    makeCosineMixtureFunction(1)
                    
)

fun <- obj.fun.list[[1]]
parameter_set <- getParamSet(fun)
design <- generateDesign(n = 10L, par.set = parameter_set, fun = lhs::randomLHS)
iters = 10L
# set Control Argument of BO 
ctrl = makeMBOControl(propose.points = 1L)
ctrl = setMBOControlTermination(ctrl, iters = iters)
ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritMeanResponse(), opt = "focussearch",
                           opt.focussearch.points = 200, opt.focussearch.maxit = 1)

# classic bo with default learner (GP)
cl_bo_res <- mbo(fun, design, ctrl, learner = NULL)


parallel_imp_bo_res <- imp_BO_univariate(fun = fun, design = design, control = ctrl, 
                  base_kernel = "powexp", imprecision_degree = 100)

batch_imp_bo_res <- imp_BO_univariate_batch(fun = fun, design = design, control = ctrl, 
                  base_kernel = "powexp", imprecision_degree = c(1,10,100,1000),
                  number_of_batches = 2, iters_per_batch = c(2,2))

cl_bo_res$y
parallel_imp_bo_res
batch_imp_bo_res$global.optimum





    