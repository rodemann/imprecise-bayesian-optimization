# This is a sensitivity analysis for Bayesian Optimization with Gaussian Processes 
# (Kriging) as surrogate model with regard to different kernel functions
# It uses various target functions from package smoof

# specify kernel functions that  will be compared
cov_types = c("matern5_2", "matern3_2", "exp", "powexp")
library(mlrMBO)
library(smoof)
library(doParallel)
# Functions to be optimized: 
obj.fun.list = list( makeBraninFunction(),
                     makeBBOBFunction(dimension = 4, fid = 2, iid = 3),
                     makeBBOBFunction(dimension = 2, fid = 3, iid = 6),
                     makeBBOBFunction(dimension = 6, fid = 20, iid = 2),
                     makeBBOBFunction(dimension = 7, fid = 17, iid = 4),
                     makeBBOBFunction(dimension = 3, fid = 5, iid = 16),
                     makeAckleyFunction(1),
                     makeAckleyFunction(2),
                     makeAckleyFunction(5),
                     makeBirdFunction(),
                     makeEngvallFunction(),
                     makeGiuntaFunction(),
                     makeZettlFunction(),
                     makePriceN4Function(),
                     makeCarromTableFunction(),
                     makeSchwefelFunction(dimensions = 3),
                     makeInvertedVincentFunction(dimensions = 4),
                     makeFreudensteinRothFunction(),
                     makeShekelFunction(m = 5),
                     makeKeaneFunction()
                     
)

# set initial design for each objective function
# so that the analysis is conditioned on one set of initial values
parameters_sets = lapply(obj.fun.list, getParamSet)
designs = lapply(parameters_sets, generateDesign, n= 10L, fun = lhs::maximinLHS) 
# set sequence of functions to loop over
sequence_functions = seq_along(obj.fun.list)

# set number of iterations per BO
iters = 10L
#set number of BO runs per function to be optimized 
BO_runs_per_fun = 1L

## parallelize on all except 2 available cores
cl <- parallel::makeForkCluster(detectCores() - 2)
doParallel::registerDoParallel(cl)

sim_results = foreach(i = seq_along(cov_types), .packages = c("mlrMBO", "mlr")) %dopar% {
 
mbo_mean_proposals = list()
covtype = cov_types[i]
  for(t in sequence_functions){
    # pick objective function
    obj.fun = obj.fun.list[[t]]
    # pick respective design. Recall that this is crucial for our analysis to 
    # be conditioned on the initialization 
    initial.design = designs[[t]]
    # set MBO 
    ctrl = makeMBOControl(propose.points = 1L)
    ctrl = setMBOControlTermination(ctrl, iters = iters)
    ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI(),
                               opt = "focussearch")
    # choose surrogate model including kernel function
    lrn = makeLearner("regr.km", predict.type = "se", covtype = covtype)

    # store results here  
    mbo_runs = list()
    
     for (j in 1:BO_runs_per_fun) {
      mbo_runs[[j]] = mbo(obj.fun, design = initial.design, control = ctrl, 
                          learner = lrn, show.info = FALSE)
     }
  proposed_y = lapply(mbo_runs, `[[`,"y")
  mbo_mean_proposals[[t]] = mean(unlist(proposed_y)) 
  }
mbo_mean_proposals
}

results.as.table = data.frame()
for (i in 1:5) {
  for (j in 1:20) {
    results.as.table[i,j] = sim_results[[i]][[j]]
    }
}

#save(results.as.table, file = "simulation-results.RData")

# mbo_results = lapply(mbo_runs[[3]], `[[`,"mbo.res")
# proposed_y = lapply(mbo_results, `[[`,"y")
# data.frame(Proposals = unlist(proposed_y), Mean = mean(unlist(proposed_y)))
# 
# mbo_runs[[3]][[1]][["mbo.res"]]$y







