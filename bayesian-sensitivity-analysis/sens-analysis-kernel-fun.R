library(tidyverse)
library(mlrMBO)
library(mlr3)
library(smoof)
library(rgenoud)
library(DiceKriging)
library(foreach)
load("bayesian-sensitivity-analysis/test-functions")
# set seed
set.seed(628496)

####
# Kernel Functions
####


# set initial design size 
initial_design_size = 10


# set number of BO iterations
iters = 20L

#set number of BO runs per function to be optimized 
BO_runs_per_fun = 40L



parameters_sets = lapply(obj.fun.list, getParamSet)
all_designs <- list()
for(i in 1:BO_runs_per_fun){
all_designs[[i]] = lapply(parameters_sets, generateDesign, n= initial_design_size, fun = lhs::maximinLHS) 
}

# set seqeunce across functions
sequence_functions = seq_along(obj.fun.list)

# kernel functions
cov_types = c("gauss","matern5_2", "matern3_2", "exp", "powexp")

# store global results here
mbo_results_paths = list()

for (i in seq_along(cov_types)) {
  
  # store local (t) results here
  opt.paths = vector("list", length(obj.fun.list))
  
  # pick kernel function
  covtype = cov_types[i]
  
  for(t in sequence_functions){
    # pick objective function
    obj.fun = obj.fun.list[[t]]
    # store results here  
    mbo_runs = list()
    
    for (j in 1:BO_runs_per_fun) {
      # pick respective design (one for each function t and each BO run j)
      # this way, in each round all five configs use the same design
      initial.design <- all_designs[[j]][[t]]
      # set MBO 
      ctrl = makeMBOControl(propose.points = 1L)
      ctrl = setMBOControlTermination(ctrl, iters = iters)
      ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI(),
                                 opt = "focussearch")
      # choose surrogate model including kernel function
      lrn = makeLearner("regr.km", predict.type = "se", 
                        covtype = covtype, optim.method = "gen", 
                        control = list(trace = FALSE))
      # ensure numerical stability in km {DiceKriging} 
      y = apply(initial.design, 1, obj.fun)
      Nuggets = 1e-8*var(y)
      lrn = setHyperPars(learner = lrn, nugget=Nuggets)
      #run BO
      mbo_runs[[j]] = mbo(obj.fun, design = initial.design, control = ctrl, 
                          learner = lrn, show.info = FALSE)
    }
    # retrieve optimization path
    opt.path = lapply(mbo_runs, `[[`,"opt.path")
    opt.path.y = lapply(opt.path, getOptPathY)
    
    # retrieve current best target value
    get_current_best <- function(proposals){
      for (o in 2:length(proposals)) {
        if(proposals[o] > proposals[o - 1]) 
          proposals[o] = proposals[o - 1] 
      }
      proposals  
    }
    opt.path.best.y = lapply(opt.path.y, get_current_best)
    
    opt.paths[[t]] = opt.path.best.y
    
  }
  mbo_results_paths[[i]] = opt.paths
}  


configs = length(cov_types)

####
# Analysis of sim results
####

BO_paths_acc_difference <- vector(length = length(obj.fun.list))

for (fun in 1:length(obj.fun.list)) {
  BO_paths_all_means <- matrix(nrow = iters, ncol = configs)
  
  for (k in 1:configs) {
    BO_paths <- mbo_results_paths[[k]][[fun]]
    #BO_paths_global_min <- lapply(BO_paths, min)
    
    BO_paths <- matrix(unlist(BO_paths), ncol = length(BO_paths))
    BO_paths <- as.data.frame(BO_paths)
    
    #remove initial design
    BO_paths_initial <- BO_paths %>% slice_head(n = initial_design_size)
    BO_paths_optim <- BO_paths %>% slice_tail(n = nrow(BO_paths) - initial_design_size)
    
    # log scale
    #BO_paths_optim <- log(BO_paths_optim) 
    
    #get lower bound/upper bound/mean per iteration:
    BO_paths_sd <- apply(BO_paths_optim, 1, sd)
    BO_paths_mean_per_iter <- apply(BO_paths_optim, 1, mean)
    BO_paths_all_means[,k] <- BO_paths_mean_per_iter
  }
  BO_paths_min <- apply(BO_paths_all_means, 1, min)
  BO_paths_max <- apply(BO_paths_all_means, 1, max)
  BO_paths_acc_difference[fun] <- sum(BO_paths_max - BO_paths_min)
}
BO_paths_acc_difference_kernel_fun <- BO_paths_acc_difference
BO_paths_acc_difference_kernel_fun
save(BO_paths_acc_difference_kernel_fun, file = "bayesian-sensitivity-analysis/acc-diff-kernel-fun")

save(BO_paths_acc_difference_kernel_fun, file = "bayesian-sensitivity-analysis/acc-diff-kernel-fun")
save(BO_paths_acc_difference_kernel_fun, file = "bayesian-sensitivity-analysis/acc-diff-kernel-fun")





