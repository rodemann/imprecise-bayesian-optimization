library(tidyverse)
library(mlrMBO)
library(mlr3)
library(smoof)
library(rgenoud)
library(DiceKriging)

load("bayesian-sensitivity-analysis/test-functions")

####
# Mean Function Pars
####
# set seed
set.seed(628496)


# set number of BO iterations
iters = 20L

#set number of BO runs per function to be optimized 
BO_runs_per_fun = 40L

initial_design_size = 10L

# set initial design for each objective function and each bo run
# and each BO run
parameters_sets = lapply(obj.fun.list, getParamSet)
all_designs <- list()
for(i in 1:BO_runs_per_fun){
  all_designs[[i]] = lapply(parameters_sets, generateDesign, n= initial_design_size, fun = lhs::maximinLHS) 
}
# ground truth to get true mean: 
ground_truths = lapply(parameters_sets, generateDesign, n = 400L, fun = lhs::maximinLHS) 

# set seqeunce across functions
sequence_functions = seq_along(obj.fun.list)


# pick kernel function
covtype = "gauss"

# store global results here
mbo_results_paths = list()

# number of mean configs
configs = 5

# mean functions
#mean_constants = c(0, 1, 5, 10, 100)

for (i in seq.int(1,configs)) {
  
  # store local (t) results here
  opt.paths = vector("list", length(obj.fun.list))
  
  for(t in sequence_functions){
    # pick objective function
    obj.fun = obj.fun.list[[t]]
    # get true y 
    ground_truth = ground_truths[[t]]
    y_true = apply(ground_truth, 1, obj.fun)
    mean_true = mean(y_true)
    # mean functions
    mean_constants = c(mean_true - abs(3*mean_true), mean_true - abs(0.2*mean_true), mean_true, mean_true + abs(0.2*mean_true), mean_true + abs(3*mean_true))
    
    
    
    # store results here  
    mbo_runs = list()
    
    for (j in 1:BO_runs_per_fun) {
      # pick respective design  
      initial.design = all_designs[[j]][[t]]
      # set MBO 
      ctrl = makeMBOControl(propose.points = 1L)
      ctrl = setMBOControlTermination(ctrl, iters = iters)
      ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI(),
                                 opt = "focussearch")
      
      
      ## choose surrogate model including mean function
      mean_fun_constant = mean_constants[i]
      
      lrn = makeLearner("regr.km", formula = ~1, 
                        coef.trend = mean_fun_constant, 
                        predict.type = "se", 
                        covtype = covtype, 
                        optim.method = "gen", 
                        control = list(trace = FALSE), 
                        config = list(on.par.without.desc = "warn"))
      
      # ensure numerical stability in km {DiceKriging} cf. github issue and recommendation by Bernd Bischl 
      y = apply(initial.design, 1, obj.fun)
      Nuggets = 1e-8*var(y)
      lrn = setHyperPars(learner = lrn, nugget=Nuggets)
      
      
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

BO_paths_acc_difference_mean_par <- BO_paths_acc_difference
BO_paths_acc_difference_mean_par
save(BO_paths_acc_difference_mean_par, file = "bayesian-sensitivity-analysis/acc-diff-mean-par")







