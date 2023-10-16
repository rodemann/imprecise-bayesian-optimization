library(tidyverse)
library(mlrMBO)
library(mlr3)
library(smoof)
library(rgenoud)
library(DiceKriging)

load("bayesian-sensitivity-analysis/test-functions")
# set seed
set.seed(628496)

####
# Mean Function Trends
####

# set seqeunce across functions
sequence_functions = seq_along(obj.fun.list)

# set number of BO iterations
iters = 20L

#set number of BO runs per function to be optimized 
BO_runs_per_fun = 40L

# set initial design for each objective function
# and each BO run
initial_design_size = 10L
parameters_sets = lapply(obj.fun.list, getParamSet)
all_designs <- list()
for(i in 1:BO_runs_per_fun){
  all_designs[[i]] = lapply(parameters_sets, generateDesign, n= initial_design_size, fun = lhs::maximinLHS) 
}


# pick kernel function
covtype = "powexp"

# store global results here
mbo_results_paths = list()

# number of different configurations
configs = 4

for (i in seq.int(1,configs)) {
  print(i)
  # store local (t) results here
  opt.paths = vector("list", length(obj.fun.list))
  
  for(t in sequence_functions){
    print(t)
    # pick objective function
    obj.fun = obj.fun.list[[t]]
   
    # pick parameter set
    parameters = parameters_sets[[t]]
    
    # specify trends according to dim of parameter space
    number_of_params <- parameters$pars$x$len
    switch (as.character(number_of_params),
            "1" = {trends = c(as.formula("y~1"), as.formula("y~x"), as.formula("y~x^2"), as.formula("y~x^3"))},
            "2" = {trends = c(as.formula("y~1"), as.formula("y~x1+x2"), as.formula("y~x1^2+x2^2"), as.formula("y~x1^3+x2^3"))},
            "3" = {trends = c(as.formula("y~1"), as.formula("y~x1+x2+x3"), as.formula("y~x1^2+x2^2+x3^2"), as.formula("y~x1^3+x2^3+x3^3"))},
            "4" = {trends = c(as.formula("y~1"), as.formula("y~x1+x2+x3+x4"), as.formula("y~x1^2+x2^2+x3^2+x4^2"), as.formula("y~x1^3+x2^3+x3^3+x4^3"))},
            "5" = {trends = c(as.formula("y~1"), as.formula("y~x1+x2+x3+x4+x5"), as.formula("y~x1^2+x2^2+x3^2+x4^2+x5^2"), as.formula("y~x1^3+x2^3+x3^3+x4^3+x5^3"))},
            "6" = {trends = c(as.formula("y~1"), as.formula("y~x1+x2+x3+x4+x5+x6"), as.formula("y~x1^2+x2^2+x3^2+x4^2+x5^2+x6^2"), as.formula("y~x1^3+x2^3+x3^3+x4^3+x5^3+x6^3"))},
            "7" = {trends = c(as.formula("y~1"), as.formula("y~x1+x2+x3+x4+x5+x6+x7"), as.formula("y~x1^2+x2^2+x3^2+x4^2+x5^2+x6^2+x7^2"), as.formula("y~x1^3+x2^3+x3^3+x4^3+x5^3+x6^3+x7^3"))}
    )
    
    #select trend
    trend = trends[[i]]
    
   
    # store results here  
    mbo_runs = list()
    
    for (j in 1:BO_runs_per_fun) {
      
      #choose initial design
      initial.design = all_designs[[j]][[t]]
      # set MBO 
      ctrl = makeMBOControl(propose.points = 1L)
      ctrl = setMBOControlTermination(ctrl, iters = iters)
      ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI(),
                                 opt = "focussearch")
      lrn = makeLearner("regr.km", formula = trend, 
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
BO_paths_acc_difference_mean_fun <- BO_paths_acc_difference
BO_paths_acc_difference_mean_fun
save(BO_paths_acc_difference_mean_fun, file = "bayesian-sensitivity-analysis/acc-diff-mean-fun-powexp-kernel")







