library(tidyr)
library(mlrMBO)
library(DiceKriging)

source("weighted-ML-GP/BO_weighted_ML.R")
# unable noisy target functions and nugget stability (because they affect the usage of loglikfun heavily)
source("weighted-ML-GP/makeMBOLearner-no-noise.R")
source("_Explore_Exploit_Measures/xplxplMBO-jr.R")
source("univariate-benchmark-functions/get-synthetic-benchmark-funs-for-minimization.R")

funs_to_test <- all_benchmark_funs[c(27)]
res_opt_paths_list_w <- list()
res_opt_paths_list_uw <- list()
global_res_weightedML <- list()
global_res_unweightedML <- list()

# hyperparams
n <- 2000
iters <- 20
init_design <- 4

for (i in 1:length(funs_to_test)) {
  print(i)
  ## BO settings
  fun <- funs_to_test[[i]]
  parameter_set <- getParamSet(fun)
  
  # set Control Argument of BO 
  # iters = Budget for classic bo
  ctrl <- makeMBOControl()
  ctrl <- setMBOControlTermination(ctrl, iters = iters)
  ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritCB(cb.lambda = 10), opt = "focussearch",
                              opt.focussearch.maxit = 1, opt.focussearch.points = 5000) # de facto random search
  ## end BO settings
  
  ## n BO runs
  for (test_run in 1:n) {
    design <- generateDesign(n = init_design, par.set = parameter_set, fun = lhs::randomLHS)
    
    tryCatch({
      res_weightedML = BO_weighted_ML(fun, design, control = ctrl, show.info = FALSE)
    }, error=function(e){cat("ERROR in weighted bo", test_run, conditionMessage(e), "\n")
      res_weightedML = NULL})
    
    tryCatch({
      res_unweightedML = mbo(fun, design, control = ctrl, show.info = FALSE)
    }, error=function(e){
      cat("ERROR in unweighted bo", test_run, conditionMessage(e), "\n")
      res_unweightedML = NULL})
    
    res_opt_paths_list_w[[test_run]] <- res_weightedML$res_mbo[["opt.path"]]
    res_opt_paths_list_uw[[test_run]] <- res_unweightedML[["opt.path"]]
    
    #weighted_ML_prior_res <- append(weighted_ML_prior_res,res_weightedML$res_mbo$y )
    #unweighted_ML_prior_res <- append(unweighted_ML_prior_res, res_unweightedML$y)
  }
  
  #remove NULLs from errors
  res_opt_paths_list_w <- res_opt_paths_list_w[!sapply(res_opt_paths_list_w, is.null)]
  res_opt_paths_list_uw <- res_opt_paths_list_uw[!sapply(res_opt_paths_list_uw, is.null)]
 
  opt_path_y_w = lapply(res_opt_paths_list_w, getOptPathY)
  opt_path_y_uw = lapply(res_opt_paths_list_uw , getOptPathY)
  
  # retrieve incumbent best target value (MINIMIZATION!)
  get_current_best <- function(proposals){
    for (o in 2:length(proposals)) {
      if(proposals[o] > proposals[o - 1]) 
        proposals[o] = proposals[o - 1] 
    }
    proposals  
  }
  
  opt_path_best_y_w = lapply(opt_path_y_w, get_current_best)
  opt_path_best_y_uw = lapply(opt_path_y_uw, get_current_best)
  
  global_res_weightedML[[i]] = opt_path_best_y_w
  global_res_unweightedML[[i]] = opt_path_best_y_uw
  
}




####
# Analysis of results
####

# Mean Optimization Paths (mop)
w_paths_mop <- list()
uw_paths_mop <- list()
w_paths_mop_sd <- list()
uw_paths_mop_sd <- list()


# weighted 
for (i in 1:length(funs_to_test)) {
  BO_paths_all_means <- matrix(nrow = iters, ncol = length(funs_to_test))
  
  BO_paths <- global_res_weightedML[[i]]
  
  BO_paths <- matrix(unlist(BO_paths), ncol = length(BO_paths))
  BO_paths <- as.data.frame(BO_paths)
  
  # remove initial design
  BO_paths_initial <- BO_paths %>% dplyr::slice_head(n = init_design)
  BO_paths_optim <- BO_paths %>% dplyr::slice_tail(n = nrow(BO_paths) - init_design)
  
  # log scale
  #BO_paths_optim <- log(BO_paths_optim) 
  
  # get lower bound/upper bound/mean per iteration:
  BO_paths_sd <- apply(BO_paths_optim, 1, sd)
  BO_paths_mean_per_iter <- apply(BO_paths_optim, 1, mean)
  
  w_paths_mop[[i]] <- BO_paths_mean_per_iter
  w_paths_mop_sd[[i]] <- BO_paths_sd
}


#unweighted 
for (i in 1:length(funs_to_test)) {
  BO_paths_all_means <- matrix(nrow = iters, ncol = length(funs_to_test))
  
  BO_paths <- global_res_unweightedML[[i]]
  
  BO_paths <- matrix(unlist(BO_paths), ncol = length(BO_paths))
  BO_paths <- as.data.frame(BO_paths)
  
  # remove initial design
  BO_paths_initial <- BO_paths %>% dplyr::slice_head(n = init_design)
  BO_paths_optim <- BO_paths %>% dplyr::slice_tail(n = nrow(BO_paths) - init_design)
  
  # log scale
  #BO_paths_optim <- log(BO_paths_optim) 
  
  # get mean per iteration:
  BO_paths_sd <- apply(BO_paths_optim, 1, sd)
  BO_paths_mean_per_iter <- apply(BO_paths_optim, 1, mean)
  
  uw_paths_mop[[i]] <- BO_paths_mean_per_iter
  uw_paths_mop_sd[[i]] <- BO_paths_sd
}

for (i in 1:length(funs_to_test)) {
  plot(w_paths_mop[[i]], type = "l")
  lines(uw_paths_mop[[i]], col = "red")
}

save(w_paths_mop, file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-10")
save(w_paths_mop_sd, file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-10")
save(uw_paths_mop, file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-10")
save(uw_paths_mop_sd, file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-10")

#save(global_res_weightedML, file = "weighted-ML-GP/results/weighted-ML-results-global-n-100-iters-20-ninit-4-funs26-27-29-34")
#save(global_res_unweightedML, file = "weighted-ML-GP/results/weighted-ML-results-global-n-100-iters-20-ninit-4-funs26-27-29-34")

#     
# compare results

#BO_paths_all_means[,k] <- BO_paths_mean_per_iter
#   
# BO_paths_min <- apply(BO_paths_all_means, 1, min)
# BO_paths_max <- apply(BO_paths_all_means, 1, max)
# BO_paths_acc_difference[fun] <- sum(BO_paths_max - BO_paths_min)
# 
# 
# BO_paths_acc_difference_mean_par <- BO_paths_acc_difference
# BO_paths_acc_difference_mean_par
# 


