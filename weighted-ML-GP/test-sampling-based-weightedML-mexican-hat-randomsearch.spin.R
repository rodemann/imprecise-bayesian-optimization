library(tidyr)
library(mlrMBO)
library(DiceKriging)

# Instructions: 
trace(logLikFun, edit=TRUE)
#and edit as follows from line 8:

# weights <- eval(weights, envir = globalenv())
# i <- eval(i, envir = globalenv())
# BO_type <- eval(BO_type, envir = globalenv())
# if (i != 1L & BO_type == "weighted")
#   if (nrow(model@X) == length(weights)) {
#     indices <- sample(1:length(weights), prob = weights,
#                       replace = TRUE)
#     model@X <- as.matrix(model@X[indices] + rnorm(indices,
#                                                   0, 0.01))
#     model@y <- as.matrix(model@y[indices] + rnorm(indices,
#                                                   0, 0.01))
#     warning("weights used")
#   }



source("weighted-ML-GP/BO_weighted_ML_sampling-based.R")
# unable noisy target functions and nugget stability (because they affect the usage of loglikfun heavily)
source("weighted-ML-GP/makeMBOLearner-no-noise.R")
source("_Explore_Exploit_Measures/xplxplMBO-jr.R")
#source("univariate-benchmark-functions/get-synthetic-benchmark-funs-for-minimization.R")

funs_to_test <- all_benchmark_funs[c(31)]
res_opt_paths_list_w <- list()
res_opt_paths_list_uw <- list()
global_res_weightedML <- list()
global_res_unweightedML <- list()

# hyperparams
n <- 40
iters <- 25
init_design <- 5

for (fun_i in 1:length(funs_to_test)) {
  print(fun_i)
  ## BO settings
  fun <- funs_to_test[[fun_i]]
  parameter_set <- getParamSet(fun)
  #parameter_set$pars$x$lower = -0.2
  #parameter_set$pars$x$upper = 0.2
  
  # set Control Argument of BO 
  # iters = Budget for classic bo
  ctrl <- makeMBOControl()
  ctrl <- setMBOControlTermination(ctrl, iters = iters)
  ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI(), opt = "focussearch",
                              opt.focussearch.maxit = 1, opt.focussearch.points = 1000) # de facto random search
  ## end BO settings
  
  ## n BO runs
  for (test_run in 1:n) {
    print(test_run)
    design <- generateDesign(n = init_design, par.set = parameter_set, fun = lhs::randomLHS)
    
    BO_type = "weighted"
    assign("BO_type", BO_type, envir = globalenv()) # needed in edited version of loglikfun
    tryCatch({
      res_weightedML = BO_weighted_ML(fun, design, control = ctrl, show.info = FALSE)
    }, error=function(e){
      cat("ERROR in weighted bo", test_run, conditionMessage(e), "\n")
      res_weightedML = NULL
      assign("res_weightedML", res_weightedML, envir = globalenv()) 
      })
    
    BO_type = "unweighted"
    assign("BO_type", BO_type, envir = globalenv()) # needed in edited version of loglikfun
    tryCatch({
      res_unweightedML = mbo(fun, design, control = ctrl, show.info = FALSE)
    }, error=function(e){
      cat("ERROR in unweighted bo", test_run, conditionMessage(e), "\n")
      res_unweightedML = NULL
      assign("res_unweightedML", res_unweightedML, envir = globalenv()) 
      })
    
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
  
  global_res_weightedML[[fun_i]] = opt_path_best_y_w
  global_res_unweightedML[[fun_i]] = opt_path_best_y_uw
  
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
for (fun_i in 1:length(funs_to_test)) {
  BO_paths_all_means <- matrix(nrow = iters, ncol = length(funs_to_test))
  
  BO_paths <- global_res_weightedML[[fun_i]]
  
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
  
  w_paths_mop[[fun_i]] <- BO_paths_mean_per_iter
  w_paths_mop_sd[[fun_i]] <- BO_paths_sd
}


#unweighted 
for (fun_i in 1:length(funs_to_test)) {
  BO_paths_all_means <- matrix(nrow = iters, ncol = length(funs_to_test))
  
  BO_paths <- global_res_unweightedML[[fun_i]]
  
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
  
  uw_paths_mop[[fun_i]] <- BO_paths_mean_per_iter
  uw_paths_mop_sd[[fun_i]] <- BO_paths_sd
}

for (fun_i in 1:length(funs_to_test)) {
  plot(w_paths_mop[[fun_i]], type = "l")
  lines(uw_paths_mop[[fun_i]], col = "red")
}

save(w_paths_mop, file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs")
save(w_paths_mop_sd, file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs")
save(uw_paths_mop, file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs")
save(uw_paths_mop_sd, file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs")

save(global_res_weightedML, file = "weighted-ML-GP/results/weighted-ML-results-global-n-100-iters-20-ninit-4-mexican-hat-rs")
save(global_res_unweightedML, file = "weighted-ML-GP/results/weighted-ML-results-global-n-100-iters-20-ninit-4-mexican-hat-rs")

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


