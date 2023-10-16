library(mlrMBO)
library(doParallel)
library(foreach)
source("imprecise-bayes-opt-parallel/igp_upper_mean.R")
source("imprecise-bayes-opt-parallel/igp_lower_mean.R")


#' Run Imprecise BO 
#' 
#' This functions starts an imprecise version of Bayesian Optimization mlrMBO::mbo() 
#'
#' @param 
#' @return 
#' @examples
#'  
#'   
#' 
imp_BO_univariate <- function(fun,
                              design = NULL,
                              control = NULL,
                              base_kernel = "powexp",
                              imprecision_degree) {
  # input checking of imp_BO specific arguments (rest is checked in mbo() below)
  assert_character(base_kernel)
  assert_numeric(imprecision_degree)

  if(length(imprecision_degree) == 1){
  # create 3 surrogate models (upper/lower Imprecise Gaussian Process (igp) and classic
  # precise Gaussian Process)
  lrn_upper <- makeLearner("igp_upper_mean", predict.type = "se", base_kernel = base_kernel, imprecision = imprecision_degree)
  lrn_lower <- makeLearner("igp_lower_mean", predict.type = "se", base_kernel = base_kernel, imprecision = imprecision_degree)
  lrn_classic <- makeLearner("regr.km", covtype = base_kernel, predict.type = "se", optim.method = "gen", 
                            control = list(trace = FALSE), config = list(on.par.without.desc = "warn"))
  # ensure numerical stability in km {DiceKriging} cf. github issue and recommendation by Bernd Bischl 
  y = apply(design, 1, fun)
  Nuggets = 1e-8*var(y)
  lrn_classic = setHyperPars(learner = lrn_classic, nugget=Nuggets)
  
  # run BOs in parallel
  surrogate_models <- list(lrn_upper, lrn_lower, lrn_classic)
  cl <- parallel::makeForkCluster(length(surrogate_models))
  doParallel::registerDoParallel(cl)
  BO_results <- foreach(i = seq_along(surrogate_models), .packages = c("mlrMBO")) %dopar% {
    mbo(fun = fun, design = design, control = control, learner = surrogate_models[[i]])
  }
  stopCluster(cl)
  
  returned_list <- return_stuff(BO_results)
  returned_list %>% return()
  }else{
  ###
  # multiple imprecision degrees:
  
  # create 2S +1  surrogate models (upper/lower Imprecise Gaussian Process (igp) and classic
  # precise Gaussian Process)
  lrns_upper <- list()
  lrns_upper <- foreach(i = seq_along(imprecision_degree)) %do% {
    makeLearner("igp_upper_mean", base_kernel = base_kernel, imprecision = imprecision_degree[i])
  }
  
  lrns_lower <- list()
  lrns_lower <- foreach(i = seq_along(imprecision_degree)) %do% {
    makeLearner("igp_lower_mean", base_kernel = base_kernel, imprecision = imprecision_degree[i])
  }
  
  lrn_classic <- makeLearner("regr.km", covtype = base_kernel, predict.type = "se", optim.method = "gen", 
                             control = list(trace = FALSE), config = list(on.par.without.desc = "warn"))
  # ensure numerical stability in km {DiceKriging} cf. github issue and recommendation by Bernd Bischl 
  y = apply(design, 1, fun)
  Nuggets = 1e-8*var(y)
  lrn_classic = setHyperPars(learner = lrn_classic, nugget = Nuggets)
  
  #combine models
  lrns_upper_lower <- c(lrns_upper, lrns_lower)
  surrogate_models <- append(lrns_upper_lower, list(lrn_classic))
  # run BOs in parallel
  registerDoParallel(length(surrogate_models))  # use multicore structure 
  BO_results <- foreach(i = seq_along(surrogate_models)) %dopar% {
    mbo(fun = fun, design = design, control = control, learner = surrogate_models[[i]])
  }
  
  returned_list <- return_stuff(BO_results)
  returned_list %>% return()  
  }
}
  
  
return_stuff <- function(BO_results){
  # list that will be returned:  
  returned_list <- list()
  
  optimal_params <- sapply(BO_results, "[[", "x") %>% unlist()  
  optimal_targets <- sapply(BO_results, "[[", "y")    
    all_optima <- data.frame("Optimal Parameters " = optimal_params, "Optimal Targets" = optimal_targets,
                           row.names = c("Upper Mean", "Lower Mean", "Precise GP"))
    if (BO_results[[1]]$control$minimize == TRUE) {
      global_optimum <- min(optimal_targets)
    }else{
      global_optimum <- max(optimal_targets)
    }
  # attach results  
  returned_list$all_optima <- all_optima
  returned_list$global_optimum <- global_optimum
  
  # global opt path
  opt_paths <- lapply(BO_results, "[[", "opt.path")
  opt_paths_y <- lapply(opt_paths, getOptPathY) 
  opt_paths_init <- opt_paths_y[[1]][1:nrow(design)]
  opt_paths_prop <- lapply(opt_paths_y, "[", -(1:nrow(design)))
  opt_paths_prop_per_iter <- do.call(rbind, opt_paths_prop) 
  if (BO_results[[1]]$control$minimize == TRUE) {
    opt_paths_prop_per_iter <- apply(opt_paths_prop_per_iter, 2, max)
    opt_path_global_prop <- lapply(opt_paths_prop_per_iter, rep, times = length(opt_paths_prop)) %>% unlist()
    opt_path_global_prop <- sort(opt_path_global_prop, decreasing = TRUE)
    opt_path_y_global <- c(opt_paths_init, opt_path_global_prop)
    for (o in 2:length(opt_path_y_global)) {
      if (opt_path_y_global[o] > opt_path_y_global[o - 1]) 
        opt_path_y_global[o] = opt_path_y_global[o - 1] 
    }
  }else{
    opt_paths_prop_per_iter <- apply(opt_paths_prop_per_iter, 2, max)
    opt_path_global_prop <- lapply(opt_paths_prop_per_iter, rep, times = length(opt_paths_prop)) %>% unlist()
    opt_path_global_prop <- sort(opt_path_global_prop, decreasing = FALSE)
    opt_path_y_global <- c(opt_paths_init, opt_path_global_prop)
    for (o in 2:length(opt_path_y_global)) {
      if (opt_path_y_global[o] < opt_path_y_global[o - 1]) 
        opt_path_y_global[o] = opt_path_y_global[o - 1] 
    }
  }


  returned_list$opt_path_y <- opt_path_y_global 
  
  #return results
  returned_list
}


  
