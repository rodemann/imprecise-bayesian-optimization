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
imp_BO_univariate_batch <- function(fun,
                              design = NULL,
                              control = NULL,
                              base_kernel = "powexp",
                              imprecision_degree,
                              #iters_per_batch,
                              number_of_batches = NULL,
                              iters_per_batch_per_model) {
  #input checking: 
  assert_numeric(imprecision_degree)
  assert_numeric(iters_per_batch_per_model)
  #assert_integer(iters_per_batch)
  assert_class(control, "MBOControl")
  # if(is.null(number_of_batches))
  #   number_of_batches <- control$iters  
  assert_number(number_of_batches)

  

  initial_design_size <- nrow(design)  
  initial_targets <- apply(design, 1, fun) # save fun values of initial design for opt.path
  # create 2S +1  surrogate models (upper/lower Imprecise Gaussian Process (igp) and classic
  # precise Gaussian Process)
  lrns_upper <- list()
  lrns_upper <- foreach(i = seq_along(imprecision_degree)) %do% {
    makeLearner("igp_upper_mean", predict.type = "se", base_kernel = base_kernel, imprecision = imprecision_degree[i])
  }
  
  lrns_lower <- list()
  lrns_lower <- foreach(i = seq_along(imprecision_degree)) %do% {
    makeLearner("igp_lower_mean", predict.type = "se", base_kernel = base_kernel, imprecision = imprecision_degree[i])
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
  
  BO_all_results <- list()
  
  evaluated_points_y <- list()  
  # TODO parallelize outer loop?
  for (b in 1:number_of_batches) {
  
  # set iterations
  control <- setMBOControlTermination(control, iters = iters_per_batch_per_model[b])  
  # make sure all surrogate models are stored for computation of surprise function
  control$store.model.at <- 1:control$iters
  # run BOs in parallel
  cl <- parallel::makeForkCluster(length(surrogate_models))
  doParallel::registerDoParallel(cl)
  # run BO
  BO_results <- foreach(i = seq_along(surrogate_models), .packages = "mlrMBO") %dopar% {
    #source("imprecise-bayes-opt-parallel/igp_upper_mean.R", local = TRUE)
    #source("imprecise-bayes-opt-parallel/igp_lower_mean.R", local = TRUE)
    mbo(fun = fun, design = design, control = control, learner = surrogate_models[[i]])
  }
  stopCluster(cl)
  # save results of this batch
  BO_all_results[[b]] <- BO_results
  
  # get opt path and proposed points as data frame
  opt.paths.raw <- lapply(BO_results, "[[", "opt.path")
  opt.paths <- lapply(opt.paths.raw, as.data.frame)
  n.initial <- nrow(design) # need to know
  opt.paths <- lapply(opt.paths, slice, n.initial + 1:control$iters)
  true_targets <- lapply(opt.paths, "[[", "y")
  prop_points <- lapply(opt.paths, "[[", "x")
  prop_points <- lapply(prop_points, function(x) data.frame("x" = x))
  
  models_per_iter <- lapply(BO_results, "[[", "models" )
  
  sm.at.prop.points <- matrix(nrow = length(surrogate_models), ncol = control$iter)
  for(i in seq_along(surrogate_models)) {
    for (j in 1:control$iter) {
      pred <- predict(models_per_iter[[i]][[j]], newdata = data.frame("x" = prop_points[[i]][j,]))
      sm.at.prop.points[i,j] <- pred$data$response
      }
  }
  # proposed points as matrix 
  true_targets_matrix <- do.call(rbind, true_targets)
  
  # surprise matrix (see definition of suprise function in paper)
  surprise <- true_targets_matrix - sm.at.prop.points
  # accumulated surprise per surrogate model
  acc_surprise <- apply(surprise, 1, sum)
  ranked_acc_suprise <- rank(-abs(acc_surprise)) # absolute deviation from true target, then rank so that the 
  # lower the value, the better -> -abs()
  
  #throw away worst (length(surrogate_models)-1)/2 models
  # how many loosers?
  no_kills <- (length(surrogate_models)+1)/2
  # where to find them?
  looser_indices <- which(ranked_acc_suprise %in% 1:no_kills)
  
  surrogate_models <- surrogate_models[-c(looser_indices)]
  ## update design for next batch
  # retrieve all evaluated points (proposed points that have been evaluated) across all batches
  evaluated_points <- do.call(rbind, opt.paths) %>% select(x) 

    # save target values for opt path 
  evaluated_points_y[[b]] <- true_targets
  # update design 
  design <- rbind(design, evaluated_points)

  
  }  
  
  # return results in a handy format
  BO_optimal_targets <- lapply(BO_all_results, function(x) lapply(x, '[[', 'y'))
  BO_optimal_x <- lapply(BO_all_results, function(x) lapply(x, '[[', 'x'))

  BO_opt_paths <- lapply(BO_all_results, function(x) lapply(x, '[[', 'opt.path'))
  BO_opt_paths_y <- lapply(BO_opt_paths, function(x) lapply(x, getOptPathY))
  # get initial design
  BO_opt_path_initial <- BO_opt_paths_y[[1]][[1]][1:initial_design_size]

    # get all proposed points (opt path)
  global_opt_path <- list()
  for (l in 1:number_of_batches) {
    global_opt_path[[l]] <- do.call(rbind, evaluated_points_y[[l]])
    if (BO_all_results[[1]][[1]]$control$minimize == TRUE){
      global_opt_path[[l]] <- apply(global_opt_path[[l]], 2, min)
      no_models <- length(evaluated_points_y[[l]])
      global_opt_path[[l]] <- rep(global_opt_path[[l]], each = no_models)
    }else{
      global_opt_path[[l]] <- apply(global_opt_path[[l]], 2, max)
      no_models <- length(evaluated_points_y[[l]])
      global_opt_path[[l]] <- rep(global_opt_path[[l]], each = no_models)
      }
  }
  # append initial design
  global_opt_path <- c(initial_targets, global_opt_path) %>% unlist()
  
  if(BO_all_results[[1]][[1]]$control$minimize == TRUE)
    opt_path_y <- global_opt_path %>% overwrite_decr()
  else
    opt_path_y <- global_opt_path %>% overwrite_incr()
  
  # remove names
  opt_path_y <- unname(opt_path_y)
  
  ## get global optimum
  # check direction
  if (BO_all_results[[1]][[1]]$control$minimize == TRUE){
    global_opt_target <- min(unlist(BO_optimal_targets))
    global_opt_index <- which.min(unlist(BO_optimal_targets))

  }else{
    global_opt_target <- max(unlist(BO_optimal_targets))
    global_opt_index <- which.max(unlist(BO_optimal_targets))
  }
  global_opt_x <- unlist(BO_optimal_x)[global_opt_index]
  BO_returned_results <- list()
  BO_returned_results$opt_path_y_global <- opt_path_y
  BO_returned_results$all.bo.results <- BO_all_results
  BO_returned_results$global.optimum <- data.frame("Optimal Target Value" = global_opt_target, 
                                                   "Optimal Covariate Value" = global_opt_x)
  BO_returned_results

}


# overwrite function (make opt-path monotonically decreasing/increasing)

overwrite_decr <- function(proposals){
  for (o in 2:length(proposals)) {
    if(proposals[o] > proposals[o - 1]) 
      proposals[o] = proposals[o - 1] 
  }
  proposals
}



overwrite_incr <- function(proposals){
  for (o in 2:length(proposals)) {
    if(proposals[o] < proposals[o - 1]) 
      proposals[o] = proposals[o - 1] 
  }
  proposals
}


  
# optimal_params <- sapply(BO_results, "[[", "x") %>% unlist()  
# optimal_targets <- sapply(BO_results, "[[", "y")    
# 
# incumbent_results <- data.frame("Optimal Parameters " = optimal_params, "Optimal Targets" = optimal_targets)
# 


  # # create surrogate models (upper/lower Imprecise Gaussian Process (igp) and classic
  # # precise Gaussian Process)
  # lrn_upper <- makeLearner("igp_upper_mean", base_kernel = base_kernel, imprecision = imprecision_degree)
  # lrn_lower <- makeLearner("igp_lower_mean", base_kernel = base_kernel, imprecision = imprecision_degree)
  # lrn_classic <- makeLearner("regr.km", covtype = base_kernel, predict.type = "se", optim.method = "gen", 
  #                            control = list(trace = FALSE), config = list(on.par.without.desc = "warn"))
  # # ensure numerical stability in km {DiceKriging} cf. github issue and recommendation by Bernd Bischl 
  # y = apply(design, 1, obj.fun)
  # Nuggets = 1e-8*var(y)
  # lrn_classic = setHyperPars(learner = lrn_classic, nugget=Nuggets)
  # 
  # # run BOs in parallel
  # surrogate_models <- list(lrn_upper, lrn_lower, lrn_classic)
  # registerDoParallel(3)  # use multicore structure 
  # BO_results <- foreach(i = seq_along(surrogate_models)) %dopar% {
  #   mbo(fun = fun, design = design, control = control, learner = surrogate_models[[i]])
  # }
  # 
  # optimal_params <- sapply(BO_results, "[[", "x") %>% unlist()  
  # optimal_targets <- sapply(BO_results, "[[", "y")    
  # 
  # data.frame("Optimal Parameters " = optimal_params, "Optimal Targets" = optimal_targets, 
  #            row.names = c("Upper Mean", "Lower Mean", "Precise GP"))
  # 
  
