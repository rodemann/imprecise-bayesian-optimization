library(tidyverse)
library(mlr)
library(BBmisc)
library(parallel)
library(foreach)
library(doParallel)


# this script creates an S3 object and integrates it in the mlr framework
# It first defines the 0bject (1), then the training function (2) and the predicition function (3)

# (1) create S3 learner object
makeRLearner.igp_upper_mean = function() {
  makeRLearnerRegr(
    cl = "igp_upper_mean",
    package = "DiceKriging",
    par.set = makeParamSet(
      makeDiscreteLearnerParam(id = "covtype", default = "matern5_2",
                               values = list("gauss", "matern5_2", "matern3_2", "exp", "powexp")),
      makeDiscreteLearnerParam(id = "base_kernel", default = "matern5_2",
                               values = list("gauss", "matern5_2", "matern3_2", "exp", "powexp")),
      makeNumericVectorLearnerParam(id = "coef.trend"),
      makeNumericVectorLearnerParam(id = "coef.cov"),
      makeNumericVectorLearnerParam(id = "coef.var"),
      makeNumericLearnerParam(id = "nugget"),
      makeLogicalLearnerParam(id = "nugget.estim", default = FALSE),
      makeNumericVectorLearnerParam(id = "noise.var"),
      makeDiscreteLearnerParam(id = "estim.method", default = "MLE",
                               values = c("MLE", "LOO")),
      makeDiscreteLearnerParam(id = "optim.method", default = "BFGS",
                               values = c("BFGS", "gen")),
      makeNumericVectorLearnerParam(id = "lower"),
      makeNumericVectorLearnerParam(id = "upper"),
      makeNumericVectorLearnerParam(id = "parinit"),
      makeIntegerLearnerParam(id = "multistart", default = 1L, lower = 1L),
      makeUntypedLearnerParam(id = "control"),
      makeLogicalLearnerParam(id = "gr", default = TRUE),
      makeLogicalLearnerParam(id = "iso", default = FALSE),
      makeLogicalLearnerParam(id = "scaling", default = FALSE),
      makeUntypedLearnerParam(id = "knots"),
      makeLogicalLearnerParam(id = "jitter", default = FALSE, when = "predict"),
      makeNumericLearnerParam(id = "nugget.stability", requires = quote(!nugget.estim && is.null(nugget))),
      makeNumericLearnerParam(id = "imprecision", default = 1L)
    ),
    par.vals = list(jitter = FALSE),
    properties = c("numerics", "se"),
    name = "Upper Posterior Mean of Imprecise Gaussian Process",
    short.name = "igp",
    note = ""
  )
}


# (2) training function of learner. this is where all the magic happens. Note that the predict functions is
# already designed here and will just be called from the "official" mlr S3 predict function's environment
trainLearner.igp_upper_mean = function(.learner, .task, .subset,
                                 .weights = NULL, ...) {
  args = list(...)
  data = getTaskData(.task, target.extra = TRUE)
  design = data$data
  # check if design is univariate 
  checkmate::assert_vector(design)
  response = data$target
  base_kernel_fun = get_base_kernel(args$base_kernel)
  
  ### estimate parameters of base kernel (via Maximum Likelihood as integrated in DiceKriging::km())
  # change design for this precise fit if data are too close to each other for numerical stability 
  # (otherwise error in kmNuggets.init is thrown)
  
  # get data points that are too close to each other and number of those
  too_close <- lapply(list(design, response), function(data){
  data <- as.matrix(dist(data, upper = FALSE))
  dist_vector <- apply(data, 2, function(x){
    n <- length(x)
    sort(x)[2]
  })
  too_close_data <- which(dist_vector < 0.0001) %>% unname()
  too_close_nr <- sum(dist_vector < 0.0001)
  list(too_close_data, too_close_nr)
  })
  too_close_indices <- c(too_close[[1]][[1]], too_close[[2]][[1]]) %>% unique()
  
  #remove close data and reponse points if there are more than 5 in design OR in data
  if(too_close[[1]][[2]] > 5 | too_close[[2]][[2]] > 5){  
     design_precise_fit <- data.frame(x = design[-c(too_close_indices),])
     response_precise_fit <- response[-c(too_close_indices)]  
  }else{
     design_precise_fit <- design
     response_precise_fit <- response
   }
  if(nrow(design_precise_fit) <= 2 | length(response_precise_fit) <= 2)
    stop("data contains (approx.) identical values")
  
  # note that a small nugget is added in order to ensure numerical stability
  precise_GP_fit <- DiceKriging::km(response = response_precise_fit, design = design_precise_fit, covtype = args$base_kernel,
                    estim.method = "MLE", optim.method = "gen",
                    nugget = 1e-4*var(response), nugget.estim = TRUE, control = list(trace = FALSE))
  # get estimated theta (range/length-scale parameter) from S4 km class
  theta <- precise_GP_fit@covariance@range.val
  nugget <- precise_GP_fit@covariance@nugget
  std_dev <- precise_GP_fit@covariance@sd2
  var <- (std_dev)^2
  if (args$base_kernel == "powexp")
    power <- precise_GP_fit@covariance@shape.val
  
  # compute Kernel Matrix and sum_rows and sum_all
  design_dist_matrix <- as.matrix(dist(design, upper = TRUE))
  if (args$base_kernel == "powexp")
    Kernel_matrix <- base_kernel_fun(design_dist_matrix, theta = theta, var = var, p = power) else
      Kernel_matrix <- base_kernel_fun(design_dist_matrix, theta = theta, var = var)
  Kernel_matrix_n <- Kernel_matrix + diag(nugget, nrow(Kernel_matrix))
  Kernel_matrix_n_inverse <- solve(Kernel_matrix_n)
  sum_rows_k <- solve(Kernel_matrix_n) %*% rep(1, nrow(Kernel_matrix_n))
  sum_all_k <- rep(1, nrow(Kernel_matrix_n)) %*% solve(Kernel_matrix_n) %*% rep(1, nrow(Kernel_matrix_n))
  sum_all_k <- as.double(sum_all_k) # seems necessary
  
  igp_upper_mean_scalar <- function(scalar) {
    # compute kernel vectors and matrices required for formulas in Theorem 3 
    # in Mangili (2015). Notation follows the theorem.
    distance <- abs(design - rep(scalar,nrow(design)))
    if(args$base_kernel == "powexp")
      kernel_vector <- base_kernel_fun(distance, theta, var = var, p = power) else
        kernel_vector <- base_kernel_fun(distance, theta, var = var)
      
      # compute upper bound of mean prediction
      if (abs(t(sum_rows_k) %*% response / sum_all_k) <= 1 + args$imprecision/sum_all_k){
        upper_mean <- (t(kernel_vector) %*% Kernel_matrix_n_inverse %*% response 
                       + (1 - t(kernel_vector) %*% sum_rows_k) %*% (t(sum_rows_k)/sum_all_k) %*% response
                       + args$imprecision * abs(1 - t(kernel_vector) %*% sum_rows_k)/sum_all_k)
      }
      
      if (abs(t(sum_rows_k) %*% response / sum_all_k) > 1 + args$imprecision/sum_all_k){
        upper_mean <- (t(kernel_vector) %*% Kernel_matrix_n_inverse %*% response 
                       + (1 - t(kernel_vector) %*% sum_rows_k) %*% (t(sum_rows_k)/sum_all_k) %*% response
                       + args$imprecision * (1 - t(kernel_vector) %*% sum_rows_k)/sum_all_k)
      } 
      upper_mean
  }
  
  igp_upper_mean_se_scalar <- function(scalar) {
    # compute kernel vectors and matrices required for definition of standard error (std. deviation) in Theorem 4 
    # in Mangili (2015). Notation follows the theorem.Distance in kernel is 0, since k(x,x) in theorem 4
    distance <- abs(design - rep(scalar,nrow(design)))
    if(args$base_kernel == "powexp")
      kernel_vector <- base_kernel_fun(distance, theta, var = var, p = power) else
        kernel_vector <- base_kernel_fun(distance, theta, var = var)
      
    if(args$base_kernel == "powexp")
      kernel_null <- base_kernel_fun(0, theta, p = power, var = var) else
        kernel_null <- base_kernel_fun(0, theta, var = var)
    #compute standard error  
      upper_mean_var <- (kernel_null - (t(kernel_vector) %*% Kernel_matrix_n_inverse %*% as.matrix(kernel_vector)) 
                       + ((1 - t(kernel_vector) %*% sum_rows_k)^2 / sum_all_k))
      # control for tiny negative variances due to numerical instability
      #if(-0.00001 <= upper_mean_var & upper_mean_var < 0) upper_mean_var <- 0
      if(upper_mean_var < 0) stop(paste("negative var estimates produced by igp upper mean of value",upper_mean_var))
      upper_mean_se <- sqrt(upper_mean_var)
      upper_mean_se
                       
  }
  
  # this is essentially the predict function (the above is the "training" with provided training data)
  # it is attached to the model object below and will then be called in the predictLearner function (also below)
  igp_upper_mean <- function(x){  
  # input checking for argument newdata (see below)
  checkmate::assert_data_frame(x, ncols = 1)
  #TODO parallelize apply fun
  pred <- list()
  if(.learner$predict.type == "response"){
    pred$mean <- apply(x, 1, igp_upper_mean_scalar)
  }
  if(.learner$predict.type == "se"){
    pred$mean <- apply(x, 1, igp_upper_mean_scalar)
    pred$se <- apply(x, 1, igp_upper_mean_se_scalar)
  }
  pred
  } 

  
  model = structure(list(x = design, y = response), class = "igp_upper_mean") 
  
  predict.igp_upper_mean = function(newdata) {
    igp_upper_mean(newdata)
  }
  model$learner.model = predict.igp_upper_mean
  return(model)
}    


# (3) mlr predict function
predictLearner.igp_upper_mean = function(.learner, .model, .newdata, ...) {
  se = (.learner$predict.type != "response")
  p =   .model$learner.model$learner.model(newdata = .newdata) 
  if (!se) {
    return(p$mean)
  } else {
    cbind(p$mean, p$se)
  }
}

# helper function to get base kernel
get_base_kernel <- function(base_kernel){
  switch(base_kernel,
          "gauss" = {base_kernel <- function(dist, theta, var){var * exp(-1/2*(dist/theta)^2)}},
          "exp" = {base_kernel <- function(dist, theta, var){var * exp(-dist/theta)}},
          "matern3_2" = {base_kernel <- function(dist, theta, var){var * (1+sqrt(3)*dist/theta)*exp(-sqrt(3)*dist/theta)}},
          "matern5_2" = {base_kernel <- function(dist, theta, var){var * (1+sqrt(5)*dist/theta+(1/3)*5*(dist/theta)^2)*exp(-sqrt(5)*dist/theta)}},
          "powexp" = {base_kernel <- function(dist, theta, p, var){var * exp(-(dist/theta)^p)}}
  )
  base_kernel
}





