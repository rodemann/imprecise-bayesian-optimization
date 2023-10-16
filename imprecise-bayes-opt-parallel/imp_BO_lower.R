
## TODO: adapt to S3 type changes in imp_BO_upper!!! 

makeRLearner.regr.igp.upper = function() {
  makeRLearnerRegr(
    cl = "regr.igp",
    package = "DiceKriging",
    par.set = makeParamSet(
      makeIntegerLearnerParam(id = "base.kernel", default = "gaussian")
    ),
    properties = c("numerics", "factors"),
    name = "Imprecise Gaussian Process",
    short.name = "igp",
    note = ""
  )
}

trainLearner.regr.igp.upper = function(.learner, .task, .subset,
                                       .weights = NULL, .imp_degree, ...) {
  data = getTaskData(.task)
  target = getTaskTargets(.task)
  base_kernel_fun = get_base_kernel(base_kernel)
  
  # estimate parameters of base kernel (via Maximum Likelihood as integrated in DiceKriging::km())
  # note that a small nugget is added in order to ensure numerical stability
  precise_GP_fit <- km(response = target, design = data, covtype = base_kernel,
                       estim.method = "MLE", optim.method = "gen", 
                       nugget = 1e-8*var(target), control = list(trace = FALSE))
  # get estimated theta (range/length-scale parameter) from S4 km class
  theta <- precise_GP_fit@covariance@range.val
  nugget <- precise_GP_fit@covariance@nugget
  
  # check if data is univariate 
  assert_vector(data)
  # compute kernel vectors and matrices required for formulas in Theorem 3 
  # in Mangili (2015). Notation follows the theorem
  distance <- abs(data - x)
  kernel_vector <- base_kernel_fun(distance, theta)
  data_dist_matrix <- as.matrix(dist(data, upper = TRUE))
  Kernel_matrix <- apply(data_dist_matrix, c(1,2), base_kernel_fun, theta = theta)
  Kernel_matrix_n <- Kernel_matrix + diag(nugget, nrow(Kernel_matrix))
  sum_rows_k <- solve(Kernel_matrix_n) %*% rep(1, nrow(Kernel_matrix_n))
  sum_all_k <- rep(1, nrow(Kernel_matrix_n)) %*% solve(Kernel_matrix_n) %*% rep(1, nrow(Kernel_matrix_n))
  
  
  #train_igp <- function(x, data, target) UseMethod("train_igp")
  #train_igp.default <- function(x, data, target){
  igp_lower_mean <- function(x){  
    
    if (abs(sum_rows_k %*% target) <= 1 + .imp_degree/sum_all_k){
      
      lower_mean <- function(x) { t(kernel_vector) %*% solve(Kernel_matrix_n) %*% target 
        + (1 - t(kernel_vector) %*% sum_rows_k) %*% t(sum_rows_k)/sum_all_k %*% target
        - .imp_degree * abs(1 - t(kernel_vector) %*% sum_rows_k)/sum_all_k
      }    
    }
    
    if (abs(sum_rows_k %*% target) > 1 + .imp_degree/sum_all_k){
      
      lower_mean <- function(x) { t(kernel_vector) %*% solve(Kernel_matrix_n) %*% target 
        + (1 - t(kernel_vector) %*% sum_rows_k) %*% t(sum_rows_k) %*% target/(.imp_degree + sum_all_k)
      }
    } else{stop("Condition for estimating lower mean posterior not fulfilled.")}
    lower_mean(x)
  }
  
}    






  
  get_base_kernel <- function(base_kernel){
    switch(base_kernel,
           "gauss" = {base_kernel <- function(dist, theta){exp(-1/2*(dist/theta)^2)}},
           "exp" = {base_kernel <- function(dist, theta){exp(-dist/theta)}},
           "matern3_2" = {base_kernel <- function(dist, theta){(1+sqrt(3)*dist/theta)*exp(-sqrt(3)*dist/theta)}},
           "matern5_2" = {base_kernel <- function(dist, theta){(1+sqrt(5)*dist/theta+(1/3)*5*(dist/theta)^2)*exp(-sqrt(5)*dist/theta)}},
           "powexp" = {base_kernel <- function(dist, theta, p){exp(-(dist/theta)^p)}}
    )
    base_kernel
  }
  
  
  
  
  
  