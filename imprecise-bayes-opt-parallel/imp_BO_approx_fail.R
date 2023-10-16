library(kergp) # use predefined kernel functions as base kernel in Imprecise Gaussian Process, cf. Mangili (2015)
library(DiceKriging)
library(kergp)
library(mlrMBO)

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
#' 
#' 
imp_BO_univariate <- function(fun,
                   design = NULL,
                   control = NULL,
                   base_kernel = "gauss",
                   imp_degree,  
                   approximate = TRUE,
                   M, 
                   h,
                   viz = FALSE,
                   iters_to_show = c(1)) {
  if(approximate == TRUE){
    # imprecision degree (notation follows Mangili (2015))
    c <- imp_degree
    
    # sample parameters that determine the GP, notation as in Mangili (2015)
    #M <- runif(1,0,100)
    #h <- sample(c(-1,1), 1)
    # define base kernel function
    base_kernel_fun <- get_base_kernel(base_kernel)
    # estimate parameters of base kernel (via Maximum Likelihood as integrated in DiceKriging::km())
    target <- apply(design, 1, fun) 
    initial_fit <- km(response = target, design = design, covtype = base_kernel,
                      estim.method = "MLE", optim.method = "gen", 
                      nugget = 1e-8*var(target), control = list(trace = FALSE))
    # get estimated theta (range/length-scale parameter)
    theta <- initial_fit@covariance@range.val
    # in case of power-exponential: get estimated power p
    if (base_kernel == "powexp")
      p <- initial_fit@covariance@shape.val
    M <- M
    h <- h
    # define kernel 
    kernel <- function(x, y) {
      dist <- abs(x-y)
      if(base_kernel != "powexp"){
        base_kernel_fun(dist, theta) + (1 + M)/c
      }else{
      base_kernel_fun(dist, theta, p) + (1 + M)/c
      }    
    }
    # define mean function
    mean <- M*h
    
    # define learner
    lrn <- makeLearner("regr.km", coef.trend = mean, 
    kernel = kernel, predict.type = "se", optim.method = "gen", 
    control = list(trace = FALSE), config = list(on.par.without.desc = "warn"))
    
     
    # run BO
    if (viz == FALSE)
      mbo(fun = fun, design = design, control = control, learner = lrn)
    else{
      mbo_run = exampleRun(fun = fun, design = design, control = control, learner = lrn)
      plotExampleRun(mbo_run, iters = iters_to_show)
      }
      
  }
  else{
  print("to be done")  
  }
  
}


get_base_kernel <- function(base_kernel){
  switch (base_kernel,
    "gauss" = {base_kernel <- function(dist, theta){exp(-1/2*(dist/theta)^2)}},
    "exp" = {base_kernel <- function(dist, theta){exp(-dist/theta)}},
    "matern3_2" = {base_kernel <- function(dist, theta){(1+sqrt(3)*dist/theta)*exp(-sqrt(3)*dist/theta)}},
    "matern5_2" = {base_kernel <- function(dist, theta){(1+sqrt(5)*dist/theta+(1/3)*5*(dist/theta)^2)*exp(-sqrt(5)*dist/theta)}},
    "powexp" = {base_kernel <- function(dist, theta, p){exp(-(dist/theta)^p)}}
  )
  base_kernel
}






