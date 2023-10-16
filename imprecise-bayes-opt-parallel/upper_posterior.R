
# function to compute upper posterior of IGP

# idea: create own deterministic learner and pass it to mbo() as surrogate model  https://mlr.mlr-org.com/articles/tutorial/create_learner.html


#' This function will compute the upper bound of the posterior Gaussian Process as in \cite[Page 190]{mangili-15-sipta}
#' 
#' 
#'
#' @param k kernel vector
#' @param y list of observations
#' @return The sum of \code{x} and \code{y}
#' @examples
#' add(1, 1)
#' add(10, 1)

upper_posterior <- function(){
  
}