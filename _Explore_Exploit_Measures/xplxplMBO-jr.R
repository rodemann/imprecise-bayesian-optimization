#library(mlrMBO)
library(checkmate)

#source modified functions from mlrMBO
source("_Explore_Exploit_Measures/proposePointsByInfillOptimization-jr.R")
source("_Explore_Exploit_Measures/makeMBOResult.OptState-jr.R")
source("_Explore_Exploit_Measures/getSupportedInfillOptFunctions-jr.R")
source("_Explore_Exploit_Measures/proposePointsByInfillOptimization-jr.R")
source("_Explore_Exploit_Measures/getInfillOptFunction-jr.R")
source("_Explore_Exploit_Measures/checkStuff-jr.R")

# multi point proposal stuff
source("_Explore_Exploit_Measures/proposePointsConstantLiarXplxpl-jr.R")
source("_Explore_Exploit_Measures/proposePoints-jr.R")
source("_Explore_Exploit_Measures/getSupportedMultipointInfillOptFunctions-jr.R")

# source new infill optimization functions "...Savepts"
source("_Explore_Exploit_Measures/infillOptFocusSavepts-jr.R")
source("_Explore_Exploit_Measures/infillOptEASavepts-jr.R")
source("_Explore_Exploit_Measures/infillOptCMAESSavepts-jr.R")

# source function that computes explore exploit measures and returns them as data frame
source("_Explore_Exploit_Measures/computeXplxpl-jr.R")

# source functions that calculate se based explore exploit measures
source("_Explore_Exploit_Measures/xplxplSESingleLocal-jr.R")
source("_Explore_Exploit_Measures/xplxplSESingleGlobal-jr.R")
source("_Explore_Exploit_Measures/xplxplSEMultiLocal-jr.R")
source("_Explore_Exploit_Measures/xplxplSEMultiGlobal-jr.R")

# source functions that calculate distance based explore exploit statistics
source("_Explore_Exploit_Measures/xplxplDistNum-jr.R")
source("_Explore_Exploit_Measures/xplxplDistCat-jr.R")


# ... or simply (final version): 
#files.sources = list.files("2_acquisition_function/2.3_Explore_Exploit_Ratio")
#sapply(paste("2_acquisition_function/2.3_Explore_Exploit_Ratio/", files.sources, sep = ""), source)



#' ExploreExploitmbo Function
#' Runs mbo (using mlrMBO::mbo())and computes a statistic for each iteration 
#' indicating whether the mbo exploits or explores. Inputs from mbo():
#' 
#' @param fun [\code{smoof_function}]\cr
#'   Fitness function to optimize.
#'   For one dimensional target functions you can obtain a \code{smoof_function} by using \code{\link[smoof]{makeSingleObjectiveFunction}}.
#'   For multi dimensional functions use \code{\link[smoof]{makeMultiObjectiveFunction}}.
#'   It is possible to return even more information which will be stored
#'   in the optimization path. To achieve this, simply append the attribute \dQuote{extras}
#'   to the return value of the target function. This has to be a named list of scalar values.
#'   Each of these values will be stored additionally in the optimization path.
#' @param design [\code{data.frame}]\cr
#'   Initial design as data frame.
#'   If the y-values are not already present in design, mbo will evaluate the points.
#'   If the parameters have corresponding trafo functions, the design must not be transformed before it is passed!
#'   Functions to generate designs are available in \code{ParamHelpers}: \code{\link[ParamHelpers]{generateDesign}}, \code{\link[ParamHelpers]{generateGridDesign}}, \code{\link[ParamHelpers]{generateRandomDesign}}.
#'   Default is \code{NULL}, which means \code{\link[ParamHelpers]{generateDesign}} is called and a design of size 4 times number of all parameters is created
#'   The points are drawn via \code{\link[lhs]{maximinLHS}} to maximize the minimal distance between design points.
#' @param learner [\code{\link[mlr]{Learner}}]\cr
#'   Regression learner from mlr, which is used as a surrogate to model our fitness function.
#'   If \code{NULL} (default), the default learner is determined as described here: \link{mbo_default_learner}.
#'   
#'   Additional parameters: 
#' @param local Determines whether the explore exploit ratio is computed locally 
#' (i.e. using the seen points from the infill optimization of each iteration) 
#' or globally (i.e. using the seen points from the infill optimization of all 
#' iterations so far, that is, untill the current iteration). Default is TRUE.
#' 
#'
#'
xplxplMBO = function(fun = NULL,
                        design = NULL,
                        learner = NULL,
                        control = NULL,
                        show.info = getOption("mlrMBO.show.info", TRUE),
                        more.args = list(),
                        scope = "local",
                        dist.based = TRUE
                     ) {
  
  # input checking of scope and weighted.ML (rest will be done in mbo())
  assert_character(scope)
  ## TODO check if infill opt matches class of initial data AND parameter set of objfun
  ## TODO test cmaes/ea for non-numeric data AND parameter space
  # Fehler in `[.data.frame`(plyr::ldply(result$population.trace), , 1:length(rep.pids)) : 
  # nicht definierte Spalten gewählt
  

  # change parsed infill optimization function so that it stores the seen points:
  control$infill.opt = editInfillOpt(control$infill.opt)
  
  # make the mbo store surrogate models of each iteration
  control$store.model.at = 1:(control$iters + 1)
  
  # in case of multi-point proposal, enable storage of all surrogate models and all seen points
  if (!is.null(control$multipoint.method)){
      if(control$multipoint.method == "cl")
        control$multipoint.method = "clXplxpl"
      else
        stop("Error: xplxplMBO currently does not support multi point proposal methods ´cb´ and ´moimbo´ ")
  }
    
  
  # run mbo
  res.mbo = mlrMBO::mbo(fun = fun, design = design, learner = learner, control = control,
                 show.info = show.info)
  
  # comute xplxpl measures and return as data frame
  computeXplxpl(res.mbo = res.mbo, scope = scope, dist.based = dist.based)

}
      




editInfillOpt = function(infill.opt) {
  switch(infill.opt,
         "focussearch" = "focussearchSavepts",
         "focussearchSavepts" = "focussearchSavepts",
         "ea" = "eaSavepts",
         "eaSavepts" = "eaSavepts",
         "cmaes" = "cmaesSavepts",
         "cmaesSavepts" = "cmaesSavepts",
         "nsga2" = stop("xplxplMBO does currently not support nsga2 infill optimization")
  )
}

# Run tests
#testthat::test_file("_Explore_Exploit_Measures/test-xplxplMBO-jr.R")








  