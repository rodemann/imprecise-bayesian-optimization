library(mlr)
library(mlrMBO)
library(DiceKriging)
library(checkmate)

source("_Explore_Exploit_Measures/xplxplMBO-jr.R")


#' Explore_Exploit Function
#' computes a statistic for each iteration indicating whether the mbo exploits or explores
#' Beware: res.mbo must contain BOTH all surrogate models, i.e. in makeMBOControl() store.model.at must have
#' been set to 1:iters+1 (default is iters+1), AND all candidate (= seen) points from
#' the infill optimization in each iteration. If the parsed res.mbo does not
#' contain those, xplxpl() will ask the user to rerun the mbo using xplxplMBO()
#'   
#' @param res.mbo MBOobject (Resulting Object from mbo())
#' @param local Determines whether the explore exploit ratio is computed locally 
#' (i.e. using the seen points from the infill optimization of each iteration) 
#' or globally (i.e. using the seen points from the infill optimization of all 
#' iterations so far, that is, untill the current iteration). Default is "local".
#'
xplxpl = function(res.mbo, scope = "local") {
  
  # some input checking:
  assert(check_class(res.mbo, "MBOSingleObjResult"), check_class(res.mbo,
  "MBOMultiObjResult"), check_class(res.mbo, "MBOResult"))
  assert_character(scope)
  
  
  # extract info from parsed mbo object
  targets = res.mbo$control$y.name
  pars = names(res.mbo$opt.path$par.set$pars)
  par.set = res.mbo$opt.path$par.set
  opt.path = as.data.frame(res.mbo$opt.path, discretes.as.factor = FALSE)
  n.init = length(opt.path$prop.type[opt.path$prop.type == "initdesign"])
  n.tot = nrow(opt.path)
  n.iter = res.mbo$control$iters * res.mbo$control$propose.points
  initial.data = opt.path[1:n.init,pars]
  objfun = res.mbo[["final.opt.state"]][["opt.problem"]][["fun"]]
  learner = res.mbo[["models"]][[as.character(n.iter)]][["learner"]]
  
  # further input checking: check if res.mbo contains all surrogate models and 
  # seen points from infill optimization
  n.iter = res.mbo$control$iters * res.mbo$control$propose.points
  if ((length(res.mbo$models) != n.iter + 1) | (is.null(res.mbo$seen.points))) {
    YesNo = askYesNo("The passed MBO object does not contain necessary information (seen points and all surrogate models) for explore exploit statistics to be computed. Do you want to rerun the MBO (same settings, different RNG seed) and store all surrogate models in order to compute the explore/exploit statistic?")
    if (YesNo == TRUE) {
      xpl.xpl = xplxplMBO(fun = objfun,
                design = initial.data,
                learner = learner,
                control = res.mbo$control,
                scope = "local"
                )
      return(xpl.xpl)
           }
    return("No analysis of exploit/explore trade-off possible. Blame yourself.")
  }
  
  # comute xplxpl measures and return as data frame
  computeXplxpl(res.mbo = res.mbo, scope = scope)
  
}






