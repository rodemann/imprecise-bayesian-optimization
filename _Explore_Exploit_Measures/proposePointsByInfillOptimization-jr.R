# modification of proposePointsByInfillOptimization that stores seen points (from infill.opt ...Savepts) in opt.state
# and surrogate models in case of multi point MBO


library(BBmisc)

# use iter.multi.pts and npoints only for multi point proposal
proposePointsByInfillOptimization <- function (opt.state, par.set = NULL, control = NULL, opt.path = NULL, 
          models = NULL, designs = NULL, iter.multi.pts = 1, npoints = 1, ...) 
{
  opt.problem = getOptStateOptProblem(opt.state)
  models = models %??% getOptStateModels(opt.state)$models
  if (inherits(models, "WrappedModel")) 
    models = list(models)
  par.set = par.set %??% getOptProblemParSet(opt.problem)
  designs = designs %??% getOptStateDesigns(opt.state)
  if (inherits(designs, "data.frame")) 
    designs = list(designs)
  control = control %??% getOptProblemControl(opt.problem)
  opt.path = opt.path %??% getOptStateOptPath(opt.state)
  iter = getOptStateLoop(opt.state)
  infill.crit.id = getMBOInfillCritId(control$infill.crit)
  progress = getOptStateProgress(opt.state)
  n = control$propose.points
  prop.type = rep(paste0("infill_", infill.crit.id), n)
  ch = checkFailedModels(models, par.set, n, control = control)
  if (!ch$ok) 
    return(ch$prop)
  infill.crit.fun = control$infill.crit$fun
  infill.opt.fun = getInfillOptFunction(control$infill.opt)
  secs = measureTime({
    result = infill.opt.fun(infill.crit.fun, models = models, 
                                 control = control, par.set = par.set, opt.path = opt.path, 
                                 designs = designs, iter = iter, progress = progress, 
                                 ...)
  })
  prop.points <- getProposedPoints(result, infill.opt = control$infill.opt)
  # get global iteration that takes multiple points per iteration into account
  if (npoints > 1)
    global.iter = (iter - 1)*npoints + iter.multi.pts
  else
    global.iter = iter
  # store seen points in opt.state according to their global iters
  opt.state$seen.points[[global.iter]] = storeSeenPoints(result, infill.opt = control$infill.opt)
  # do the same thing for the surrogate model in case of multipoint proposal
  if (npoints > 1)
    opt.state$all.models[[global.iter]] = models[[1]]
  
  prop.points.converted = convertDataFrameCols(prop.points, 
                                               ints.as.num = TRUE, logicals.as.factor = TRUE)
  crit.vals = infill.crit.fun(prop.points.converted, models, 
                              control, par.set, designs, iter, progress = progress, 
                              attributes = TRUE, ...)
  crit.components = attr(crit.vals, "crit.components")
  crit.vals = matrix(crit.vals, ncol = 1L)
  makeProposal(control = control, prop.points = prop.points, 
               propose.time = secs, prop.type = prop.type, crit.vals = crit.vals, 
               crit.components = crit.components)
}


getProposedPoints <- function(result, infill.opt){
  if (grepl("Savepts", infill.opt)){
    return(result$prop.points)
  }
  return(result)
  
}

storeSeenPoints <- function(result, infill.opt){
  if (grepl("Savepts", infill.opt)){
    return(result$seen.points)
  }
}


# assure that the function will be able to call other hidden functions from the package.
environment(proposePointsByInfillOptimization) <- asNamespace('mlrMBO')

# assure that the function will be able to call other hidden functions from the package. (credits: TMS, stackoverflow)
assignInNamespace("proposePointsByInfillOptimization", proposePointsByInfillOptimization, ns = "mlrMBO") 



