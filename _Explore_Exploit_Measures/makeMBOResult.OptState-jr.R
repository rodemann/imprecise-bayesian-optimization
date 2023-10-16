# modification of makeMBOResult.OptState that allows attaching seen points via getOptStateSeenPoints 
# to MBOResult and all surrogate models (in case of multipoint MBO) via getOptStateAllModels to MBOResult



  
makeMBOResult.OptState <- function (opt.state) 
{
  opt.problem = getOptStateOptProblem(opt.state)
  control = getOptProblemControl(opt.problem)
  final.points = getOptStateFinalPoints(opt.state)
  opt.result = getOptStateOptResult(opt.state)
  x.df = final.points$x
  
  # get seen points and all surrogate models in case of multi-point MBO
  seen.points = getOptStateSeenPoints(opt.state)
  all.models = getOptStateAllModels(opt.state)
  
    
  if (!is.null(x.df) && nrow(x.df) == 1) {
    if (getOptProblemControl(opt.problem)$final.evals > 
        0) {
      ys = evalFinalPoint(opt.state, x.df)
      final.points$y = mean(ys)
    }
    makeS3Obj(c("MBOSingleObjResult", "MBOResult"), x = dfRowToList(x.df, 
                                                                    par.set = getOptProblemParSet(opt.problem), i = 1), 
              y = final.points$y, best.ind = final.points$best.ind, 
              opt.path = getOptStateOptPath(opt.state), resample.results = getOptResultResampleResults(opt.result), 
              final.state = getOptStateState(opt.state), models = getOptResultStoredModels(opt.result), 
              final.opt.state = opt.state, control = control, seen.points = seen.points, all.models = all.models)
  }
  else {
    makeS3Obj(c("MBOMultiObjResult", "MBOResult"), pareto.front = final.points$pareto.front, 
              pareto.set = final.points$pareto.set, pareto.inds = final.points$inds, 
              opt.path = getOptStateOptPath(opt.state), final.state = getOptStateState(opt.state), 
              models = getOptResultStoredModels(opt.result), final.opt.state = opt.state, 
              control = control, seen.points = seen.points)
  }


}


getOptStateSeenPoints <- function(opt.state){
  opt.state$seen.points
}

getOptStateAllModels <- function(opt.state){
  opt.state$all.models
}

environment(makeMBOResult.OptState) <- asNamespace('mlrMBO')
# assures that the function will be able to call other hidden functions from the package.

assignInNamespace("makeMBOResult.OptState", makeMBOResult.OptState, ns = "mlrMBO") 
# assures that the function will be able to call other hidden functions from the package. (credits: TMS, stackoverflow)






