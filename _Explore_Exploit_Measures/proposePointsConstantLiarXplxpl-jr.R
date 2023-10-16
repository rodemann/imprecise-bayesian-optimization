# modification of proposePointsConstantLiar that allows a) storage of seen points in proposePointsByInfillOptimization
# and b) the storage of "intermediate" models inside one global iteration für multi point proposal in proposePointsByInfillOptimization

proposePointsConstantLiarXplxpl = function(opt.state) {
  
  opt.problem = mlrMBO:::getOptStateOptProblem(opt.state)
  model = mlrMBO:::getOptStateModels(opt.state)$models[[1L]]
  par.set = mlrMBO:::getOptProblemParSet(opt.problem)
  control = mlrMBO:::getOptProblemControl(opt.problem)
  opt.path = mlrMBO:::getOptStateOptPath(opt.state)
  
  npoints = control$propose.points
  liar = control$multipoint.cl.lie
  # copy control, and propose 1 point each
  control2 = control
  control2$propose.points = 1L
  # copy opt.path to store lies
  opt.path2 = mlrMBO:::deepCopyOptPath(opt.path)
  lie = liar(getOptPathY(opt.path, control$y.name))
  dob = max(getOptPathDOB(opt.path)) + 1L
  props = list()
  for (i in seq_len(npoints)) {
    # propose point, add to opt.path2 with y = lie, then update model
    props[[i]] = proposePointsByInfillOptimization(opt.state, control = control2, 
                                                   opt.path = opt.path2, models = list(model), 
                                                   iter.multi.pts = i, npoints = npoints)
    if (i==npoints) break # we don't need to update the model when we aleady have the n-th proposal
    x = dfRowToList(props[[i]]$prop.points, par.set, 1)
    addOptPathEl(opt.path2, x = x, y = lie, dob = dob)
    rt = mlrMBO:::makeTaskSingleObj(opt.path2, control)
    model = train(model$learner, rt)
  }
  mlrMBO:::joinProposedPoints(props)
}




