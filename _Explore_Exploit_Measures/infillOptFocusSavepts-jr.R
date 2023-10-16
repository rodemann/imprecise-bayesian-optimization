# Focus search function for infill optimization that additionaly stores those points
# that have been proposed/evaluated/seen during the focus search porcess
# 
# infillOptFocus:
# Random search, where we shrink the region of interest after restarts
# around the currently best point. only numeric / ints are currently "shrunken"
# works for ALL parameter sets
#
#
#
# @param infill.crit [\code{function}]\cr
#   Infill criterion function.
# @param models [\code{\link{WrappedModel}}]\cr
#   Model fitted on designs.
# @param control [\code{\link{MBOControl}}]\cr
#   Control object.
# @param par.set [\code{\link[ParamHelpers]{ParamSet}}]\cr
#   Parameter set.
# @param opt.path [\code{\link[ParamHelpers{OptPath}}]\cr
#   Optimization path / archive.
# @param designs [list of \code{data.frame}]\cr
#   List with usually one element containing the Design of already visited points.
# @param iter [\integer(1)]
#   Current iteration
# @param ... [\code{ANY}]\cr
#   Additional arguments passed to \code{infill.crit}.
# @return [\code{data.frame}]. One proposed point that should be evaluated.

library(lhs) #required for argument "fun" in generateDesign() 

infillOptFocusSavepts = function(infill.crit, models, control, par.set, opt.path, 
                                  designs, iter, ...) {
  
   
  global.y = Inf
  
  # empty matrix to store evaluated points
  vars = colnames(designs[[1]]) # this way works for both vectorized and scalarized variables
  y.names = control$y.name
  x.names = vars[which(vars != y.names)] 
  nrow = control$infill.opt.restarts * control$infill.opt.focussearch.maxit * control$infill.opt.focussearch.points
  seen.points = generateDesign(nrow, par.set, randomLHS)
  seen.points.dummies = matrix(NA, nrow = nrow, ncol = length(x.names))
  seen.points[1:nrow,] = seen.points.dummies


 
  
  # restart the whole crap some times
  for (restart.iter in seq_len(control$infill.opt.restarts)) {
    # copy parset so we can shrink it
    ps.local = par.set
    
    
    
    # do iterations where we focus the region-of-interest around the current best point
    for (local.iter in seq_len(control$infill.opt.focussearch.maxit)) {
      # predict on design where NAs were imputed, but return proposed points with NAs
      newdesign = generateDesign(control$infill.opt.focussearch.points, ps.local, randomLHS)
      
      ## SAVEPTS: store points that have been evaluated by focus search
      global.iter = control$infill.opt.focussearch.maxit * nrow(newdesign) * (restart.iter - 1) +  nrow(newdesign) * (local.iter - 1) 
      seen.points[(global.iter + 1):(global.iter + nrow(newdesign)),] = as.data.frame(newdesign)
      
      # convert to param encoding our model was trained on and can use
      newdesign = convertDataFrameCols(newdesign, ints.as.num = TRUE, logicals.as.factor = TRUE, factors.as.char = TRUE)
      
      
      y = infill.crit(newdesign, models, control, ps.local, designs, iter, ...)
      
      
      # get current best value
      local.index = getMinIndex(y, ties.method = "random")
      local.y = y[local.index]
      local.x.df = newdesign[local.index, , drop = FALSE]
      local.x.list = dfRowToList(recodeTypes(local.x.df, ps.local), ps.local, 1)
      
      # if we found a new best value, store it
      if (local.y < global.y) {
        global.x.df = local.x.df
        global.y = local.y
      }
      
      # now shrink ps.local object so we search more locally
      ps.local$pars = lapply(ps.local$pars, function(par) {
        # only shrink when there is a value
        val = local.x.list[[par$id]]
        if (!isScalarNA(val)) {
          if (isNumeric(par)) {
            # shrink to range / 2, centered at val
            range = par$upper - par$lower
            par$lower = pmax(par$lower, val - (range / 4))
            par$upper = pmin(par$upper, val + (range / 4))
            if (isInteger(par)) {
              par$lower = floor(par$lower)
              par$upper = ceiling(par$upper)
            }
          } else if (isDiscrete(par)) {
            # randomly drop a level, which is not val
            if (length(par$values) > 1L) {
              val.names = names(par$values)
              to.del = sample(which(val.names != val), 1L)
              par$values = par$values[-to.del]
            }
          }
        }
        return(par)
      })
    }
  }
  # return prosed points (and convert back, see below) and seen points
  recodeTypes_return(global.x.df, par.set, seen.points)

}

# as we operate on other types for the learner (ints are nums, logs are factors),
# we have to convert back sometimes for dfRowsToList to work
recodeTypes = function(df, par.set) {
  types = unlist(lapply(par.set$pars, function(p) rep(p$type, p$len)))
  for (i in seq_col(df)) {
    if (types[i] %in% c("integer", "integervector"))
      df[, i] = as.integer(df[,i])
    else if (types[i] %in% c("logical", "logicalvector"))
      df[, i] = as.logical(as.character(df[,i]))
  }
  return(df)
}

# modify recodeTypes for final return of function (which should include seen points) 
recodeTypes_return = function(df, par.set, seen.points) {
  df = recodeTypes(df, par.set)
  result = list("prop.points" = df, "seen.points" = seen.points)
  return(result)
}


