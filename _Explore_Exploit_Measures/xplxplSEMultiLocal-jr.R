

#' computes explore exploit ratio locally (i.e. using the seen points from the 
#' infill optimization of each iteration) for multi-crit optimization
#' @param n.iter number of iterations in the mbo process
#' @param res.mbo MBOobject (Resulting Object from mbo())
#' @param n.init number of initial training points in mbo
#' @param opt.path Optimization path of mbo, obtained via res.mbo
#' 
xplxplSEMultiLocal = function(n.iter, res.mbo, n.init, opt.path, pars) {
  xpl.xpl.rat = vector(length = n.iter) #note: would be nice to vectorize this loop 
  distr.val.se.prop = vector(length = n.iter) 
  for (i in 1:(n.iter)) { 
    
    # extract model object from mbo object from iteration i
    sm = res.mbo$models[i][[1]] 
    # extract seen points from mbo object from iteration i
    seen.points = res.mbo$seen.points[[i]] 
    # compute standard error of seen points for all objective functions
    sm.seen = lapply(sm, FUN = predict, newdata = seen.points)
    se.seen = lapply(sm.seen, function(x) x[["data"]][["se"]])
    # average over those
    se.seen.mean = mean(unlist(se.seen))
    
    # compute standard error of actually proposed point
    prop.points = opt.path[pars]
    prop.point = prop.points[n.init + i, pars]
    prop.point.df = data.frame(prop.point)
    names(prop.point.df) = pars
    sm.prop = lapply(sm, FUN = predict, newdata = prop.point.df)
    se.prop = lapply(sm.prop, function(x) x[["data"]][["se"]])
      
    # average over those
    se.prop.mean = mean(unlist(se.prop))
    
    # compute ratio of mean(se.prop) to mean(se.seen)
    xpl.xpl.rat[i] = se.prop.mean/se.seen.mean
    # compute empirical cumulative distribution of mean standard errors
    emp.cum.distr = ecdf(se.seen.mean)
    # compute (estimate) distr.val of se.prop 
    distr.val.se.prop[i] = emp.cum.distr(se.prop.mean)
    
  }
  xpl.xpl = data.frame("Local Standard Error Ratio" = xpl.xpl.rat, "Local Standard Error Distribution Value" = distr.val.se.prop)
  
  return(xpl.xpl)
}

