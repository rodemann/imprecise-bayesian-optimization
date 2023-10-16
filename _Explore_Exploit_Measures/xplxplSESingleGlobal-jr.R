

#' computes explore exploit ratio globally ((i.e. using the seen points from the
#' infill optimization of all iterations so far, that is, untill the current iteration)
#' @param n.iter number of iterations in the mbo process
#' @param res.mbo MBOobject (Resulting Object from mbo())
#' @param n.init number of initial training points in mbo
#' @param opt.path Optimization path of mbo, obtained via res.mbo
#' 
xplxplSESingleGlobal = function(n.iter, res.mbo, n.init, opt.path, pars, n.prop.points) {
  require(purrr)
  xpl.xpl.rat = vector(length = n.iter) #note: would be nice to vectorize this loop 
  distr.val.se.prop = vector(length = n.iter) 
  for (i in 1:(n.iter)) { 
    
    # extract model object from mbo object that was used in iteration i
    if(n.prop.points == 1)
      sm = res.mbo$models[i][[1]] 
    # in case of multi-point MBO use res.mbo$all.models to get all surrogate models (also from within one iter)
    else 
      sm = res.mbo$all.models[i][[1]] 
    
    # make sure predict type is "se" (not given if infill.crit = crit.mr)
    sm$learner$predict.type = "se"
    # extract seen points from mbo object from iteration 1 untill iteration i
    all.seen.points = reduce(res.mbo$seen.points[1:i], rbind)
    # compute standard error of seen points
    sm.seen = predict(object = sm, newdata = all.seen.points) 
    se.seen = sm.seen[["data"]][["se"]]
    # compute standard error of actually proposed point
    if (res.mbo$control$infill.crit$id == "mean"  | res.mbo$control$infill.crit$id == "se") # infill criterion mean response does not estimate se
    {
      prop.points = opt.path[pars]
      prop.point = prop.points[n.init + i, pars]
      prop.point.df = data.frame(prop.point)
      names(prop.point.df) = pars
      se.prop = predict(object = sm, newdata = prop.point.df)[["data"]][["se"]]
    }
      else  
      se.prop = opt.path[["se"]][[n.init + i]]
    # compute ratio of se.prop to mean of seen points' se
    xpl.xpl.rat[i] = se.prop/mean(se.seen)
    # compute empirical cumulative distribution of standard errors
    emp.cum.distr = ecdf(se.seen)
    # compute (estimate) distr.val of se.prop 
    distr.val.se.prop[i] = emp.cum.distr(se.prop)
    
    
  }
  xpl.xpl = data.frame("Global Standard Error Ratio" = xpl.xpl.rat, "Global Standard Error Distribution Value" = distr.val.se.prop)
  
  return(xpl.xpl)
}

