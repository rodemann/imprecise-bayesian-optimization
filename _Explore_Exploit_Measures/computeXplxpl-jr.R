
# 
#' Computes xplxpl statistics (is called from xplxpl() and xplxplMBO())
#'  
#' @param res.mbo parsed mbo object from mbo()
#' @param scope parsed scope argument from xplxpl() and xplxplMBO(), either "local" or "global"
#' 
computeXplxpl <- function(res.mbo, scope = "local", dist.based = TRUE) {
  
  # extract info from parsed mbo object
  vars.total = colnames(res.mbo$seen.points[[1]]) # this way works for both vectorized and scalarized variables
  targets = res.mbo$control$y.name
  pars = vars.total[which(vars.total != targets)] 
  par.set = res.mbo$opt.path$par.set
  opt.path = as.data.frame(res.mbo$opt.path, discretes.as.factor = FALSE)
  n.init = length(opt.path$prop.type[opt.path$prop.type == "initdesign"])
  n.tot = nrow(opt.path)
  n.iter = res.mbo$control$iters * res.mbo$control$propose.points
  initial.data = opt.path[1:n.init,pars]
  infill.crit = res.mbo$control$infill.crit
  infill.crit.id = res.mbo$control$infill.crit$id
  n.objectives = res.mbo$control$n.objectives
  n.prop.points = res.mbo$control$propose.points
  
  
  ## Compute Explore Exploit Statistics
  
  
  # 1. compute explore/exploit ratio and quantile based on variance (standard error)
  # we only need n.prop.points for single objective case, as multi point proposal and multi objective are not possible
  
  if (n.objectives > 1) {
    xpl.xpl.se = switch(scope,
                        "local" = xplxplSEMultiLocal(n.iter = n.iter, res.mbo = res.mbo, opt.path = opt.path, 
                                                     n.init = n.init, pars = pars),
                        "global" = xplxplSEMultiGlobal(n.iter = n.iter, res.mbo = res.mbo, opt.path = opt.path, 
                                                       n.init = n.init, pars = pars)    
    )
  }
  else {
    xpl.xpl.se = switch(scope,
                        "local" = xplxplSESingleLocal(n.iter = n.iter, res.mbo = res.mbo, opt.path = opt.path, 
                                                      n.init = n.init, pars = pars, n.prop.points = n.prop.points),
                        "global" = xplxplSESingleGlobal(n.iter = n.iter, res.mbo = res.mbo, opt.path = opt.path, 
                                                        n.init = n.init, pars = pars, n.prop.points = n.prop.points)
    )
  }   
  
  # 2. Compute distance-based statistics
  if(dist.based == TRUE){
  if (isNumeric(par.set)) 
    xpl.xpl.distance =  xplxplDistNum(opt.path = opt.path, pars = pars, n.init = n.init, n.tot = n.tot)
  else # in case of non-numeric (i.e. categorical) data, we need to use gower distance
    xpl.xpl.distance =  xplxplDistCat(opt.path = opt.path, pars = pars, n.init = n.init, n.tot = n.tot)
  }
  
  
  # 3. finalize data.frame so that it can be returned and return it
  returned_list <- list()
  if(dist.based == TRUE)
    returned_list$xplxpl_df <- cbind(xpl.xpl.se, xpl.xpl.distance, "Iteration" = 1:n.iter)
  else
    returned_list$xplxpl_df <- cbind(xpl.xpl.se, "Iteration" = 1:n.iter)

  returned_list$res_mbo <- res.mbo
  returned_list
  }







