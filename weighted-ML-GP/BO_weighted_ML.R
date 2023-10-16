BO_weighted_ML <- function(fun = NULL,
                            design = NULL,
                            learner = NULL,
                            control = NULL,
                            show.info = getOption("mlrMBO.show.info", TRUE),
                            more.args = list(),
                            scope = "local") {
 
   # reset km-estimate-function and logLikFun in {DiceKriging} (in case modified version was sourced before)
  back_to_normal_kmEstimate()

  # reset km-estimate-function and logLikFun in {DiceKriging} (in case modified version was sourced before)
  #back_to_normal_kmEstimate()
  #back_to_normal_LogLikFun()
  ### NOTE: If alter_kmEstimate gives troubles, you may want to switch back to altering loglikFun
  
  global_iters <- control$iters
  control <- setMBOControlTermination(control, iters = 1L)  
  # initialize xplxpl values
  n_init <- nrow(design)
  xplxpl_init_design <- rep(1, n_init)
  xplxpl_prop_total <- vector()
  
for (i in 1:global_iters) {
  res_xplxpl <- xplxplMBO(fun, design, learner, control, show.info, more.args, scope, dist.based = FALSE)  
  # retrieve objects from result object 
  opt_path <- as.data.frame(res_xplxpl$res_mbo$opt.path)
  inc_design_size <- length(opt_path$prop.type[opt_path$prop.type == "initdesign"])
  proposed_point <- opt_path[inc_design_size+1,]
  proposed_x <- data.frame(x = proposed_point$x)
  proposed_y <- data.frame(y = proposed_point$y)
  xplxpl_prop_i <- res_xplxpl$xplxpl_df$Local.Standard.Error.Distribution.Value
  xplxpl_prop_total <- append(xplxpl_prop_total, xplxpl_prop_i)
  xplxpl <- c(xplxpl_init_design, xplxpl_prop_total)
  # compute weights for ML estimation of prior hyperparams in Gaussian Process
  weights <- xplxpl/sum(xplxpl)
  # update design for next BO proposal
  design <- dplyr::bind_rows(design, proposed_x)

    # alter kmEstimate Function in DiceKriging such that it takes xplxpl-weights into account
  #weights_envir <- environment()
  #env_contents <- list(ls(env))
  # put weights in global environment so that it can be assessed by loglikFun
  assign("weights", weights, envir = globalenv())
  alter_LogLikFun()
} 
  # reset log-likelihood-function in {DiceKriging}
  #back_to_normal_kmEstimate()
  back_to_normal_LogLikFun()
  # attach xplxpl-values and return the results object
  res_xplxpl$xplxpl_prop_total <- xplxpl_prop_total
  res_xplxpl
  
}
  

alter_kmEstimate <- function(weights, env){
  #env_name <- environmentName(env)
  kmEstimate <- function(model, envir) {
      
      ## begin edit
      #force(weights_envir)
      weights <- eval(weights, envir = weights_envir)
      if(nrow(X) == length(weights)){
      # compute weighted x values
      model@X <- model@X * weights
      # and weighted y values
      model@y <- model@y * weights
      warning("weights used")
      
      } #else --> in this case, mlrMBO::mbo() estimates final model that will be returned
      #warning(list(weights))  
      #warning(list(X,model@X))  
       
      ## end edit
      
      X <- model@X
      y <- model@y
      F <- model@F
      
      # these 3 lines for retro-compatibility
      if (model@case=="NoNugget") model@case <- "LLconcentration_beta_sigma2"
      if (model@case=="1Nugget") model@case <- "LLconcentration_beta_v_alpha"
      if (model@case=="Nuggets") model@case <- "LLconcentration_beta"
      
      if (model@case=="LLconcentration_beta_sigma2") case <- "Default"
      if (model@case=="LLconcentration_beta_v_alpha") case <- "Nugget"
      if (model@case=="LLconcentration_beta") case <- "Noisy"
      
      if (model@case=="LLconcentration_beta") {
        if (model@covariance@nugget.flag & !model@covariance@nugget.estim) nugget <- model@covariance@nugget
        if (model@noise.flag) nugget <- model@noise.var
      }
      

      if (is.element(model@method, c("MLE", "PMLE"))){
        fn <- logLikFun
        fnscale <- -1
      } else if (model@method=="LOO") {
        fn <- leaveOneOutFun
        fnscale <- 1
      } 
      if (model@gr==TRUE) {
        if (is.element(model@method, c("MLE", "PMLE"))){
          gr <- logLikGrad
        } else if (model@method=="LOO") {
          gr <- leaveOneOutGrad
        } 
      } else {
        gr <- NULL
      }

      # initialization: starting points and boundaries
      
      multistart <- model@control$multistart
      # below: useful if someone calls directly kmEstimate
      if (length(multistart)==0){
        model@control$multistart <- multistart <- 1
      }
      
      initList <- switch(case,
                         Default = kmNoNugget.init(model, fn, fnscale),
                         Nugget = km1Nugget.init(model),
                         Noisy = kmNuggets.init(model))
      
      lower <- model@lower <- as.numeric(initList$lower)
      upper <- model@upper <- as.numeric(initList$upper)
      parinit <- initList$par
      lp <- nrow(parinit)
      
      # printing
      control <- model@control
      if (control$trace!=0) {
        cat("\n")
        cat("optimisation start\n")
        cat("------------------\n")
        cat("* estimation method   :", model@method, "\n")
        cat("* optimisation method :", model@optim.method, "\n")
        cat("* analytical gradient :")
        if (is.null(gr)) {
          cat(" not used\n")}
        else cat(" used\n")
        cat("* trend model : "); print(model@trend.formula, showEnv=FALSE)
        cat("* covariance model : \n")
        cat("  - type : ",  model@covariance@name, "\n") 
        if (case=="Default") { 
          cat("  - nugget : NO\n")
          cat("  - parameters lower bounds : ", lower, "\n")
          cat("  - parameters upper bounds : ", upper, "\n")
        } else if (case=="Nugget") { 
          cat("  - nugget : unknown homogenous nugget effect \n")
          #if (!is.null(nugget)) cat("with initial value : ", nugget)
          cat("  - parameters lower bounds : ", lower[1:(lp-1)], "\n")
          cat("  - parameters upper bounds : ", upper[1:(lp-1)], "\n")
          cat("  - upper bound for alpha   : ", upper[lp],"\n")
        } else if (case=="Noisy") {  # technically: includes the "known nugget" case
          if (model@covariance@nugget.flag) {
            cat("  - nugget :", model@covariance@nugget, "\n")
          } else {
            cat("  - noise variances :\n")
            print(model@noise.var)
          }
          cat("  - parameters lower bounds : ", lower[1:(lp-1)], "\n")
          cat("  - parameters upper bounds : ", upper[1:(lp-1)], "\n")
          cat("  - variance bounds : ", c(lower[lp], upper[lp]), "\n")
        }
        cat("  - best initial criterion value(s) : ", initList$value, "\n")
        if (model@optim.method=="BFGS") cat("\n")     
      } # end printing
      
      # optimization
      
      if (model@optim.method=="BFGS") {
        BFGSargs <- c("trace", "parscale", "ndeps", "maxit", "abstol", "reltol", "REPORT", "lnm", "factr", "pgtol")
        commonNames <- intersect(BFGSargs, names(control))
        controlChecked <- control[commonNames]
        if (length(control$REPORT)==0) {
          controlChecked$REPORT <- 1
        }
        
        forced <- list(fnscale = fnscale)
        controlChecked[names(forced)] <- forced
        
        # multistart in parallel with foreach
        multistart <- control$multistart
        
        if (multistart==1){
          model@parinit <- parinit <- as.numeric(parinit[, 1])
          model@covariance <- initList$cov[[1]]
          o <- optim(par = parinit, fn = fn, gr = gr,
                     method = "L-BFGS-B", lower = lower, upper = upper,
                     control = controlChecked, hessian = FALSE, model, envir=envir)
          model@control$convergence <- o$convergence
        } else {
          # multistart with foreach
          if (!requireNamespace("foreach", quietly = TRUE)){
            stop("Package \"foreach\" not found")
          }
          olist <- foreach::"%dopar%"(foreach::foreach(i=1:multistart, 
                                                       .errorhandling='remove'), {
                                                         model@covariance <- initList$cov[[i]]
                                                         optim(par = parinit[, i], fn = fn, gr = gr,
                                                               method = "L-BFGS-B", lower = lower, upper = upper,
                                                               control = controlChecked, hessian = FALSE, model, envir=envir)
                                                       })
          
          
          # get the best result
          nOK <- length(olist)
          if (nOK == 0) stop("All model optimizations returned an error") 
          bestValue <- Inf
          bestIndex <- NA
          vecValue <- c()
          for (i in 1:nOK){
            currentValue <- fnscale * olist[[i]]$value
            vecValue <- c(vecValue, currentValue)
            if (currentValue < bestValue) {
              o <- olist[[i]]
              bestValue <- currentValue
              bestIndex <- i
            }
          } # end multistart
          
          model@covariance <- initList$cov[[bestIndex]]
          parinit <- parinit[, bestIndex]
          model@parinit <- as.numeric(parinit)
          model@control$convergence <- o$convergence
          
          # we need to initiate a final optimization from the best point
          # in order to recompute the intermediate variables to be stored in environment 'envir'
          # (during the 'foreach' loop, they were stored in a copy of 'envir')
          
          controlChecked$maxit <- 0
          controlChecked$trace <- 0
          o <- optim(par = o$par, fn = fn, gr = gr,
                     method = "L-BFGS-B", lower = lower, upper = upper,
                     control = controlChecked, hessian = FALSE, model, envir=envir)
          
          if (control$trace!=0){
            cat("\n")
            cat("* The", multistart, "best values (multistart) obtained are:\n", vecValue, "\n")
            cat("* The model corresponding to the best one (", bestValue, ") is stored. \n", sep="")
          }
        } # end multistart loop
      } # end BFGS loop
      
      model@control$multistart <- multistart
      
      if (model@optim.method=="gen") {
        if (!requireNamespace("rgenoud", quietly = TRUE)) {
          stop("Package \"rgenoud\" not found")
        }
        genoudArgs <- formals(rgenoud::genoud)
        commonNames <- intersect(names(genoudArgs), names(control))
        genoudArgs[commonNames] <- control[commonNames]
        if (length(control$print.level)==0) {
          genoudArgs$print.level <- control$trace
        } else {
          genoudArgs$print.level <- control$print.level
        }
        
        forced <- list(fn=fn, nvars=length(parinit), max=(fnscale < 0), starting.values=parinit, 
                       Domains=cbind(lower, upper), gr=gr, gradient.check=FALSE, boundary.enforcement=2, 
                       hessian=TRUE, optim.method="L-BFGS-B", model=model, envir=envir)
        
        genoudArgs[names(forced)] <- forced
        genoudArgs$... <- NULL
        
        o <- do.call(rgenoud::genoud, genoudArgs)   
        
      }
      
      model@logLik <- as.numeric(o$value)
      
      
      if (model@method=="LOO"){
        
        model@covariance <- vect2covparam(model@covariance, o$par)
        
        errorsLOO <- envir$errorsLOO
        sigma2LOO <- envir$sigma2LOO
        sigma2.hat <- mean(errorsLOO^2/sigma2LOO)
        sigma.hat <- sqrt(sigma2.hat)
        model@covariance@sd2 <- as.numeric(sigma2.hat)
        
        R <- envir$R
        T <- chol(R)
        x <- backsolve(t(T), y, upper.tri=FALSE)		# x:=(T')^(-1)*y
        M <- backsolve(t(T), F, upper.tri=FALSE)		# M:=(T')^(-1)*F
        
        if (identical(model@known.param, "Trend")) {
          z <- x - M %*% model@trend.coef
        } else {
          l <- lm(x ~ M-1)
          beta.hat <- as.numeric(l$coef)
          #       z <- backsolve(t(T), y - F%*%beta.hat, upper.tri=FALSE)
          model@trend.coef <- beta.hat
          Q <- qr.Q(qr(M))
          H <- Q %*% t(Q)
          z <- x - H %*% x
        }
        
        model@T <- T*sigma.hat
        model@z <- as.numeric(z/sigma.hat)
        model@M <- M/sigma.hat
        
        return(model)   # Fin cas LOO
      } else {
        
        # beta.hat ?
        T <- envir$T
        z <- envir$z
        
        x <- backsolve(t(T), y, upper.tri=FALSE)		# x:=(T')^(-1)*y
        M <- backsolve(t(T), F, upper.tri=FALSE)		# M:=(T')^(-1)*F
        
        if (!identical(model@known.param, "Trend")) {
          l <- lm(x ~ M-1)
          beta.hat <- as.numeric(l$coef)
          model@trend.coef <-beta.hat
        }
        
        if (case=="Default") {
          sigma2.hat <- envir$sigma2.hat
          sigma2.hat <- as.numeric(sigma2.hat)
          sigma.hat <- sqrt(sigma2.hat)
          model@T <- T*sigma.hat
          model@z <- as.numeric(z/sigma.hat)
          model@M <- M/sigma.hat
          model@covariance@sd2 <- as.numeric(sigma2.hat)		
          model@covariance <- vect2covparam(model@covariance, o$par)  
        } else if (case=="Nugget") {
          v <- envir$v 
          v <- as.numeric(v)
          s <- sqrt(v)
          model@T <- T * s
          model@z <- as.numeric(z / s)
          model@M <- M / s
          
          param <- as.numeric(o$par)
          lp <- length(param)
          alpha <- param[lp]
          model@covariance@sd2 <- alpha*v
          model@covariance@nugget <- (1-alpha)*v
          model@covariance <- vect2covparam(model@covariance, param[1:(lp-1)])
          
          model@lower <- model@lower[1:(lp-1)]    
          model@upper <- model@upper[1:(lp-1)]
          model@parinit <- model@parinit[1:(lp-1)]    
        }
        if (case=="Noisy") {
          model@T <- T
          model@z <- as.numeric(z)
          model@M <- M
          
          param <- as.numeric(o$par)
          lp <- length(param)
          model@covariance@sd2 <- param[lp]
          model@covariance <- vect2covparam(model@covariance, param[1:(lp-1)])
          
          model@lower <- model@lower[1:(lp-1)]    
          model@upper <- model@upper[1:(lp-1)]
          model@parinit <- model@parinit[1:(lp-1)]
        }
        
        return(model)   # fin cas MLE, PMLE
      }
    }
  
  #R.utils::reassignInPackage("kmEstimate", pkgName = "DiceKriging", kmEstimate)
  utils::assignInNamespace("kmEstimate", kmEstimate, ns = "DiceKriging") 
}





back_to_normal_kmEstimate <- function(){
  kmEstimate <- function(model, envir) {
      
      X <- model@X
      y <- model@y
      F <- model@F
      
      # these 3 lines for retro-compatibility
      if (model@case=="NoNugget") model@case <- "LLconcentration_beta_sigma2"
      if (model@case=="1Nugget") model@case <- "LLconcentration_beta_v_alpha"
      if (model@case=="Nuggets") model@case <- "LLconcentration_beta"
      
      if (model@case=="LLconcentration_beta_sigma2") case <- "Default"
      if (model@case=="LLconcentration_beta_v_alpha") case <- "Nugget"
      if (model@case=="LLconcentration_beta") case <- "Noisy"
      
      if (model@case=="LLconcentration_beta") {
        if (model@covariance@nugget.flag & !model@covariance@nugget.estim) nugget <- model@covariance@nugget
        if (model@noise.flag) nugget <- model@noise.var
      }
      
      
      if (is.element(model@method, c("MLE", "PMLE"))){
        fn <- logLikFun
        fnscale <- -1
      } else if (model@method=="LOO") {
        fn <- leaveOneOutFun
        fnscale <- 1
      } 
      if (model@gr==TRUE) {
        if (is.element(model@method, c("MLE", "PMLE"))){
          gr <- logLikGrad
        } else if (model@method=="LOO") {
          gr <- leaveOneOutGrad
        } 
      } else {
        gr <- NULL
      }
      
      # initialization: starting points and boundaries
      
      multistart <- model@control$multistart
      # below: useful if someone calls directly kmEstimate
      if (length(multistart)==0){
        model@control$multistart <- multistart <- 1
      }
      
      initList <- switch(case,
                         Default = kmNoNugget.init(model, fn, fnscale),
                         Nugget = km1Nugget.init(model),
                         Noisy = kmNuggets.init(model))
      
      lower <- model@lower <- as.numeric(initList$lower)
      upper <- model@upper <- as.numeric(initList$upper)
      parinit <- initList$par
      lp <- nrow(parinit)
      
      # printing
      control <- model@control
      if (control$trace!=0) {
        cat("\n")
        cat("optimisation start\n")
        cat("------------------\n")
        cat("* estimation method   :", model@method, "\n")
        cat("* optimisation method :", model@optim.method, "\n")
        cat("* analytical gradient :")
        if (is.null(gr)) {
          cat(" not used\n")}
        else cat(" used\n")
        cat("* trend model : "); print(model@trend.formula, showEnv=FALSE)
        cat("* covariance model : \n")
        cat("  - type : ",  model@covariance@name, "\n") 
        if (case=="Default") { 
          cat("  - nugget : NO\n")
          cat("  - parameters lower bounds : ", lower, "\n")
          cat("  - parameters upper bounds : ", upper, "\n")
        } else if (case=="Nugget") { 
          cat("  - nugget : unknown homogenous nugget effect \n")
          #if (!is.null(nugget)) cat("with initial value : ", nugget)
          cat("  - parameters lower bounds : ", lower[1:(lp-1)], "\n")
          cat("  - parameters upper bounds : ", upper[1:(lp-1)], "\n")
          cat("  - upper bound for alpha   : ", upper[lp],"\n")
        } else if (case=="Noisy") {  # technically: includes the "known nugget" case
          if (model@covariance@nugget.flag) {
            cat("  - nugget :", model@covariance@nugget, "\n")
          } else {
            cat("  - noise variances :\n")
            print(model@noise.var)
          }
          cat("  - parameters lower bounds : ", lower[1:(lp-1)], "\n")
          cat("  - parameters upper bounds : ", upper[1:(lp-1)], "\n")
          cat("  - variance bounds : ", c(lower[lp], upper[lp]), "\n")
        }
        cat("  - best initial criterion value(s) : ", initList$value, "\n")
        if (model@optim.method=="BFGS") cat("\n")     
      } # end printing
      
      # optimization
      
      if (model@optim.method=="BFGS") {
        BFGSargs <- c("trace", "parscale", "ndeps", "maxit", "abstol", "reltol", "REPORT", "lnm", "factr", "pgtol")
        commonNames <- intersect(BFGSargs, names(control))
        controlChecked <- control[commonNames]
        if (length(control$REPORT)==0) {
          controlChecked$REPORT <- 1
        }
        
        forced <- list(fnscale = fnscale)
        controlChecked[names(forced)] <- forced
        
        # multistart in parallel with foreach
        multistart <- control$multistart
        
        if (multistart==1){
          model@parinit <- parinit <- as.numeric(parinit[, 1])
          model@covariance <- initList$cov[[1]]
          o <- optim(par = parinit, fn = fn, gr = gr,
                     method = "L-BFGS-B", lower = lower, upper = upper,
                     control = controlChecked, hessian = FALSE, model, envir=envir)
          model@control$convergence <- o$convergence
        } else {
          # multistart with foreach
          if (!requireNamespace("foreach", quietly = TRUE)){
            stop("Package \"foreach\" not found")
          }
          olist <- foreach::"%dopar%"(foreach::foreach(i=1:multistart, 
                                                       .errorhandling='remove'), {
                                                         model@covariance <- initList$cov[[i]]
                                                         optim(par = parinit[, i], fn = fn, gr = gr,
                                                               method = "L-BFGS-B", lower = lower, upper = upper,
                                                               control = controlChecked, hessian = FALSE, model, envir=envir)
                                                       })
          
          
          # get the best result
          nOK <- length(olist)
          if (nOK == 0) stop("All model optimizations returned an error") 
          bestValue <- Inf
          bestIndex <- NA
          vecValue <- c()
          for (i in 1:nOK){
            currentValue <- fnscale * olist[[i]]$value
            vecValue <- c(vecValue, currentValue)
            if (currentValue < bestValue) {
              o <- olist[[i]]
              bestValue <- currentValue
              bestIndex <- i
            }
          } # end multistart
          
          model@covariance <- initList$cov[[bestIndex]]
          parinit <- parinit[, bestIndex]
          model@parinit <- as.numeric(parinit)
          model@control$convergence <- o$convergence
          
          # we need to initiate a final optimization from the best point
          # in order to recompute the intermediate variables to be stored in environment 'envir'
          # (during the 'foreach' loop, they were stored in a copy of 'envir')
          
          controlChecked$maxit <- 0
          controlChecked$trace <- 0
          o <- optim(par = o$par, fn = fn, gr = gr,
                     method = "L-BFGS-B", lower = lower, upper = upper,
                     control = controlChecked, hessian = FALSE, model, envir=envir)
          
          if (control$trace!=0){
            cat("\n")
            cat("* The", multistart, "best values (multistart) obtained are:\n", vecValue, "\n")
            cat("* The model corresponding to the best one (", bestValue, ") is stored. \n", sep="")
          }
        } # end multistart loop
      } # end BFGS loop
      
      model@control$multistart <- multistart
      
      if (model@optim.method=="gen") {
        if (!requireNamespace("rgenoud", quietly = TRUE)) {
          stop("Package \"rgenoud\" not found")
        }
        genoudArgs <- formals(rgenoud::genoud)
        commonNames <- intersect(names(genoudArgs), names(control))
        genoudArgs[commonNames] <- control[commonNames]
        if (length(control$print.level)==0) {
          genoudArgs$print.level <- control$trace
        } else {
          genoudArgs$print.level <- control$print.level
        }
        
        forced <- list(fn=fn, nvars=length(parinit), max=(fnscale < 0), starting.values=parinit, 
                       Domains=cbind(lower, upper), gr=gr, gradient.check=FALSE, boundary.enforcement=2, 
                       hessian=TRUE, optim.method="L-BFGS-B", model=model, envir=envir)
        
        genoudArgs[names(forced)] <- forced
        genoudArgs$... <- NULL
        
        o <- do.call(rgenoud::genoud, genoudArgs)   
        
      }
      
      model@logLik <- as.numeric(o$value)
      
      
      if (model@method=="LOO"){
        
        model@covariance <- vect2covparam(model@covariance, o$par)
        
        errorsLOO <- envir$errorsLOO
        sigma2LOO <- envir$sigma2LOO
        sigma2.hat <- mean(errorsLOO^2/sigma2LOO)
        sigma.hat <- sqrt(sigma2.hat)
        model@covariance@sd2 <- as.numeric(sigma2.hat)
        
        R <- envir$R
        T <- chol(R)
        x <- backsolve(t(T), y, upper.tri=FALSE)		# x:=(T')^(-1)*y
        M <- backsolve(t(T), F, upper.tri=FALSE)		# M:=(T')^(-1)*F
        
        if (identical(model@known.param, "Trend")) {
          z <- x - M %*% model@trend.coef
        } else {
          l <- lm(x ~ M-1)
          beta.hat <- as.numeric(l$coef)
          #       z <- backsolve(t(T), y - F%*%beta.hat, upper.tri=FALSE)
          model@trend.coef <- beta.hat
          Q <- qr.Q(qr(M))
          H <- Q %*% t(Q)
          z <- x - H %*% x
        }
        
        model@T <- T*sigma.hat
        model@z <- as.numeric(z/sigma.hat)
        model@M <- M/sigma.hat
        
        return(model)   # Fin cas LOO
      } else {
        
        # beta.hat ?
        T <- envir$T
        z <- envir$z
        
        x <- backsolve(t(T), y, upper.tri=FALSE)		# x:=(T')^(-1)*y
        M <- backsolve(t(T), F, upper.tri=FALSE)		# M:=(T')^(-1)*F
        
        if (!identical(model@known.param, "Trend")) {
          l <- lm(x ~ M-1)
          beta.hat <- as.numeric(l$coef)
          model@trend.coef <-beta.hat
        }
        
        if (case=="Default") {
          sigma2.hat <- envir$sigma2.hat
          sigma2.hat <- as.numeric(sigma2.hat)
          sigma.hat <- sqrt(sigma2.hat)
          model@T <- T*sigma.hat
          model@z <- as.numeric(z/sigma.hat)
          model@M <- M/sigma.hat
          model@covariance@sd2 <- as.numeric(sigma2.hat)		
          model@covariance <- vect2covparam(model@covariance, o$par)  
        } else if (case=="Nugget") {
          v <- envir$v 
          v <- as.numeric(v)
          s <- sqrt(v)
          model@T <- T * s
          model@z <- as.numeric(z / s)
          model@M <- M / s
          
          param <- as.numeric(o$par)
          lp <- length(param)
          alpha <- param[lp]
          model@covariance@sd2 <- alpha*v
          model@covariance@nugget <- (1-alpha)*v
          model@covariance <- vect2covparam(model@covariance, param[1:(lp-1)])
          
          model@lower <- model@lower[1:(lp-1)]    
          model@upper <- model@upper[1:(lp-1)]
          model@parinit <- model@parinit[1:(lp-1)]    
        }
        if (case=="Noisy") {
          model@T <- T
          model@z <- as.numeric(z)
          model@M <- M
          
          param <- as.numeric(o$par)
          lp <- length(param)
          model@covariance@sd2 <- param[lp]
          model@covariance <- vect2covparam(model@covariance, param[1:(lp-1)])
          
          model@lower <- model@lower[1:(lp-1)]    
          model@upper <- model@upper[1:(lp-1)]
          model@parinit <- model@parinit[1:(lp-1)]
        }
        
        return(model)   # fin cas MLE, PMLE
      }
    }
  
  #R.utils::reassignInPackage("kmEstimate", pkgName = "DiceKriging", kmEstimate)
  utils::assignInNamespace("kmEstimate", kmEstimate, ns = "DiceKriging") 
}







alter_LogLikFun <- function(){
  
  # this file overwries the function km() in DiceKriging
  # WARNING: Only source for specific purposes

  logLikFun <-
    function(param, model, envir=NULL) {
      if (model@known.param=="Trend") {
        beta <- model@trend.coef
      } else {
        beta <- NULL
      }
      #browser()
      ## begin edit
      weights <- eval(weights, envir = globalenv())
      if(nrow(model@X) == length(weights)){
       # compute weighted x values
      model@X <- model@X * weights
      # and weighted y values
      model@y <- model@y * weights
      }
      
      #warning(list(weights))
      ## end edit
      
      if (identical(model@case, "LLconcentration_beta_sigma2")) {
        
        model@covariance <- vect2covparam(model@covariance, param)
        model@covariance@sd2 <- 1		# to get the correlation matrix
        
        aux <- covMatrix(model@covariance, model@X)
        
        R <- aux[[1]]
        T <- chol(R)
        
        x <- backsolve(t(T), model@y, upper.tri = FALSE)
        M <- backsolve(t(T), model@F, upper.tri = FALSE)
        z <- compute.z(x=x, M=M, beta=beta)
        sigma2.hat <- compute.sigma2.hat(z)
        logLik <- -0.5*(model@n * log(2*pi*sigma2.hat) + 2*sum(log(diag(T))) + model@n)
        
        if (!is.null(envir)) { 
          envir$T <- T
          envir$R <- R
          envir$z <- z
          envir$sigma2.hat <- sigma2.hat
        }
        
      } else if (identical(model@case, "LLconcentration_beta")) {
        
        nparam <- length(param)
        if (class(model@covariance) != "covAdditive0") {
          model@covariance <- vect2covparam(model@covariance, param[1:(nparam-1)])
          model@covariance@sd2 <- param[nparam]
        } else {
          model@covariance <- vect2covparam(model@covariance, param)
        }
        
        # 		if (model@covariance@nugget.estim) {
        #       coef(model@covariance, type="all") <- param
        # 		} else {
        #       coef(model@covariance, type="all-nugget") <- param
        # 		}
        # 		# pb ici : reconnaitre les variables a estimer ? 
        aux <- covMatrix(model@covariance, model@X, noise.var=model@noise.var)
        
        C <- aux[[1]]
        vn <- aux[[2]]
        
        T <- chol(C)
        x <- backsolve(t(T), model@y, upper.tri = FALSE)
        M <- backsolve(t(T), model@F, upper.tri = FALSE)
        z <- compute.z(x=x, M=M, beta=beta)
        
        logLik <-  -0.5*(model@n * log(2*pi) + 2*sum(log(diag(T))) + t(z)%*%z)     
        
        if (!is.null(envir)) {
          envir$T <- T
          envir$C <- C
          envir$vn <- vn
          envir$z <- z
        }
        
      } else if (identical(model@case, "LLconcentration_beta_v_alpha")) {
        
        nparam <- length(param)
        
        model@covariance <- vect2covparam(model@covariance, param[1:(nparam-1)])
        model@covariance@sd2 <- 1
        model@covariance@nugget <- 0
        alpha <- param[nparam]
        
        aux <- covMatrix(model@covariance, model@X)
        R0 <- aux[[1]]-diag(aux[[2]])
        R <- alpha*R0 + (1-alpha)*diag(model@n)
        
        T <- chol(R)
        x <- backsolve(t(T), model@y, upper.tri = FALSE)
        M <- backsolve(t(T), model@F, upper.tri = FALSE)
        z <- compute.z(x=x, M=M, beta=beta)
        v <- compute.sigma2.hat(z)
        
        logLik <- -0.5*(model@n * log(2*pi*v) + 2*sum(log(diag(T))) + model@n)
        
        if (!is.null(envir)) {
          envir$T <- T
          envir$R0 <- R0
          envir$v <- v
          envir$z <- z
        }	
        
      }
      
      
      if (model@method=="PMLE") {
        fun <- match.fun(model@penalty$fun)
        penalty <- - model@n * sum(fun(1/model@covariance@range.val^2, model@penalty$value))
        logLik <- logLik + penalty
      }
      
      return(logLik)
      
    }
  
  utils::assignInNamespace("logLikFun", logLikFun, ns = "DiceKriging")   
  # environment(logLikFun) <- asNamespace('DiceKriging')
  # # assures that the function will be able to call other hidden functions from the package.
  # 
  # assignInNamespace("logLikFun", logLikFun, ns = "DiceKriging") 
  # # assures that the function will be able to call other hidden functions from the package. (credits: TMS, stackoverflow)
 
}

back_to_normal_LogLikFun <- function() {
  
  # this file overwries the function km() in DiceKriging
  # WARNING: Only source for specific purposes
  
  
  logLikFun <-
    function(param, model, envir=NULL) {
      if (model@known.param=="Trend") {
        beta <- model@trend.coef
      } else {
        beta <- NULL
      }

      if (identical(model@case, "LLconcentration_beta_sigma2")) {
        
        model@covariance <- vect2covparam(model@covariance, param)
        model@covariance@sd2 <- 1		# to get the correlation matrix
        
        aux <- covMatrix(model@covariance, model@X)
        
        R <- aux[[1]]
        T <- chol(R)
        
        x <- backsolve(t(T), model@y, upper.tri = FALSE)
        M <- backsolve(t(T), model@F, upper.tri = FALSE)
        z <- compute.z(x=x, M=M, beta=beta)
        sigma2.hat <- compute.sigma2.hat(z)
        logLik <- -0.5*(model@n * log(2*pi*sigma2.hat) + 2*sum(log(diag(T))) + model@n)
        
        if (!is.null(envir)) { 
          envir$T <- T
          envir$R <- R
          envir$z <- z
          envir$sigma2.hat <- sigma2.hat
        }
        
      } else if (identical(model@case, "LLconcentration_beta")) {
        
        nparam <- length(param)
        if (class(model@covariance) != "covAdditive0") {
          model@covariance <- vect2covparam(model@covariance, param[1:(nparam-1)])
          model@covariance@sd2 <- param[nparam]
        } else {
          model@covariance <- vect2covparam(model@covariance, param)
        }
        
        # 		if (model@covariance@nugget.estim) {
        #       coef(model@covariance, type="all") <- param
        # 		} else {
        #       coef(model@covariance, type="all-nugget") <- param
        # 		}
        # 		# pb ici : reconnaitre les variables a estimer ? 
        aux <- covMatrix(model@covariance, model@X, noise.var=model@noise.var)
        
        C <- aux[[1]]
        vn <- aux[[2]]
        
        T <- chol(C)
        x <- backsolve(t(T), model@y, upper.tri = FALSE)
        M <- backsolve(t(T), model@F, upper.tri = FALSE)
        z <- compute.z(x=x, M=M, beta=beta)
        
        logLik <-  -0.5*(model@n * log(2*pi) + 2*sum(log(diag(T))) + t(z)%*%z)     
        
        if (!is.null(envir)) {
          envir$T <- T
          envir$C <- C
          envir$vn <- vn
          envir$z <- z
        }
        
      } else if (identical(model@case, "LLconcentration_beta_v_alpha")) {
        
        nparam <- length(param)
        
        model@covariance <- vect2covparam(model@covariance, param[1:(nparam-1)])
        model@covariance@sd2 <- 1
        model@covariance@nugget <- 0
        alpha <- param[nparam]
        
        aux <- covMatrix(model@covariance, model@X)
        R0 <- aux[[1]]-diag(aux[[2]])
        R <- alpha*R0 + (1-alpha)*diag(model@n)
        
        T <- chol(R)
        x <- backsolve(t(T), model@y, upper.tri = FALSE)
        M <- backsolve(t(T), model@F, upper.tri = FALSE)
        z <- compute.z(x=x, M=M, beta=beta)
        v <- compute.sigma2.hat(z)
        
        logLik <- -0.5*(model@n * log(2*pi*v) + 2*sum(log(diag(T))) + model@n)
        
        if (!is.null(envir)) {
          envir$T <- T
          envir$R0 <- R0
          envir$v <- v
          envir$z <- z
        }	
        
      }
      
      
      if (model@method=="PMLE") {
        fun <- match.fun(model@penalty$fun)
        penalty <- - model@n * sum(fun(1/model@covariance@range.val^2, model@penalty$value))
        logLik <- logLik + penalty
      }
      
      return(logLik)
      
    }
  
  utils::assignInNamespace("logLikFun", logLikFun, ns = "DiceKriging")   
  # environment(logLikFun) <- asNamespace('DiceKriging')
  # # assures that the function will be able to call other hidden functions from the package.
  # 
  # assignInNamespace("logLikFun", logLikFun, ns = "DiceKriging") 
  # # assures that the function will be able to call other hidden functions from the package. (credits: TMS, stackoverflow)
  
}