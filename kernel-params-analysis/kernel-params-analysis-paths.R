# This is a sensitivity analysis for Bayesian Optimization with Gaussian Processes 
# (Kriging) as surrogate model with regard to different MEAN functions
#
# It uses various target functions from package smoof


library(mlrMBO)
library(mlr3)
library(smoof)
library(rgenoud)
library(DiceKriging)
library(foreach)
# Functions to be optimized: 
obj.fun.list = list(makeBraninFunction(),
                    makeBBOBFunction(dimension = 3, fid = 2, iid = 3),
                    makeBBOBFunction(dimension = 2, fid = 3, iid = 6),
                    makeBBOBFunction(dimension = 6, fid = 20, iid = 2),
                    makeBBOBFunction(dimension = 7, fid = 17, iid = 4),
                    makeBBOBFunction(dimension = 3, fid = 5, iid = 16),
                    makeAckleyFunction(1),
                    makeAckleyFunction(2),
                    makeAckleyFunction(5),
                    makeBirdFunction(),
                    makeEngvallFunction(),
                    makeGiuntaFunction(),
                    makeZettlFunction(),
                    makePriceN4Function(),
                    makeCarromTableFunction(),
                    makeSchwefelFunction(dimensions = 3),
                    makeFreudensteinRothFunction(),
                    makeShekelFunction(m = 5)
                    
)

# set initial design for each objective function
# so that the analysis is conditioned on one set of initial values
parameters_sets = lapply(obj.fun.list, getParamSet)
designs = lapply(parameters_sets, generateDesign, n = 10L, fun = lhs::maximinLHS) 

# ground truth to get true mean: 
ground_truths = lapply(parameters_sets, generateDesign, n = 400L, fun = lhs::maximinLHS) 

# set seqeunce across functions
sequence_functions = seq_along(obj.fun.list)

#length-scale parameters/kernel bandwith/theta 
#kernel_bandwith_vec = c(1, 5, 10, 100, 500) 

# store global results here
mbo_results_paths = list()

# set number of configurations
configs = 5


## create vector of kernel bandwidth parameters per function
kernel_bandwith_per_fun = list()
for (fun in sequence_functions) {
  # pick objective function
  obj.fun = obj.fun.list[[fun]]  
  ## choose kernel bandwidth parameter based on ground truth for given function
  ground_truth = ground_truths[[fun]]
  y_true = apply(ground_truth, 1, obj.fun)
  # numerical stability:
  nugget_true = 1e-8*var(y_true)
  baseline_km <- km(design = ground_truth, response = y_true, 
                    nugget = nugget_true, covtype = "matern5_2",  control = list(trace = FALSE) )
  kernel_bandwith_true = baseline_km@covariance@range.val  
  kernel_bandwith_vec <- foreach::foreach(i = 1:length(kernel_bandwith_true)) %do% (c(kernel_bandwith_true[i] - abs(3*kernel_bandwith_true[i]), kernel_bandwith_true[i] - abs(0.2*kernel_bandwith_true[i]), kernel_bandwith_true[i], kernel_bandwith_true[i] + abs(0.2*kernel_bandwith_true[i]), kernel_bandwith_true[i] + abs(3*kernel_bandwith_true[i])))
  
  kernel_bandwith_per_fun[[fun]] = kernel_bandwith_vec 
}


# set number of BO iterations
iters = 20L

#set number of BO runs per function to be optimized 
BO_runs_per_fun = 40L

# pick kernel function
covtype = "matern5_2"

for (i in 1:configs) {
  
  # store local (t) results here
  opt.paths = vector("list", length(obj.fun.list))
  
  for(t in sequence_functions){
    # pick objective function
    obj.fun = obj.fun.list[[t]]
    # pick respective design; recall that this is crucial for
    # our analysis to be conditioned on the intialization 
    initial.design = designs[[t]]
    
    
    # set MBO 
    ctrl = makeMBOControl(propose.points = 1L)
    ctrl = setMBOControlTermination(ctrl, iters = iters)
    ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI(),
                               opt = "focussearch")
    
    
    #choose kernel parameter
    kernel_bandwith_vec = kernel_bandwith_per_fun[[t]]
    kernel_bandwith = lapply(kernel_bandwith_vec, "[", i)
    kernel_bandwith = as.vector(unlist(kernel_bandwith)) # Admittedly, this is ugly. Yet, it works.
    lrn = makeLearner("regr.km", 
                      predict.type = "se", 
                      covtype = covtype,
                      coef.cov = kernel_bandwith,
                      optim.method = "gen", 
                      control = list(trace = FALSE), 
                      config = list(on.par.without.desc = "warn"))
    
    
    # ensure numerical stability in km {DiceKriging} cf. github issue and recommendation by Bernd Bischl 
    y = apply(initial.design, 1, obj.fun)
    Nuggets = 1e-8*var(y)
    lrn = setHyperPars(learner = lrn, nugget=Nuggets)
    
    
    # store results here  
    mbo_runs = list()
    
    for (j in 1:BO_runs_per_fun) {
      mbo_runs[[j]] = mbo(obj.fun, design = initial.design, control = ctrl, 
                          learner = lrn, show.info = FALSE)
    }
    
    # retrieve optimization path
    opt.path = lapply(mbo_runs, `[[`,"opt.path")
    opt.path.y = lapply(opt.path, getOptPathY)
    
    # retrieve current best target value
    get_current_best <- function(proposals){
      for (o in 2:length(proposals)) {
        if(proposals[o] > proposals[o - 1]) 
          proposals[o] = proposals[o - 1] 
      }
      proposals  
    }
    opt.path.best.y = lapply(opt.path.y, get_current_best)
    
    opt.paths[[t]] = opt.path.best.y
    
  }
  mbo_results_paths[[i]] = opt.paths
}  



mbo_results_paths
# svae/load results list as R object    
save(mbo_results_paths, file = "kernel-params-analysis/BO_results_40_20_kernel_bandwith_mat52_kernel")

#load(file = paste(getwd(),"BO_results"))





