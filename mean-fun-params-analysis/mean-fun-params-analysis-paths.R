# This is a sensitivity analysis for Bayesian Optimization with Gaussian Processes 
# (Kriging) as surrogate model with regard to different MEAN functions
#
# It uses various target functions from package smoof


library(mlrMBO)
library(mlr3)
library(smoof)
library(rgenoud)
library(DiceKriging)
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
ground_truths = lapply(parameters_sets, generateDesign, n = 4000L, fun = lhs::maximinLHS) 

# set seqeunce across functions
sequence_functions = seq_along(obj.fun.list)

# set number of BO iterations
iters = 20L

#set number of BO runs per function to be optimized 
BO_runs_per_fun = 40L

# pick kernel function
covtype = "gauss"

# store global results here
mbo_results_paths = list()

# mean functions
#mean_constants = c(0, 1, 5, 10, 100)

for (i in seq.int(1,5)) {
  
  # store local (t) results here
  opt.paths = vector("list", length(obj.fun.list))

  for(t in sequence_functions){
    # pick objective function
    obj.fun = obj.fun.list[[t]]
    # pick respective design; recall that this is crucial for
    # our analysis to be conditioned on the intialization 
    initial.design = designs[[t]]
    # get true y 
    ground_truth = ground_truths[[t]]
    y_true = apply(ground_truth, 1, obj.fun)
    mean_true = mean(y_true)
    # mean functions
    mean_constants = c(mean_true - abs(3*mean_true), mean_true - abs(0.2*mean_true), mean_true, mean_true + abs(0.2*mean_true), mean_true + abs(3*mean_true))
    
    
    # set MBO 
    ctrl = makeMBOControl(propose.points = 1L)
    ctrl = setMBOControlTermination(ctrl, iters = iters)
    ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI(),
                               opt = "focussearch")
    
    
    ## choose surrogate model including mean function
    mean_fun_constant = mean_constants[i]
    
    lrn = makeLearner("regr.km", formula = ~1, 
                      coef.trend = mean_fun_constant, 
                      predict.type = "se", 
                      covtype = covtype, 
                      optim.method = "gen", 
                      control = list(trace = FALSE), 
                      config = list(on.par.without.desc = "warn"))
    
    # ensure numerical stability in km {DiceKriging} cf. github issue and recommendation by Bernd Bischl 
    y = apply(initial.design, 1, obj.fun)
    Nuggets = 1e-8*var(y)
    lrn = setHyperPars(learner = lrn, nugget=Nuggets)
    
    
    ## just test something:
    # base_model = km(response = y, design = initial.design,
    #                 formula = ~1, #coef.trend = mean_fun_constant, 
    #                 covtype = covtype, 
    #                 optim.method = "gen", 
    #                 control = list(trace = FALSE),
    #                 nugget = Nuggets)
    # 
    
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
save(mbo_results_paths, file = "BO_results_40_20_meanfun")

#load(file = paste(getwd(),"BO_results"))





