source("imprecise-bayes-opt-plug-in/makeMBOInfillCritGLCB.R")
source("imprecise-bayes-opt-plug-in/initCrit.InfillCritGLCB.R")
source("read-data-kapton.R")
library(dplyr)
library(mlrMBO)

kapton_data <- data
kapton_data <- select(kapton_data, c(ratio, time))

ground_truth = train(makeLearner("regr.randomForest"), 
                     makeRegrTask(data = kapton_data, target = "ratio"))

fun = function(x) {
  df = as.data.frame(x)
  # df$gas = factor(df$gas, levels = levels(data$gas))
  return(getPredictionResponse(predict(ground_truth, newdata = df)))
}

parameter_set = makeParamSet(
  #makeIntegerParam("power", lower = 10, upper = 5555)#,
  makeNumericParam("time", lower = 500, upper = 21000)#,
  #makeDiscreteParam("gas", values = c("Nitrogen", "Air", "Argon")),
  #makeIntegerParam("pressure", lower = 0, upper = 1000)
)

obj_fun = makeSingleObjectiveFunction(
  name = "Kapton",
  fn = fun,
  par.set = parameter_set,
  has.simple.signature = FALSE,
  minimize = FALSE
)
# visualize ground truth
grid <- data.frame("time" = seq(0, 21000, 1))
plot(fun(grid), type = "l")

# unbiased sample
design <- generateDesign(n = 40L, par.set = parameter_set, fun = lhs::randomLHS)
# biased sample
design_biased <- design %>% filter(time < 10000)

# set up BO 
iters = 15L
# set Control Argument of BO 
ctrl = makeMBOControl(propose.points = 1L)
ctrl = setMBOControlTermination(ctrl, iters = iters)
# plug in GLCB
infill_crit = makeMBOInfillCritGLCB(cb.lambda = 1, cb.rho = 100)
infill_crit = makeMBOInfillCritCB(cb.lambda = 1)
ctrl = setMBOControlInfill(ctrl, crit = infill_crit, opt = "focussearch")#,
#opt.focussearch.points = 200, opt.focussearch.maxit = 1)

lrn = makeLearner("regr.km", covtype = "powexp", predict.type = "se", optim.method = "gen", 
                  control = list(trace = FALSE), config = list(on.par.without.desc = "warn"))
# ensure numerical stability in km {DiceKriging} cf. github issue and recommendation by Bernd Bischl 
y = obj_fun(design_biased)
Nuggets = 1e-8*var(y)
lrn = setHyperPars(learner = lrn, nugget=Nuggets)
initial_opt = max(y)


res_mbo = mbo(fun = obj_fun, design = design_biased, control = ctrl)

as.data.frame(res_mbo$opt.path)
c( res_mbo$y, res_mbo$best.ind)
initial_opt







