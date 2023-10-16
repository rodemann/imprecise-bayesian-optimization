library(testthat)
library(mlrMBO)

#gobal setting
iters = 3 

# 1. Test on numerical parameter set (material design)

# read material design data
source("_Explore_Exploit_Measures/analyze-material-design-xplxpl/read-data-material-design-jr.R")


# set hyperparameter
model = train(makeLearner("regr.randomForest"), makeRegrTask(data = data, target = "target"))
fun = function(x) {
  df = as.data.frame(x)
  return(getPredictionResponse(predict(model, newdata = df)))
}
ps = makeParamSet(
  makeNumericParam("f", lower = 0, upper = 0.25),
  makeNumericParam("k", lower = 0, upper = 0.1),
  makeNumericParam("du", lower = 0, upper = 2e-5),
  makeNumericParam("dv", lower = 0, upper = 2e-5),
  makeIntegerParam("x", lower = 0, upper = 200),
  makeIntegerParam("y", lower = 0, upper = 200)
)
objfun = makeSingleObjectiveFunction(
  name = "Synthesis",
  fn = fun,
  par.set = ps,
  has.simple.signature = FALSE,
  minimize = FALSE
)
# sample 10 points for the initial surrogate model
initial.data = generateDesign(10, attributes(objfun)$par.set) 

# set global control settings
ctrl = makeMBOControl(y.name = "target")
ctrl = setMBOControlTermination(ctrl, iters = iters)


test_that("basically works for infill.opt = focussearch with various control 
          settings (all except DIB, which is for multi-objective only)", {
  ctrl = setMBOControlInfill(ctrl, opt = "focussearch", crit = makeMBOInfillCritAdaCB())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "focussearch", crit = makeMBOInfillCritAEI())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "focussearch", crit = makeMBOInfillCritCB())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "focussearch", crit = makeMBOInfillCritEI())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "focussearch", crit = makeMBOInfillCritEQI())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "focussearch", crit = makeMBOInfillCritMeanResponse())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "focussearch", crit = makeMBOInfillCritStandardError())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  }
)
  




test_that("basically works for infill.opt = ea with various control settings", {
  ctrl = setMBOControlInfill(ctrl, opt = "ea", crit = makeMBOInfillCritAdaCB())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "ea", crit = makeMBOInfillCritAEI())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "ea", crit = makeMBOInfillCritCB())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "ea", crit = makeMBOInfillCritEI())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "ea", crit = makeMBOInfillCritEQI())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "ea", crit = makeMBOInfillCritMeanResponse())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "ea", crit = makeMBOInfillCritStandardError())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
}
)


# cmaes only allows numeric data, no integers

# kill integer parameters for cmaes:
data = data[,-(5:6)]

# initialize surrogate surrogate model
model = train(makeLearner("regr.randomForest"), makeRegrTask(data = data, target = "target"))
ps = makeParamSet(
  makeNumericParam("f", lower = 0, upper = 0.25),
  makeNumericParam("k", lower = 0, upper = 0.1),
  makeNumericParam("du", lower = 0, upper = 2e-5),
  makeNumericParam("dv", lower = 0, upper = 2e-5)
  # kill integer parameters for cmaes:
  #makeIntegerParam("x", lower = 0, upper = 200),
  #makeIntegerParam("y", lower = 0, upper = 200)
)

fun = function(x) {
  df = as.data.frame(x)
  return(getPredictionResponse(predict(model, newdata = df)))
}

objfun = makeSingleObjectiveFunction(
  name = "Synthesis",
  fn = fun,
  par.set = ps,
  has.simple.signature = FALSE,
  minimize = FALSE
)

infill.opt = "cmaes"

initial.data = generateDesign(20, attributes(objfun)$par.set) #from par space


test_that("basically works for infill.opt = cmaes with various control settings", {
  ctrl = setMBOControlInfill(ctrl, opt = "cmaes", crit = makeMBOInfillCritAdaCB())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "cmaes", crit = makeMBOInfillCritAEI())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "cmaes", crit = makeMBOInfillCritCB())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "cmaes", crit = makeMBOInfillCritEI())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "cmaes", crit = makeMBOInfillCritEQI())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "cmaes", crit = makeMBOInfillCritMeanResponse())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
  ctrl = setMBOControlInfill(ctrl, opt = "cmaes", crit = makeMBOInfillCritStandardError())
  expect_data_frame(xplxplMBO(objfun, design = initial.data, control = ctrl, show.info = FALSE),
                    min.rows = 1, max.rows = iters,
                    ncols = 7)
}
)



# 2. Test on mixed parameter space (including categorical pars)

source("_Explore_Exploit_Measures/analyze-kapton-xplxpl/read-data-kapton-jr.R")


model = train(makeLearner("regr.randomForest"), makeRegrTask(data = data, target = "ratio"))

# Set hyperparameters

fun = function(x) {
  df = as.data.frame(x)
  df$gas = factor(df$gas, levels = levels(data$gas))
  return(getPredictionResponse(predict(model, newdata = df)))
}

ps = makeParamSet(
  makeIntegerParam("power", lower = 10, upper = 5555),
  makeIntegerParam("time", lower = 500, upper = 20210),
  makeDiscreteParam("gas", values = c("Nitrogen", "Air", "Argon")),
  makeIntegerParam("pressure", lower = 0, upper = 1000)
)

objfun = makeSingleObjectiveFunction(
  name = "Kapton",
  fn = fun,
  par.set = ps,
  has.simple.signature = FALSE,
  minimize = FALSE
)
'samples.argon = sample(rownames(data[data$gas == "Argon", ]), 3)'
# sample 9 points for the initial surrogate model, stratified across gases
samples.argon = sample(rownames(data[data$gas == "Argon", ]), 3)
samples.nitro = sample(rownames(data[data$gas == "Nitrogen", ]), 3)
samples.air = sample(rownames(data[data$gas == "Air", ]), 3)
initial.data = data[c(samples.argon, samples.nitro, samples.air), ]
cat(paste("Best training ratio: ", max(initial.data$ratio), "\n", sep = ""))

ctrl = makeMBOControl(y.name = "ratio")
ctrl = setMBOControlInfill(ctrl, opt = "focussearch", opt.focussearch.maxit = 20, opt.focussearch.points = 5, crit = makeMBOInfillCritEI())
ctrl = setMBOControlTermination(ctrl, iters = iters)




