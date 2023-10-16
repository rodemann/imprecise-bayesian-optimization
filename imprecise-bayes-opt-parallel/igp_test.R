library(mlrMBO)
source("imprecise-bayes-opt-parallel/igp_upper_mean.R")
source("imprecise-bayes-opt-parallel/igp_lower_mean.R")
#source("data/make-kapton-rf.R")
## playground: 
fun = smoof::makeAlpine01Function(1)
#fun = kapton_fun_rf
set.seed(152020)
x = runif(39,-10,8) %>% sort() %>% as.vector() 
#x = rep(1,40) %>% sort() %>% as.vector() 
#x = data.frame(cl_bo_res$opt.path)$x

y = sapply(x, FUN = fun) %>% as.vector()
data = data.frame(x=x, y=y)
.task = makeRegrTask(data = data, target = "y")

base_kernel = "powexp"
lrn = makeLearner("igp_upper_mean", base_kernel = base_kernel, imprecision = 1, predict.type = "se")
mod_upper = train(lrn, .task)
lrn = makeLearner("igp_lower_mean", base_kernel = base_kernel, imprecision = 1, predict.type = "se")
mod_lower = train(lrn, .task)

#test_data = subset(data, select = x)
#test_data = data.frame(x = rep(1,50), row.names = NULL)

# test_data = data.frame(x = c(-3,0,3))
#predict(mod, newdata = data.frame(4))
#predict(mod, newdata = test_data)

grid = seq(-5,5,0.1)
test_data = data.frame(x = grid, row.names = NULL)

##viz results
#lower
preds = predict(mod_lower, newdata = test_data)
preds_raw = preds$data$response
par(bg = 'lightgrey')
plot(x = grid , y = preds_raw, type = "l", col = "steelblue", ylim = c(-10,10))
# 95 % CIs
preds_up <- preds$data$response + qnorm(0.975) * preds$data$se
lines(x = grid, y = preds_up, lty = 3, col = "steelblue")
preds_low <- preds$data$response - qnorm(0.975) * preds$data$se
lines(x = grid, y = preds_low, lty = 3, col = "steelblue")

#upper
preds = predict(mod_upper, newdata = test_data)
preds_raw = preds$data$response
lines(x = grid , y = preds_raw, type = "l", col = "red")
# 95 % CIs
preds_up <- preds$data$response + qnorm(0.975) * preds$data$se
lines(x = grid, y = preds_up, lty = 3, col = "red")
preds_low <- preds$data$response - qnorm(0.975) * preds$data$se
lines(x = grid, y = preds_low, lty = 3, col = "red")



#add training data points
truth = sapply(x, fun) %>% as.vector()
points(x = x, y = truth)
# compare with precise GP model
data = getTaskData(.task, target.extra = TRUE)
design = data$data
response = data$target
precise_GP_fit <- DiceKriging::km(response = response, design = design, covtype = base_kernel,
                                  estim.method = "MLE", optim.method = "gen", 
                                  nugget = 1e-8*var(response), control = list(trace = FALSE))
## add precise GP model to plot
test_points = grid
precise_GP <- predict(precise_GP_fit, newdata = test_points, type = "SK", se.compute = TRUE)
pred_grid_mean = cbind(test_points, precise_GP$mean)
lines(pred_grid_mean)
# CIs from classic GP
pred_grid_lower = cbind(test_points, precise_GP$lower95)
pred_grid_upper = cbind(test_points, precise_GP$upper95)
lines(pred_grid_upper, lty = 3)
lines(pred_grid_lower, , lty = 3)


grid(col = "gray", lty = "dashed")




# 
# ###
# ## BO with igp:
# fun <- smoof::makeAckleyFunction(1)
# parameter_set <- getParamSet(fun)
# design <- generateDesign(n = 10L, par.set = parameter_set, fun = lhs::randomLHS)
# iters = 4L
# # set Control Argument of BO 
# ctrl = makeMBOControl(propose.points = 1L, store.model.at = 1:iters)
# ctrl = setMBOControlTermination(ctrl, iters = iters)
# ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritMeanResponse(), opt = "focussearch",
#                            opt.focussearch.points = 200, opt.focussearch.maxit = 1)
# lrn <- makeLearner("igp_upper_mean", base_kernel = "powexp", imprecision = 100)
# lrn <- makeLearner("igp_lower_mean", base_kernel = "powexp", imprecision = 100)
# 
# res_mbo <- mbo(fun = fun, design = design, control = ctrl, learner = lrn)
# as.data.frame(res_mbo$opt.path)
# res_mbo$models$`5`$learner.model
# 
# #lineprof(mbo(fun = fun, design = design, control = ctrl, learner = lrn))
# 
# 
# 
# # classic BO:
# ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())
# lrn_classic = makeLearner("regr.km", covtype = "powexp", predict.type = "se", optim.method = "gen", 
#                           control = list(trace = FALSE), config = list(on.par.without.desc = "warn"))
# # ensure numerical stability in km {DiceKriging} cf. github issue and recommendation by Bernd Bischl 
# y = apply(design, 1, fun)
# Nuggets = 1e-8*var(y)
# lrn_classic = setHyperPars(learner = lrn_classic, nugget=Nuggets)
# 
# res_mbo_classic = mbo(fun = fun, design = design, control = ctrl, learner = lrn_classic)
# 
# res_mbo_classic$models
# 
# 
# 
# 
# 





# true_grid = cbind(design, response)
# plot(pred_grid, type = "l")
# points(true_grid)


# 
# as.vector(preds)
# 
# plot(preds)
# select(unlist(preds))
# 
# predict(mod, newdata = data.frame(1))
# ackley.fun(1)
# 
# task.pred
# task.pred = predict(mod, newdata = data.frame(x = c(0,1,2,3)))
# task.pred
# 
# predict(mod, task = .task)$data

## TODO
# vectorize predict function, can only handle doubles as of now 
# predictions do not make any sense! go through it!
