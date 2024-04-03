source("imprecise-bayes-opt-plug-in/makeMBOInfillCritGLCB.R")
source("imprecise-bayes-opt-plug-in/initCrit.InfillCritGLCB.R")
source("data/read-data-kapton.R")
library(dplyr)
library(mlrMBO)
library(ggplot2)
kapton_data <- select(data, -ratio)

# Convert 'gas' to a factor
kapton_data$gas <- as.factor(kapton_data$gas)

# Create a model matrix to get dummy variables for 'gas'
kapton_data_numeric <- model.matrix(~ . + 0, data = kapton_data)

# Perform PCA on the dataset including the dummy variables for 'gas'
pca_result_with_gas <- prcomp(kapton_data_numeric, center = TRUE, scale. = TRUE)

# Print summary of PCA results
summary(pca_result_with_gas)

# Print the rotation (loadings) of the principal components
#print(pca_result_with_gas$rotation)

# Extracting the scores (principal component scores) for the first component
scores_first_component <- pca_result_with_gas$x[,1]


kapton_data_embed = cbind(data$ratio,scores_first_component) %>% as.data.frame()
colnames(kapton_data_embed) <- c("y","x")

ground_truth = train(makeLearner("regr.randomForest"), 
                     makeRegrTask(data = kapton_data_embed, target = "y"))

fun = function(x) {
  df = data.frame("x" = x)
  # df$gas = factor(df$gas, levels = levels(data$gas))
  return(getPredictionResponse(predict(ground_truth, newdata = df)))
}

parameter_set = makeParamSet(
  #makeIntegerParam("power", lower = 10, upper = 5555)#,
  makeNumericParam("x", lower = 500, upper = 21000)#,
  #makeDiscreteParam("gas", values = c("Nitrogen", "Air", "Argon")),
  #makeIntegerParam("pressure", lower = 0, upper = 1000)
)

kapton_fun_rf = makeSingleObjectiveFunction(
  name = "Kapton",
  fn = fun,
  par.set = parameter_set,
  has.simple.signature = FALSE,
  minimize = FALSE
)
# visualize ground truth
grid <- data.frame("x" = seq(min(scores_first_component), max(scores_first_component), 0.005))
# # langsam: plot(apply(grid, 1, fun), type = "l")
# #schneller:
plot(fun(grid), type = "l")

time_df = cbind(grid, fun(grid))
names(time_df) = c("embedding", "graphene quality")
graph_time_plot =  qplot(data = time_df, x=time_df$"embedding", y=time_df$"graphene quality", 
                         geom = "line", xlab = "embedding", ylab = "graphene quality",
                         main = "Graphene-embedding target function (Random Forest)" )
graph_time_plot
#combine with power plot
#ggarrange(graph_time_plot, graph_power_plot, ncol = 2, nrow = 1)




# # set up BO 
# iters = 15L
# # set Control Argument of BO 
# ctrl = makeMBOControl(propose.points = 1L)
# ctrl = setMBOControlTermination(ctrl, iters = iters)
# # plug in GLCB
# infill_crit = makeMBOInfillCritGLCB(cb.lambda = 1, cb.rho = 100)
# infill_crit = makeMBOInfillCritCB(cb.lambda = 1)
# ctrl = setMBOControlInfill(ctrl, crit = infill_crit, opt = "focussearch")#,
# #opt.focussearch.points = 200, opt.focussearch.maxit = 1)
# 
# lrn = makeLearner("regr.km", covtype = "powexp", predict.type = "se", optim.method = "gen", 
#                   control = list(trace = FALSE), config = list(on.par.without.desc = "warn"))
# # ensure numerical stability in km {DiceKriging} cf. github issue and recommendation by Bernd Bischl 
# y = obj_fun(design_biased)
# Nuggets = 1e-8*var(y)
# lrn = setHyperPars(learner = lrn, nugget=Nuggets)
# initial_opt = max(y)
# 
# 
# res_mbo = mbo(fun = obj_fun, design = design_biased, control = ctrl)
# 
# as.data.frame(res_mbo$opt.path)
# c( res_mbo$y, res_mbo$best.ind)
# initial_opt







