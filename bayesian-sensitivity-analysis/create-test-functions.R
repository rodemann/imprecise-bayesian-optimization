library(smoof)
## randomly create test set of functions
## WARNING: some functions require specification of "m" as argument in makeFunctionsByName(). 
## They are excluded, i.e. in case such a function is randomly selected, the sample is withdrawn.

# 2d test set
all_funs <- filterFunctionsByTags("single-objective")
random_ids <- sample(1:length(all_funs),10)
test_functions <-  makeFunctionsByName(all_funs[random_ids], dimensions = 2L)

# 3d test set
all_funs <- filterFunctionsByTags("single-objective")
random_ids <- sample(1:length(all_funs),10)
test_functions <-   c(test_functions, makeFunctionsByName(all_funs[random_ids], dimensions = 3L))

# 4d test set
all_funs <- filterFunctionsByTags("single-objective")
random_ids <- sample(1:length(all_funs),10)
test_functions <- c(test_functions, makeFunctionsByName(all_funs[random_ids], dimensions = 4L))

# 7d test set
all_funs <- filterFunctionsByTags("single-objective")
random_ids <- sample(1:length(all_funs),10)
test_functions <- c(test_functions, makeFunctionsByName(all_funs[random_ids], dimensions = 7L))

# ## there is some weird bug in michalewiz function, hence:
# for (i in 1:length(test_functions)) {
#   if(getName(test_functions[[i]]) == "1-d Michalewicz Function (m = 10)")
#     print("Warning! michalewicz fun included")
#     
# }

obj.fun.list <- test_functions
save(obj.fun.list, file = "bayesian-sensitivity-analysis/test-functions" )
