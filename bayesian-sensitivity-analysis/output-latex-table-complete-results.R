library(xtable)

## combine results from all analyses and output a latex ready table

load("bayesian-sensitivity-analysis/univariate-acc-diff-mean-fun")
load("bayesian-sensitivity-analysis/univariate-acc-diff-mean-par")
load("bayesian-sensitivity-analysis/univariate-acc-diff-kernel-fun")
load("bayesian-sensitivity-analysis/univariate-acc-diff-kernel-par")
source("univariate-benchmark-functions/univariate-test-functions-smoof.R")
source("univariate-benchmark-functions/univariate-test-functions-soobench.R")
source("univariate-benchmark-functions/get-synthetic-benchmark-funs-for-minimization.R")
obj.fun.list <- all_benchmark_funs[1:10]


fun_names <- lapply(obj.fun.list, getName)
unlist(fun_names)
acc_diff_df_1d <- data.frame( "Test Function" = unlist(fun_names),
                              "Mean Functional Form" = BO_paths_acc_difference_mean_fun,
                              "Mean Parameters" = BO_paths_acc_difference_mean_par,
                              "Kernel Functional Form" = BO_paths_acc_difference_kernel_fun,
                              "Kernel Parameters" = BO_paths_acc_difference_kernel_par)

tex_table <- xtable(acc_diff_df_1d, display = c("g","g","g","g","g","g") )
print(tex_table, include.rownames=FALSE)

load("bayesian-sensitivity-analysis/acc-diff-mean-fun")
load("bayesian-sensitivity-analysis/acc-diff-mean-par")
load("bayesian-sensitivity-analysis/acc-diff-kernel-fun")
load("bayesian-sensitivity-analysis/acc-diff-kernel-par")
load("bayesian-sensitivity-analysis/test-functions")


fun_names <- lapply(obj.fun.list, getName)
unlist(fun_names)
acc_diff_df <- data.frame( "Test Function" = unlist(fun_names),
                           "Mean Functional Form" = BO_paths_acc_difference_mean_fun,
                           "Mean Parameters" = BO_paths_acc_difference_mean_par,
                           "Kernel Functional Form" = BO_paths_acc_difference_kernel_fun,
                           "Kernel Parameters" = BO_paths_acc_difference_kernel_par)
tex_table <- xtable(acc_diff_df, display = c("g","g","g","g","g","g") )
print(tex_table, include.rownames=FALSE)

acc_diff_df <- rbind(acc_diff_df_1d, acc_diff_df)

tex_table <- xtable(acc_diff_df, display = c("g","g","g","g","g","g") )
print(tex_table, include.rownames=FALSE)

sum(acc_diff_df$Mean.Functional.Form)
sum(acc_diff_df$Mean.Parameters)
sum(acc_diff_df$Kernel.Functional.Form)
sum(acc_diff_df$Kernel.Parameters)


## standardization
acc_diff_df_plain <-acc_diff_df[,-1]
mean_acc_diff <- apply(acc_diff_df_plain, 1, mean)
acc_diff_df_standardized <- acc_diff_df %>% mutate("Mean Functional Form stand." = Mean.Functional.Form/mean_acc_diff ) %>% mutate("Mean Parameters stand." = Mean.Parameters/mean_acc_diff ) %>% mutate("Kernel Functional Form stand." = Kernel.Functional.Form/mean_acc_diff ) %>% mutate("Kernel Parameters stand." = Kernel.Parameters/mean_acc_diff )
acc_diff_df_standardized_named <- acc_diff_df_standardized[,-c(2,3,4,5)]
acc_diff_df_standardized_plain  <- acc_diff_df_standardized[,-c(1,2,3,4,5)]
#sums: 
apply(acc_diff_df_standardized_plain, 2, sum)

tex_table <- xtable(acc_diff_df_standardized_named, display = c("g","g","g","g","g","g") )
print(tex_table, include.rownames=FALSE)




