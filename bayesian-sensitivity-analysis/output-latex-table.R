library(xtable)
library(smoof)
library(dplyr)
## combine results from all analyses and output a latex ready table

load("bayesian-sensitivity-analysis/acc-diff-mean-fun-powexp-kernel")
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

sum(BO_paths_acc_difference_mean_fun)
sum(BO_paths_acc_difference_mean_par)
sum(BO_paths_acc_difference_kernel_fun)
sum(BO_paths_acc_difference_kernel_par)

## standardization
acc_diff_df_plain <-acc_diff_df[,-1]
mean_acc_diff <- apply(acc_diff_df_plain, 1, mean)
acc_diff_df_standardized <- acc_diff_df %>% mutate("Mean Functional Form stand." = Mean.Functional.Form/mean_acc_diff ) %>% mutate("Mean Parameters stand." = Mean.Parameters/mean_acc_diff ) %>% mutate("Kernel Functional Form stand." = Kernel.Functional.Form/mean_acc_diff ) %>% mutate("Kernel Parameters stand." = Kernel.Parameters/mean_acc_diff )
acc_diff_df_standardized_named <- acc_diff_df_standardized[,-c(2,3,4,5)]
acc_diff_df_standardized_plain  <- acc_diff_df_standardized[,-c(1,2,3,4,5)]
#sums: 
apply(acc_diff_df_standardized_plain, 2, mean)

tex_table <- xtable(acc_diff_df_standardized_named, display = c("g","g","g","g","g","g") )
print(tex_table, include.rownames=FALSE)



