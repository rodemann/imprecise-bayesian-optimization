#source("univariate-benchmark-functions/get-synthetic-benchmark-funs-for-minimization.R")
library(ggplot2)
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-funs9-12-16-20")
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-funs9-12-16-20")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-funs9-12-16-20")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-funs9-12-16-20")
# load(file = "weighted-ML-GP/results/weighted-ML-results-global-n-100-iters-20-ninit-4-funs9-12-16-20")

# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-funs26-27-29-34")
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-funs26-27-29-34")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-funs26-27-29-34")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-funs26-27-29-34")
# load(file = "weighted-ML-GP/results/weighted-ML-results-global-n-100-iters-20-ninit-4-funs26-27-29-34")

# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-funs19-1-3-4-6-8")
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-funs19-1-3-4-6-8")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-funs19-1-3-4-6-8")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-funs19-1-3-4-6-8")
#load(file = "weighted-ML-GP/results/weighted-ML-results-global-n-100-iters-20-ninit-4-funs26-27-29-34")

# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat")
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat")
# load(file = "weighted-ML-GP/results/weighted-ML-results-global-n-100-iters-20-ninit-4-mexican-hat")
# 
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs")
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs")
# load(file = "weighted-ML-GP/results/weighted-ML-results-global-n-100-iters-20-ninit-4-mexican-hat-rs")
# 
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-0.5")
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-0.5")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-0.5")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-0.5")
# 
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-1")
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-1")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-1")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-1")
# 
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-10")
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-10")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-10")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-10")
# 
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-funs1-3-4-6-8-lcb-10")
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-funs1-3-4-6-8-lcb-10")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-funs1-3-4-6-8-lcb-10")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-funs1-3-4-6-8-lcb-10")

# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-4")
# load(file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-4")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-4")
# load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-mexican-hat-rs-lcb-4")

load(file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-funs-1-3-4-6-8-27-31-sampling-based")
load(file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-funs-1-3-4-6-8-27-31-sampling-based")
load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-funs-1-3-4-6-8-27-31-sampling-based")
load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-funs-1-3-4-6-8-27-31-sampling-based")
# note: fun 1 is not included


load(file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-20-ninit-4-mhat-lcb-sampling-based")
load(file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-20-ninit-4-funs-mhat-lcb-sampling-based")
load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-20-ninit-4-funs-mhat-lcb-sampling-based")
load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-20-ninit-4-funs-mhat-lcb-sampling-based")

load(file = "weighted-ML-GP/results/weighted-ML-results-mop-n-100-iters-25-ninit-4-alpine-2-lcb-sampling-based-low-vals")
load(file = "weighted-ML-GP/results/weighted-ML-results-mop-sd-n-100-iters-25-ninit-4-funs-alpine-2-lcb-sampling-based-low-vals")
load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-n-100-iters-25-ninit-4-funs-alpine-2-lcb-sampling-based-low-vals")
load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-n-100-iters-25-ninit-4-funs-alpine-2-lcb-sampling-based-low-vals")


#load(file = "weighted-ML-GP/results/weighted-ML-results-mop-funs21-30-n-60-iters-25-ninit-5")
#load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-funs21-30-n-60-iters-25-ninit-5")

#load(file = "weighted-ML-GP/results/unweighted-ML-results-mop-sd-funs1-10-n-100-iters-25-ninit-5")

#funs_to_test <- all_benchmark_funs[1:10]
# 
# for (i in 1:4) {
#   plot(w_paths_mop[[i]], type = "p")
#   points(uw_paths_mop[[i]], col = "red")
# }

initial_design_size = 4
fun_nr = 1
fun = fun_nr # same same
n = 200

# #BO_paths_global_min <- lapply(BO_paths, min)
# BO_paths_matrix <- matrix(unlist(global_res_unweightedML[[fun_nr]]), ncol = n) %>% as.data.frame()
# #remove initial design
# BO_paths_initial <- BO_paths_matrix %>% slice_head(n = initial_design_size)
# BO_paths_optim <- BO_paths_matrix %>% slice_tail(n = nrow(BO_paths_matrix) - initial_design_size)

# check for normality 
#apply(BO_paths_optim, 1, hist, breaks = 500)
#plot(1:100, BO_paths_optim[5,])
# --> aprrox. normal (only few outlier in some iterations)

# compute CIs in case of normal distribution of mops (X approx. normal, \sigma^2 unknown -> t-quantiles) 
alpha = 0.05
compute_CI_upper <- function(mean, sd){
  mean + qt(p = 1 - alpha/2, df = n - 1) * sd/sqrt(n)
}
CI_upper_bound_w <- mapply(compute_CI_upper, w_paths_mop, w_paths_mop_sd)
CI_upper_bound_uw <- mapply(compute_CI_upper, uw_paths_mop, uw_paths_mop_sd)

compute_CI_lower <- function(mean, sd){
  mean - qt(p = 1 - alpha/2, df = n - 1) * sd/sqrt(n)
}
CI_lower_bound_w <- mapply(compute_CI_lower, w_paths_mop, w_paths_mop_sd)
CI_lower_bound_uw <- mapply(compute_CI_lower, uw_paths_mop, uw_paths_mop_sd)

# compute CIs in case of non-normal distribution
#qnorm(1 - alpha/2)

alpha = 0.05
compute_CI_upper <- function(mean, sd){
  mean + qnorm(1 - alpha/2) * sd/sqrt(n)
}
CI_upper_bound_w <- mapply(compute_CI_upper, w_paths_mop, w_paths_mop_sd)
CI_upper_bound_uw <- mapply(compute_CI_upper, uw_paths_mop, uw_paths_mop_sd)

compute_CI_lower <- function(mean, sd){
  mean - qnorm(1 - alpha/2) * sd/sqrt(n)
}
CI_lower_bound_w <- mapply(compute_CI_lower, w_paths_mop, w_paths_mop_sd)
CI_lower_bound_uw <- mapply(compute_CI_lower, uw_paths_mop, uw_paths_mop_sd)




# BO_paths_ub_per_iter <- apply(BO_paths_optim, 1, function(optima){
#   boot_means <- boot(data=optima, statistic=Bmean, R=1000)
#   upper <- boot.ci(boot_means, type="basic")[["basic"]][5]
#   upper
# })
# BO_paths_lb_per_iter <- apply(BO_paths_optim, 1, function(optima){
#   boot_means <- boot(data=optima, statistic=Bmean, R=1000)
#   lower <- boot.ci(boot_means, type="basic")[["basic"]][4]
#   lower
# })
# 

fun = fun_nr
BO_paths_plot = data.frame("mean_w" = w_paths_mop[[fun]], "upper_w" = CI_upper_bound_w[,fun], "lower_w" = CI_lower_bound_w[,fun], 
                           "mean_uw" = uw_paths_mop[[fun]], "upper_uw" = CI_upper_bound_uw[,fun], "lower_uw" = CI_lower_bound_uw[,fun], 
                           "iter" = 1:20)

ggplot(data = BO_paths_plot) +
  geom_point(aes(x = iter, y = mean_w), color = "blue")  + 
  geom_line(aes(x = iter, y = mean_w), color = "blue", alpha = 0.5) +
  geom_errorbar(aes(x = iter, ymin = lower_w, ymax = upper_w), color = "blue", alpha = 0.5)  +
  geom_point(aes(x = iter, y = mean_uw), color = "magenta")  + 
  geom_line(aes(x = iter, y = mean_uw), color = "magenta") +
  geom_errorbar(aes(x = iter, ymin = lower_uw, ymax = upper_uw), color = "magenta", alpha = 0.5)  




# make results figure
funs_to_test <- all_benchmark_funs[c(3,4,6,8,27,31)]
plots_res = list()
for (f in 1:length(lambdas)) {
fun = f
BO_paths_plot = data.frame("mean_w" = w_paths_mop[[fun]], "upper_w" = CI_upper_bound_w[,fun], "lower_w" = CI_lower_bound_w[,fun], 
                           "mean_uw" = uw_paths_mop[[fun]], "upper_uw" = CI_upper_bound_uw[,fun], "lower_uw" = CI_lower_bound_uw[,fun], 
                           "iter" = 1:20)

plots_res[[f]] = ggplot(data = BO_paths_plot) +
  geom_point(aes(x = iter, y = mean_w), color = "blue")  + 
  geom_line(aes(x = iter, y = mean_w), color = "blue", alpha = 0.5) +
  geom_errorbar(aes(x = iter, ymin = lower_w, ymax = upper_w), color = "blue", alpha = 0.5)  +
  geom_point(aes(x = iter, y = mean_uw), color = "magenta")  + 
  geom_line(aes(x = iter, y = mean_uw), color = "magenta") +
  geom_errorbar(aes(x = iter, ymin = lower_uw, ymax = upper_uw), color = "magenta", alpha = 0.5) + 
  labs(x = "Evaluations", y = "Mean Best Target Value", title = paste("tau = ", as.character(lambdas[f])))

  }

plots_res <- plots_res[!sapply(plots_res,is.null)]
plots_res <- plots_res[1:6]
ggpubr::ggarrange(plotlist = plots_res)

# # weighted performs better: 
# plot(all_benchmark_funs[[27]])
# plot(all_benchmark_funs[[9]]) # minimal
# 
# # uw performs better 
# plot(all_benchmark_funs[[6]])
# plot(all_benchmark_funs[[19]]) # minimal
