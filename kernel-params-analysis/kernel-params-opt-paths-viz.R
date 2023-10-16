library(ggplot2)
library(tidyverse)
library(wesanderson)
library(RColorBrewer)
library(ggsci)
library(smoof)
library(viridis)
source("kernel-params-analysis/draw_plot.R")
source("kernel-params-analysis/draw_legend_kernel-params.R")

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


# select palette
#pal <- wes_palette("Zissou1", 5)
pal <- c("steelblue", "magenta", "forestgreen", "darkred", "orange1")
pal <- brewer.pal(5, "Oranges")
#al <- pal_ucscgb("default")(5)
#pal <- pal_tron("legacy")(5)
pal <- viridis(5)

optimum_col = "#ff00f7"

#  create leegend
names <- c("l = 0.01", "l = 0.1", "l = 1", "l = 2", "l = 5")
draw_legend_kernel_param(pal, names)

# load results to environment
load("kernel-params-analysis/BO_results_40_20_kernel_bandwith")
# select function
#f = 4
plots = list()
for (f in 1:length(obj.fun.list)) {

obj.fun <- obj.fun.list[f]
global.opt <- getGlobalOptimum(obj.fun[[1]])[["value"]]
fun_name <- getName(obj.fun[[1]])
  
# plot all opt paths for selected function
plot =draw_plot(fun = f, mbo_results_paths = mbo_results_paths,
          initial_design_size = 10, pal, configs = 5) +
labs(x = "Iteration", y = "Mean Best Target Value") +
geom_hline(yintercept = global.opt, linetype = "dashed",  color = optimum_col) + 
  labs(title = paste("Bayesian Optimization of", fun_name), 
       subtitle = "40 BO runs per Kernel Bandwith with 20 iterations each. 
Dotted pink line: Global Optimum.
Errorbars show 0.95-CI of best target value.")  
plots[[f]] = plot

}
plots

#obj.fun
#plot(obj.fun[[1]])
#plot3D(obj.fun[[1]])
#plot3D(obj.fun[[1]], package = "plotly")








