library(ggplot2)
library(tidyverse)
library(wesanderson)
library(RColorBrewer)
library(ggsci)
library(smoof)

source("kernel-fun-analysis/draw_plot.R")
source("kernel-fun-analysis/draw_legend.R")

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


# select palette
#pal <- wes_palette("Zissou1", 5)
#pal <- c("steelblue", "magenta", "forestgreen", "darkred", "orange1")
#pal <- display.brewer.pal(5, "Set 1")
pal <- pal_ucscgb("default")(5)
#pal <- pal_tron("legacy")(5)
optimum_col = "#ff00f7"

#  create leegend
draw_legend(pal)
# load results to environment
load("kernel-fun-analysis/BO_results_20_60_kernel_fun")
# select function


plots = list()

for (f in 1:length(obj.fun.list)) {

obj.fun <- obj.fun.list[f]
global.opt <- getGlobalOptimum(obj.fun[[1]])[["value"]]
fun_name <- getName(obj.fun[[1]])
if(f == 8)
  fun_name <- "Bivariate Ackley Function"
  
# plot all opt paths for selected function
plots[[f]] <- draw_plot(fun = f, mbo_results_paths = mbo_results_paths,
          initial_design_size = 10, pal, configs = 5) +
labs(x = "Iteration", y = "Mean Best Target Value") +
geom_hline(yintercept = global.opt, linetype = "dashed",  color = optimum_col) + 
  labs(title = paste("Bayesian Optimization of", fun_name), 
       subtitle = "40 BO runs per Kernel with 20 iterations each. 
Dotted pink line: Global Optimum.
Errorbars represent 0.95-CI of best target value.")  

}

plots

#obj.fun
#plot(obj.fun[[1]])
#plot3D(obj.fun[[1]])


