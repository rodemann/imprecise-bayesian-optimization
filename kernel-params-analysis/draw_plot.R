library(ggplot2)
library(tidyverse)
library(wesanderson)
library(RColorBrewer)
library(ggsci)
library(ggshadow)


# get BO path of Kernel k and function f
draw_plot <- function(fun, mbo_results_paths, initial_design_size, pal, configs){
  plot <- ggplot()

  for (k in 1:configs) {
    BO_paths <- mbo_results_paths[[k]][[fun]]
    #BO_paths_global_min <- lapply(BO_paths, min)
    
    BO_paths <- matrix(unlist(BO_paths), ncol = length(BO_paths))
    BO_paths <- as.data.frame(BO_paths)
    
    #remove initial design
    BO_paths_initial <- BO_paths %>% slice_head(n = initial_design_size)
    BO_paths_optim <- BO_paths %>% slice_tail(n = nrow(BO_paths) - initial_design_size)
    
    # log scale
    #BO_paths_optim <- log(BO_paths_optim) 
    
    #get lower bound/upper bound/mean per iteration:
    BO_paths_sd <- apply(BO_paths_optim, 1, sd)
    BO_paths_mean_per_iter <- apply(BO_paths_optim, 1, mean)
    BO_paths_min_per_iter <- apply(BO_paths_optim, 1, min)
    BO_paths_max_per_iter <- apply(BO_paths_optim, 1, max)
    BO_paths_ub_per_iter <- BO_paths_mean_per_iter + qnorm(.975) * BO_paths_sd/sqrt(ncol(BO_paths))
    BO_paths_lb_per_iter <- BO_paths_mean_per_iter - qnorm(.975) * BO_paths_sd/sqrt(ncol(BO_paths))
    
    paths <- data.frame(iter = 1:nrow(BO_paths_optim), "Upper CB" = BO_paths_ub_per_iter, 
                        "Lower CB" = BO_paths_lb_per_iter,
                        "Mean Target Value" = BO_paths_mean_per_iter)
    
    plot <- plot + geom_point(data = paths, aes(x = iter, y = Mean.Target.Value),
                                              color = pal[k], show.legend = FALSE) + 
      geom_glowline(data = paths, aes(x = iter, y = Mean.Target.Value), color = pal[k]) +
      geom_errorbar(data = paths, aes(x = iter, ymin = Lower.CB, 
                                    ymax = Upper.CB), colour = pal[k]) 

}
  #plot <- plot + scale_color_discrete(name = "Kernel", labels = kernel_names)
  plot + theme(panel.background = element_rect(fill = "#575656"))
               
}
