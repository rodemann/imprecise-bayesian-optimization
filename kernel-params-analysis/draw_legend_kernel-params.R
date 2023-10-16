draw_legend_kernel_param <- function(pal, names){ 
  kernel_names <- names
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("top", legend = kernel_names, pch=16, pt.cex=3, cex=1.5, bty='n',
         col = pal)
  mtext("Length-Scale Parameter", at=0.5, cex=2)
}