draw_legend <- function(pal){ 
  kernel_names <- c("Gaussian","Matern (5,2)", "Matern (3,2)", 
                    "Exponential", "Power Exponential")
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("top", legend =kernel_names, pch=16, pt.cex=3, cex=1.5, bty='n',
         col = pal)
  mtext("Kernel", at=0.2, cex=2)
}