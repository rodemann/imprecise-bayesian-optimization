



# we analyse the correlations of our xplxpl measures and print a latex-ready table
library(xtable)



infill.crits = list(makeMBOInfillCritEI(),makeMBOInfillCritAdaCB(),
                    makeMBOInfillCritStandardError(), makeMBOInfillCritMeanResponse(),
                    makeMBOInfillCritAEI(), makeMBOInfillCritEQI(),
                    makeMBOInfillCritCB(cb.lambda = 0.1), 
                    makeMBOInfillCritCB(cb.lambda = 1),
                    makeMBOInfillCritCB(cb.lambda = 2),
                    makeMBOInfillCritCB(cb.lambda = 3),
                    makeMBOInfillCritCB(cb.lambda = 4),
                    makeMBOInfillCritCB(cb.lambda = 5),
                    makeMBOInfillCritCB(cb.lambda = 10),
                    makeMBOInfillCritCB(cb.lambda = 50),
                    makeMBOInfillCritCB(cb.lambda = 100),
                    makeMBOInfillCritCB(cb.lambda = 1000))

infill.opts = c("focussearch", "ea")
scope = c("local", "global")


scope = scope[1]
cor.matrix.list = matrix(NA, nrow = 2*16, ncol = 6*6)

path = paste(getwd(),"/_Explore_Exploit_Measures/results-kapton-xplxpl/", sep = "")
path = paste(getwd(),"/_Explore_Exploit_Measures/results-material-design-xplxpl/", sep = "")


for (j in seq_along(infill.opts)) {
  for (i in seq_along(infill.crits)) {
    infill.opt = infill.opts[j]  
    infill.criterion = infill.crits[[i]]
    ctrl = setMBOControlInfill(ctrl, opt = infill.opt, crit = infill.criterion)

    lambda = infill.criterion$params$cb.lambda
    if (is.null(lambda)) {
      lambda = "NA"
    }
    
    # save stored data frame as xpl.xpl
    load(paste(path, "xpl.xpl.", scope, ".", infill.opt, ".", infill.criterion$id,".lambda=", lambda, ".RData", sep = ""))
    xpl.xpl = xpl.xpl[,-7]
    print(paste(scope, ".", infill.opt, ".", infill.criterion$id,".lambda=", lambda, ".RData", sep = ""))
    print(cor(xpl.xpl))
    global.iter = (j - 1)*length(infill.crits) + i
    cor.matrix.list[global.iter,] = as.vector(cor(xpl.xpl))
   
  }
}

mean.cor = apply(cor.matrix.list, 2, mean, na.rm = TRUE)
mean.cor = matrix(mean.cor, nrow = 6)
colnames(mean.cor) = c("Standard Error Ratio", "Standard Error Distribution Value", "Pairwise Distance", "Minimal Distance",
                       "Maximal Distance", "Mean Distance")
rownames(mean.cor) = c("Standard Error Ratio", "Standard Error Distribution Value", "Pairwise Distance", "Minimal Distance",
                       "Maximal Distance", "Mean Distance")
mean.cor = mean.cor * lower.tri(mean.cor)
xtable(mean.cor)

