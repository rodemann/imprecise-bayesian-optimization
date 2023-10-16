



# we compare infill optimizers with regard to exploration/exploitation

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

infill.opt = "ea"
scope = "local"

saved.xplxpl = matrix(NA, nrow = 16, ncol = 6)

#path = paste(getwd(),"/_Explore_Exploit_Measures/results-kapton-xplxpl/", sep = "")
#path = paste(getwd(),"/_Explore_Exploit_Measures/results-material-design-xplxpl/", sep = "")




  for (i in seq_along(infill.crits)) {
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
    mean.xpl.xpl = data.frame(mean(xpl.xpl$Local.Standard.Error.Ratio),
                             mean(xpl.xpl$Local.Standard.Error.Distribution.Value),
                             mean(xpl.xpl$Pairwise.Euclidean.Distance),
                             mean(xpl.xpl$Mean.Euclidean.Distance),
                             mean(xpl.xpl$Minimal.Euclidean.Distance),
                             mean(xpl.xpl$Maximal.Euclidean.Distance)
                             )
    print(mean.xpl.xpl)
   saved.xplxpl[i,] = unlist(mean.xpl.xpl)
  
}

# focus search:
#mean.xpl.xpl.fs = apply(saved.xplxpl, 2, mean)

# ea:
mean.xpl.xpl.ea = apply(saved.xplxpl, 2, mean)





colnames(mean.cor) = c("Standard Error Ratio", "Standard Error Distribution Value", "Pairwise Distance", "Mean Distance",
                       "Minimal Distance", "Maximal Distance")
rownames(mean.cor) = c("Standard Error Ratio", "Standard Error Distribution Value", "Pairwise Distance", "Mean Distance",
                       "Minimal Distance", "Maximal Distance")
mean.cor = mean.cor * lower.tri(mean.cor)
xtable(mean.cor)

