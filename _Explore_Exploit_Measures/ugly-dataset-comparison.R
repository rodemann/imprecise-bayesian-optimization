



# we compare data sets kapton(matdesign) with regard to exploration/exploitation

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

infill.opt = "focussearch"
scope = "local"

saved.xplxpl = matrix(NA, nrow = 16, ncol = 6)

path = paste(getwd(),"/_Explore_Exploit_Measures/results-kapton-xplxpl-effect of initial data/", sep = "")
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
                            mean(xpl.xpl$Pairwise.Gower.Distance),
                            mean(xpl.xpl$Mean.Gower.Distance),
                            mean(xpl.xpl$Minimal.Gower.Distance),
                            mean(xpl.xpl$Maximal.Gower.Distance)
  )
  print(mean.xpl.xpl)
  saved.xplxpl[i,] = unlist(mean.xpl.xpl)
  
}

# material design
#mean.xpl.xpl.md = apply(saved.xplxpl, 2, mean)

# kapton:
mean.xpl.xpl.k = apply(saved.xplxpl, 2, mean)

mean.xpl.xpl.datasets = rbind(mean.xpl.xpl.md, mean.xpl.xpl.k)
colnames(mean.xpl.xpl.datasets) = c("Standard Error Ratio", "Standard Error Distribution Value", "Pairwise Distance", "Mean Distance",
                       "Minimal Distance", "Maximal Distance")
rownames(mean.xpl.xpl.datasets) = c("Material Design", "Kapton")
xtable(mean.xpl.xpl.datasets)

