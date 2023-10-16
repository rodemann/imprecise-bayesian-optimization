# Trivial dispatcher for infill optimizers.
#
# @note Keep in mind to update getSupportedInfillOptFunctions too,
# if a new method is implemented.
#
# @param infill.opt [\code{character(1)}]\cr
#   String key for infill optimizer.
# @return [\code{function}]


getInfillOptFunction <- function(infill.opt) {
  switch(infill.opt,
         cmaes = infillOptCMAES,
         cmaesSavepts = infillOptCMAESSavepts,
         focussearch = infillOptFocus,
         focussearchSavepts = infillOptFocusSavepts,
         ea = infillOptEA,
         eaSavepts = infillOptEASavepts,
         nsga2 = infillOptMultiObjNSGA2,
         # default: try to match the fun which is given as string
         match.fun(infill.opt)
  )
}


environment(getInfillOptFunction) <- asNamespace('mlrMBO')
# assures that the function will be able to call other hidden functions from the package.

assignInNamespace("getInfillOptFunction", getInfillOptFunction, ns = "mlrMBO") 
# assures that the function will be able to call other hidden functions from the package. (credits: TMS, stackoverflow)









