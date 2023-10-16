getSupportedMultipointInfillOptFunctions = function() {
  c("cl", "clXplxpl", "cb", "moimbo")
}

# assure that the function will be able to call other hidden functions from the package.
environment(getSupportedMultipointInfillOptFunctions) <- asNamespace('mlrMBO')

# assure that the function will be able to call other hidden functions from the package. (credits: TMS, stackoverflow)
assignInNamespace("getSupportedMultipointInfillOptFunctions", getSupportedMultipointInfillOptFunctions, ns = "mlrMBO") 

