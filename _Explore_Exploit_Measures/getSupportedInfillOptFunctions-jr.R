#' @title Get names of supported infill-criteria optimizers.
#' @description
#' None.
#' @return [\code{character}]
#' @export
getSupportedInfillOptFunctions = function() {
  c("cmaes", "cmaesSavepts", "focussearch", "focussearchSavepts", "ea", "eaSavepts", "nsga2")
}


# assure that the function will be able to call other hidden functions from the package.
environment(getSupportedInfillOptFunctions) <- asNamespace('mlrMBO')

# assure that the function will be able to call other hidden functions from the package. (credits: TMS, stackoverflow)
assignInNamespace("getSupportedInfillOptFunctions", getSupportedInfillOptFunctions, ns = "mlrMBO") 


