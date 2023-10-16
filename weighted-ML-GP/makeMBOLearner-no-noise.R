

makeMBOLearner = function(control, fun, config = list(), ...) {
  assertClass(control, "MBOControl")
  assertClass(fun, "smoof_function")
  
  ps = getParamSet(fun)
  if (isSimpleNumeric(ps)) {
    lrn = makeLearner("regr.km", covtype = "matern3_2", optim.method = "gen", config = config)
    if (isNoisy(fun)){
      stop("BO_weighted_ML does not support noisy target functions")}
  } else {
    lrn = makeLearner("regr.randomForest", se.method = "jackknife", keep.inbag = TRUE, config = config)
    if (hasRequires(ps))
      lrn = makeImputeWrapper(lrn, classes = list(numeric = imputeMax(2), factor = imputeConstant("__miss__")))
  }
  
  if (control$infill.crit$requires.se || (!is.null(control$multipoint.method) && control$multipoint.method == "moimbo"))
    lrn = setPredictType(lrn, "se")
  
  lrn = setHyperPars(lrn, ...)
  return(lrn)
}

R.utils::reassignInPackage("makeMBOLearner", pkgName = "mlrMBO", makeMBOLearner)
