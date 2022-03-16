checkParametersDE_time <- function(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, step_method, direction_method,
                                   preprocess_method, tol1, tol2, nfolds, nsimulations, heatStep, heatIter, search)
{
  ################################################## Parameters Check ##################################################
  if (is.null(data))
    stop("'data' required; is NULL.")
  else{
    if(any(is.na(data)))
      stop("Missing values are not admitted in 'data'.")
  }

  if (is.null(data_time))
    stop("'data_time' required; is NULL.")
  else{
    if(any(is.na(data_time)))
      stop("Missing values are not admitted in 'data_time'.")
  }

  if (is.null(FEMbasis))
    stop("'FEMbasis' required; is NULL.")
  if(class(FEMbasis)!="FEMbasis")
    stop("'FEMbasis' is not class 'FEMbasis'.")

  if (is.null(mesh_time))
    stop("'mesh_time' required; is NULL.")
  else{
    for(i in 1:(length(mesh_time)-1)){
      if(mesh_time[i]>=mesh_time[i+1])
        stop("'mesh_time' nodes must be distinct and sorted in increasing order.")
    }
  }

  if (is.null(lambda))
    stop("'lambda' required; is NULL.")
  else{
    for(i in 1:length(lambda)){
      if(lambda[i]<=0)
        stop("'lambda' must have positive members.")
    }
  }

  if (is.null(lambda_time))
    stop("'lambda_time' requires; is NULL.")
  else{
    for(i in 1:length(lambda_time)){
      if(lambda_time[i]<=0)
        stop("'lambda_time' must have positive members.")
    }
  }

  if (is.null(step_method))
    stop("'step_method' is required; is NULL.")
  else{
    if(step_method!="Fixed_Step" && step_method!="Backtracking_Method" && step_method!="Wolfe_Method")
      stop("'step_method' needs to be either 'Fixed_Step' or 'Backtarcking_Method' or 'Wolfe_Method'.")
  }

  if (is.null(direction_method))
    stop("'direction_method' is required; is NULL.")
  else{
    if(direction_method!="Gradient" && direction_method!="ConjugateGradientFR" && direction_method!="ConjugateGradientPRP" && direction_method!="ConjugateGradientHS" && direction_method!="ConjugateGradientDY" && direction_method!="ConjugateGradientCD" && direction_method!="ConjugateGradientLS" && direction_method!="BFGS" && direction_method!="L-BFGS5" && direction_method!="L-BFGS10")
      stop("'direction_method' needs to be 'Gradient', 'ConjugateGradientFR', 'ConjugateGradientPRP', 'ConjugateGradientHS', 'ConjugateGradientDY', 'ConjugateGradientCD', 'ConjugateGradientLS', 'BFGS', 'L-BFGS5' or 'L-BFGS10'.")
  }

  if((length(lambda)>1 || length(lambda_time)>1) && preprocess_method!="RightCV" && preprocess_method!="SimplifiedCV")
    stop("'preprocess_method' needs to be either 'RightCV' or 'SimplifiedCV' if there is more than one possible pair of smoothing parameters 'lambda' and 'lambda_time'.")

  if((length(lambda)==1 && length(lambda_time)==1) && preprocess_method!="NoCrossValidation")
    stop("'preprocess_method' needs to be 'NoCrossValidation' if there is only one possible pair of smoothing parameters 'lambda' and 'lambda_time'.")

  if(preprocess_method=="SimplifiedCV" && (length(lambda)*length(lambda_time))!=nfolds)
    stop("'SimplifiedCV' requires the number of possible pairs of smoothing parameters in space and in time to be equal to the number of folds.")

  if(tol1 < 0 || tol2 < 0)
    stop("Tolerances 'tol1' and 'tol2' needs to be non negative numbers")

  if(length(lambda) > 1 && (!is.numeric(nfolds) || floor(nfolds)<=1))
    stop("'nfolds' needs to be an integer greater or equal than 2.")

  if(!is.numeric(nsimulations) || nsimulations<1)
    stop("'nsimulations' needs to be a positive integer.")

  if(!is.numeric(heatStep) || heatStep<0 || heatStep>1)
    stop("'heatStep' needs to be a positive real number not greater than 1.")

  if(!is.numeric(heatIter) || heatIter<1)
    stop("'heatIter' needs to be a positive integer.")

  if(!is.numeric(search))
    stop("'search' needs to be an integer.")


}


checkParametersSizeDE_time <- function(data, data_time, FEMbasis, mesh_time, ndim, fvec, preprocess_method, nfolds)
{
  if(nrow(data) < 1)
    stop("'data' must contain at least one element.")
  if(preprocess_method!="NoCrossValidation" && nrow(data) < floor(nfolds))
    stop("The number of folds needs to be less than the number of data.")
  if(ncol(data) != ndim)
    stop("'data' and the mesh points have incompatible size.")

  if(length(data_time) < 1)
    stop("'data_time' must contain at least one element.")
  if(nrow(data) != length(data_time))
    stop("'data' and 'data_time' must have the same number of rows (equal to the number of observations).")

  if(!is.null(fvec)){
    SPLINE_DEGREE <- 3
    if(length(fvec) != nrow(FEMbasis$mesh$nodes) * (length(mesh_time)+SPLINE_DEGREE-1))
      stop("The length of fvec has to be equal to the product of the number of mesh nodes and the number of B-spline basis functions.")
  }

  if(preprocess_method!="NoCrossValidation")
    if (nrow(data)*(floor(nfolds)-1)/floor(nfolds) < 30)
      stop("The training set needs to have at least 30 data: increase the number of folds.")

}
