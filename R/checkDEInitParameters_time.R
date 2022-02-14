checkParametersDE_init_time <- function(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, heatStep, heatIter, init, search)
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
  if(class(FEMbasis)!= "FEMbasis")
    stop("'FEMbasis' is not class 'FEMbasis'.")

  if (is.null(mesh_time))
    stop("'mesh_time' required; is NULL.")
  else{
    for(i in 1:(length(mesh_time)-1)){
      if(mesh_time[i]>=mesh_time[i+1])
        stop("'mesh_time' nodes must be distinct and sorted in increasing order.")
    }
  }

  if(init=="Heat"){
    if (is.null(lambda) || is.null(lambda_time))
      stop("'lambda' and 'lambda_time' required if init='Heat'; at least one of them is NULL.")
    else{
      for(i in 1:length(lambda)){
        if(lambda[i]<=0)
          stop("'lambda' must have positive members.")
      }
      for(i in 1:length(lambda_time)){
        if(lambda_time[i]<=0)
          stop("'lambda_time' must have positive members.")
      }
    }
  }

  if(!is.numeric(search))
    stop("'search' needs to be an integer.")

  if (is.null(init))
    stop("'init' is required;  is NULL.")
  else{
    if(init!="Heat" && init!="CV")
      stop("'init' needs to be either 'Heat' or 'CV'.")
  }

  if(init=="CV" && (length(lambda)>1 || length(lambda_time)>1))
    stop("The initialization procedure via cross-validation is only for a given pair of (lambda, lambda_time).")
}


checkParametersSizeDE_init_time <- function(data, data_time, ndim)
{
  if(nrow(data) < 1)
    stop("'data' must contain at least one element.")
  if(ncol(data) != ndim)
    stop("'data' and the mesh points have incompatible size.")

  if(length(data_time) < 1)
    stop("'data_time' must contain at least one element.")
  if(nrow(data) != length(data_time))
    stop("'data' and 'data_time' must have the same number of rows (equal to the number of observations).")
}
