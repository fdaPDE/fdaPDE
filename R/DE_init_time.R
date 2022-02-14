#' Spatio-temporal density initialization
#'
#' @param data A matrix of dimensions #observations-by-ndim. Data are locations: each row corresponds to one point,
#' the first column corresponds to the \code{x}-coordinates, the second column corresponds to the \code{y}-coordinates
#' and, if ndim=3, the third column corresponds to the \code{z}-coordinates.
#' @param data_time A vector of dimensions #observations. The i-th datum is the time instant during which the i-th location is observed
#' (according to the order in which data are provided).
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param mesh_time A vector containing the b-splines knots for separable smoothing. It is the time mesh of the considered time domain
#' (interval). Its nodes are in increasing order.
#' @param lambda A scalar or vector of smoothing parameters in space. Default is NULL. It is useful only if \code{init='Heat'}.
#' @param lambda_time A scalar or vector of smoothing parameters in time. Default is NULL. It is useful only if \code{init='Heat'}.
#' @param heatStep A real specifying the time step for the discretized heat diffusionn process.
#' @param heatIter An integer specifying the number of iteriations to perform the discretized heat diffusion process.
#' @param init A string specifying the initialization procedure. It can be either 'Heat' or 'CV'.
#' @param nFolds An integer specifying the number of folds used in the cross validation techinque. It is useful only
#' for the case \code{init = 'CV'}.
#' @param search A flag to decide the search algorithm type (tree or naive or walking search algorithm).
#' @param isTimeDiscrete A boolean specifying the time data type: \code{TRUE} for discrete (with many duplicates) time data;
#' \code{FALSE} for continuous time data. Default is \code{FALSE}.
#' @param flagMass A boolean specifying whether to consider full mass matrices (\code{TRUE}) or identity mass matrices
#' (\code{FALSE}) for the computation of space and time penalty matrices. Default is \code{FALSE}.
#' @param flagLumped A boolean specifying whether to perform mass lumping. This numerical technique presents computational
#' advantages during the procedure involving a mass matrix inversion for the computation of the space penalty matrix.
#' Default is \code{FALSE}.
#' N.B. We suggest to put it \code{TRUE} in case of a large spatial domain or in case of a dense/refined spatial mesh.
#' @return If \code{init = 'Heat'} it returns a matrix in which each column contains the initial vector
#' for each possible pair (\code{lambda}, \code{lambda_time}). If \code{init = 'CV'} it returns the initial vector associated
#' to the unique pair (\code{lambda}, \code{lambda_time}) given.
#' @description This function implements two methods for the density initialization procedure.
#' @usage DE.heat.FEM(data, data_time, FEMbasis, mesh_time, lambda=NULL, lambda_time=NULL, heatStep=0.1, heatIter=500,
#'                    init="Heat", nFolds=5, search="tree", isTimeDiscrete=0, flagMass=0, flagLumped=0)
#' @export
#' @examples
#'
#' ...................................
#'
#'

DE.heat.FEM_time <- function(data, data_time, FEMbasis, mesh_time, lambda=NULL, lambda_time=NULL, heatStep=0.1,
                             heatIter=500, init="Heat", nFolds=5, search="tree", isTimeDiscrete=0, flagMass=0, flagLumped=0)
{
  if(class(FEMbasis$mesh) == "mesh.2D"){
    ndim = 2
    mydim = 2
  }else if(class(FEMbasis$mesh) == "mesh.2.5D"){
    ndim = 3
    mydim = 2
  }else if(class(FEMbasis$mesh) == "mesh.3D"){
    ndim = 3
    mydim = 3
  }else{
    stop('Unknown mesh class')
  }

  fvec=NULL
  stepProposals=NULL
  tol1=NULL
  tol2=NULL
  print=NULL
  nfolds=NULL
  nsimulations=NULL
  step_method=NULL
  direction_method=NULL
  preprocess_method=NULL

  # Search algorithm
  if(search=="naive"){
    search=1
  }else if(search=="tree"){
    search=2
  }else if(search=="walking" & class(FEMbasis$mesh) == "mesh.2.5D"){
    stop("walking search is not available for mesh class mesh.2.5D.")
  }else if(search=="walking" & class(FEMbasis$mesh) != "mesh.2.5D"){
    search=3
  }else{
    stop("'search' must must belong to the following list: 'naive', 'tree' or 'walking'.")
  }


  ###################################### Checking parameters, sizes and conversion #####################################
  checkParametersDE_init(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, heatStep, heatIter, init, search)

  ## Converting to format for internal usage
  data = as.matrix(data)
  data_time = as.vector(data_time)
  lambda = as.vector(lambda)
  lambda_time = as.vector(lambda_time)
  mesh_time = as.vector(mesh_time)

  checkParametersSizeDE_init(data, data_time, ndim)
  #################################### End checking parameters, sizes and conversion ###################################


  ################################################# C++ Code Execution #################################################
  bigsol = NULL
  if(class(FEMbasis$mesh) == 'mesh.2D'){

    bigsol = CPP_FEM.DE_init_time(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, fvec, heatStep, heatIter, ndim,
                                  mydim, step_method, direction_method, preprocess_method, stepProposals, tol1, tol2, print,
                                  nfolds, nsimulations, search, isTimeDiscrete, flagMass, flagLumped, init, nFolds)

  } else if(class(FEMbasis$mesh) == 'mesh.2.5D'){

    bigsol = CPP_FEM.manifold.DE_init_time(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, fvec, heatStep, heatIter, ndim,
                                           mydim, step_method, direction_method, preprocess_method, stepProposals, tol1, tol2, print,
                                           nfolds, nsimulations, search, isTimeDiscrete, flagMass, flagLumped, init, nFolds)

  } else if(class(FEMbasis$mesh) == 'mesh.3D'){
    bigsol = CPP_FEM.volume.DE_init_time(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, fvec, heatStep, heatIter, ndim,
                                         mydim, step_method, direction_method, preprocess_method, stepProposals, tol1, tol2, print,
                                         nfolds, nsimulations, search, isTimeDiscrete, flagMass, flagLumped, init, nFolds)
  }

  ################################################### Collect Results ##################################################

  f_init = bigsol[[1]]

  reslist = list(f_init = f_init)
  return(reslist)
}
