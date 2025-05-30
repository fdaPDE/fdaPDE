#' Nonparametric spatio-temporal density estimation with differential regularization
#'
#' @param data A matrix of dimensions #observations-by-ndim. Data are locations: each row corresponds to one point,
#' the first column corresponds to the \code{x}-coordinates, the second column corresponds to the \code{y}-coordinates
#' and, if ndim=3, the third column corresponds to the \code{z}-coordinates.
#' @param data_time A vector of length #observations. The i-th datum is the time instant during which the i-th location is
#' observed (according to the order in which locations are provided in data).
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param mesh_time A vector containing the b-splines knots for separable smoothing. It is the time mesh of the considered time domain
#' (interval). Its nodes are in increasing order.
#' @param lambda A scalar or vector of smoothing parameters in space. If it is a vector, the optimal smoothing parameter in space
#' is chosen, together with the optimal smoothing parameter in time, with a \code{k}-fold cross-validation procedure based on the L2 norm.
#' @param lambda_time A scalar or vector of smoothing parameters in time. If it is a vector, the optimal smoothing parameter in time
#' is chosen, together with the optimal smoothing parameter in space, with a \code{k}-fold cross-validation procedure based on the L2 norm.
#' @param scaling A positive factor needed to scale the smoothing parameters in the construction of confidence intervals.
#' If the scaling is not specified, it is automatically set as the square root of the number of observations.
#' @param fvec A vector of length #\code{nodes} of the spatial mesh times #\code{B-spline} temporal functional basis. It corresponds to the
#' node values of the initial density function. If this is \code{NULL} the initial density is estimated thanks to a discretized heat diffusion
#' process that starts from the empirical density of the data. Default is \code{NULL}.
#' N.B. This vector cannot be the constant vector of zeros since the algorithm works with the log(f).
#' @param heatStep A real specifying the time step for the discretized heat diffusion process. Default is \code{0.1}.
#' @param heatIter An integer specifying the number of iterations to perform the discretized heat diffusion process. Default is \code{10}.
#' @param stepProposals A scalar or a vector containing the step parameters useful for the descent algorithm. If there is a
#' vector of parameters, the biggest one such that the functional decreases at each iteration is chosen. If it is \code{NULL}
#' the following vector \code{c(0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 1e-7, 1e-8, 1e-9)} is proposed. Default is \code{NULL}.
#' N.B. If the program does not receive a right parameter, it aborts the R session. Try a smaller parameter.
#' @param tol1 A scalar specifying the tolerance to use for the termination criterion based on the percentage difference
#' between two consecutive iterations of the minimization algorithm of the loss function, the log-likelihood and the
#' penalizations. Default is 1e-5.
#' @param tol2 A scalar specifying the tolerance to use for the termination criterion based on the norm of the gradient
#' of the functional to be minimized (the true minimum is such that this norm is zero). The default version does not use this
#' criterion. Default is 0.
#' @param print A boolean that is \code{TRUE} if the user wants the value of the functional, of the loglikelihood and of the
#' penalization terms printed on console at each iteration of the descent algorithm (plus some other information/warnings). Default is \code{FALSE}.
#' N.B. We suggest to let it \code{FALSE} if \code{preprocess_method} is 'RightCV' or 'SimplifiedCV'.
#' @param nfolds An integer specifying the number of folds used in cross validation technique to find the best pair of
#' (\code{lambda}, \code{lambda_time}) smoothing parameters.
#' If there is only one pair of (\code{lambda}, \code{lambda_time}) it can be \code{0}. Default is \code{0}.
#' @param nsimulations An integer specifying the number of iterations used in the optimization algorithms. Default value is 500.
#' @param step_method A string specifying which step method to use in the descent algorithm.
#' If it is \code{Fixed_Step}, the step is constant during the algorithm and it is chosen according to \code{stepProposals};
#' if it is \code{Backtracking_Method}, the step is computed at each iteration according to the backtracking method; finally
#' if it is \code{Wolfe_Method}, the step is computed at each iteration according to the Wolfe method. Default is \code{Fixed_Step}.
#' @param direction_method A string specifying which descent direction to use in the descent algorithm.
#' If it is \code{Gradient}, the direction is the one given by the gradient descent method (the opposite to the gradient of
#' the functional); if instead it is \code{BFGS} the direction is the one given by the BFGS method
#' (Broyden-Fletcher-Goldfarb-Shanno, a Quasi-Newton method). Default is \code{BFGS}. Other possible choices:
#' Conjugate Gradient direction with Fletcher-Reeves formula (\code{ConjugateGradientFR}), Conjugate Gradient direction with
#' Polak-Ribiére-Polyak formula (\code{ConjugateGradientPRP}), Conjugate Gradient direction with Hestenes-Stiefel formula
#' (\code{ConjugateGradientHS}), Conjugate Gradient direction with Dai-Yuan formula (\code{ConjugateGradientDY}), Conjugate
#' Gradient direction with Conjugate-Descent formula (\code{ConjugateGradientCD}), Conjugate Gradient direction with Liu-Storey
#' formula (\code{ConjugateGradientLS}), L-BFGS direction with 5 correction vectors (\code{L-BFGS5}), L-BFGS direction with 10
#' correction vectors (\code{L-BFGS10}).
#' @param preprocess_method A string specifying the k fold cross validation technique to use, if there is more than one pair of
#' smoothing parameters in space and in time (\code{lambda}, \code{lambda_time}); otherwise it should be \code{NULL}.
#' If it is \code{RightCV} the usual k fold cross validation method is performed. If it is \code{SimplifiedCV} a simplified
#' version is performed. In the latter case the number of possible pairs of smoothing parameters in space and in time
#' (\code{lambda}, \code{lambda_time}) must be equal to the number of folds \code{nfolds}. Default is \code{NULL}.
#' @param search A flag to decide the search algorithm type (tree or naive or walking search algorithm). Default is \code{tree}.
#' @param isTimeDiscrete A boolean specifying the time data type: \code{TRUE} for discrete (with many duplicates) time data;
#' \code{FALSE} for continuous time data. Default is \code{FALSE}.
#' @param flagMass A boolean specifying whether to consider full mass matrices (\code{TRUE}) or identity mass matrices
#' (\code{FALSE}) for the computation of space and time penalty matrices. Default is \code{FALSE}.
#' @param flagLumped A boolean specifying whether to perform mass lumping. This numerical technique presents computational
#' advantages during the procedure involving a mass matrix inversion for the computation of the space penalty matrix.
#' Default is \code{FALSE}.
#' N.B. We suggest to use it as \code{TRUE} in case of a large spatial domain or in case of a dense/refined spatial mesh.
#' @param inference A boolean that is \code{TRUE} if the user wants to estimate confidence intervals. Default is \code{FALSE}.
#' @return A list with the following variables:
#' \item{\code{FEMbasis}}{Given FEMbasis with tree information.}
#' \item{\code{g}}{A vector of length #\code{nodes} times #\code{B-splines} that represents the value of the g-function estimated for
#' each \code{node} of the spatial mesh and at each time instant of the time mesh. The density is the exponential of this function.}
#' \item{\code{f_init}}{A #\code{nodes}-by-#\code{lambda}x#\code{lambda_time} parameters matrix. Each column contains the node values of the initial
#' density used for the pair (\code{lambda}, \code{lambda_time}) given by the column.}
#' \item{\code{lambda}}{A scalar representing the optimal smoothing parameter in space selected, together with \code{lambda_time},
#' via k fold cross validation, if in the input there is a vector of parameters (in space and/or in time); the scalar given in input otherwise.}
#' \item{\code{lambda_time}}{A scalar representing the optimal smoothing parameter in time selected, together with \code{lambda},
#' via k fold cross validation, if in the input there is a vector of parameters (in space and/or in time); the scalar given in input otherwise.}
#' \item{\code{data}}{A matrix of dimensions #observations-by-ndim containing the spatial data used in the algorithm. They are the
#' same given in input if the domain is 2D pr 3D; they are the original data projected on the mesh if the domain is 2.5D. Data lying
#' outside the spatial domain, defined through its mesh, are not considered.}
#' \item{\code{data_time}}{A vector of length #observations containing the time data used in the algorithm. Data lying
#' outside the temporal domain, defined through its mesh, are not considered.}
#' \item{\code{CV_err}}{A vector of length \code{nfolds} containing the cross validation errors obtained in each fold, if
#' \code{preprocess_method} is either \code{RightCV} or \code{SimplifiedCV}.}
#' @description This function implements a nonparametric spatio-temporal density estimation method with differential regularization
#' (given by the sum of the square of the L2 norm of the laplacian of the density function and the square of the L2 norm of the second-
#' order time-derivative), when points are located over a planar mesh. The computation relies only on the C++ implementation of the algorithm.
#' @usage DE.FEM.time(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, scaling=NULL,
#'                    fvec=NULL, heatStep=0.1, heatIter=10, stepProposals=NULL, tol1=1e-4,
#'                    tol2=0, print=FALSE, nfolds=0, nsimulations=500,
#'                    step_method="Fixed_Step", direction_method="BFGS",
#'                    preprocess_method="NoCrossValidation", search="tree",
#'                    isTimeDiscrete=FALSE, flagMass=FALSE, flagLumped=FALSE,
#'                    inference = FALSE)
#' @export
#' @examples
#' library(fdaPDE)
#'
#' ## Create a 2D mesh over a squared domain
#' Xbound <- seq(-3, 3, length.out = 10)
#' Ybound <- seq(-3, 3, length.out = 10)
#' grid_XY <- expand.grid(Xbound, Ybound)
#' Bounds <- grid_XY[(grid_XY$Var1 %in% c(-3, 3)) | (grid_XY$Var2 %in% c(-3, 3)), ]
#' mesh <- create.mesh.2D(nodes = Bounds, order = 1)
#' mesh <- refine.mesh.2D(mesh, maximum_area = 0.25, minimum_angle = 20)
#' FEMbasis <- create.FEM.basis(mesh)
#'
#' ## Create a 1D time mesh over a (non-negative) interval
#' mesh_time <- seq(0, 1, length.out=3)
#'
#' ## Generate data
#' n <- 50
#' set.seed(10)
#' x <- rnorm(n,0,2)
#' y <- rnorm(n,0,2)
#' locations <- cbind(x,y)
#' times <- runif(n,0,1)
#' data <- cbind(locations, times)
#'
#' t <- 0.5 # time instant in which to evaluate the solution
#'
#' plot(mesh)
#' sample <- data[abs(data[,3]-t)<0.05,1:2]
#' points(sample, col="red", pch=19, cex=1, main=paste('Sample | ', t-0.05,' < t < ', t+0.05))
#'
#' ## Density Estimation
#' lambda <- 0.1
#' lambda_time <- 0.001
#' sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis, mesh_time = mesh_time,
#'                    lambda = lambda, lambda_time = lambda_time, nsimulations=300, inference=TRUE)
#'
#' ## Visualization
#' n = 100
#' X <- seq(-3, 3, length.out = n)
#' Y <- seq(-3, 3, length.out = n)
#' grid <- expand.grid(X, Y)
#'
#' FEMfunction = FEM.time(sol$g, mesh_time, FEMbasis, FLAG_PARABOLIC = FALSE)
#' evaluation <- eval.FEM.time(FEM.time = FEMfunction, locations = grid, time.instants = t)
#' FEMfunction_L = FEM.time(sol$g_CI_L, mesh_time, FEMbasis, FLAG_PARABOLIC = FALSE)
#' evaluation_L <- eval.FEM.time(FEM.time = FEMfunction_L, locations = grid, time.instants = t)
#' FEMfunction_U = FEM.time(sol$g_CI_U, mesh_time, FEMbasis, FLAG_PARABOLIC = FALSE)
#' evaluation_U <- eval.FEM.time(FEM.time = FEMfunction_U, locations = grid, time.instants = t)
#'
#' image2D(x = X, y = Y, z = matrix(exp(evaluation_L), n, n), col = heat.colors(100),
#'         xlab = "x", ylab = "y", contour = list(drawlabels = FALSE),
#'         main = paste("Estimated CI lower bound at t = ", t), zlim=c(0,0.3), asp = 1)
#' image2D(x = X, y = Y, z = matrix(exp(evaluation), n, n), col = heat.colors(100),
#'         xlab = "x", ylab = "y", contour = list(drawlabels = FALSE),
#'         main = paste("Estimated density at t = ", t), zlim=c(0,0.3), asp = 1)
#' image2D(x = X, y = Y, z = matrix(exp(evaluation_U), n, n), col = heat.colors(100),
#'         xlab = "x", ylab = "y", contour = list(drawlabels = FALSE),
#'         main = paste("Estimated CI upper bound at t = ", t), zlim=c(0,0.3), asp = 1)


DE.FEM.time <- function(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, scaling=NULL, fvec=NULL, heatStep=0.1, heatIter=10,
                        stepProposals=NULL, tol1=1e-4, tol2=0, print=FALSE, nfolds=0, nsimulations=500,
                        step_method="Fixed_Step", direction_method="BFGS", preprocess_method="NoCrossValidation",
                        search="tree", isTimeDiscrete=FALSE, flagMass=FALSE, flagLumped=FALSE, inference=FALSE)
{
  if(is(FEMbasis$mesh, "mesh.1.5D")){
    ndim = 2
    mydim = 1
  }else if(is(FEMbasis$mesh, "mesh.2D")){
    ndim = 2
    mydim = 2
  }else if(is(FEMbasis$mesh, "mesh.2.5D")){
    ndim = 3
    mydim = 2
  }else if(is(FEMbasis$mesh, "mesh.3D")){
    ndim = 3
    mydim = 3
  }else{
    stop('Unknown mesh class')
  }

  # Search algorithm
  if(search=="naive"){
    search=1
  }else if(search=="tree"){
    search=2
  }else if(search=="walking" & is(FEMbasis$mesh, "mesh.1.5D")){
    stop("walking search is not available for mesh class mesh.1.5D.")
  }else if(search=="walking" & is(FEMbasis$mesh, "mesh.2.5D")){
    stop("walking search is not available for mesh class mesh.2.5D.")
  }else if(search=="walking" & !is(FEMbasis$mesh, "mesh.1.5D") &
           !is(FEMbasis$mesh, "mesh.2.5D")){
    search=3
  }else{
    stop("'search' must must belong to the following list: 'naive', 'tree' or 'walking'.")
  }

  if(is.null(scaling))
    scaling=sqrt(dim(data)[1])

  ###################################### Checking parameters, sizes and conversion #####################################
  checkParametersDE_time(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, scaling, step_method, direction_method,
                         preprocess_method, tol1, tol2, nfolds, nsimulations, heatStep, heatIter, search)

  ## Converting to format for internal usage
  data = as.matrix(data)
  data_time = as.vector(data_time)
  lambda = as.vector(lambda)
  lambda_time = as.vector(lambda_time)
  mesh_time = as.vector(mesh_time)
  if(!is.null(fvec))
    fvec = as.vector(fvec)
  if(!is.null(stepProposals))
    stepProposals = as.vector(stepProposals)

  checkParametersSizeDE_time(data, data_time, FEMbasis, mesh_time, ndim, fvec, preprocess_method, nfolds)
  #################################### End checking parameters, sizes and conversion ###################################


  ################################################# C++ Code Execution #################################################
  bigsol = NULL
  if(is(FEMbasis$mesh, "mesh.2D")){

    bigsol = CPP_FEM.DE_time(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, scaling, fvec, heatStep, heatIter, ndim,
                             mydim, step_method, direction_method, preprocess_method, stepProposals, tol1, tol2, print,
                             nfolds, nsimulations, search, isTimeDiscrete, flagMass, flagLumped, inference)

  } else if(is(FEMbasis$mesh, "mesh.2.5D")){

    bigsol = CPP_FEM.manifold.DE_time(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, scaling, fvec, heatStep, heatIter, ndim,
                                      mydim, step_method, direction_method, preprocess_method, stepProposals, tol1, tol2, print,
                                      nfolds, nsimulations, search, isTimeDiscrete, flagMass, flagLumped, inference)

  } else if(is(FEMbasis$mesh, "mesh.3D")){
    bigsol = CPP_FEM.volume.DE_time(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, scaling, fvec, heatStep, heatIter, ndim,
                                    mydim, step_method, direction_method, preprocess_method, stepProposals, tol1, tol2, print,
                                    nfolds, nsimulations, search, isTimeDiscrete, flagMass, flagLumped, inference)
    
  } else if(is(FEMbasis$mesh, "mesh.1.5D")){
    bigsol = CPP_FEM.graph.DE_time(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, scaling, fvec, heatStep, heatIter, ndim,
                                   mydim, step_method, direction_method, preprocess_method, stepProposals, tol1, tol2, print,
                                   nfolds, nsimulations, search, isTimeDiscrete, flagMass, flagLumped, inference)
    
  }

  ################################################### Collect Results ##################################################

  g = bigsol[[1]]
  f_init = bigsol[[2]]
  lambda = bigsol[[3]]
  lambda_time = bigsol[[4]]
  data = bigsol[[5]]
  data_time = bigsol[[6]]
  CV_err = bigsol[[7]]
  g_CI_L = bigsol[[8]]
  g_CI_U = bigsol[[9]]

  # Save information of Tree Mesh
  tree_mesh = list(
    treelev = bigsol[[10]][1],
    header_orig= bigsol[[11]],
    header_scale = bigsol[[12]],
    node_id = bigsol[[13]][,1],
    node_left_child = bigsol[[13]][,2],
    node_right_child = bigsol[[13]][,3],
    node_box= bigsol[[14]])

  # Reconstruct FEMbasis with tree mesh
  mesh.class = class(FEMbasis$mesh)
  if (is.null(FEMbasis$mesh$treelev)) { # If does not exist the tree information
    FEMbasis$mesh = append(FEMbasis$mesh, tree_mesh)
  } # If already exists the tree information, do not append
  class(FEMbasis$mesh) = mesh.class

  reslist = list(FEMbasis = FEMbasis, g = g, f_init = f_init, lambda = lambda, lambda_time = lambda_time, data = data,
                 data_time = data_time, CV_err = CV_err, g_CI_L = g_CI_L, g_CI_U = g_CI_U)
  return(reslist)
}

