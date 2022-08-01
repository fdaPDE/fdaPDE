CPP_smooth.volume.FEM.basis<-function(locations, observations, FEMbasis, covariates = NULL, ndim, mydim, BC = NULL, incidence_matrix = NULL, areal.data.avg = TRUE, search, bary.locations, optim, lambda = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1, lambda.optimization.tolerance = 0.05, inference.data.object)
{

  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(DOF.matrix))
  {
    DOF.matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }

  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  {
    BC$BC_indices<-as.vector(BC$BC_indices)-1

  }

  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }
  
  if(is.null(lambda))
  {
    lambda<-vector(length=0)
  }else
  {
    lambda<-as.vector(lambda)
  }
  
  ## Extract the parameters for inference from inference.data.object to prepare them for c++ reading
  test_Type<-as.vector(inference.data.object@test)
  interval_Type<-as.vector(inference.data.object@interval)
  implementation_Type<-as.vector(inference.data.object@type)
  component_Type<-as.vector(inference.data.object@component)
  exact_Inference<-inference.data.object@exact
  locs_Inference<-as.matrix(inference.data.object@locations)
  locs_index_Inference<-as.vector(inference.data.object@locations_indices - 1) #converting indices from R to c++ ones
  locs_are_nodes_Inference <- inference.data.object@locations_are_nodes
  coeff_Inference<-as.matrix(inference.data.object@coeff)
  beta_0<-as.vector(inference.data.object@beta0)
  f_0_eval<-as.vector(inference.data.object@f0_eval)
  f_var_Inference<-inference.data.object@f_var
  inference_Quantile<-as.vector(inference.data.object@quantile)
  inference_Alpha<-as.vector(inference.data.object@alpha)
  inference_N_Flip<-inference.data.object@n_flip
  inference_Tol_Fspai<-inference.data.object@tol_fspai
  inference_Defined<-inference.data.object@definition
  

  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  data <- as.vector(observations)
  storage.mode(data) <- "double"
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$faces) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"  
  storage.mode(lambda) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  
  ## Set proper type for correct C++ reading for inference parameters
  storage.mode(test_Type) <- "integer"
  storage.mode(interval_Type) <- "integer"
  storage.mode(implementation_Type) <- "integer"
  storage.mode(component_Type) <- "integer"
  storage.mode(exact_Inference) <- "integer"
  storage.mode(locs_Inference) <- "double"
  storage.mode(locs_index_Inference) <- "integer"
  storage.mode(locs_are_nodes_Inference) <- "integer"
  storage.mode(coeff_Inference) <- "double"
  storage.mode(beta_0) <- "double"
  storage.mode(f_0_eval) <- "double"
  storage.mode(f_var_Inference) <- "integer"
  storage.mode(inference_Quantile) <- "double"
  storage.mode(inference_Alpha) <- "double"
  storage.mode(inference_N_Flip) <- "integer"
  storage.mode(inference_Tol_Fspai) <- "double"
  storage.mode(inference_Defined) <- "integer"
  
  ## Call C++ function
  bigsol <- .Call("regression_Laplace", locations, bary.locations, data, FEMbasis$mesh, FEMbasis$mesh$order, mydim, ndim, covariates,
                  BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg, search,
                  optim, lambda, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance,
                  test_Type,interval_Type,implementation_Type,component_Type,exact_Inference,locs_Inference,locs_index_Inference,locs_are_nodes_Inference,coeff_Inference,beta_0,
                  f_0_eval,f_var_Inference,inference_Quantile,inference_Alpha,inference_N_Flip,inference_Tol_Fspai,inference_Defined,
                  PACKAGE = "fdaPDE")

  return(bigsol)
}

CPP_smooth.volume.FEM.PDE.basis<-function(locations, observations, FEMbasis, covariates = NULL, PDE_parameters, ndim, mydim, BC = NULL, incidence_matrix = NULL, areal.data.avg = TRUE, search, bary.locations, optim, lambda = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1, lambda.optimization.tolerance = 0.05, inference.data.object)
{

  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  #
  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(DOF.matrix))
  {
    DOF.matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }

  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  {
    BC$BC_indices<-as.vector(BC$BC_indices)-1
  }

  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }
  
  ## Extract the parameters for inference from inference.data.object to prepare them for c++ reading
  test_Type<-as.vector(inference.data.object@test)
  interval_Type<-as.vector(inference.data.object@interval)
  implementation_Type<-as.vector(inference.data.object@type)
  component_Type<-as.vector(inference.data.object@component)
  exact_Inference<-inference.data.object@exact
  locs_Inference<-as.matrix(inference.data.object@locations)
  locs_index_Inference<-as.vector(inference.data.object@locations_indices - 1) #converting indices from R to c++ ones
  locs_are_nodes_Inference <- inference.data.object@locations_are_nodes
  coeff_Inference<-as.matrix(inference.data.object@coeff)
  beta_0<-as.vector(inference.data.object@beta0)
  f_0_eval<-as.vector(inference.data.object@f0_eval)
  f_var_Inference<-inference.data.object@f_var
  inference_Quantile<-as.vector(inference.data.object@quantile)
  inference_Alpha<-as.vector(inference.data.object@alpha)
  inference_N_Flip<-inference.data.object@n_flip
  inference_Tol_Fspai<-inference.data.object@tol_fspai
  inference_Defined<-inference.data.object@definition

  ## Set propr type for correct C++ reading

  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  data <- as.vector(observations)
  storage.mode(data) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(FEMbasis$mesh$faces) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"

  storage.mode(PDE_parameters$K) <- "double"
  storage.mode(PDE_parameters$b) <- "double"
  storage.mode(PDE_parameters$c) <- "double"

  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"  
  storage.mode(lambda) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  
  ## Set proper type for correct C++ reading for inference parameters
  storage.mode(test_Type) <- "integer"
  storage.mode(interval_Type) <- "integer"
  storage.mode(implementation_Type) <- "integer"
  storage.mode(component_Type) <- "integer"
  storage.mode(exact_Inference) <- "integer"
  storage.mode(locs_Inference) <- "double"
  storage.mode(locs_index_Inference) <- "integer"
  storage.mode(locs_are_nodes_Inference) <- "integer"
  storage.mode(coeff_Inference) <- "double"
  storage.mode(beta_0) <- "double"
  storage.mode(f_0_eval) <- "double"
  storage.mode(f_var_Inference) <- "integer"
  storage.mode(inference_Quantile) <- "double"
  storage.mode(inference_Alpha) <- "double"
  storage.mode(inference_N_Flip) <- "integer"
  storage.mode(inference_Tol_Fspai) <- "double"
  storage.mode(inference_Defined) <- "integer"

  ## Call C++ function
  bigsol <- .Call("regression_PDE", locations, bary.locations, data, FEMbasis$mesh, FEMbasis$order, mydim, ndim, PDE_parameters$K, PDE_parameters$b, PDE_parameters$c, covariates,
                  BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg, search,
                  optim, lambda, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance,
                  test_Type,interval_Type,implementation_Type,component_Type,exact_Inference,locs_Inference,locs_index_Inference,locs_are_nodes_Inference,coeff_Inference,beta_0,
                  f_0_eval,f_var_Inference,inference_Quantile,inference_Alpha,inference_N_Flip,inference_Tol_Fspai,inference_Defined,
                  PACKAGE = "fdaPDE")

  return(bigsol)
}

CPP_smooth.volume.FEM.PDE.sv.basis<-function(locations, observations, FEMbasis, covariates = NULL, PDE_parameters, ndim, mydim, BC = NULL, incidence_matrix = NULL, areal.data.avg = TRUE, search, bary.locations, optim, lambda = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1, lambda.optimization.tolerance = 0.05, inference.data.object)
{

  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(DOF.matrix))
  {
    DOF.matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }

  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  {
    BC$BC_indices<-as.vector(BC$BC_indices)-1
  }

  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }


  PDE_param_eval = NULL
  points_eval = matrix(CPP_get_evaluations_points(mesh = FEMbasis$mesh, order = FEMbasis$order),ncol = 2)
  PDE_param_eval$K = (PDE_parameters$K)(points_eval)
  PDE_param_eval$b = (PDE_parameters$b)(points_eval)
  PDE_param_eval$c = (PDE_parameters$c)(points_eval)
  PDE_param_eval$u = (PDE_parameters$u)(points_eval)
  
  if(inference.data.object@definition==1 & mean(PDE_param_eval$u != rep(0, nrow(points_eval)))!=0){
    warning("Inference for linear estimators is implemented only if reaction term is zero, \nInference Data are ignored")
    inference.data.object=new("inferenceDataObject", test = as.integer(0), interval = as.integer(0), type = as.integer(0), component = as.integer(0), exact = as.integer(0), dim = as.integer(0), n_cov = as.integer(0),
                                locations = matrix(data=0, nrow = 1 ,ncol = 1), locations_indices = as.integer(0), locations_are_nodes = as.integer(0), coeff = matrix(data=0, nrow = 1 ,ncol = 1), beta0 = -1, f0 = function(){}, 
                                f0_eval = -1, f_var = as.integer(0), quantile = -1, n_flip = as.integer(1000), tol_fspai = -1, definition=as.integer(0))
  }
  

  ## Extract the parameters for inference from inference.data.object to prepare them for c++ reading
  test_Type<-as.vector(inference.data.object@test)
  interval_Type<-as.vector(inference.data.object@interval)
  implementation_Type<-as.vector(inference.data.object@type)
  component_Type<-as.vector(inference.data.object@component)
  exact_Inference<-inference.data.object@exact
  locs_Inference<-as.matrix(inference.data.object@locations)
  locs_index_Inference<-as.vector(inference.data.object@locations_indices - 1) #converting indices from R to c++ ones
  locs_are_nodes_Inference <- inference.data.object@locations_are_nodes
  coeff_Inference<-as.matrix(inference.data.object@coeff)
  beta_0<-as.vector(inference.data.object@beta0)
  f_0_eval<-as.vector(inference.data.object@f0_eval)
  f_var_Inference<-inference.data.object@f_var
  inference_Quantile<-as.vector(inference.data.object@quantile)
  inference_Alpha<-as.vector(inference.data.object@alpha)
  inference_N_Flip<-inference.data.object@n_flip
  inference_Tol_Fspai<-inference.data.object@tol_fspai
  inference_Defined<-inference.data.object@definition
  
  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  data <- as.vector(observations)
  storage.mode(data) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(FEMbasis$mesh$faces) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"

  storage.mode(PDE_param_eval$K) <- "double"
  storage.mode(PDE_param_eval$b) <- "double"
  storage.mode(PDE_param_eval$c) <- "double"
  storage.mode(PDE_param_eval$u) <- "double"

  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"  
  storage.mode(lambda) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  
  ## Set proper type for correct C++ reading for inference parameters
  storage.mode(test_Type) <- "integer"
  storage.mode(interval_Type) <- "integer"
  storage.mode(implementation_Type) <- "integer"
  storage.mode(component_Type) <- "integer"
  storage.mode(exact_Inference) <- "integer"
  storage.mode(locs_Inference) <- "double"
  storage.mode(locs_index_Inference) <- "integer"
  storage.mode(locs_are_nodes_Inference) <- "integer"
  storage.mode(coeff_Inference) <- "double"
  storage.mode(beta_0) <- "double"
  storage.mode(f_0_eval) <- "double"
  storage.mode(f_var_Inference) <- "integer"
  storage.mode(inference_Quantile) <- "double"
  storage.mode(inference_Alpha) <- "double"
  storage.mode(inference_N_Flip) <- "integer"
  storage.mode(inference_Tol_Fspai) <- "double"
  storage.mode(inference_Defined) <- "integer"

  ## Call C++ function
  bigsol <- .Call("regression_PDE_space_varying", locations, bary.locations, data, FEMbasis$mesh, FEMbasis$order, mydim, ndim, PDE_param_eval$K, PDE_param_eval$b, PDE_param_eval$c, PDE_param_eval$u, covariates,
                  BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg, search,
                  optim, lambda, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance,
                  test_Type,interval_Type,implementation_Type,component_Type,exact_Inference,locs_Inference,locs_index_Inference,locs_are_nodes_Inference,coeff_Inference,beta_0,
                  f_0_eval,f_var_Inference,inference_Quantile,inference_Alpha,inference_N_Flip,inference_Tol_Fspai,inference_Defined,
                  PACKAGE = "fdaPDE")

  return(bigsol)
}

CPP_eval.volume.FEM = function(FEM, locations, incidence_matrix, redundancy, ndim, mydim, search, bary.locations)
{
  FEMbasis = FEM$FEMbasis

  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  # Imposing types, this is necessary for correct reading from C++
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(FEMbasis$mesh$faces) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  coeff <- as.matrix(FEM$coeff)
  storage.mode(coeff) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(locations) <- "double"
  storage.mode(redundancy) <- "integer"
  storage.mode(search) <- "integer"

  if(!is.null(bary.locations))
  {
    storage.mode(bary.locations$element_ids) <- "integer"
    element_ids <- as.matrix(bary.locations$element_ids)
    storage.mode(bary.locations$barycenters) <- "double"
    barycenters <- as.matrix(bary.locations$barycenters)
  }

  # if (search == 1) { #use Naive search
  #   print('This is Naive Search')
  # } else if (search == 2)  { #use Tree search (default)
  #   print('This is Tree Search')
  # }

  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  evalmat = matrix(0,max(nrow(locations),nrow(incidence_matrix)),ncol(coeff))
  for (i in 1:ncol(coeff)){
    evalmat[,i] <- .Call("eval_FEM_fd", FEMbasis$mesh, locations, incidence_matrix, coeff[,i],
                         FEMbasis$order, redundancy, mydim, ndim, search, bary.locations, PACKAGE = "fdaPDE")
  }

  #Returning the evaluation matrix
  evalmat
}

