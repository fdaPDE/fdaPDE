CPP_smooth.graph.FEM.time<-function(locations, time_locations, observations, FEMbasis, time_mesh,
                                       covariates = NULL, ndim, mydim, BC = NULL,
                                       incidence_matrix = NULL, areal.data.avg = TRUE,
                                       FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, threshold = 10^(-4), max.steps = 50, IC,
                                       search, bary.locations, optim , lambdaS = NULL, lambdaT = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1, lambda.optimization.tolerance = 0.05, inference.data.object)
{
  
  # Indexes in C++ starts from 0, in R from 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  
  num_sides = 2*dim(FEMbasis$mesh$edges)[1] 
  for(i in 1:num_sides){
    if( dim(FEMbasis$mesh$neighbors[[i]] )[1] > 0)
      FEMbasis$mesh$neighbors[[i]] = FEMbasis$mesh$neighbors[[i]] - 1
  }
  
  max.steps = max.steps - 1
  
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
  
  if(is.null(IC))
  {
    IC<-matrix(nrow = 0, ncol = 1)
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
  if(is.null(lambdaS))
  {
    lambdaS<-vector(length=0)
  }else
  {
    lambdaS<-as.vector(lambdaS)
  }
  
  if(is.null(lambdaT))
  {
    lambdaT<-vector(length=0)
  }else
  {
    lambdaT<-as.vector(lambdaT)
  }
  
  # Create a null inference object for preliminary computations 
  inference.data.object.null=new("inferenceDataObject", test = as.integer(0), interval =as.integer(0), type = as.integer(0), component = as.integer(0), exact = as.integer(0), dim = as.integer(0), n_cov = as.integer(0), 
                                 locations = matrix(data=0, nrow = 1 ,ncol = 1), locations_indices = as.integer(0), locations_are_nodes = as.integer(0), coeff = matrix(data=0, nrow = 1 ,ncol = 1), beta0 = -1, f0 = function(){}, 
                                 f0_eval = -1, f_var = as.integer(0), quantile = -1, alpha = 0, n_flip = as.integer(1000), tol_fspai = -1, definition=as.integer(0))
  
  ## Extract the parameters for inference from inference.data.object to prepare them for c++ reading
  test_Type<-as.vector(inference.data.object@test)
  interval_Type<-as.vector(inference.data.object@interval)
  implementation_Type<-as.vector(inference.data.object@type)
  component_Type<-as.vector(inference.data.object@component)
  exact_Inference<-inference.data.object@exact
  coeff_Inference=as.matrix(inference.data.object@coeff)
  beta_0=as.vector(inference.data.object@beta0)
  f_var_Inference<-inference.data.object@f_var
  inference_Quantile=as.vector(inference.data.object@quantile)
  inference_Alpha=inference.data.object@alpha
  inference_N_Flip=inference.data.object@n_flip
  inference_Tol_Fspai=inference.data.object@tol_fspai
  inference_Defined=inference.data.object@definition
  
  ## Extract the parameters for preliminary computations from inference.data.object.null to prepare them for c++ reading
  test_Type_Null<-as.vector(inference.data.object.null@test)
  interval_Type_Null<-as.vector(inference.data.object.null@interval)
  implementation_Type_Null<-as.vector(inference.data.object.null@type)
  component_Type_Null<-as.vector(inference.data.object.null@component)
  exact_Inference_Null<-inference.data.object.null@exact
  locs_Inference_Null<-as.matrix(inference.data.object.null@locations)
  locs_index_Inference_Null<-as.vector(inference.data.object.null@locations_indices - 1) #converting the indices from R to c++ ones 
  locs_are_nodes_Inference_Null<-inference.data.object.null@locations_are_nodes
  coeff_Inference_Null=as.matrix(inference.data.object.null@coeff)
  beta_0_Null=as.vector(inference.data.object.null@beta0)
  f_0_eval_Null<-as.vector(inference.data.object.null@f0_eval)
  f_var_Inference_Null<-inference.data.object.null@f_var
  inference_Quantile_Null=as.vector(inference.data.object.null@quantile)
  inference_Alpha_Null=inference.data.object.null@alpha
  inference_N_Flip_Null=inference.data.object.null@n_flip
  inference_Tol_Fspai_Null=inference.data.object.null@tol_fspai
  inference_Defined_Null=inference.data.object.null@definition
  
  
  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  observations <- as.vector(observations)
  storage.mode(observations) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
 
  for(i in 1:num_sides)
    storage.mode(FEMbasis$mesh$neighbors[[i]]) <- "integer" 
  
  storage.mode(FEMbasis$order) <- "integer"
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
  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <-"integer"
  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <-"integer"
  FLAG_ITERATIVE <- as.integer(FLAG_ITERATIVE)
  storage.mode(FLAG_ITERATIVE) <-"integer"
  storage.mode(max.steps) <- "integer"
  storage.mode(threshold) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
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
  storage.mode(coeff_Inference) <- "double"
  storage.mode(beta_0) <- "double"
  storage.mode(f_var_Inference) <- "integer"
  storage.mode(inference_Quantile) <- "double"
  storage.mode(inference_Alpha) <- "double"
  storage.mode(inference_N_Flip) <- "integer"
  storage.mode(inference_Tol_Fspai) <- "double"
  storage.mode(inference_Defined) <- "integer"
  
  ## Set proper type for correct C++ reading for preliminary computations inference parameters
  storage.mode(test_Type_Null) <- "integer"
  storage.mode(interval_Type_Null) <- "integer"
  storage.mode(implementation_Type_Null) <- "integer"
  storage.mode(component_Type_Null) <- "integer"
  storage.mode(exact_Inference_Null) <- "integer"
  storage.mode(locs_Inference_Null) <- "double"
  storage.mode(locs_index_Inference_Null) <- "integer"
  storage.mode(locs_are_nodes_Inference_Null) <- "integer"
  storage.mode(coeff_Inference_Null) <- "double"
  storage.mode(beta_0_Null) <- "double"
  storage.mode(f_0_eval_Null) <- "double"
  storage.mode(f_var_Inference_Null) <- "integer"
  storage.mode(inference_Quantile_Null) <- "double"
  storage.mode(inference_Alpha_Null) <- "double"
  storage.mode(inference_N_Flip_Null) <- "integer"
  storage.mode(inference_Tol_Fspai_Null) <- "double"
  storage.mode(inference_Defined_Null) <- "integer"
  
  ## Call C++ function
  ICsol=NA
  #empty dof matrix
  DOF.matrix_IC<-matrix(nrow = 0, ncol = 1)
  if(nrow(IC)==0 && FLAG_PARABOLIC)
  {
    NobsIC = length(observations)%/%nrow(time_locations)
    
    if(nrow(covariates)==0)
      covariatesIC = covariates
    else
    {
      covariatesIC = covariates[1:NobsIC,]
      covariatesIC = as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }
    
    ## set of lambdas for GCV in IC estimation
    lambdaSIC <- 10^seq(-7,3,0.1)
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ## call the smoothing function with initial observations to estimates the IC
    ICsol <- .Call("regression_Laplace", locations, bary.locations, observations[1:NobsIC],
                   FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC,
                   BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
                   search, as.integer(c(0,2,1)), lambdaSIC, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor, lambda.optimization.tolerance, 
                   test_Type_Null,interval_Type_Null,implementation_Type_Null,component_Type_Null,exact_Inference_Null,locs_Inference_Null,locs_index_Inference_Null,locs_are_nodes_Inference_Null,coeff_Inference_Null,
                   beta_0_Null,f_0_eval_Null,f_var_Inference_Null,inference_Quantile_Null,inference_Alpha_Null,inference_N_Flip_Null, inference_Tol_Fspai_Null, inference_Defined_Null,
                   PACKAGE = "fdaPDE")
    
    ## shifting the lambdas interval if the best lambda is the smaller one and retry smoothing
    if(ICsol[[6]]==1)
    {
      lambdaSIC <- 10^seq(-9,-7,0.1)
      lambdaSIC <- as.matrix(lambdaSIC)
      storage.mode(lambdaSIC) <- "double"
      ICsol <- .Call("regression_Laplace", locations, bary.locations, observations[1:NobsIC],
                     FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC,
                     BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
                     search, as.integer(c(0,2,1)), lambdaSIC, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor, lambda.optimization.tolerance, 
                     test_Type_Null,interval_Type_Null,implementation_Type_Null,component_Type_Null,exact_Inference_Null,locs_Inference_Null,locs_index_Inference_Null,locs_are_nodes_Inference_Null,coeff_Inference_Null,
                     beta_0_Null,f_0_eval_Null,f_var_Inference_Null,inference_Quantile_Null,inference_Alpha_Null,inference_N_Flip_Null, inference_Tol_Fspai_Null, inference_Defined_Null,
                     PACKAGE = "fdaPDE")
    }
    else
    {
      ## shifting the lambdas interval if the best lambda is the higher one and retry smoothing
      if(ICsol[[6]]==length(lambdaSIC))
      {
        lambdaSIC <- 10^seq(3,5,0.1)
        lambdaSIC <- as.matrix(lambdaSIC)
        storage.mode(lambdaSIC) <- "double"
        ICsol <- .Call("regression_Laplace", locations, bary.locations, observations[1:NobsIC],
                       FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC,
                       BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
                       search, as.integer(c(0,2,1)), lambdaSIC, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor,lambda.optimization.tolerance, 
                       test_Type_Null,interval_Type_Null,implementation_Type_Null,component_Type_Null,exact_Inference_Null,locs_Inference_Null,locs_index_Inference_Null,locs_are_nodes_Inference_Null,coeff_Inference_Null,
                       beta_0_Null,f_0_eval_Null,f_var_Inference_Null,inference_Quantile_Null,inference_Alpha_Null,inference_N_Flip_Null, inference_Tol_Fspai_Null, inference_Defined_Null,
                       PACKAGE = "fdaPDE")
      }
    }
    
    if(nrow(covariates)!=0)
    {
      betaIC = ICsol[[15]]
      IC = ICsol[[1]][1:nrow(FEMbasis$mesh$nodes),] ## best IC estimation
      covariates=covariates[(NobsIC+1):nrow(covariates),]
      covariates <- as.matrix(covariates)
    }
    else
    {
      IC = ICsol[[1]][1:nrow(FEMbasis$mesh$nodes),] ## best IC estimation
      betaIC = NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol = list(IC.FEM=FEM(ICsol[[1]][1:nrow(FEMbasis$mesh$nodes),],FEMbasis),bestlambdaindex=ICsol[[6]],bestlambda=ICsol[[5]],beta=betaIC)
    time_locations=time_locations[2:nrow(time_locations)]
    observations = observations[(NobsIC+1):length(observations)]
  }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  
  M = ifelse(FLAG_PARABOLIC,length(time_mesh)-1,length(time_mesh) + 2);
  BC$BC_indices = rep((0:(M-1))*nrow(FEMbasis$mesh$nodes),each=length(BC$BC_indices)) + rep(BC$BC_indices,M)
  BC$BC_values = rep(BC$BC_values,M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  
  bigsol <- .Call("regression_Laplace_time", locations, bary.locations, time_locations, observations, FEMbasis$mesh, time_mesh, FEMbasis$order,
                  mydim, ndim, covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold,
                  IC, search, optim, lambdaS, lambdaT, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance, 
                  test_Type,interval_Type,implementation_Type,component_Type,exact_Inference,coeff_Inference,beta_0,f_var_Inference,inference_Quantile,
                  inference_Alpha,inference_N_Flip,inference_Tol_Fspai, inference_Defined,
                  PACKAGE = "fdaPDE")
  
  return(c(bigsol,ICsol))
}

CPP_eval.graph.FEM.time = function(FEM.time, locations, time_locations, incidence_matrix, FLAG_PARABOLIC, redundancy, ndim, mydim, search, bary.locations)
{
  FEMbasis = FEM.time$FEMbasis
  # C++ function for manifold works with vectors not with matrices
  
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  
  num_sides = 2*dim(FEMbasis$mesh$edges)[1] 
  for(i in 1:num_sides){
    if( dim(FEMbasis$mesh$neighbors[[i]] )[1] > 0)
      FEMbasis$mesh$neighbors[[i]] = FEMbasis$mesh$neighbors[[i]] - 1
  }
  
  if(is.null(time_locations)){
    time_locations<-matrix(ncol=0, nrow=0)
  }else
  {
    time_locations <- as.matrix(time_locations)
  }
  
  # Imposing types, this is necessary for correct reading from C++
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(time_locations) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
 
  for(i in 1:num_sides)
    storage.mode(FEMbasis$mesh$neighbors[[i]]) <- "integer" 
  
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(FEM.time$mesh_time) <- "double"
  coeff <- as.matrix(FEM.time$coeff)
  storage.mode(coeff) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(locations) <- "double"
  storage.mode(redundancy) <- "integer"
  storage.mode(FLAG_PARABOLIC) <- "integer"
  
  storage.mode(search) <- "integer"
  
  if(!is.null(bary.locations))
  {
    storage.mode(bary.locations$element_ids) <- "integer"
    element_ids <- as.matrix(bary.locations$element_ids)
    storage.mode(bary.locations$barycenters) <- "double"
    barycenters <- as.matrix(bary.locations$barycenters)
  }else{
    bary.locations = list(locations=matrix(nrow=0,ncol=ndim), element_ids=matrix(nrow=0,ncol=1), barycenters=matrix(nrow=0,ncol=2))
    storage.mode(bary.locations$locations) <- "double"
    storage.mode(bary.locations$element_ids) <- "integer"
    storage.mode(bary.locations$barycenters) <- "double"
  }
  
  #Calling the C++ function "eval_FEM_time" in FEM_Eval.cpp
  evalmat = matrix(0,max(nrow(locations),nrow(incidence_matrix)),ncol(coeff))
  for (i in 1:ncol(coeff)){
    evalmat[,i] <- .Call("eval_FEM_time", FEMbasis$mesh, FEM.time$mesh_time, locations, time_locations, incidence_matrix, as.matrix(coeff[,i]),
                         FEMbasis$order, redundancy, FLAG_PARABOLIC, mydim, ndim, search, bary.locations, PACKAGE = "fdaPDE")
  }
  
  #Returning the evaluation matrix
  return(evalmat)
}
