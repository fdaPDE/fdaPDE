CPP_smooth.FEM.time<-function(locations, time_locations, observations, FEMbasis, time_mesh,
                              covariates = NULL, ndim, mydim, BC = NULL,
                              incidence_matrix = NULL, areal.data.avg = TRUE,
                              FLAG_MASS, FLAG_PARABOLIC, IC,
                              search, bary.locations, optim , lambdaS = NULL, lambdaT = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1, lambda.optimization.tolerance = 0.05, R_Inference_Data_Object)
{
  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
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
    locations<-matrix(nrow = 0, ncol = 2)
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
  R_Inference_Data_Object_Null=new("inferenceDataObject", test = as.integer(0), interval =as.integer(0), type = as.integer(0), exact = as.integer(0), dim = as.integer(0), 
                                   coeff = matrix(data=0, nrow = 1 ,ncol = 1), beta0 = -1, f_var = as.integer(0), quantile = -1, n_flip = as.integer(1000), tol_fspai = -1, definition=as.integer(0))
  
  ## Extract the parameters for inference from R_Inference_Data_Object to prepare them for c++ reading
  test_Type<-as.vector(R_Inference_Data_Object@test)
  interval_Type<-as.vector(R_Inference_Data_Object@interval)
  implementation_Type<-as.vector(R_Inference_Data_Object@type)
  exact_Inference<-R_Inference_Data_Object@exact
  coeff_Inference=as.matrix(R_Inference_Data_Object@coeff)
  beta_0=as.vector(R_Inference_Data_Object@beta0)
  f_var_Inference<-R_Inference_Data_Object@f_var
  inference_Quantile=as.vector(R_Inference_Data_Object@quantile)
  inference_n_flip=R_Inference_Data_Object@n_flip
  inference_tol_fspai=R_Inference_Data_Object@tol_fspai
  inference_Defined=R_Inference_Data_Object@definition
  
  ## Extract the parameters for preliminary computations from R_Inference_Data_Object_Null to prepare them for c++ reading
  test_Type_Null<-as.vector(R_Inference_Data_Object_Null@test)
  interval_Type_Null<-as.vector(R_Inference_Data_Object_Null@interval)
  implementation_Type_Null<-as.vector(R_Inference_Data_Object_Null@type)
  exact_Inference_Null<-R_Inference_Data_Object_Null@exact
  coeff_Inference_Null=as.matrix(R_Inference_Data_Object_Null@coeff)
  beta_0_Null=as.vector(R_Inference_Data_Object_Null@beta0)
  f_var_Inference_Null<-R_Inference_Data_Object_Null@f_var
  inference_Quantile_Null=as.vector(R_Inference_Data_Object_Null@quantile)
  inference_n_flip_Null=R_Inference_Data_Object_Null@n_flip
  inference_tol_fspai_Null=R_Inference_Data_Object_Null@tol_fspai
  inference_Defined_Null=R_Inference_Data_Object_Null@definition

  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
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
  storage.mode(exact_Inference) <- "integer"
  storage.mode(coeff_Inference) <- "double"
  storage.mode(beta_0) <- "double"
  storage.mode(f_var_Inference) <- "integer"
  storage.mode(inference_Quantile) <- "double"
  storage.mode(inference_n_flip) <- "integer"
  storage.mode(inference_tol_fspai) <- "double"
  storage.mode(inference_Defined) <- "integer"
  
  ## Set proper type for correct C++ reading for preliminary computations inference parameters
  storage.mode(test_Type_Null) <- "integer"
  storage.mode(interval_Type_Null) <- "integer"
  storage.mode(implementation_Type_Null) <- "integer"
  storage.mode(exact_Inference_Null) <- "integer"
  storage.mode(coeff_Inference_Null) <- "double"
  storage.mode(beta_0_Null) <- "double"
  storage.mode(f_var_Inference_Null) <- "integer"
  storage.mode(inference_Quantile_Null) <- "double"
  storage.mode(inference_n_flip_Null) <- "integer"
  storage.mode(inference_tol_fspai_Null) <- "double"
  storage.mode(inference_Defined_Null) <- "integer"

  ## IC estimation for parabolic smoothing from the first column of observations
  ICsol = NA
  #empty dof matrix
  DOF.matrix_IC<-matrix(nrow = 0, ncol = 1)
  if(nrow(IC)==0 && FLAG_PARABOLIC)
  {
    NobsIC = length(observations)%/%nrow(time_locations)
    notNAIC = which(!is.na(observations[1:NobsIC]))
    observationsIC = observations[notNAIC]

    if(nrow(locations)==0)
      locationsIC=locations
    else
    {
      locationsIC = as.matrix(locations[notNAIC,])
      storage.mode(locationsIC) <- "double"
    }

    if(nrow(covariates)==0)
      covariatesIC = covariates
    else
    {
      covariatesIC = covariates[notNAIC,]
      covariatesIC = as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }

    ## set of lambdas for GCV in IC estimation
    lambdaSIC <- 10^seq(-7,3,0.1)
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ## call the smoothing function with initial observations to estimates the IC

    ICsol <- .Call("regression_Laplace", locationsIC, bary.locations, observationsIC,
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC,
      BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg, search,
      as.integer(c(0,1,1)), lambdaSIC, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor, lambda.optimization.tolerance, 
      test_Type_Null,interval_Type_Null,implementation_Type_Null,exact_Inference_Null,coeff_Inference_Null,beta_0_Null,f_var_Inference_Null,inference_Quantile_Null,inference_n_flip_Null, inference_tol_fspai_Null, inference_Defined_Null,
      PACKAGE = "fdaPDE")

    ## shifting the lambdas interval if the best lambda is the smaller one and retry smoothing
    if(ICsol[[6]]==1)
    {
      lambdaSIC <- 10^seq(-9,-7,0.1)
      lambdaSIC <- as.matrix(lambdaSIC)
      storage.mode(lambdaSIC) <- "double"
      ICsol <- .Call("regression_Laplace", locationsIC, bary.locations, observationsIC,
       FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC,
       BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg, search,
       as.integer(c(0,1,1)), lambdaSIC, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor, lambda.optimization.tolerance,
       test_Type_Null,interval_Type_Null,implementation_Type_Null,exact_Inference_Null,coeff_Inference_Null,beta_0_Null,f_var_Inference_Null,inference_Quantile_Null,inference_n_flip_Null, inference_tol_fspai_Null, inference_Defined_Null,
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
        ICsol <- .Call("regression_Laplace", locationsIC, bary.locations, observationsIC,
         FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC,
         BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg, search,
         as.integer(c(0,1,1)), lambdaSIC, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor, lambda.optimization.tolerance, 
         test_Type_Null,interval_Type_Null,implementation_Type_Null,exact_Inference_Null,coeff_Inference_Null,beta_0_Null,f_var_Inference_Null,inference_Quantile_Null,inference_n_flip_Null, inference_tol_fspai_Null, inference_Defined_Null,
         PACKAGE = "fdaPDE")
      }
    }

    if(nrow(covariates)!=0)
    {
      betaIC = ICsol[[15]]
      IC = ICsol[[1]][1:nrow(FEMbasis$mesh$nodes),] ## best IC estimation
      covariates=covariates[(NobsIC+1):nrow(covariates)]
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


  ## Call C++ function
  bigsol <- .Call("regression_Laplace_time", locations, bary.locations, time_locations, observations, FEMbasis$mesh, time_mesh, FEMbasis$order,
    mydim, ndim, covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg, FLAG_MASS, FLAG_PARABOLIC,
    IC, search, optim, lambdaS, lambdaT, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance, 
    test_Type,interval_Type,implementation_Type,exact_Inference,coeff_Inference,beta_0,f_var_Inference,inference_Quantile,inference_n_flip,inference_tol_fspai, inference_Defined,
    PACKAGE = "fdaPDE")
  return(c(bigsol,ICsol))
}

CPP_smooth.FEM.PDE.time<-function(locations, time_locations, observations, FEMbasis, time_mesh,
                                  covariates = NULL, PDE_parameters, ndim, mydim, BC = NULL,
                                  incidence_matrix = NULL, areal.data.avg = TRUE,
                                  FLAG_MASS, FLAG_PARABOLIC, IC,
                                  search, bary.locations, optim , lambdaS = NULL, lambdaT = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1, lambda.optimization.tolerance = 0.05, R_Inference_Data_Object)
{

  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
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
    locations<-matrix(nrow = 0, ncol = 2)
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

  if(is.null(time_locations)){
    time_locations<-matrix(ncol=0, nrow=0)
  }else
  {
    time_locations <- as.matrix(time_locations)
  }

  if(is.null(time_mesh)){
    time_mesh<-matrix(ncol=0, nrow=0)
  }else
  {
    time_mesh <- as.matrix(time_mesh)
  }
  
  # Create a null inference object for preliminary computations 
  R_Inference_Data_Object_Null=new("inferenceDataObject", test = as.integer(0), interval =as.integer(0), type = as.integer(0), exact = as.integer(0), dim = as.integer(0), 
                                   coeff = matrix(data=0, nrow = 1 ,ncol = 1), beta0 = -1, f_var = as.integer(0), quantile = -1, n_flip = as.integer(1000), tol_fspai = -1, definition=as.integer(0))
  
  ## Extract the parameters for inference from R_Inference_Data_Object to prepare them for c++ reading
  test_Type<-as.vector(R_Inference_Data_Object@test)
  interval_Type<-as.vector(R_Inference_Data_Object@interval)
  implementation_Type<-as.vector(R_Inference_Data_Object@type)
  exact_Inference<-R_Inference_Data_Object@exact
  coeff_Inference=as.matrix(R_Inference_Data_Object@coeff)
  beta_0=as.vector(R_Inference_Data_Object@beta0)
  f_var_Inference<-R_Inference_Data_Object@f_var
  inference_Quantile=as.vector(R_Inference_Data_Object@quantile)
  inference_n_flip=R_Inference_Data_Object@n_flip
  inference_tol_fspai=R_Inference_Data_Object@tol_fspai
  inference_Defined=R_Inference_Data_Object@definition
  
  ## Extract the parameters for preliminary computations from R_Inference_Data_Object_Null to prepare them for c++ reading
  test_Type_Null<-as.vector(R_Inference_Data_Object_Null@test)
  interval_Type_Null<-as.vector(R_Inference_Data_Object_Null@interval)
  implementation_Type_Null<-as.vector(R_Inference_Data_Object_Null@type)
  exact_Inference_Null<-R_Inference_Data_Object_Null@exact
  coeff_Inference_Null=as.matrix(R_Inference_Data_Object_Null@coeff)
  beta_0_Null=as.vector(R_Inference_Data_Object_Null@beta0)
  f_var_Inference_Null<-R_Inference_Data_Object_Null@f_var
  inference_Quantile_Null=as.vector(R_Inference_Data_Object_Null@quantile)
  inference_n_flip_Null=R_Inference_Data_Object_Null@n_flip
  inference_tol_fspai_Null=R_Inference_Data_Object_Null@tol_fspai
  inference_Defined_Null=R_Inference_Data_Object_Null@definition

  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(time_locations) <- "double"
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(PDE_parameters$K) <- "double"
  storage.mode(PDE_parameters$b) <- "double"
  storage.mode(PDE_parameters$c) <- "double"
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
  storage.mode(exact_Inference) <- "integer"
  storage.mode(coeff_Inference) <- "double"
  storage.mode(beta_0) <- "double"
  storage.mode(f_var_Inference) <- "integer"
  storage.mode(inference_Quantile) <- "double"
  storage.mode(inference_n_flip) <- "integer"
  storage.mode(inference_tol_fspai) <- "double"
  storage.mode(inference_Defined) <- "integer"
  
  ## Set proper type for correct C++ reading for preliminary computations inference parameters
  storage.mode(test_Type_Null) <- "integer"
  storage.mode(interval_Type_Null) <- "integer"
  storage.mode(implementation_Type_Null) <- "integer"
  storage.mode(exact_Inference_Null) <- "integer"
  storage.mode(coeff_Inference_Null) <- "double"
  storage.mode(beta_0_Null) <- "double"
  storage.mode(f_var_Inference_Null) <- "integer"
  storage.mode(inference_Quantile_Null) <- "double"
  storage.mode(inference_n_flip_Null) <- "integer"
  storage.mode(inference_tol_fspai_Null) <- "double"
  storage.mode(inference_Defined_Null) <- "integer"

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

    lambdaSIC <- 10^seq(-7,3,0.1)
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ICsol <- .Call("regression_PDE", locations, bary.locations, observations[1:NobsIC], FEMbasis$mesh, FEMbasis$order,
      mydim, ndim, PDE_parameters$K, PDE_parameters$b, PDE_parameters$c, covariatesIC, BC$BC_indices, BC$BC_values,
      incidence_matrix, areal.data.avg, search, as.integer(c(0,1,1)), lambdaSIC, DOF.stochastic.realizations, DOF.stochastic.seed,  DOF.matrix_IC, GCV.inflation.factor, lambda.optimization.tolerance, 
      test_Type_Null,interval_Type_Null,implementation_Type_Null,exact_Inference_Null,coeff_Inference_Null,beta_0_Null,f_var_Inference_Null,inference_Quantile_Null,inference_n_flip_Null, inference_tol_fspai_Null, inference_Defined_Null,
      PACKAGE = "fdaPDE")

    if(ICsol[[6]]==1)
    {
      lambdaSIC <- 10^seq(-9,-7,0.1)
      lambdaSIC <- as.matrix(lambdaSIC)
      storage.mode(lambdaSIC) <- "double"
      ICsol <- .Call("regression_PDE", locations, bary.locations, observations[1:NobsIC], FEMbasis$mesh, FEMbasis$order,
       mydim, ndim, PDE_parameters$K, PDE_parameters$b, PDE_parameters$c, covariatesIC, BC$BC_indices, BC$BC_values,
       incidence_matrix, areal.data.avg, search, as.integer(c(0,1,1)), lambdaSIC, DOF.stochastic.realizations, DOF.stochastic.seed,  DOF.matrix_IC, GCV.inflation.factor, lambda.optimization.tolerance, 
       test_Type_Null,interval_Type_Null,implementation_Type_Null,exact_Inference_Null,coeff_Inference_Null,beta_0_Null,f_var_Inference_Null,inference_Quantile_Null,inference_n_flip_Null, inference_tol_fspai_Null, inference_Defined_Null,
       PACKAGE = "fdaPDE")
    }
    else
    {
      if(ICsol[[6]]==length(lambdaSIC))
      {
        lambdaSIC <- 10^seq(3,5,0.1)
        lambdaSIC <- as.matrix(lambdaSIC)
        storage.mode(lambdaSIC) <- "double"
        ICsol <- .Call("regression_PDE", locations, bary.locations, observations[1:NobsIC], FEMbasis$mesh, FEMbasis$order,
         mydim, ndim, PDE_parameters$K, PDE_parameters$b, PDE_parameters$c, covariatesIC, BC$BC_indices, BC$BC_values,
         incidence_matrix, areal.data.avg, search, as.integer(c(0,1,1)), lambdaSIC, DOF.stochastic.realizations, DOF.stochastic.seed,  DOF.matrix_IC, GCV.inflation.factor, lambda.optimization.tolerance, 
         test_Type_Null,interval_Type_Null,implementation_Type_Null,exact_Inference_Null,coeff_Inference_Null,beta_0_Null,f_var_Inference_Null,inference_Quantile_Null,inference_n_flip_Null, inference_tol_fspai_Null, inference_Defined_Null,
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

  ## Call C++ function
  bigsol <- .Call("regression_PDE_time", locations, bary.locations, time_locations, observations, FEMbasis$mesh, time_mesh, FEMbasis$order,
                  mydim, ndim, PDE_parameters$K, PDE_parameters$b, PDE_parameters$c, covariates,
                  BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg, FLAG_MASS, FLAG_PARABOLIC,
                  IC, search, optim, lambdaS, lambdaT, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance,
                  test_Type,interval_Type,implementation_Type,exact_Inference,coeff_Inference,beta_0,f_var_Inference,inference_Quantile,inference_n_flip,inference_tol_fspai, inference_Defined,
                  PACKAGE = "fdaPDE")
  return(c(bigsol,ICsol))
}

CPP_smooth.FEM.PDE.sv.time<-function(locations, time_locations, observations, FEMbasis, time_mesh,
                                     covariates = NULL, PDE_parameters, ndim, mydim, BC = NULL,
                                     incidence_matrix = NULL, areal.data.avg = TRUE,
                                     FLAG_MASS, FLAG_PARABOLIC, IC,
                                     search, bary.locations, optim , lambdaS = NULL, lambdaT = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1, lambda.optimization.tolerance = 0.05, R_Inference_Data_Object)
{

  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
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
    locations<-matrix(nrow = 0, ncol = 2)
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

  if(is.null(time_locations)){
    time_locations<-matrix(ncol=0, nrow=0)
  }else
  {
    time_locations <- as.matrix(time_locations)
  }

  if(is.null(time_mesh)){
    time_mesh<-matrix(ncol=0, nrow=0)
  }else
  {
    time_mesh <- as.matrix(time_mesh)
  }


  PDE_param_eval = NULL
  points_eval = matrix(CPP_get_evaluations_points(mesh = FEMbasis$mesh, order = FEMbasis$order),ncol = 2)
  PDE_param_eval$K = (PDE_parameters$K)(points_eval)
  PDE_param_eval$b = (PDE_parameters$b)(points_eval)
  PDE_param_eval$c = (PDE_parameters$c)(points_eval)
  PDE_param_eval$u = (PDE_parameters$u)(points_eval)
  
  if(R_Inference_Data_Object@definition==1 && mean(PDE_param_eval$u != rep(0, nrow(points_eval)))!=0){
    warning("Inference for linear estimators is implemented only if reaction term is zero, \nInference Data are ignored")
    R_Inference_Data_Object=new("inferenceDataObject", test = as.integer(0), interval =as.integer(0), type = as.integer(0), exact = as.integer(0), dim = as.integer(0), 
                                coeff = matrix(data=0, nrow = 1 ,ncol = 1), beta0 = -1, f_var = as.integer(0), quantile = -1, n_flip = as.integer(1000), tol_fspai = -1, definition=as.integer(0))
  }
  
  # Create a null inference object for preliminary computations 
  R_Inference_Data_Object_Null=new("inferenceDataObject", test = as.integer(0), interval =as.integer(0), type = as.integer(0), exact = as.integer(0), dim = as.integer(0), 
                                   coeff = matrix(data=0, nrow = 1 ,ncol = 1), beta0 = -1, f_var = as.integer(0), quantile = -1, n_flip = as.integer(1000), tol_fspai = -1, definition=as.integer(0))
  
  ## Extract the parameters for inference from R_Inference_Data_Object to prepare them for c++ reading
  test_Type<-as.vector(R_Inference_Data_Object@test)
  interval_Type<-as.vector(R_Inference_Data_Object@interval)
  implementation_Type<-as.vector(R_Inference_Data_Object@type)
  exact_Inference<-R_Inference_Data_Object@exact
  coeff_Inference=as.matrix(R_Inference_Data_Object@coeff)
  beta_0=as.vector(R_Inference_Data_Object@beta0)
  f_var_Inference<-R_Inference_Data_Object@f_var
  inference_Quantile=as.vector(R_Inference_Data_Object@quantile)
  inference_n_flip=R_Inference_Data_Object@n_flip
  inference_tol_fspai=R_Inference_Data_Object@tol_fspai
  inference_Defined=R_Inference_Data_Object@definition
  
  ## Extract the parameters for preliminary computations from R_Inference_Data_Object_Null to prepare them for c++ reading
  test_Type_Null<-as.vector(R_Inference_Data_Object_Null@test)
  interval_Type_Null<-as.vector(R_Inference_Data_Object_Null@interval)
  implementation_Type_Null<-as.vector(R_Inference_Data_Object_Null@type)
  exact_Inference_Null<-R_Inference_Data_Object_Null@exact
  coeff_Inference_Null=as.matrix(R_Inference_Data_Object_Null@coeff)
  beta_0_Null=as.vector(R_Inference_Data_Object_Null@beta0)
  f_var_Inference_Null<-R_Inference_Data_Object_Null@f_var
  inference_Quantile_Null=as.vector(R_Inference_Data_Object_Null@quantile)
  inference_n_flip_Null=R_Inference_Data_Object_Null@n_flip
  inference_tol_fspai_Null=R_Inference_Data_Object_Null@tol_fspai
  inference_Defined_Null=R_Inference_Data_Object_Null@definition

  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(time_locations) <- "double"
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(PDE_param_eval$K) <- "double"
  storage.mode(PDE_param_eval$b) <- "double"
  storage.mode(PDE_param_eval$c) <- "double"
  storage.mode(PDE_param_eval$u) <- "double"
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
  storage.mode(exact_Inference) <- "integer"
  storage.mode(coeff_Inference) <- "double"
  storage.mode(beta_0) <- "double"
  storage.mode(f_var_Inference) <- "integer"
  storage.mode(inference_Quantile) <- "double"
  storage.mode(inference_n_flip) <- "integer"
  storage.mode(inference_tol_fspai) <- "double"
  storage.mode(inference_Defined) <- "integer"
  
  ## Set proper type for correct C++ reading for preliminary computations inference parameters
  storage.mode(test_Type_Null) <- "integer"
  storage.mode(interval_Type_Null) <- "integer"
  storage.mode(implementation_Type_Null) <- "integer"
  storage.mode(exact_Inference_Null) <- "integer"
  storage.mode(coeff_Inference_Null) <- "double"
  storage.mode(beta_0_Null) <- "double"
  storage.mode(f_var_Inference_Null) <- "integer"
  storage.mode(inference_Quantile_Null) <- "double"
  storage.mode(inference_n_flip_Null) <- "integer"
  storage.mode(inference_tol_fspai_Null) <- "double"
  storage.mode(inference_Defined_Null) <- "integer"

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

    lambdaSIC <- 10^seq(-7,3,0.1)
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"

    ICsol <- .Call("regression_PDE_space_varying", locations, bary.locations, observations[1:NobsIC],
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, PDE_param_eval$K, PDE_param_eval$b, PDE_param_eval$c, PDE_param_eval$u,
      covariatesIC, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
      search, as.integer(c(0,1,1)), lambdaSIC, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor, lambda.optimization.tolerance, 
      test_Type_Null,interval_Type_Null,implementation_Type_Null,exact_Inference_Null,coeff_Inference_Null,beta_0_Null,f_var_Inference_Null,inference_Quantile_Null,inference_n_flip_Null, inference_tol_fspai_Null, inference_Defined_Null,
      PACKAGE = "fdaPDE")

    if(ICsol[[6]]==1)
    {
      lambdaSIC <- 10^seq(-9,-7,0.1)
      lambdaSIC <- as.matrix(lambdaSIC)
      storage.mode(lambdaSIC) <- "double"
      ICsol <- .Call("regression_PDE_space_varying", locations, bary.locations, observations[1:NobsIC],
       FEMbasis$mesh, FEMbasis$order, mydim, ndim, PDE_param_eval$K, PDE_param_eval$b, PDE_param_eval$c, PDE_param_eval$u,
       covariatesIC, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
       search, as.integer(c(0,1,1)), lambdaSIC, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor, lambda.optimization.tolerance, 
       test_Type_Null,interval_Type_Null,implementation_Type_Null,exact_Inference_Null,coeff_Inference_Null,beta_0_Null,f_var_Inference_Null,inference_Quantile_Null,inference_n_flip_Null, inference_tol_fspai_Null, inference_Defined_Null,
       PACKAGE = "fdaPDE")
    }
    else
    {
      if(ICsol[[6]]==length(lambdaSIC))
      {
        lambdaSIC <- 10^seq(3,5,0.1)
        lambdaSIC <- as.matrix(lambdaSIC)
        storage.mode(lambdaSIC) <- "double"
        ICsol <- .Call("regression_PDE_space_varying", locations, bary.locations, observations[1:NobsIC],
         FEMbasis$mesh, FEMbasis$order, mydim, ndim, PDE_param_eval$K, PDE_param_eval$b, PDE_param_eval$c, PDE_param_eval$u,
         covariatesIC, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
         search, as.integer(c(0,1,1)), lambdaSIC, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor, lambda.optimization.tolerance, 
         test_Type_Null, interval_Type_Null, implementation_Type_Null, exact_Inference_Null, coeff_Inference_Null, beta_0_Null, f_var_Inference_Null, inference_Quantile_Null, inference_n_flip_Null, inference_tol_fspai_Null, inference_Defined_Null,
         PACKAGE = "fdaPDE")
      }
    }

    if(nrow(covariates)!=0)
    {
      betaIC = ICsol[[15]]
      IC = ICsol[[1]][1:nrow(FEMbasis$mesh$nodes),]# best IC estimation
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


  ## Call C++ function
  bigsol <- .Call("regression_PDE_space_varying_time", locations, bary.locations, time_locations, observations, FEMbasis$mesh, time_mesh, FEMbasis$order,
    mydim, ndim,PDE_param_eval$K, PDE_param_eval$b, PDE_param_eval$c, PDE_param_eval$u, covariates,
    BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg, FLAG_MASS, FLAG_PARABOLIC,
    IC, search,  optim, lambdaS, lambdaT, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance,
    test_Type, interval_Type, implementation_Type, exact_Inference, coeff_Inference, beta_0, f_var_Inference, inference_Quantile, inference_n_flip, inference_tol_fspai, inference_Defined,
    PACKAGE = "fdaPDE")
  return(c(bigsol,ICsol))
}

CPP_eval.FEM.time <- function(FEM.time, locations, time_locations, incidence_matrix, FLAG_PARABOLIC, redundancy, ndim, mydim, search, bary.locations)
{

  # EVAL_FEM_FD evaluates the FEM fd object at points (X,Y)
  #
  #        arguments:
  # X         an array of x-coordinates.
  # Y         an array of y-coordinates.
  # FELSPLOBJ a FELspline object
  # FAST      a boolean indicating if the walking algorithm should be apply
  #        output:
  # EVALMAT   an array of the same size as X and Y containing the value of
  #           FELSPLOBJ at (X,Y).

  FEMbasis = FEM.time$FEMbasis
  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1


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
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
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
  }

  # if (search == 1) { #use Naive search
  #   print('This is Naive Search')
  # } else if (search == 2)  { #use Tree search (default)
  #   print('This is Tree Search')
  # } else if (search == 3) { #use Walking search
  #     print('This is Walking Search')
  # }

  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  evalmat = matrix(0,max(nrow(locations),nrow(incidence_matrix)),ncol(coeff))
  for (i in 1:ncol(coeff)){
    evalmat[,i] <- .Call("eval_FEM_time", FEMbasis$mesh, FEM.time$mesh_time, locations, time_locations, incidence_matrix, coeff[,i],
                         FEMbasis$order, redundancy, FLAG_PARABOLIC, mydim, ndim, search, bary.locations, PACKAGE = "fdaPDE")
  }

  #Returning the evaluation matrix
  evalmat
}
