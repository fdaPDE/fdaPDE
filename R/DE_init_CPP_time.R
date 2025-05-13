CPP_FEM.DE_init_time <- function(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, scaling, fvec, heatStep, heatIter, ndim,
                                 mydim, step_method, direction_method, preprocess_method, stepProposals, tol1, tol2, print,
                                 nfolds, nsimulations, search, isTimeDiscrete, flagMass, flagLumped, init, nFolds, inference)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  if(is.null(stepProposals))
    stepProposals = c(0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 1e-7, 1e-8, 1e-9)

  if(is.null(preprocess_method))
    preprocess_method = ""
    
  ## Set proper type for correct C++ reading
  storage.mode(data) <- "double"
  storage.mode(data_time) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(mesh_time) <- "double"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(lambda_time) <- "double"
  storage.mode(scaling) <- "double"
  storage.mode(fvec) <- "double"
  storage.mode(heatStep) <- "double"
  heatIter <- as.integer(heatIter)
  storage.mode(heatIter) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(stepProposals) <- "double"
  storage.mode(tol1) <- "double"
  storage.mode(tol2) <- "double"
  storage.mode(print) <- "logical"
  nfolds <- as.integer(nfolds)
  storage.mode(nfolds) <- "integer"
  nsimulations <- as.integer(nsimulations)
  storage.mode(nsimulations) <- "integer"
  storage.mode(search) <- "integer"
  storage.mode(isTimeDiscrete) <- "logical"
  storage.mode(flagMass) <- "logical"
  storage.mode(flagLumped) <- "logical"
  init <- as.character(init)
  storage.mode(init) <- "character"
  nFolds <- as.integer(nFolds)
  storage.mode(nFolds) <- "integer"
  storage.mode(inference) <- "logical"


  ## Call C++ function
  bigsol <- .Call("Density_Initialization_time", data, data_time, FEMbasis$mesh, mesh_time, FEMbasis$order, mydim, ndim,
                  scaling, fvec, heatStep, heatIter, lambda, lambda_time, nfolds, nsimulations, stepProposals, tol1, tol2, print,
                  search, isTimeDiscrete, flagMass, flagLumped, init, nFolds, inference, PACKAGE = "fdaPDE")

  return(bigsol)
}


CPP_FEM.manifold.DE_init_time <- function(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, scaling, fvec, heatStep, heatIter, ndim,
                                     mydim, step_method, direction_method, preprocess_method, stepProposals, tol1, tol2, print,
                                     nfolds, nsimulations, search, isTimeDiscrete, flagMass, flagLumped, init, nFolds, inference)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  ## Set proper type for correct C++ reading
  storage.mode(data) <- "double"
  storage.mode(data_time) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(mesh_time) <- "double"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(lambda_time) <- "double"
  storage.mode(scaling) <- "double"
  storage.mode(fvec) <- "double"
  storage.mode(heatStep) <- "double"
  heatIter <- as.integer(heatIter)
  storage.mode(heatIter) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(stepProposals) <- "double"
  storage.mode(tol1) <- "double"
  storage.mode(tol2) <- "double"
  storage.mode(print) <- "logical"
  nfolds <- as.integer(nfolds)
  storage.mode(nfolds) <- "integer"
  nsimulations <- as.integer(nsimulations)
  storage.mode(nsimulations) <- "integer"
  storage.mode(search) <- "integer"
  storage.mode(isTimeDiscrete) <- "logical"
  storage.mode(flagMass) <- "logical"
  storage.mode(flagLumped) <- "logical"
  init <- as.character(init)
  storage.mode(init) <- "character"
  nFolds <- as.integer(nFolds)
  storage.mode(nFolds) <- "integer"
  storage.mode(inference) <- "logical"


  ## Call C++ function
  bigsol <- .Call("Density_Initialization_time", data, data_time, FEMbasis$mesh, mesh_time, FEMbasis$order, mydim, ndim,
                  scaling, fvec, heatStep, heatIter, lambda, lambda_time, nfolds, nsimulations, stepProposals, tol1, tol2, print,
                  search, isTimeDiscrete, flagMass, flagLumped, init, nFolds, inference, PACKAGE = "fdaPDE")

  return(bigsol)
}


CPP_FEM.volume.DE_init_time <- function(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, scaling, fvec, heatStep, heatIter, ndim,
                                   mydim, step_method, direction_method, preprocess_method, stepProposals, tol1, tol2, print,
                                   nfolds, nsimulations, search, isTimeDiscrete, flagMass, flagLumped, init, nFolds, inference)
{

  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1


  ## Set proper type for correct C++ reading
  storage.mode(data) <- "double"
  storage.mode(data_time) <- "double"
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(FEMbasis$mesh$faces) <- "integer"
  storage.mode(mesh_time) <- "double"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(lambda_time) <- "double"
  storage.mode(scaling) <- "double"
  storage.mode(fvec) <- "double"
  storage.mode(heatStep) <- "double"
  heatIter <- as.integer(heatIter)
  storage.mode(heatIter) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(stepProposals) <- "double"
  storage.mode(tol1) <- "double"
  storage.mode(tol2) <- "double"
  storage.mode(print) <- "logical"
  storage.mode(nfolds) <- "integer"
  storage.mode(nsimulations) <- "integer"
  storage.mode(search) <- "integer"
  storage.mode(isTimeDiscrete) <- "logical"
  storage.mode(flagMass) <- "logical"
  storage.mode(flagLumped) <- "logical"
  init <- as.character(init)
  storage.mode(init) <- "character"
  nFolds <- as.integer(nFolds)
  storage.mode(nFolds) <- "integer"
  storage.mode(inference) <- "logical"


  bigsol <- .Call("Density_Initialization_time", data, data_time, FEMbasis$mesh, mesh_time, FEMbasis$order, mydim, ndim,
                  scaling, fvec, heatStep, heatIter, lambda, lambda_time, nfolds, nsimulations, stepProposals, tol1, tol2, print,
                  search, isTimeDiscrete, flagMass, flagLumped, init, nFolds, inference, PACKAGE = "fdaPDE")

  return(bigsol)
}

CPP_FEM.graph.DE_init_time <- function(data, data_time, FEMbasis, mesh_time, lambda, lambda_time, scaling, fvec, heatStep, heatIter, ndim,
                                       mydim, step_method, direction_method, preprocess_method, stepProposals, tol1, tol2, print,
                                       nfolds, nsimulations, search, isTimeDiscrete, flagMass, flagLumped, init, nFolds, inference)
{
  
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  num_sides = 2*dim(FEMbasis$mesh$edges)[1] 
  for(i in 1:num_sides){
    if( dim(FEMbasis$mesh$neighbors[[i]] )[1] > 0)
      FEMbasis$mesh$neighbors[[i]] = FEMbasis$mesh$neighbors[[i]] - 1
  }
  
  ## Set proper type for correct C++ reading
  data <- as.matrix(data)
  storage.mode(data) <- "double"
  storage.mode(data_time) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  for(i in 1:num_sides)
    storage.mode(FEMbasis$mesh$neighbors[[i]]) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(mesh_time) <- "double"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(lambda_time) <- "double"
  storage.mode(scaling) <- "double"
  storage.mode(fvec) <- "double"
  storage.mode(heatStep) <- "double"
  heatIter <- as.integer(heatIter)
  storage.mode(heatIter) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(stepProposals) <- "double"
  storage.mode(tol1) <- "double"
  storage.mode(tol2) <- "double"
  storage.mode(print) <- "logical"
  storage.mode(nfolds) <- "integer"
  storage.mode(nsimulations) <- "integer"
  storage.mode(search) <- "integer"
  storage.mode(isTimeDiscrete) <- "logical"
  storage.mode(flagMass) <- "logical"
  storage.mode(flagLumped) <- "logical"
  init <- as.character(init)
  storage.mode(init) <- "character"
  nFolds <- as.integer(nFolds)
  storage.mode(nFolds) <- "integer"
  storage.mode(inference) <- "logical"
  
  
  bigsol <- .Call("Density_Initialization_time", data, data_time, FEMbasis$mesh, mesh_time, FEMbasis$order, mydim, ndim,
                  scaling, fvec, heatStep, heatIter, lambda, lambda_time, nfolds, nsimulations, stepProposals, tol1, tol2, print,
                  search, isTimeDiscrete, flagMass, flagLumped, init, nFolds, inference, PACKAGE = "fdaPDE")
  
  return(bigsol)
}