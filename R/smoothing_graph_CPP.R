
CPP_eval.graph.FEM <-function(FEM, locations, incidence_matrix, redundancy, ndim, mydim, search, bary.locations){
  
  FEMbasis = FEM$FEMbasis
  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation
  
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  
  # Imposing types, this is necessary for correct reading from C++
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  
  storage.mode(FEMbasis$order) <- "integer"
  coeff <- as.matrix(FEM$coeff)
  storage.mode(coeff) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(locations) <- "double"
  storage.mode(redundancy) <- "integer"
  storage.mode(search) <- "integer"
  
  if(!is.null(bary.locations)){
    storage.mode(bary.locations$element_ids) <- "integer"
    element_ids <- as.matrix(bary.locations$element_ids)
    storage.mode(bary.locations$barycenters) <- "double"
    barycenters <- as.matrix(bary.locations$barycenters)
  }else{   
    bary.locations = list(locations=matrix(nrow=0,ncol=ndim), element_ids=matrix(nrow=0,ncol=1), barycenters=matrix(nrow=0,ncol=2))
    storage.mode(bary.locations$element_ids) <- "integer"
    storage.mode(bary.locations$barycenters) <- "double"
  }
  
  #Calling the C++ function "eval_FEM_fd" in FEM_Eval.cpp
  evalmat = matrix(0,max(nrow(locations),nrow(incidence_matrix)),ncol(coeff))
  for (i in 1:ncol(coeff)){
    evalmat[,i] <- .Call("eval_FEM_fd", FEMbasis$mesh, locations, incidence_matrix, as.matrix(coeff[,i]),
                         FEMbasis$order, redundancy, mydim, ndim, search, bary.locations, PACKAGE = "fdaPDE")
  }
  
  #Returning the evaluation matrix
  return(evalmat)
}
