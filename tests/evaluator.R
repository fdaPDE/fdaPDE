#NUOVA VERSIONE 
eval.FEM_new <- function(FEM, locations = NULL, incidence_matrix = NULL, search = "tree", bary.locations = NULL)
{
  ##################### Checking parameters, sizes and conversion ################################
  if (is.null(FEM))
    stop("FEM required;  is NULL.")
  if(class(FEM) != "FEM")
    stop("'FEM' is not of class 'FEM'")
  #if locations is null but bary.locations is not null, use the locations in bary.locations
  if(is.null(locations) && !is.null(bary.locations) && is.null(incidence_matrix)) {
    # print("'locations' and 'incidence_matrix' are NULL, evalutation performed in bary.locations")
    locations = bary.locations$locations
    locations = as.matrix(locations)
  }
  if (is.null(locations) && is.null(incidence_matrix))
    stop("'locations' NOR 'incidence_matrix' required;  both are NULL")
  if (!is.null(locations) && !is.null(incidence_matrix))
    stop("'locations' NOR 'incidence_matrix' required; both are given")
  
  # if(!is.null(locations))
  #  if(dim(locations)[1]==dim(FEM$FEMbasis$mesh$nodes)[1] & dim(locations)[2]==dim(FEM$FEMbasis$mesh$nodes)[2])
  #   warning("The locations matrix has the same dimensions as the mesh nodes. If you want to get the FEM object evaluation
  #           at the mesh nodes, use FEM$coeff instead")
  
  if(search == "naive" || search == 1){
    search=1
  }else if(search == "tree" || search == 2){
    search=2
  }else if(search == "walking" || search == 3){
    search=3
  }
  
  if(class(FEM$FEMbasis$mesh)=='mesh.2.5D' && search ==3){
    stop("2.5D search must be either 'tree' or 'naive'")
  }else if(class(FEM$FEMbasis$mesh)=='mesh.1D' && search != 1){
    stop("search must be 'naive'")
  }else if (search != 1 && search != 2 && search != 3){
    stop("search must be either 'tree' or 'naive' or 'walking'")
  }  
  #Check the locations in 'bary.locations' and 'locations' are the same
  if(!is.null(bary.locations) && !is.null(locations)){
    flag=TRUE
    for (i in 1:nrow(locations)) {
      if (!(locations[i,1]==bary.locations$locations[i,1] & locations[i,2] == bary.locations$locations[i,2])) {
        flag = FALSE
        break
      }
    }
    if (flag == FALSE) {
      stop("Locations are not same as the one in barycenter information.")
    }
  }  # end of bary.locations
  
  #MOLTO BENE (era già presente)
  if (is.null(locations)){
    locations <- matrix(nrow = 0, ncol = 2)
  }else{
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  ################## End checking parameters, sizes and conversion #############################
  res <- NULL
  
  if(class(FEM$FEMbasis$mesh)=='mesh.2D'){
    ndim = 2
    mydim = 2
    res = fdaPDE::CPP_eval.FEM(FEM, locations, incidence_matrix, TRUE, ndim, mydim, search, bary.locations)
    
  }else if(class(FEM$FEMbasis$mesh)=='mesh.2.5D'){
    ndim = 3
    mydim = 2
    res = fdaPDE::CPP_eval.manifold.FEM(FEM, locations, incidence_matrix, TRUE, ndim, mydim, search, bary.locations)
  }else if(class(FEM$FEMbasis$mesh)=='mesh.3D'){
    ndim = 3
    mydim = 3
    res = fdaPDE::CPP_eval.volume.FEM(FEM, locations, incidence_matrix, TRUE, ndim, mydim, search, bary.locations)
  }else if(class(FEM$FEMbasis$mesh)=='mesh.1D'){
    ndim = 2
    mydim = 1
    res = CPP_eval.graph.FEM(FEM, locations, incidence_matrix, TRUE, ndim, mydim, search, bary.locations)
  }
    
  return(as.matrix(res))
}


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
  }else{   # to use RObjects! 
    bary.locations = list(locations=locations, element_ids=matrix(nrow=0,ncol=1), barycenters=matrix(nrow=0,ncol=2))
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
  for (i in 1:ncol(coeff)){# as.matrix(coeff[,i]) Altrimenti non riesco ad utilizzare RObject!
    evalmat[,i] <- .Call("eval_FEM_fd_Auxiliary", FEMbasis$mesh, locations, incidence_matrix, as.matrix(coeff[,i]),
                         FEMbasis$order, redundancy, mydim, ndim, search, bary.locations, PACKAGE = "fdaPDE")
  }
  
  #Returning the evaluation matrix
  return(evalmat)
}

