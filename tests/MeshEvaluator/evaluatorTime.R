eval.FEM.time.new <- function(FEM.time, locations = NULL, time.instants = NULL, space.time.locations = NULL, 
                          incidence_matrix = NULL, lambdaS = 1, lambdaT = 1, search = "tree", bary.locations = NULL)
{
  if (is.null(FEM.time))
    stop("FEM.time required;  is NULL.")
  if(class(FEM.time) != "FEM.time")
    stop("'FEM.time' is not of class 'FEM.time'")
  
  # save flag for areal evaluation
  FLAG_AREAL_EVALUATION = FALSE
  if (!is.null(incidence_matrix)) 
    FLAG_AREAL_EVALUATION = TRUE
  
  # save flag for tensorization: if TRUE tensorization has to be performed, 
  # if FALSE the 'space.time.locations' has already have been constructed
  FLAG_TENSORIZE = TRUE
  if (!is.null(space.time.locations)) 
    FLAG_TENSORIZE = FALSE
  
  # check that when space.time.locations is provided the other fields are NULL
  if (!FLAG_TENSORIZE && FLAG_AREAL_EVALUATION) 
    stop("'incidence_matrix' and 'space.time.locations' both provided")
  if (!FLAG_TENSORIZE && !is.null(locations)) 
    stop("'locations' and 'space.time.locations' both provided")
  if (!FLAG_TENSORIZE && !is.null(time.instants)) 
    stop("'time.instants' and 'space.time.locations' both provided")
  
  if (FLAG_TENSORIZE) {
    if(is.null(time.instants))
      stop("'time.instants' must be given if space.time.locations is NULL")
    if (!is.vector(time.instants) && nrow(time.instants)!=1 && ncol(time.instants)!=1) 
      stop("'time.instants' must be a vector")
    time.instants = as.vector(time.instants)
    
    if (!FLAG_AREAL_EVALUATION) {
      incidence_matrix <- matrix(nrow = 0, ncol = 1)
      
      if(is.null(locations) && !is.null(bary.locations)) {
        # print("space.time.locations', 'locations' and 'incidence_matrix' are NULL, evalutation performed in bary.locations")
        locations = bary.locations$locations
        locations = as.matrix(locations)
      }
      if (is.null(locations))
        stop("'space.time.locations', 'locations' and 'incidence_matrix' are NULL")
      if (ncol(locations)!=2 && ncol(locations)!=3)
        stop("'locations' and 'bary.locations$locations' when provided must have 2 or 3 columns")
      
      #Check the locations in 'bary.locations' and 'locations' are the same
      if(!is.null(bary.locations) && !is.null(locations))
      {
        if (dim(locations)!=dim(bary.locations$locations))
          stop("'locations' and 'bary.locations$locations' are different")
        if (sum(abs(locations - bary.locations$locations))!=0)
          stop("'locations' and 'bary.locations$locations' are different")
        
        #Repeat 'bary.locations' to have the same size as 'space.time.locations' (needed to convert size for C++)
        time_len = length(FEM.time$mesh_time)
        bary.locations$locations = matrix( rep( t(bary.locations$locations) , time_len ) , ncol =  ncol(bary.locations$locations) , byrow = TRUE )
        bary.locations$element_ids = rep(bary.locations$element_ids, time_len)
        bary.locations$barycenters = matrix( rep( t(bary.locations$barycenters) , time_len ) , ncol =  ncol(bary.locations$barycenters) , byrow = TRUE ) 
      }  # end of bary.locations
      
      time_locations = rep(time.instants,each=nrow(locations))
      space.time.locations = cbind(rep(locations[,1],length(time.instants)),rep(locations[,2],length(time.instants)))
      if(ncol(locations)==3)
        space.time.locations=cbind(space.time.locations,rep(locations[,3],length(time.instants)))
    }else{
      space.time.locations <- matrix(nrow=0, ncol=1)
      
      check2D =   class(FEM.time$FEMbasis$mesh)   == 'mesh.2D'   && ncol(incidence_matrix)!=nrow(FEM.time$FEMbasis$mesh$triangles)
      check2.5D = class(FEM.time$FEMbasis$mesh)   == 'mesh.2.5D' && ncol(incidence_matrix)!=nrow(FEM.time$FEMbasis$mesh$triangles)
      check3D =   class(FEM.time$FEMbasis$mesh)   == 'mesh.3D'   && ncol(incidence_matrix)!=nrow(FEM.time$FEMbasis$mesh$tetrahedrons)
      
      if( check2D || check2.5D || check3D )
        stop("incidence_matrix has wrong number of columns")
      time_locations = rep(time.instants,each=nrow(incidence_matrix))
      incidence_matrix = matrix(rep(incidence_matrix,length(time.instants)), nrow = nrow(incidence_matrix)*length(time.instants),
                                                                             ncol = ncol(incidence_matrix),
                                                                             byrow = TRUE)
    }
  }else{
    if(dim(space.time.locations)[2]<3)
      stop("'space.time.locations' requires at least t,X,Y")
    # if(dim(space.time.locations)[1]==dim(FEM.time$FEMbasis$mesh$nodes)[1] & (dim(space.time.locations)[2]-1)==dim(FEM.time$FEMbasis$mesh$nodes)[2])
    #   warning("The space.time.locations matrix has the same dimensions as the mesh nodes. If you want to get the FEM.time object evaluation
    #       at the mesh nodes, use FEM.time$coeff instead")
    time_locations <- space.time.locations[,1]
    space.time.locations <- space.time.locations[,2:dim(space.time.locations)[2]]
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  #conversion of search type
  if(search == "naive" || search == 1)
    search=1
  else if(search == "tree" || search == 2)
    search=2
  else if(search == "walking" || search == 3)
    search=3
  if (search != 1 & search != 2 & search != 3)
    stop("search must be either tree or naive or walking.")
  
  if(class(FEM.time$FEMbasis$mesh)=='mesh.2.5D' & search ==3)
    stop("2.5D search must be either tree or naive.")
  
  
  if(dim(FEM.time$coeff)[2]>1||dim(FEM.time$coeff)[3]>1)
  {
    if(dim(FEM.time$coeff)[2]>1 && lambdaS==1)
      warning("the first value of lambdaS is being used")
    if(dim(FEM.time$coeff)[3]>1 && lambdaT==1)
      warning("the first value of lambdaT is being used")
    f = FEM.time(coeff=array(FEM.time$coeff[,lambdaS,lambdaT]),time_mesh=FEM.time$mesh_time,FEMbasis=FEM.time$FEMbasis,FLAG_PARABOLIC=FEM.time$FLAG_PARABOLIC)
  }
  else
    f = FEM.time
  
  res <- NULL
  if(class(FEM.time$FEMbasis$mesh)=='mesh.2D'){
    ndim = 2
    mydim = 2
    res = CPP_eval.FEM.time.new(f, space.time.locations, time_locations, incidence_matrix, FEM.time$FLAG_PARABOLIC, TRUE, ndim, mydim, search, bary.locations)
  }else if(class(FEM.time$FEMbasis$mesh)=='mesh.2.5D'){
    ndim = 3
     mydim = 2
   res = CPP_eval.manifold.FEM.time.new(f, space.time.locations, time_locations, incidence_matrix, FEM.time$FLAG_PARABOLIC, TRUE, ndim, mydim, search, bary.locations)
  }else if(class(FEM.time$FEMbasis$mesh)=='mesh.3D'){
    ndim = 3
    mydim = 3
    res = CPP_eval.volume.FEM.time.new(f, space.time.locations, time_locations, incidence_matrix, FEM.time$FLAG_PARABOLIC, TRUE, ndim, mydim, search, bary.locations)
  }
  
  return(as.matrix(res))
}


CPP_eval.FEM.time.new <- function(FEM.time, locations, time_locations, incidence_matrix, FLAG_PARABOLIC, redundancy, ndim, mydim, search, bary.locations)
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
  
  #same trick of CPP_eval.FEM.new in eval.R "else part"
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
  
  # same trick in evaluator.R eval_FEM_fd_Auxiliary_new as.matrix(coeff[,i])
  evalmat = matrix(0,max(nrow(locations),nrow(incidence_matrix)),ncol(coeff)) #nb) locations sono ripetute x time_instant ma non mi importa
  for (i in 1:ncol(coeff)){
    evalmat[,i] <- .Call("eval_FEM_time_Auxiliary_new", FEMbasis$mesh, FEM.time$mesh_time, locations, time_locations, incidence_matrix, as.matrix(coeff[,i]),
                         FEMbasis$order, redundancy, FLAG_PARABOLIC, mydim, ndim, search, bary.locations, PACKAGE = "fdaPDE")
  }
  
  #Returning the evaluation matrix
  return(evalmat)
}

CPP_eval.manifold.FEM.time.new = function(FEM.time, locations, time_locations, incidence_matrix, FLAG_PARABOLIC, redundancy, ndim, mydim, search, bary.locations)
{
  FEMbasis = FEM.time$FEMbasis
  # C++ function for manifold works with vectors not with matrices
  
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
  }else{
    bary.locations = list(locations=matrix(nrow=0,ncol=ndim), element_ids=matrix(nrow=0,ncol=1), barycenters=matrix(nrow=0,ncol=2))
    storage.mode(bary.locations$locations) <- "double"
    storage.mode(bary.locations$element_ids) <- "integer"
    storage.mode(bary.locations$barycenters) <- "double"
  }
  
  # if (search == 1) { #use Naive search
  #    print('This is Naive Search')
  #  } else if (search == 2)  { #use Tree search (default)
  #    print('This is Tree Search')
  #  }
  
  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  evalmat = matrix(0,max(nrow(locations),nrow(incidence_matrix)),ncol(coeff))
  for (i in 1:ncol(coeff)){
    evalmat[,i] <- .Call("eval_FEM_time_Auxiliary_new", FEMbasis$mesh, FEM.time$mesh_time, locations, time_locations, incidence_matrix, as.matrix(coeff[,i]),
                         FEMbasis$order, redundancy, FLAG_PARABOLIC, mydim, ndim, search, bary.locations, PACKAGE = "fdaPDE")
  }
  
  #Returning the evaluation matrix
  return(evalmat)
}

CPP_eval.volume.FEM.time.new = function(FEM.time, locations, time_locations, incidence_matrix, FLAG_PARABOLIC, redundancy, ndim, mydim, search, bary.locations)
{
  FEMbasis = FEM.time$FEMbasis
  
  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
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
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$faces) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
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
  
  
  # if (search == 1) { #use Naive search
  #   print('This is Naive Search')
  # } else if (search == 2)  { #use Tree search (default)
  #   print('This is Tree Search')
  # }
  
  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  evalmat = matrix(0,max(nrow(locations),nrow(incidence_matrix)),ncol(coeff))
  for (i in 1:ncol(coeff)){
    evalmat[,i] <- .Call("eval_FEM_time_Auxiliary_new", FEMbasis$mesh, FEM.time$mesh_time, locations, time_locations, incidence_matrix, as.matrix(coeff[,i]),
                         FEMbasis$order, redundancy, FLAG_PARABOLIC, mydim, ndim, search, bary.locations, PACKAGE = "fdaPDE")
  }
  
  #Returning the evaluation matrix
 return(evalmat)
}
