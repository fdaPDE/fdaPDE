
# it returns the ID of the element to which the point belong or -1 if the point is not 
# on the graph

isInside <- function(mesh,points, search="tree", redundancy=TRUE){
  
  if(is.null(mesh))
    stop("mesh is null")
  if(is.null(points))
    stop("points is null")
  
  if(class(mesh)=="mesh.2D"){
    mydim=2
    ndim=2
  }else if(class(mesh)=="mesh.2.5D"){
    mydim=2
    ndim=3
  }else if(class(mesh)=="mesh.3D"){
    mydim=3
    ndim=3
    }else if(class(mesh)=="mesh.1.5D"){
    mydim=1
    ndim=2
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
  
  if((class(mesh)=='mesh.2.5D' || class(mesh)=="mesh.1.5D") && search ==3)
    stop("For manifold mesh search must be either tree or naive.")
  
  if(class(mesh)=="mesh.2D" || class(mesh)=="mesh.2.5D"){
    mesh$triangles = mesh$triangles - 1
    mesh$edges = mesh$edges -1
    mesh$neighbors[mesh$neighbors != -1] = mesh$neighbors[mesh$neighbors != -1] - 1
    
    storage.mode(mesh$triangles) <-"integer"
    storage.mode(mesh$edges) <- "integer"
    storage.mode(mesh$neighbors) <-"integer"
    
  }else if(class(mesh)=="mesh.3D"){
    mesh$tetrahedrons = mesh$tetrahedrons - 1
    mesh$faces = mesh$faces - 1
    mesh$neighbors[mesh$neighbors != -1] = mesh$neighbors[mesh$neighbors != -1] - 1
    
    storage.mode(mesh$tetrahedrons) <-"integer"
    storage.mode(mesh$faces) <- "integer"
    storage.mode(mesh$neighbors) <-"integer"
    
  }else if(class(mesh)=="mesh.1.5D"){
    mesh$edges = mesh$edges - 1
    storage.mode(mesh$edges) <-"integer"
    
    # NB neighbors is a #edges x 2 matrix of "vector"
    N = 2*dim(mesh$edges)[1] 
    for(i in 1:N){
      if( dim(mesh$neighbors[[i]] )[1] > 0)
        mesh$neighbors[[i]] = mesh$neighbors[[i]] - 1
      
      storage.mode(mesh$neighbors[[i]]) <- "integer" 
    }
  
  }
  
  storage.mode(mesh$order) <- "integer"
  storage.mode(redundancy) <- "integer"
  storage.mode(search) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(ndim) <-"integer"
  
  res = .Call("isInside",mesh,points,mesh$order,mydim,ndim,search,redundancy)
  
  return(res)
}
  
  
  