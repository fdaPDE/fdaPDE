# compute integral of function
integrate_f <- function(FEM){
  mesh = FEM$FEMbasis$mesh
  if( class(mesh) != "mesh.1.5D")
    stop("Not implemented yet")
  search=1
  mesh$edges=mesh$edges-1

  storage.mode(mesh$edges)<-"integer"
  storage.mode(mesh$order)<-"integer"
  storage.mode(search) <-"integer"
  storage.mode(FEM$coeff) <- "double"
  
  res <- .Call("CPP_integrate_f",mesh,search,as.matrix(FEM$coeff))
  return(res)
}