test_that("Spatial Regression (no covariates) - Linear Network", {
  
  foldername <- test_path("../data/SR-PDE/test_4/")
  
  eps = 1 / 2
  x = c(0., 1)
  y = c(0.,eps)
  vertices = expand.grid(x,y)
  vertices = cbind(vertices[,1], vertices[,2])
  edges = matrix(c(1,2,1,3,3,4), nrow=3,ncol=2, byrow=T)
  
  mesh = create.mesh.1.5D(vertices, edges)
  mesh = refine.mesh.1.5D(mesh,delta=0.0125)
  
  nnodes=dim(mesh$nodes)[1]
  FEMbasis=create.FEM.basis(mesh)
  
  # Exact solution (pointwise at nodes)
  aux.4 = function(x,y){
    h = 1
    source = 4 
    points_ = cbind(x,y)
    idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
    coef = vector(mode="numeric", length=length(x))
    
    for(i in idx.ok){
      delta = abs(points_[i,1] - mesh$nodes[source,1])
      if(delta < h ){
        coef[i] = 3*exp(1/(delta^2 / h^2 -1)+1 ) - 1
        
      }
    }
    
    return(coef)
  }
  aux.3 = function(x,y){
    
    h = eps
    source = 3 
    points_ = cbind(x,y)
    idx.ok = which( abs(points_[,1]-mesh$nodes[source,1]) < 10 * .Machine$double.eps )
    coef = vector(mode="numeric", length=length(x))
    
    for(i in idx.ok){
      delta = abs(points_[i,2] - mesh$nodes[source,2])
      if(delta < h ){
        coef[i] = -1 - 1/h*delta
      }
      
      
    }
    return(coef)
  }
  aux.1 = function(x,y){
    
    h = 1
    source = 1 
    points_ = cbind(x,y)
    idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
    coef = vector(mode="numeric", length=length(x))
    
    for(i in idx.ok){
      delta = abs(points_[i,1] - mesh$nodes[source,1])
      if(delta <= h ){
        coef[i] = -2 - 1/h*delta
      }
    }
    return(coef)
    
  }  
  AUX = function(x,y){
    
    res = aux.1(x,y) + aux.3(x,y) + aux.4(x,y)
    return(res)
  }
  
  sol_exact=AUX(mesh$nodes[,1],mesh$nodes[,2])
  
  # Add error to simulate data
  set.seed(32)
  ran=range(sol_exact)
  data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  # Set smoothing parameter
  lambda = 10^seq(-4,-2,length.out=10)
  
  #### Test 4.1: grid with exact GCV
  invisible(capture.output(sol<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                         lambda.selection.criterion='grid', DOF.evaluation='exact', 
                         lambda.selection.lossfunction='GCV')))
  
  load(file=paste0(foldername,"/test_4_1.RData"))
  expect_equal( max(abs((sol$fit.FEM$coeff-sol_ref$fit.FEM$coeff))) < tol, TRUE);
  
  #### Test 4.5: Inference on f, hypothesis testing and CI, Wald and SF implementations
  inf_obj <- inferenceDataObjectBuilder(test = "sim", interval = "oat", type = c("w", "sf"), component = "nonparametric", dim = 2, n_cov = 0, f0 = AUX, locations_by_nodes = T)
  
  invisible(capture.output(sol<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                         lambda.selection.criterion='grid', DOF.evaluation='exact', 
                         lambda.selection.lossfunction='GCV', inference.data.object = inf_obj)))
  
  load(file=paste0(foldername,"/test_4_2.RData"))
  expect_equal( max(abs((sol$fit.FEM$coeff-sol_ref$fit.FEM$coeff))) < tol, TRUE);
  expect_equal( max(abs(sol$inference$f$p_values$wald[[1]] - 
                   sol_ref$inference$f$p_values$wald[[1]]))< tol, TRUE); 
  
})

