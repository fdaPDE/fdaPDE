test_that("tGSR-PDE Network",{
  options(warn=-1)
  foldername <- test_path("../data/tGSR-PDE/test_24/")
  
  FAMILY = "poisson"
  l<-make.link("log")
  link<-l$linkfun
  inv.link<-l$linkinv
  
  set.seed(0)
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
  
  # Locations at nodes
  nodesLocations=mesh$nodes
  
  nnodes = nrow(mesh$nodes)
  time_locations = seq(0,pi,length.out = 6)
  space_time_locations = cbind(rep(time_locations,each=nnodes),rep(nodesLocations[,1],length(time_locations)),rep(nodesLocations[,2],length(time_locations)) )
  
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
  func = function(x,y,t){
    
    res = (aux.1(x,y) + aux.3(x,y) + aux.4(x,y)) * cos(t)
    return(res)
  }
  
  func_evaluation = func(space_time_locations[,2], 
                         space_time_locations[,3], 
                         space_time_locations[,1])
  
  ran = range(func_evaluation)
  data = func_evaluation +rnorm(nnodes*length(time_locations),
                                mean=0,sd=0.05*(ran[2]-ran[1]))
  mu = inv.link(data)
  response = rpois(nnodes*length(time_locations), lambda = mu)
  
  data = matrix(data,nrow(mesh$nodes),length(time_locations))
  
  storage.mode(response)<-"numeric"
  response <- matrix(response, nrow(mesh$nodes), length(time_locations))
  
  # Set smoothing parameters
  lambdaS = 10^seq(1, 2, by=0.25)
  lambdaT = 10^seq(0, 1, by=0.25)
  
  #### Test 24.1: with exact GCV (PARABOLIC)
  invisible(capture.output(sol <- smooth.FEM.time(observations=response,
                                                      FEMbasis = FEMbasis, time_mesh = time_locations, 
                                                      time_locations = time_locations,
                                                      lambdaS = lambdaS, lambdaT = lambdaT, 
                                                      family = FAMILY, 
                                                      lambda.selection.criterion = "grid",
                                                      DOF.evaluation = "exact", 
                                                      lambda.selection.lossfunction = "GCV",
                                                      FLAG_PARABOLIC = TRUE)))
  load(file=paste0(foldername,"/test_24_1.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  #### Test 24.2: with exact GCV (SEPARABLE)
  invisible(capture.output(sol <- smooth.FEM.time(observations=response,
                                                      FEMbasis = FEMbasis, time_mesh = time_locations, 
                                                      time_locations = time_locations,
                                                      lambdaS = lambdaS, lambdaT = lambdaT, 
                                                      family = FAMILY, 
                                                      lambda.selection.criterion = "grid",
                                                      DOF.evaluation = "exact", 
                                                      lambda.selection.lossfunction = "GCV",
                                                      FLAG_PARABOLIC = FALSE)))
  load(file=paste0(foldername,"/test_24_2.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
})