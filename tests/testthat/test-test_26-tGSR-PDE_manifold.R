test_that("tGSR-PDE Manifold",{
  options(warn=-1)
  foldername <- test_path("../data/tGSR-PDE/test_26/")
  
  FAMILY = "poisson"
  
  l<-make.link("log")
  link<-l$linkfun
  inv.link<-l$linkinv
  
  data("hub2.5D")
  mesh <- create.mesh.2.5D(nodes = hub2.5D$hub2.5D.nodes,
                           triangles = hub2.5D$hub2.5D.triangles)
  FEMbasis <- create.FEM.basis(mesh)
  
  ### generating data ###
  nodesLocations=mesh$nodes
  set.seed(32)
  x = runif(400, min = -0.6, max = 0.6)
  y = runif(400, min = -0.6, max = 0.6)
  z = runif(400, min =  0.0, max = 1.0)
  points_ = cbind(x,y,z)
  loc = projection.points.2.5D(mesh, points_)
  
  nloc = nrow(loc) 
  a1 = rnorm(1,mean = 1, sd = 1)
  a2 = rnorm(1,mean = 1, sd = 1)
  a3 = rnorm(1,mean = 1, sd = 1)
  time_mesh = as.numeric(0:4)
  
  time_locations = time_mesh
  space_time_locations = cbind(rep(time_locations,each=nloc),
                               rep(loc[,1],length(time_locations)),
                               rep(loc[,2],length(time_locations)),
                               rep(loc[,3],length(time_locations)))
  
  func = function(x)
  {
    (a1*sin(2*pi*x[,2])+a2*sin(2*pi*x[,3])+a3*sin(2*pi*x[,4])+1)*cos(x[,1])
  }
  
  func_evaluation = func(space_time_locations)
  
  lambdaS=10^seq(from=-3, to=0, by=1)
  lambdaT=10^seq(from=-2, to=1, by=1)
  
  cov1=4*sin(2*pi*space_time_locations[,2])*cos(2*pi*space_time_locations[,3])
  cov2=rnorm(nloc*length(time_locations), mean=3, sd=0.1)*rep(exp(-time_locations/length(time_locations)),each=nloc)
  W=cbind(cov1,cov2)
  
  # Fix betas
  beta_exact=c(0.45,0.3)
  
  ran=range(W%*%beta_exact + func_evaluation)
  ran=range(func_evaluation)
  
  ran = range(func_evaluation+ W%*%beta_exact)
  datacov=func_evaluation+ W%*%beta_exact +rnorm(nloc*length(time_locations),mean=0,sd=0.05*(ran[2]-ran[1]))
  
  mu<-inv.link(datacov)
  range(mu)
  response <- rpois(length(mu), lambda = mu)
  
  response <- matrix(response, nloc, length(time_locations))
  storage.mode(response) <- "numeric"
  
  #### Test 26.1: without GCV (PARABOLIC)
  invisible(capture.output(sol <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                                                      locations = loc, observations = response, FEMbasis = FEMbasis, 
                                                      covariates = W, 
                                                      lambdaS = lambdaS[1], lambdaT =  lambdaT[1],
                                                      family = FAMILY, FLAG_PARABOLIC = T)))
  load(file=paste0(foldername,"/test_26_1.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  #### Test 26.2: with exact GCV (PARABOLIC)
  invisible(capture.output(sol <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                                                      locations = loc, observations = response, FEMbasis = FEMbasis, 
                                                      covariates = W,
                                                      lambda.selection.lossfunction = "GCV",
                                                      DOF.evaluation = "exact",
                                                      lambdaS = lambdaS, lambdaT =  lambdaT,
                                                      family = FAMILY, FLAG_PARABOLIC = T)))
  load(file=paste0(foldername,"/test_26_2.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  ### Test 26.3: without GCV (SEPARABLE)
  invisible(capture.output(sol <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                                                      locations = loc, observations = response, FEMbasis = FEMbasis, 
                                                      covariates = W, 
                                                      lambdaS = lambdaS[1], lambdaT =  lambdaT[1],
                                                      family = FAMILY, FLAG_PARABOLIC = F)))
  load(file=paste0(foldername,"/test_26_3.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  ### Test 26.4; with exact GCV (SEPARABLE)
  invisible(capture.output(sol <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                                                      locations = loc, observations = response, FEMbasis = FEMbasis, 
                                                      covariates = W, 
                                                      lambda.selection.lossfunction = "GCV",
                                                      DOF.evaluation = "exact",
                                                      lambdaS = lambdaS, lambdaT =  lambdaT,
                                                      family = FAMILY, FLAG_PARABOLIC = F)))
  load(file=paste0(foldername,"/test_26_4.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
})
