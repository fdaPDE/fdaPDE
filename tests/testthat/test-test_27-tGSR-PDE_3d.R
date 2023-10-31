test_that("tGSR-PDE 3d",{
  options(warn=-1)
  foldername <- test_path("../data/tGSR-PDE/test_27/")
  
  set.seed(0)
  FAMILY = "poisson"
  l<-make.link("log")
  link<-l$linkfun
  inv.link<-l$linkinv
  
  data(sphere3Ddata)
  sphere3D<-create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)
  
  FEMbasis <- create.FEM.basis(sphere3D)
  
  nodesLocations=sphere3D$nodes
  nnodes = nrow(sphere3D$nodes)
  time_locations = seq(0,1,length.out = 3)
  space_time_locations = cbind(rep(time_locations,each=nnodes),rep(nodesLocations[,1],length(time_locations)),rep(nodesLocations[,2],length(time_locations)),rep(nodesLocations[,3],length(time_locations)))
  
  # Exact test function
  a1 = -1
  a2 = -2
  a3 = -3
  
  func = function(x){
    (a1* sin(x[,2]) +  a2* sin(x[,3]) +  a3*sin(x[,4])) * (sin(x[,1]))
  }
  
  func_evaluation = func(space_time_locations)
  ran=range(func_evaluation)
  
  # generating data
  data = func_evaluation +rnorm(length(func_evaluation),mean=0,sd=0.05*(ran[2]-ran[1]))
  mu = inv.link(data)
  response  <- rpois(length(mu), lambda = mu)
  range(response)
  
  response = matrix(response, nnodes, length(time_locations))
  storage.mode(response) <- "numeric"
  
  # Set smoothing parameter
  lambdaS=10^seq(0, 1, 0.5)
  lambdaT=10^seq(0, 1, 0.5)
  
  #### Test 27.1: without GCV (SEPARABLE)
  invisible(capture.output(sol <- smooth.FEM.time(locations = NULL, observations = response, FEMbasis =FEMbasis, 
                                                      covariates = NULL, time_mesh = time_locations, time_locations = time_locations,
                                                      max.steps=15, family=FAMILY,
                                                      FLAG_PARABOLIC = F,
                                                      lambdaS = lambdaS[1], lambdaT = lambdaT[1])))
  load(file=paste0(foldername,"/test_27_1.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  #### Test 27.2: with exact GCV (SEPARABLE)
  invisible(capture.output(sol <- smooth.FEM.time(locations = NULL, observations = response, FEMbasis =FEMbasis, 
                                                      covariates = NULL, time_mesh = time_locations, time_locations = time_locations,
                                                      max.steps=15, family=FAMILY,
                                                      lambda.selection.lossfunction = "GCV",
                                                      DOF.evaluation = "exact",
                                                      FLAG_PARABOLIC = F,
                                                      lambdaS = lambdaS, lambdaT = lambdaT)))
  load(file=paste0(foldername,"/test_27_2.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
})