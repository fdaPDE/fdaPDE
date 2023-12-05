test_that("tSR-PDE Manifold",{
  options(warn=-1)
  foldername <- test_path("../data/tSR-PDE/test_16/")
  
  data(hub2.5D)
  set.seed(0)
  mesh <- create.mesh.2.5D(nodes = hub2.5D$hub2.5D.nodes,triangles = hub2.5D$hub2.5D.triangles)
  FEMbasis <- create.FEM.basis(mesh)
  
  # Locations at nodes
  nodesLocations=mesh$nodes
  
  # Exact data - Locations at nodes
  nnodes = nrow(mesh$nodes)
  a1 = rnorm(1,mean = 1, sd = 1)
  a2 = rnorm(1,mean = 1, sd = 1)
  a3 = rnorm(1,mean = 1, sd = 1)
  TimeNodes = 0:4
  
  locations = cbind(rep(TimeNodes,each=nnodes),rep(nodesLocations[,1],length(TimeNodes)),rep(nodesLocations[,2],length(TimeNodes)),rep(nodesLocations[,3],length(TimeNodes)))
  
  func = function(x)
  {
    (a1*sin(2*pi*x[,2])+a2*sin(2*pi*x[,3])+a3*sin(2*pi*x[,4])+1)*cos(x[,1])
  }
  
  func_evaluation = func(locations)
  
  lambdaS=10^seq(-9, -7, 0.5)
  lambdaT=10^seq(-6, -4, 0.5)
  
  lambdaS_par=10^seq(-4, -3, 0.25)
  lambdaT_par=10^seq(1, 1.8, 0.2)
  
  cov1=4*sin(2*pi*locations[,2])*cos(2*pi*locations[,3])
  cov2=rnorm(nnodes*length(TimeNodes), mean=3, sd=0.1)*rep(exp(-TimeNodes/length(TimeNodes)),each=nnodes)
  W=cbind(cov1,cov2)
  
  # Fix betas
  beta_exact=c(0.45,0.3)
  
  ran=range(W%*%beta_exact + func_evaluation)
  ran=range(func_evaluation)

  ran = range(func_evaluation)
  data = func_evaluation +rnorm(nnodes*length(TimeNodes),mean=0,sd=0.05*(ran[2]-ran[1]))
  
  ran = range(func_evaluation+ W%*%beta_exact)
  datacov=func_evaluation+ W%*%beta_exact +rnorm(nnodes*length(TimeNodes),mean=0,sd=0.05*(ran[2]-ran[1]))
  
  data = matrix(data,nrow(mesh$nodes),length(TimeNodes))
  datacov = matrix(datacov,nrow(mesh$nodes),length(TimeNodes))
  
  #########################################SEPARABLE####################################################
  #### Test 16.1
  invisible(capture.output(sol <-  smooth.FEM.time(observations=data,
                           FEMbasis = FEMbasis, time_mesh = TimeNodes, time_locations = TimeNodes,
                           lambdaS = lambdaS, lambdaT = lambdaT,
                           FLAG_PARABOLIC = FALSE)))
  #save(sol_ref, file=paste0(foldername,"/test_16_1.RData"))
  
  load(file=paste0(foldername,"/test_16_1.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  #### Test 16.2
  invisible(capture.output(sol <- smooth.FEM.time(observations=datacov, covariates = W,
                              FEMbasis = FEMbasis, time_mesh = TimeNodes,
                              lambdaS = lambdaS, lambdaT = lambdaT,
                              FLAG_PARABOLIC = FALSE)))

  load(file=paste0(foldername,"/test_16_2.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  expect_equal( max(abs((sol$solution$beta-sol_ref$solution$beta))) < tol, TRUE);
  
  ##########################################PARABOLIC####################################################
  ### MONOLITIC METHOD
  #### Test 16.3
  invisible(capture.output(sol <- smooth.FEM.time(observations=data,
                           FEMbasis = FEMbasis, time_mesh = TimeNodes, time_locations = TimeNodes,
                           lambdaS = lambdaS, lambdaT = lambdaT,
                           FLAG_PARABOLIC = TRUE)))
  
  load(file=paste0(foldername,"/test_16_3.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  #### Test 16.4
  invisible(capture.output(sol <- smooth.FEM.time(observations=datacov[,2:length(TimeNodes)], 
                                                         covariates = W[(1+nrow(mesh$nodes)):(length(TimeNodes)*nrow(mesh$nodes)),],
                              FEMbasis = FEMbasis, time_mesh = TimeNodes,
                              lambdaS = lambdaS, lambdaT = lambdaT,
                              IC=func_evaluation[1:nrow(mesh$nodes)],
                              FLAG_PARABOLIC = TRUE)))
  
  load(file=paste0(foldername,"/test_16_4.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  expect_equal( max(abs((sol$solution$beta-sol_ref$solution$beta))) < tol, TRUE);
  
  ### Inference test: 
  inf_obj <- inferenceDataObjectTimeBuilder(test = 'oat', interval = 'oat', component = 'parametric', type=c('w', 's', 'esf', 'enh-esf'), beta0 = beta_exact, dim=3, n_cov=2)
  
  # #### Test 16.5: overall inference on beta, serparable case
  # invisible(capture.output(sol <- smooth.FEM.time(observations=datacov, covariates = W,
  #                           FEMbasis = FEMbasis, time_mesh = TimeNodes,
  #                           lambdaS = lambdaS, lambdaT = lambdaT,
  #                           lambda.selection.lossfunction = 'GCV', lambda.selection.criterion = 'grid',
  #                           DOF.evaluation = 'stochastic',
  #                           FLAG_PARABOLIC = FALSE,
  #                           inference.data.object = inf_obj)))
  # 
  # load(file=paste0(foldername,"/test_16_5.RData"))
  # expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  # expect_equal( max(abs((sol$solution$beta-sol_ref$solution$beta))) < tol, TRUE);
  # expect_equal( max(abs((sol$inference$beta$p_values$wald[[1]]-
  #                          sol_ref$inference$beta$p_values$wald[[1]]))) < tol, TRUE);
  # expect_equal( max(abs((sol$inference$beta$p_values$speckman[[1]]-
  #                          sol_ref$inference$beta$p_values$speckman[[1]]))) < tol, TRUE);
 
  #### Test 16.6: overall inference on beta, parabolic case
  invisible(capture.output(sol <- smooth.FEM.time(observations=datacov[,2:length(TimeNodes)], covariates = W[(1+nrow(mesh$nodes)):(length(TimeNodes)*nrow(mesh$nodes)),],
                            FEMbasis = FEMbasis, time_mesh = TimeNodes,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            IC=func_evaluation[1:nrow(mesh$nodes)],
                            FLAG_PARABOLIC = TRUE, 
                            inference.data.object = inf_obj)))
  
  load(file=paste0(foldername,"/test_16_6.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  expect_equal( max(abs((sol$solution$beta-sol_ref$solution$beta))) < tol, TRUE);
  expect_equal( max(abs((sol$inference$beta$p_values$wald[[1]]-
                           sol_ref$inference$beta$p_values$wald[[1]]))) < tol, TRUE);
  expect_equal( max(abs((sol$inference$beta$p_values$speckman[[1]]-
                           sol_ref$inference$beta$p_values$speckman[[1]]))) < tol, TRUE);
})