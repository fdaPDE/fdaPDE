test_that("tSR-PDE Cshaped - areal data - 2D",{
  options(warn=-1)
  foldername <- test_path("../data/tSR-PDE/test_14/")
  data(quasicircle2Dareal)
  mesh = quasicircle2Dareal$mesh
  incidence_matrix = quasicircle2Dareal$incidence_matrix
  DatiEsatti = quasicircle2Dareal$data
  
  FEMbasis = create.FEM.basis(mesh)
  time_mesh=seq(0,4,length.out = 11)
  
  DatiEsatti=rep(DatiEsatti,length(time_mesh))*rep(exp(time_mesh),each=length(DatiEsatti))
  
  # Add error to simulate data
  set.seed(5839745)
  data = DatiEsatti + rnorm(length(DatiEsatti), sd = 0.1)
  observations=matrix(data,nrow(incidence_matrix),length(time_mesh))
  
  # Set smoothing parameter
  lambdaS = 10^-6
  lambdaT = 10^-6
  
  # Set BC (equal to zero)
  BC = NULL
  BC$BC_indices = which(mesh$nodesmarkers == 1)
  BC$BC_values = rep(0,length(BC$BC_indices))
  
  # Set sv-PDE parameters
  R = 2.8
  K1 = 0.1
  K2 = 0.2
  beta = 0.5
  
  K_func<-function(points)
  {
    output = array(0, c(2, 2, nrow(points)))
    for (i in 1:nrow(points))
      output[,,i] = 10*rbind(c(points[i,2]^2 + K1*points[i,1]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2),
                               (K1-1)*points[i,1]*points[i,2]),
                             c((K1-1)*points[i,1]*points[i,2],
                               points[i,1]^2 + K1*points[i,2]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2)))
    output
  }
  
  b_func<-function(points)
  {
    output = array(0, c(2, nrow(points)))
    for (i in 1:nrow(points))
      output[,i] = 10*beta*c(points[i,1],points[i,2])
    output
  }
  
  c_func<-function(points)
  {
    rep(c(0), nrow(points))
  }
  
  u_func<-function(points)
  {
    rep(c(0), nrow(points))
  }
  PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)
  
  
  #### Test 14.1: Forcing term = 0, estimated IC, without GCV, monolitic method
  invisible(capture.output(sol<-smooth.FEM.time(observations=observations,
                              incidence_matrix = incidence_matrix,
                              FEMbasis = FEMbasis, time_mesh = time_mesh,
                              lambdaS = lambdaS, lambdaT = lambdaT,
                              BC = BC,
                              PDE_parameters = PDE_parameters,
                              FLAG_PARABOLIC = TRUE,
                              DOF.evaluation = NULL)))
  
  load(file=paste0(foldername,"/test_14_1.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  #### Test 14.2: Forcing term = 0 exact  GCV,  monolitic method
  invisible(capture.output(sol<-smooth.FEM.time(observations=observations,
                              incidence_matrix = incidence_matrix, areal.data.avg = TRUE,
                              FEMbasis = FEMbasis, time_mesh = time_mesh,
                              lambdaS = lambdaS, lambdaT = lambdaT,
                              BC = BC,
                              PDE_parameters = PDE_parameters,
                              FLAG_PARABOLIC = TRUE,
                              lambda.selection.criterion='grid', DOF.evaluation='exact', 
                              lambda.selection.lossfunction='GCV')))
  
  load(file=paste0(foldername,"/test_14_2.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  #### Test 14.3: Forcing term != 0 without GCV
  # forcing function != 0
  u_func<-function(points)
  {
    output = array(0, c(1, nrow(points)))
    for (i in 1:nrow(points))
      output[,i] = -ifelse((points[i,1]^2+points[i,2]^2)<1,100,0)
    output
  }
  
  PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)
  
  invisible(capture.output(sol<-smooth.FEM.time(observations=observations,
                              incidence_matrix = incidence_matrix,
                              FEMbasis = FEMbasis, time_mesh = time_mesh,
                              lambdaS = lambdaS, lambdaT = lambdaT,
                              BC = BC,
                              PDE_parameters = PDE_parameters,
                              FLAG_PARABOLIC = TRUE)))
  
  load(file=paste0(foldername,"/test_14_3.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  #### Test 14.4: Forcing term != 0 exact GCV, monolitic method
  invisible(capture.output(sol<-smooth.FEM.time(observations=observations,
                              incidence_matrix = incidence_matrix, areal.data.avg = TRUE,
                              FEMbasis = FEMbasis, time_mesh = time_mesh,
                              lambdaS = lambdaS, lambdaT = lambdaT,
                              BC = BC,
                              PDE_parameters = PDE_parameters,
                              FLAG_PARABOLIC = TRUE,
                              lambda.selection.criterion='grid', 
                              DOF.evaluation='exact', lambda.selection.lossfunction='GCV')))
  
  load(file=paste0(foldername,"/test_14_4.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  #### Test 14.5: BC != 0  without GCV
  u_func<-function(points)
  {
    rep(c(0), nrow(points))
  }
  PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)
  
  # Add a constat to the data to change true BC
  obs_backup=observations #save a copy of original data
  observations = observations + 5
  
  # Set new value for the BC
  BC$BC_values = rep(5,length(BC$BC_indices))
  
  ### monolitic method
  invisible(capture.output(sol<-smooth.FEM.time(observations=observations,
                              incidence_matrix = incidence_matrix, areal.data.avg = TRUE,
                              FEMbasis = FEMbasis, time_mesh = time_mesh,
                              lambdaS = lambdaS, lambdaT = lambdaT,
                              BC = BC,
                              PDE_parameters = PDE_parameters,
                              FLAG_PARABOLIC = TRUE)))

  load(file=paste0(foldername,"/test_14_5.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  #### Test 14.6: BC != 0 exact GCV
  ### monolitic method
  invisible(capture.output(sol<-smooth.FEM.time(observations=observations,
                              incidence_matrix = incidence_matrix, areal.data.avg = TRUE,
                              FEMbasis = FEMbasis, time_mesh = time_mesh,
                              lambdaS = lambdaS, lambdaT = lambdaT,
                              BC = BC,
                              PDE_parameters = PDE_parameters,
                              FLAG_PARABOLIC = TRUE,
                              lambda.selection.criterion='grid', 
                              DOF.evaluation='exact', lambda.selection.lossfunction='GCV')))
  
  load(file=paste0(foldername,"/test_14_6.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  ### Inference
  
  # Create covariates
  set.seed(509875)
  cov1 = rnorm(length(DatiEsatti), mean = 2, sd = 1)
  
  # Add error to simulate data
  set.seed(5839745)
  data = DatiEsatti + 3*cov1 
  data = data + rnorm(length(DatiEsatti), sd = 0.1)
  observations=matrix(data,nrow(incidence_matrix),length(time_mesh))
  
  # Inference obj: inference test, adding cov, PDE, separable case
  inf_obj <- inferenceDataObjectTimeBuilder (test='oat', interval='oat',  dim=2, n_cov=1, type=c('w', 's', 'esf', 'enh-esf'), component='parametric', n_flip=10000, f_var=T)
  
  #### Test 14.7: overall inference on beta, forcing term = 0 exact  GCV,  monolitic method
  invisible(capture.output(sol<-smooth.FEM.time(observations=observations,
                              covariates = cov1,
                              incidence_matrix = incidence_matrix, areal.data.avg = TRUE,
                              FEMbasis = FEMbasis, time_mesh = time_mesh,
                              lambdaS = lambdaS, lambdaT = lambdaT,
                              BC = BC,
                              PDE_parameters = PDE_parameters,
                              FLAG_PARABOLIC = TRUE,
                              lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV',
                              inference.data.object = inf_obj)))
  
  load(file=paste0(foldername,"/test_14_7.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  expect_equal( max(abs((sol$solution$beta-sol_ref$solution$beta))) < tol, TRUE);
  expect_equal( max(abs((sol$inference$beta$p_values$wald[[1]]-
                           sol_ref$inference$beta$p_values$wald[[1]]))) < tol, TRUE);
  expect_equal( max(abs((sol$inference$beta$p_values$speckman[[1]]-
                           sol_ref$inference$beta$p_values$speckman[[1]]))) < tol, TRUE);
})