test_that("GSR-PDE Manifold", {
  options(warn=-1)
  foldername <- test_path("../data/GSR-PDE/test_18/")
  
  # family
  FAMILY = "poisson"
  
  l<-make.link("log")
  link<-l$linkfun
  inv.link<-l$linkinv
  
  # beta
  beta1 = 2
  beta2 = -3
  betas_truth = c(beta1,beta2)
  
  # lambda 
  lambda = c(0.0001,0.001,0.005,0.006,0.007,0.008,0.009,0.0095,
             0.00975,0.01,0.0125,0.015,0.02,0.03,0.04,0.05,0.1,
             1,10,100,1000)
  
  # scale param
  scale.param = 1
  
  # mesh
  data(hub2.5D)
  mesh <- create.mesh.2.5D(nodes = hub2.5D$hub2.5D.nodes,triangles = hub2.5D$hub2.5D.triangles)
  FEMbasis <- create.FEM.basis(mesh)
  
  # locations 
  loc = mesh$nodes
  nloc = dim(loc)[1]
  
  set.seed(0)
  
  # 2.5D random field (function f)
  a1 = rnorm(1,mean = 0, sd = 1)
  a2 = rnorm(1,mean = 0, sd = 1)
  a3 = rnorm(1,mean = 0, sd = 1)
  
  sol_exact = numeric(nloc)
  for (i in 0:(nloc-1)){
    sol_exact[i+1] = a1* sin(2*pi*loc[i+1,1]) +  a2* sin(2*pi*loc[i+1,2]) +  a3*sin(2*pi*loc[i+1,3]) + 7
  }
  
  # covariates
  set.seed(42)
  
  desmat=matrix(0,nrow=nloc, ncol=2)
  desmat[,1]= sin(2*pi*loc[,1])*cos(2*pi*loc[,2])
  desmat[,2]= rnorm(nloc, mean=2, sd=0.1)
  
  # samopling response
  ran=range(desmat%*%betas_truth + sol_exact) 
  param=sol_exact+beta1*desmat[,1]+beta2*desmat[,2]
  
  mu<-inv.link(param)
  response <- rpois(nloc, lambda = mu)
  
  #### Test 18.1: Without GCV
  invisible(capture.output(sol <- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                                                 max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                                                 lambda = lambda)))
  load(file=paste0(foldername,"/test_18_1.RData"))
  expect_equal( max(abs((sol$fit.FEM$coeff-sol_ref$fit.FEM$coeff))) < tol, TRUE);
  expect_equal( max(abs((sol$solution$beta-sol_ref$solution$beta))) < tol, TRUE);
  
  
  #### Test 18.2: grid with exact GCV
  invisible(capture.output(sol <- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                                                 max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                                                 lambda = lambda, lambda.selection.criterion = 'grid', 
                                                 DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')))
  load(file=paste0(foldername,"/test_18_2.RData"))
  expect_equal( max(abs((sol$fit.FEM$coeff-sol_ref$fit.FEM$coeff))) < tol, TRUE);
  expect_equal( max(abs((sol$solution$beta-sol_ref$solution$beta))) < tol, TRUE);
  
})