test_that("GSR-PDE binomial - 2d",{
  options(warn=-1)
  foldername <- test_path("../data/GSR-PDE/test_19/")
  
  FAMILY = "binomial"
  
  logit <- function(x){qlogis(x)}
  inv.logit <- function(x){plogis(x)}
  link = logit
  inv.link = inv.logit
  
  # mesh
  x = seq(0,1, length.out = 20)
  y = x
  locations = expand.grid(x,y)
  
  mesh = create.mesh.2D(locations)
  
  nnodes = dim(mesh$nodes)[1]
  
  FEMbasis = create.FEM.basis(mesh)
  
  set.seed(42)
  # locations
  nloc = 800
  xobs=runif(min=0,max=1,n=nloc)
  yobs=runif(min=0,max=1,n=nloc)
  loc=cbind(xobs,yobs)
  
  # 2D random field (function f) 
  a1=-2.5
  a2=0.8
  
  z<-function(p)
  {
    a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1])
  }
  
  # exact solution
  sol_exact=rep(0,length(loc[,1]))
  for(i in 1:length(loc[,1])){
    sol_exact[i] <- z(loc[i,])
  }
  
  nnodes = dim(mesh$nodes)[1]
  sol_nodes = numeric(nnodes)
  for(i in 1:nnodes){
    sol_nodes[i] = z(mesh$nodes[i,])
  }
  
  param = sol_exact
  mu<-inv.link(param)
  
  # sampling response:
  response <- rbernoulli(length(loc[,1]),p = mu)
  
  # Set smoothing parameter
  lambda = 10^seq(-5,0,length.out = 20)
  
  #### Test 19.1: Without GCV
  invisible(capture.output(sol <- smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis = FEMbasis, covariates = NULL,
                                                 max.steps=15, fam=FAMILY, mu0=NULL, 
                                                 scale.param=NULL,
                                                 lambda = lambda)))
  load(file=paste0(foldername,"/test_19_1.RData"))
  expect_equal( max(abs((sol$fit.FEM$coeff-sol_ref$fit.FEM$coeff))) < tol, TRUE);
  
  #### Test 19.2: With exact GCV
  invisible(capture.output(sol <- smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                                         max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                                         lambda = lambda, lambda.selection.criterion = 'grid', 
                                                         DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')))
  
  load(file=paste0(foldername,"/test_19_2.RData"))
  expect_equal( max(abs((sol$fit.FEM$coeff-sol_ref$fit.FEM$coeff))) < tol, TRUE);
})