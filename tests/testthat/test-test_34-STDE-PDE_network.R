test_that("Space-Time Density Estimation - Linear Network", {
  
  foldername <- test_path("../data/STDE-PDE/test_34/")
  
  ## Create a linear network
  set.seed(0)
  eps = 1 / 2
  x = c(0., 1)
  y = c(0.,eps)
  vertices = expand.grid(x,y)
  vertices = cbind(vertices[,1], vertices[,2])
  edges = matrix(c(1,2,1,3,3,4), nrow=3,ncol=2, byrow=T)
  
  ## Create a discretization
  mesh = create.mesh.1.5D(vertices, edges)
  mesh = refine.mesh.1.5D(mesh,delta=0.0125)
  nnodes=dim(mesh$nodes)[1]
  FEMbasis=create.FEM.basis(mesh)
  
  mesh_time = seq(0, 1, length.out=3)
  
  ## Generate data
  n = 100
  
  data_1 = cbind(rnorm(n*0.5, mean=0.5, sd=0.1), rep(eps,times=n*0.5))
  data_2 = cbind(rnorm(n*0.5, mean=0.25, sd=0.1), rep(0.,times=n*0.5))
  data_t = runif(n,0,1)
  
  locations <- rbind(data_1,data_2)
  times <- data_t
  
  ## Density Estimation:
  # Lambda fixed
  lambda = 1e-1
  lambda_time = 1e-2
  invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                              mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                              heatStep=0.1, heatIter=10, nsimulations=300,
                                              preprocess_method="NoCrossValidation")))
  
  load(file=paste0(foldername, "/test_34.RData"))
  expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);
  
})