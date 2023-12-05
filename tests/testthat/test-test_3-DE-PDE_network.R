test_that("Density Estimation - Linear Network", {
  foldername <- test_path("../data/DE-PDE/test_3/")
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
  
  ## Generate data
  n = 100
  
  data_1 = cbind(rnorm(n*0.7, mean=0.5, sd=0.25), rep(eps,times=n*0.7))
  data_2 = cbind(rnorm(n*0.3, mean=0.25, sd=0.05), rep(0. ,times=n*0.3))
  
  data <- rbind(data_1,data_2)
  
  ## Density Estimation:
  lambda = 1e-3
  invisible(capture.output(sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda)))
  load(file=paste0(foldername, "/test_3.RData"))
  expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);
})
