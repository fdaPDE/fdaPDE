test_that("Density Estimation - Unit Sphere", {
  
  foldername <- test_path("../data/DE-PDE/test_2/")
  ## Create a 3D mesh 
  set.seed(0)
  data(sphere3Ddata)
  sphere3D<-create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)
  FEMbasis <- create.FEM.basis(sphere3D)
  
  ## Generate data
  fact = 0.01
  n = 100
  data_x <- fact*rnorm(n)
  data_y <- fact*rnorm(n)
  data_z <- fact*rnorm(n)
  data <- cbind(data_x, data_y, data_z)
  
  ## Density Estimation:
  lambda = 1e-5
  invisible(capture.output(sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda)))
  load(file=paste0(foldername, "/test_2.RData"))
  expect_equal( mean((sol_ex$g-sol$g)^2) < 10*.Machine$double.eps, TRUE);
})
