
if(!dir.exists(test_path("../data/STDE-PDE")))
  dir.create(test_path("../data/STDE-PDE"))

if(!dir.exists(test_path("../data/STDE-PDE/test_33"))){
  dir.create(test_path("../data/STDE-PDE/test_33"))

  foldername <- test_path("../data/STDE-PDE/test_33/")
  ## Create a 3D mesh
  set.seed(0)
  data(sphere3Ddata)
  sphere3D <- create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)
  FEMbasis <- create.FEM.basis(sphere3D)

  mesh_time <- seq(0, 1, length.out=3)

  ## Generate data
  fact = 0.01
  n = 100
  data_x <- fact*rnorm(n)
  data_y <- fact*rnorm(n)
  data_z <- fact*rnorm(n)
  times <- runif(n,0,1)
  locations <- cbind(data_x, data_y, data_z)

  ## Density Estimation:
  lambda = 1e-5
  lambda_time = 1e-6
  invisible(capture.output(sol_ref <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                  mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                  heatStep=0.1, heatIter=10, nsimulations=300,
                                                  preprocess_method="NoCrossValidation")))

  save(sol_ref, file=paste0(foldername, "/test_33.RData"))
}