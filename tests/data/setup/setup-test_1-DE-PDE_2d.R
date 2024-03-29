
if(!dir.exists(test_path("../data/DE-PDE")))
  dir.create(test_path("../data/DE-PDE"))

if(!dir.exists(test_path("../data/DE-PDE/test_1"))){
  dir.create(test_path("../data/DE-PDE/test_1"))

foldername <- test_path("../data/DE-PDE/test_1")
  
## Create a 2D mesh over a squared domain
Xbound <- seq(-3, 3, length.out = 10)
Ybound <- seq(-3, 3, length.out = 10)
grid_XY <- expand.grid(Xbound, Ybound)
Bounds <- grid_XY[(grid_XY$Var1 %in% c(-3, 3)) | (grid_XY$Var2 %in% c(-3, 3)), ]
mesh <- create.mesh.2D(nodes = Bounds, order = 1)
mesh <- refine.mesh.2D(mesh, maximum_area = 0.2)
FEMbasis <- create.FEM.basis(mesh)

set.seed(10)

## Generate data
n <- 50
data_x <- rnorm(n)
data_y <- rnorm(n)
data <- cbind(data_x, data_y)

## Density Estimation:

# 1) Cross-validation 
lambda = c(0.001, 0.01, 0.1, 1)
nfolds = 5
invisible(capture.output(sol_ref <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, nfolds=nfolds, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="RightCV")))

save(sol_ref, file = paste0(foldername, "/test_1_1.RData"))

# 2) Lambda fixed
lambda = 0.1
invisible(capture.output(sol_ref <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")))

save(sol_ref, file = paste0(foldername, "/test_1_2.RData"))

# 3) Cross-validation simplified version
lambda = c(0.0001, 0.001, 0.01, 0.1, 1)
nfolds = 5
invisible(capture.output(sol_ref <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, nfolds=nfolds, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="SimplifiedCV")))

save(sol_ref, file = paste0(foldername, "/test_1_3.RData"))

# 4) Initialization given
lambda = 0.1
initF = rep(1, nrow(mesh$nodes))
invisible(capture.output(sol_ref <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, fvec = initF,
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")))
save(sol_ref, file = paste0(foldername, "/test_1_4.RData"))

# 5) step_method = Backtracking method
lambda = 0.1
invisible(capture.output(sol_ref <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              step_method = "Backtracking_Method", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")))
save(sol_ref, file = paste0(foldername, "/test_1_5.RData"))

# 6) step_method = Wolfe method
lambda = 0.1
invisible(capture.output(sol_ref <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              step_method = "Wolfe_Method", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")))
save(sol_ref, file = paste0(foldername, "/test_1_6.RData"))

# 7) direction_method = Gradient
lambda = 0.1
invisible(capture.output(sol_ref <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              step_method = "Fixed_Step", direction_method = "Gradient",
              preprocess_method="NoCrossValidation")))
save(sol_ref,file = paste0(foldername, "/test_1_7.RData"))

# 8) Naive search algorithm
lambda = 0.1
invisible(capture.output(sol_ref <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation", search = "naive")))

save(sol_ref, file = paste0(foldername, "/test_1_8.RData"))

# 9) Confidence Intervals with default scaling
lambda = 0.1
invisible(capture.output(sol_ref <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation", inference = TRUE)))

save(sol_ref, file = paste0(foldername, "/test_1_9.RData"))

# 10) Confidence Intervals with given scaling
lambda = 0.1
scaling_factor = 10
invisible(capture.output(sol_ref <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, scaling = scaling_factor,
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation", inference = TRUE)))

save(sol_ref, file = paste0(foldername, "/test_1_10.RData"))

}
