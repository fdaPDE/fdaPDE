
test_that("Space-Time Density Estimation - Square domain", {

foldername <- test_path("../data/STDE-PDE/test_32")

## Create a 2D mesh over a squared domain
Xbound <- seq(-3, 3, length.out = 10)
Ybound <- seq(-3, 3, length.out = 10)
grid_XY <- expand.grid(Xbound, Ybound)
Bounds <- grid_XY[(grid_XY$Var1 %in% c(-3, 3)) | (grid_XY$Var2 %in% c(-3, 3)), ]
mesh <- create.mesh.2D(nodes = Bounds, order = 1)
mesh <- refine.mesh.2D(mesh, maximum_area = 0.2)
FEMbasis <- create.FEM.basis(mesh)

mesh_time <- seq(0, 1, length.out=3)

set.seed(10)

## Generate data
n <- 50
data_x <- rnorm(n)
data_y <- rnorm(n)
times <- runif(n,0,1)
locations <- cbind(data_x, data_y)

## Density Estimation:
# 1) Cross-validation
lambda = c(0.001, 0.01, 0.1)
lambda_time = c(0.001,0.0001,0.00001)
nfolds = 9
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nfolds=nfolds,
                                                nsimulations=300, preprocess_method="RightCV")))
load(file = paste0(foldername, "/test_32_1.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 2) Lambda fixed
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_2.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 3) Cross-validation simplified version
lambda = c(0.001, 0.01, 0.1)
lambda_time = c(0.001,0.0001,0.00001)
nfolds = 9
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nfolds=nfolds,
                                                nsimulations=300, preprocess_method="SimplifiedCV")))

load(file = paste0(foldername, "/test_32_3.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 4) Initialization given
SPLINE_DEGREE <- 3
lambda = 0.1
lambda_time = 0.0001
initF = rep(1, nrow(mesh$nodes)* (length(mesh_time)+SPLINE_DEGREE-1))
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                fvec = initF, heatStep=0.1, heatIter=10, nsimulations=300,
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_4.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 5) step_method = Backtracking method
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                step_method = "Backtracking_Method",
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_5.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 6) step_method = Wolfe method
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                step_method = "Wolfe_Method",
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_6.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 7) direction_method = Gradient
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                step_method = "Fixed_Step", direction_method = "Gradient",
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_7.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 8) direction_method = BFGS
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                step_method = "Fixed_Step", direction_method = "BFGS",
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_8.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 9) direction_method = L-BFGS
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                step_method = "Fixed_Step", direction_method = "L-BFGS5",
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_9.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 10) direction_method = L-BFGS10
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                step_method = "Fixed_Step", direction_method = "L-BFGS10",
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_10.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 11) direction_method = ConjugateGradientFR
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                step_method = "Fixed_Step", direction_method = "ConjugateGradientFR",
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_11.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 12) direction_method = ConjugateGradientFR
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                step_method = "Fixed_Step", direction_method = "ConjugateGradientPRP",
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_12.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 13) direction_method = ConjugateGradientHS
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                step_method = "Fixed_Step", direction_method = "ConjugateGradientHS",
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_13.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 14) direction_method = ConjugateGradientDY
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                step_method = "Fixed_Step", direction_method = "ConjugateGradientDY",
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_14.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 15) direction_method = ConjugateGradientCD
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                step_method = "Fixed_Step", direction_method = "ConjugateGradientCD",
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_15.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 16) direction_method = ConjugateGradientLS
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                step_method = "Fixed_Step", direction_method = "ConjugateGradientLS",
                                                preprocess_method="NoCrossValidation")))

load(file = paste0(foldername, "/test_32_16.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);

# 17) Naive search algorithm
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                preprocess_method="NoCrossValidation", search = "naive")))

load(file = paste0(foldername, "/test_32_17.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);


# 18) Confidence Intervals with default scaling
lambda = 0.1
lambda_time = 0.0001
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                heatStep=0.1, heatIter=10, nsimulations=300,
                                                preprocess_method="NoCrossValidation", inference = TRUE)))

load(file = paste0(foldername, "/test_32_18.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);
expect_equal( max(abs((sol_ref$g_CI_L-sol$g_CI_L))) < tol, TRUE);
expect_equal( max(abs((sol_ref$g_CI_U-sol$g_CI_U))) < tol, TRUE);

# 19) Confidence Intervals with given scaling
lambda = 0.1
lambda_time = 0.0001
scaling_factor = 10
invisible(capture.output(sol <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                mesh_time = mesh_time, lambda = lambda, lambda_time = lambda_time,
                                                scaling = scaling_factor, heatStep=0.1, heatIter=10, nsimulations=300,
                                                preprocess_method="NoCrossValidation", inference = TRUE)))

load(file = paste0(foldername, "/test_32_19.RData"))
expect_equal( max(abs((sol_ref$g-sol$g))) < tol, TRUE);
expect_equal( max(abs((sol_ref$g_CI_L-sol$g_CI_L))) < tol, TRUE);
expect_equal( max(abs((sol_ref$g_CI_U-sol$g_CI_U))) < tol, TRUE);
})