# 1.5D  ---------------
rm(list=ls())
graphics.off()

nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh = create.mesh.1.5D(nodes,edges,order=1)
mesh = refine.mesh.1.5D(mesh=mesh,0.01)
plot(mesh)

FEMbasis <- create.FEM.basis(mesh)

set.seed(10)

## Generate data
n <- 40
data_x <- rnorm(n, mean = 0.5,sd = 0.25*0.5)
data_y <- rnorm(n, mean = 0.5,sd = 0.25*0.25)
data <- cbind(data_x, data_y)
data <- projection.points.1.5D(mesh = mesh,locations = data)

# 1) Cross-validation 
lambda = c(0.001, 0.01, 0.1, 1)
nfolds = 5
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, nfolds=nfolds,
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="RightCV")
sol$lambda
R_plot_graph(FEM(exp(sol$g),FEMbasis))

# 2) Lambda fixed
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")
R_plot_graph(FEM(exp(sol$g),FEMbasis))

# 3) Cross-validation simplified version
lambda = c(0.0001, 0.001, 0.01, 0.1, 1)
nfolds = 5
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, nfolds=nfolds, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="SimplifiedCV")
sol$lambda

R_plot_graph(FEM(exp(sol$g),FEMbasis))

# 4) Initialization given
lambda = 0.1
initF = rep(1, nrow(mesh$nodes))
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, fvec = initF,
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")

R_plot_graph(FEM(exp(sol$g),FEMbasis))

# 5) step_method = Backtracking method
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              step_method = "Backtracking_Method", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")

R_plot_graph(FEM(exp(sol$g),FEMbasis))
# 6) step_method = Wolfe method
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              step_method = "Wolfe_Method", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")

R_plot_graph(FEM(exp(sol$g),FEMbasis))

# 7) direction_method = Gradient
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              step_method = "Fixed_Step", direction_method = "Gradient",
              preprocess_method="NoCrossValidation")

R_plot_graph(FEM(exp(sol$g),FEMbasis))
# 8) Naive search algorithm
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation", search = "naive")

R_plot_graph(FEM(exp(sol$g),FEMbasis))

# 2D ------------------
library(fdaPDE)
rm(list=ls())
graphics.off()

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
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, nfolds=nfolds, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="RightCV")
sol$lambda

## Visualization
image(FEM(exp(sol$g), FEMbasis))

## Visualization on a fine grid
n = 100
X <- seq(-3, 3, length.out = n)
Y<- seq(-3, 3, length.out = n)
grid <- expand.grid(X, Y)
evaluation <- eval.FEM(FEM(FEMbasis, coeff = sol$g), locations = grid)
evaluation <- exp(evaluation)
eval <- matrix(evaluation, n, n)
image2D(x=X, y=Y, z=eval, col=heat.colors(100), xlab="x",ylab="y", contour=list(drawlabels = FALSE), main = "Estimated density")

# 2) Lambda fixed
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")

# ## Visualization 
# image(FEM(exp(sol$g), FEMbasis))

# 3) Cross-validation simplified version
lambda = c(0.0001, 0.001, 0.01, 0.1, 1)
nfolds = 5
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, nfolds=nfolds, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="SimplifiedCV")
sol$lambda

# ## Visualization
# image(FEM(exp(sol$g), FEMbasis))

# 4) Initialization given
lambda = 0.1
initF = rep(1, nrow(mesh$nodes))
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, fvec = initF,
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")

# ## Visualization 
# image(FEM(exp(sol$g), FEMbasis))

# 5) step_method = Backtracking method
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              step_method = "Backtracking_Method", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")

# ## Visualization 
# image(FEM(exp(sol$g), FEMbasis))

# 6) step_method = Wolfe method
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              step_method = "Wolfe_Method", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")

# ## Visualization 
# image(FEM(exp(sol$g), FEMbasis))

# 7) direction_method = Gradient
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              step_method = "Fixed_Step", direction_method = "Gradient",
              preprocess_method="NoCrossValidation")

# ## Visualization 
# image(FEM(exp(sol$g), FEMbasis))

# 8) Naive search algorithm
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation", search = "naive")

# ## Visualization 
# image(FEM(exp(sol$g), FEMbasis))
