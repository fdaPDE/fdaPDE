##########################################
############## TEST SCRIPT ###############
######## for DENSITY ESTIMATION ##########
##########################################

library(fdaPDE)

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


# 2.5D -------------
library(fdaPDE)
rm(list=ls())

## Create a 2.5D mesh
data(sphere2.5D)
sphere<-create.mesh.2.5D(sphere2.5D$nodes, sphere2.5D$triangles)
# Normalize the sphere
sum_square = sphere$nodes[,1]^2 + sphere$nodes[,2]^2 + sphere$nodes[,3]^2
sphere$nodes[,1] = sphere$nodes[,1]  / sqrt(sum_square)
sphere$nodes[,2] = sphere$nodes[,2]  / sqrt(sum_square)
sphere$nodes[,3] = sphere$nodes[,3]  / sqrt(sum_square)
FEMbasis <- create.FEM.basis(sphere)

## Generating data
library(Directional)
mu <- c(0, -1, 0)
mu <- mu / sqrt( sum(mu^2) )
k=10
beta=0
data <- rkent(100, k, mu, beta)

plot(sphere)
pch3d(data, pch=19, cex=0.1, col="red")

## Density Estimation:
lambda = 0.01
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda)

## Visualization
plot(FEM(exp(sol$g), FEMbasis))

# plot(sphere)
# pch3d(data, pch=19, cex=0.1, col="red")
# pch3d(sol$data, pch=19, cex=0.1, col="blue")


# 3D -------------------------------------------
library(fdaPDE)
rm(list=ls())

## Create a 3D mesh 
data(sphere3Ddata)
sphere3D<-create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)
plot(sphere3D)
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
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda)

## Visualization
plot(FEM(exp(sol$g), FEMbasis))

# 1.5D -------------------------------------------
library(fdaPDE)
rm(list=ls())

## Create a 1.5D mesh 
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
plot(mesh, asp=1)

## Generate data
n = 100

data_1 = cbind(rnorm(n*0.7, mean=0.5, sd=0.25), rep(eps,times=n*0.7))
data_2 = cbind(rnorm(n*0.3, mean=0.25, sd=0.05), rep(0. ,times=n*0.3))

data <- rbind(data_1,data_2)

## Density Estimation:
lambda = 1e-3
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda)

## Visualization
plot(FEM(exp(sol$g), FEMbasis))
