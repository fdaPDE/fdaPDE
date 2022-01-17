# 1.5D  ---------------
rm(list=ls())
graphics.off()

data("simplenet.1.5D")

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
nfolds = 10

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

# 2D ------------------ ORDER 2
library(fdaPDE)
rm(list=ls())
graphics.off()

## Create a 2D mesh over a squared domain
Xbound <- seq(-3, 3, length.out = 10)
Ybound <- seq(-3, 3, length.out = 10)
grid_XY <- expand.grid(Xbound, Ybound)
Bounds <- grid_XY[(grid_XY$Var1 %in% c(-3, 3)) | (grid_XY$Var2 %in% c(-3, 3)), ]
mesh <- create.mesh.2D(nodes = Bounds, order = 2)

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

# 2) Lambda fixed
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")

# 3) Cross-validation simplified version
lambda = c(0.0001, 0.001, 0.01, 0.1, 1)
nfolds = 5
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, nfolds=nfolds, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="SimplifiedCV")

# 4) Initialization given
lambda = 0.1
initF = rep(1, nrow(mesh$nodes))
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, fvec = initF,
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")

# 5) step_method = Backtracking method
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              step_method = "Backtracking_Method", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")

# 6) step_method = Wolfe method
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              step_method = "Wolfe_Method", direction_method = "BFGS",
              preprocess_method="NoCrossValidation")

# 7) direction_method = Gradient
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              step_method = "Fixed_Step", direction_method = "Gradient",
              preprocess_method="NoCrossValidation")

# 8) Naive search algorithm
lambda = 0.1
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="NoCrossValidation", search = "naive")

####
source("~/Scrivania/fdaPDE/tests/integrate_f.R")
library(spatstat)

data("simplenet.1.5D")

plot(mesh)
mesh.fdaPDE = mesh
FEMbasis.fdaPDE = create.FEM.basis(mesh.fdaPDE)

mu1 = mesh.fdaPDE$nodes[6,]
mu2 = mesh.fdaPDE$nodes[8,]
sigmax = 0.1
sigmay = 0.1
aux = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], sigma_x=sigmax, sigma_y=sigmay) 
{ ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))}


integral_ex = integrate_f( FEM(aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2]),FEMbasis.fdaPDE) )

aux = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], sigma_x=sigmax, sigma_y=sigmay) 
{ 1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))}
my_dens <- linfun(aux, L) # FUNZIONE SU UN NETWORK

coef.ex = aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])

nfolds=10
niter=500
nobs = c(100,200,300,400,500)
nobs = c(200,500)
nnodes = dim(mesh.fdaPDE$nodes)[1]
#nobs = c(50)
N = length(nobs)
#cvec=matrix(0,nrow=mesh$nnodes,ncol=1)

lambda = c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1) 

M=5

mise.fdaPDE = matrix(0,nrow=M,ncol=N)
err.L2 = matrix(0,nrow=M,ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)
integrals.fdaPDE = matrix(0,nrow=M,ncol=N)
sols.fdaPDE = array(0,dim=c(M,N,nnodes) )

lambda.opt = matrix(0,nrow=M,ncol=N)

tot.start = Sys.time()
for(i in 1:N){
  for(j in 1:M){
    PPP = rlpp(nobs[i], my_dens)
    points_x = PPP$data$x
    points_y = PPP$data$y
    points = NULL 
    points = cbind(points_x, points_y)
    points = as.matrix(points)
    colnames(points) = NULL
    
    #nb. nel tempo di fdaPDE è anche compresa la crossvalidazione.
    start.fdaPDE <- Sys.time()
    sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=lambda, 
                        nsimulations=niter, 
                        nfolds=nfolds, 
                        step_method="Fixed_Step", 
                        direction_method="Gradient",
                        preprocess_method = "RightCV")
    end.fdaPDE <- Sys.time()
    
    coef.fdaPDE = exp(sol.fdaPDE$g)
    
    sols.fdaPDE[j,i,] = coef.fdaPDE
    
    mise.fdaPDE[j,i] = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
    
    times.fdaPDE[j,i] = end.fdaPDE - start.fdaPDE
    
    integrals.fdaPDE[j,i] =   integrate_f(FEM(coef.fdaPDE,FEMbasis.fdaPDE))
    
    lambda.opt[j,i] = sol.fdaPDE$lambda
  }
}
tot.end = Sys.time()
tot.time = tot.end-tot.start

