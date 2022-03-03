library(fdaPDE)
rm(list=ls())

## Create a 1.5D mesh 
delta = 1 / 2
x = c(0., 1)
y = c(0., delta)
vertices = expand.grid(x,y)
vertices = cbind(vertices[,1], vertices[,2])
edges = matrix(c(1,2,1,3,3,4), nrow=3,ncol=2, byrow=T)

mesh = create.mesh.1.5D(vertices, edges)
mesh = refine.mesh.1.5D(mesh,delta=0.0125)

nnodes=dim(mesh$nodes)[1]
FEMbasis=create.FEM.basis(mesh)
plot(mesh, asp=1, pch=16)

## Generate data
n = 100

data_1 = cbind(rnorm(n*0.7, mean=0.5, sd=0.15), rep(delta,times=n*0.7))
data_2 = cbind(rnorm(n*0.3, mean=0.25, sd=0.05), rep(0. ,times=n*0.3))

data <- rbind(data_1,data_2)
plot(mesh, type="n")
points(data[,1], data[,2], pch=16, col="red3")

## DE without GCV:
lambda = 1e-3
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda)

## DE with exact GCV:
lambda = 10^seq(from=-5, to=-2, by=0.25)
nfolds = 5
sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
              nfolds=nfolds, 
              step_method = "Fixed_Step", direction_method = "BFGS",
              preprocess_method="RightCV")

plot(log10(lambda), sol$CV_err, ylab="CV" )

