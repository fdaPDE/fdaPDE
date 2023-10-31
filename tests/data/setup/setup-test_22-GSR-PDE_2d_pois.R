if(!dir.exists(test_path("../data/GSR-PDE")))
  dir.create(test_path("../data/GSR-PDE"))

if(!dir.exists(test_path("../data/GSR-PDE/test_22"))){
  dir.create(test_path("../data/GSR-PDE/test_22"))  

options(warn=-1)
foldername <- test_path("../data/GSR-PDE/test_22/")

FAMILY = "poisson"

l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

# mesh
x = seq(0,1, length.out = 20)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)

nnodes = dim(mesh$nodes)[1]

FEMbasis = create.FEM.basis(mesh)

set.seed(42)
# locations
nloc = 800
xobs=runif(min=0,max=1,n=nloc)
yobs=runif(min=0,max=1,n=nloc)
loc=cbind(xobs,yobs)

# 2D random field (function f) 
a1=-1.5
a2=0.4

z<-function(p)
{
  a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1]) + 2
  
}

# exact solution
sol_exact=rep(0,length(loc[,1]))
for(i in 1:length(loc[,1])){
  sol_exact[i] <- z(loc[i,])
}

nnodes = dim(mesh$nodes)[1]
sol_nodes = numeric(nnodes)
for(i in 1:nnodes){
  sol_nodes[i] = z(mesh$nodes[i,])
}

param = sol_exact
mu<-inv.link(param)
# sampling response:
set.seed(95)

response <- rpois(length(loc[,1]), lambda = mu)

# Set smoothing parameter
lambda = 10^seq(-5,0,length.out = 20)

#### Test 22.1: Without GCV
invisible(capture.output(sol_ref <- smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda=lambda)))
save(sol_ref, file=paste0(foldername,"/test_22_1.RData"))
#### Test 22.2: With exact GCV
invisible(capture.output(sol_ref <- smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', 
                                 DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')))
save(sol_ref, file=paste0(foldername,"/test_22_2.RData"))
}
