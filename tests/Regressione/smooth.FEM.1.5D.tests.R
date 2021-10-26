##########################################
############## TEST SCRIPT ###############
##########################################

library(fdaPDE)

####### 1.5D ########

#### Test 1: Bifurcation domain ####
#            locations = nodes 
#            no covariates
#            no BC
#            order FE = 1

nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh = create.mesh.1.5D(nodes,edges,order=1)
mesh = refine.mesh.1.5D(mesh=mesh,0.1)
plot(mesh)

nnodes=dim(mesh$nodes)[1]
FEMbasis=create.FEM.basis(mesh)

# Exact solution (pointwise at nodes)
sol_exact=fs.test(mesh$nodes[,1],mesh$nodes[,2])
#plot(FEM(sol_exact, FEMbasis))

# Add error to simulate data
set.seed(7893475)
ran=range(sol_exact)
data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^seq(-2,0.5,by=0.25)

#### Test 1.1: Without GCV
output_CPP.1<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda)

#### Test 1.2: grid with exact GCV

output_CPP.2<-smooth.FEM.new(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

#### Test 1.3: grid with stochastic GCV
output_CPP.3<-smooth.FEM.new(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')

#### Test 1.4: Newton_fd method with stochastic GCV, default initial lambda and tolerance
output_CPP.4<-smooth.FEM.new(observations=data, FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')


### Test 2: Bifurcation domain ####
#            locations != nodes 
#            with covariates
#            no BC
#            order FE = 1
library(fdaPDE)
rm(list=ls())
graphics.off()

nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh = create.mesh.1.5D(nodes,edges,order=1)
mesh = refine.mesh.1.5D(mesh=mesh,0.1)
plot(mesh)

FEMbasis=create.FEM.basis(mesh)

# Test function 
f = function(x, y,z=1){
  phi = (asin(y/sqrt(x^2+y^2)))
  theta = acos(z/sqrt(x^2+y^2+z^2))
  # rho = 1
  
  sin(4*(1/2*sin(theta)*exp(-sin(theta)^2)+1)*theta)*cos(2*(1/2*cos(phi)*exp(-cos(phi)^2)+1)*phi)
}

# Exact solution (pointwise at nodes)
sol_exact=f(mesh$nodes[,1], mesh$nodes[,2])
#plot(FEM(sol_exact, FEMbasis))

# Generate data locations on the sphere
set.seed(598944)
ndata = 500
locations = matrix(rnorm(ndata*2, mean=0.5,sd=0.25 ), ncol = 2)

# Generate covariate and data
cov1 = runif(ndata, min = -0.5, max = 0.5)
DatiEsatti = f(locations[,1],locations[,2]) + cov1

# Add error to simulate data
set.seed(7893475)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(ndata, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Project locations on the mesh
projected_locations = projection.points.1.5D(mesh, locations)

# Set smoothing parameter
lambda = 10^seq(-4,-2,by=0.25)

#### Test 2.1: Without GCV
output_CPP.1<-smooth.FEM(observations=data, 
                       locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, 
                       lambda=lambda[1])
output_CPP$solution$beta

#### Test 2.2: grid with exact GCV
output_CPP.2<-smooth.FEM.new(observations=data, locations = projected_locations,
                       covariates = cov1,
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP.2$optimization$GCV_vector)
output_CPP.2$solution$beta

#### Test 2.3: grid with stochastic GCV
output_CPP.3<-smooth.FEM.new(observations=data, locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP.3$optimization$GCV_vector)

output_CPP.3$solution$beta

#### Test 2.4: Newton method with exact GCV, default initial lambda and tolerance
output_CPP.4<-smooth.FEM.new(observations=data, locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')


output_CPP.4$solution$beta

#### Test 2.5: Newton_fd method with exact GCV, default initial lambda and tolerance
output_CPP.5<-smooth.FEM.new(observations=data, locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton_fd', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

output_CPP.5$solution$beta

#### Test 2.6: Newton_fd method with stochastic GCV, default initial lambda and tolerance
output_CPP.6<-smooth.FEM.new(observations=data, locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')

output_CPP.6$solution$beta


### Test 3: Bifurcation domain ####
#            locations != nodes 
#            with covariates
#            no BC
#            order FE = 1
library(fdaPDE)
rm(list=ls())
graphics.off()

nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh = create.mesh.1.5D(nodes,edges,order=1)
mesh = refine.mesh.1.5D(mesh=mesh,0.1)
plot(mesh)

FEMbasis=create.FEM.basis(mesh)

# Test function 
f = function(x, y,z=1){
  phi = (asin(y/sqrt(x^2+y^2)))
  theta = acos(z/sqrt(x^2+y^2+z^2))
  # rho = 1
  
  sin(4*(1/2*sin(theta)*exp(-sin(theta)^2)+1)*theta)*cos(2*(1/2*cos(phi)*exp(-cos(phi)^2)+1)*phi)
}

# Exact solution (pointwise at nodes)
sol_exact=f(mesh$nodes[,1], mesh$nodes[,2])
#plot(FEM(sol_exact, FEMbasis))

# Generate data locations on the sphere
set.seed(598944)
ndata = 500
locations = matrix(rnorm(ndata*2, mean=0.5,sd=0.25 ), ncol = 2)

# Generate covariate and data
cov1 = rnorm(ndata, mean = 1, sd = 2)
cov2 = sin(locations[,1])
beta_ex = c(2.0,1.0)
DatiEsatti =  f(locations[,1],locations[,2]) + beta_ex[1]*cov1 + beta_ex[2]*cov2

# Add error to simulate data
set.seed(7893475)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(ndata, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Project locations on the mesh
projected_locations = projection.points.1D(mesh, locations)

# Set smoothing parameter
lambda = 10^seq(-4,-2,by=0.25)

#### Test 2.1: Without GCV
output_CPP.1<-smooth.FEM.new(observations=data, 
                             locations = projected_locations, 
                             covariates = cbind(cov1,cov2),
                             FEMbasis=FEMbasis, 
                             lambda=lambda[1])
output_CPP.1$solution$beta

#### Test 2.2: grid with exact GCV
output_CPP.2<-smooth.FEM.new(observations=data, locations = projected_locations,
                             covariates = cbind(cov1,co2),
                             FEMbasis=FEMbasis, lambda=lambda,
                             lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP.2$optimization$GCV_vector)
output_CPP.2$solution$beta

#### Test 2.3: grid with stochastic GCV
output_CPP.3<-smooth.FEM(observations=data, locations = projected_locations, 
                             covariates = cbind(cov1,cov2),
                             FEMbasis=FEMbasis, lambda=lambda,
                             lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP.3$optimization$GCV_vector)

output_CPP.3$solution$beta

#### Test 2.4: Newton method with exact GCV, default initial lambda and tolerance
output_CPP.4<-smooth.FEM(observations=data, locations = projected_locations, 
                             covariates = cbind(cov1,cov2),
                             FEMbasis=FEMbasis, 
                             lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')


output_CPP.4$solution$beta

#### Test 2.5: Newton_fd method with exact GCV, default initial lambda and tolerance
output_CPP.5<-smooth.FEM.new(observations=data, locations = projected_locations, 
                             covariates = cbind(cov1,cov2),
                             FEMbasis=FEMbasis, 
                             lambda.selection.criterion='newton_fd', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

output_CPP.5$solution$beta

#### Test 2.6: Newton_fd method with stochastic GCV, default initial lambda and tolerance
output_CPP.6<-smooth.FEM(observations=data, locations = projected_locations, 
                             covariates = cbind(cov1,cov2),
                             FEMbasis=FEMbasis, 
                             lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')

output_CPP.6$solution$beta

