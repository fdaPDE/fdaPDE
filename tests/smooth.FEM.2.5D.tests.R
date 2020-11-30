##########################################
############## TEST SCRIPT ###############
##########################################

library(fdaPDE)

####### 2.5D ########

#### Test 1: C-shaped domain ####
#            locations = nodes 
#            no covariates
#            no BC
#            order FE = 1
library(fdaPDE)
rm(list=ls())
graphics.off()

data(horseshoe2.5D)
mesh = horseshoe2.5D

nnodes=dim(mesh$nodes)[1]

FEMbasis=create.FEM.basis(mesh)

# Exact solution (pointwise at nodes)
sol_exact=fs.test.3D(mesh$nodes[,1],mesh$nodes[,3],mesh$nodes[,2])
plot(FEM(sol_exact, FEMbasis))

# Add error to simulate data
set.seed(7893475)
ran=range(sol_exact)
data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^seq(-2,0.5,by=0.25)

#### Test 1.1: Without GCV
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda)
plot(output_CPP$fit.FEM)

#### Test 1.2: grid with exact GCV
#it takes a lot of time
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

#### Test 1.3: grid with stochastic GCV
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

#### Test 1.4: Newton_fd method with stochastic GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')

plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))



#### Test 2: sphere domain ####
#            locations != nodes 
#            with covariates
#            no BC
#            order FE = 1
library(fdaPDE)
rm(list=ls())
graphics.off()

data(sphere2.5D)
mesh = sphere2.5D

FEMbasis=create.FEM.basis(mesh)

# Test function 
f = function(x, y, z){
  phi = (asin(y/sqrt(x^2+y^2)))
  theta = acos(z/sqrt(x^2+y^2+z^2))
  # rho = 1
  
  sin(4*(1/2*sin(theta)*exp(-sin(theta)^2)+1)*theta)*cos(2*(1/2*cos(phi)*exp(-cos(phi)^2)+1)*phi)
}

# Exact solution (pointwise at nodes)
sol_exact=f(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3])
plot(FEM(sol_exact, FEMbasis))

# Generate data locations on the sphere
set.seed(598944)
ndata = 500
locations = matrix(rnorm(ndata*3), ncol = 3)
locations = locations/sqrt(locations[,1]^2 + locations[,2]^2 + locations[,3]^2)

# Generate covariate and data
cov1 = runif(ndata, min = -1, max = 1)
DatiEsatti = f(locations[,1],locations[,2], locations[,3]) + cov1

# Add error to simulate data
set.seed(7893475)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(ndata, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Project locations on the mesh
projected_locations = projection.points.2.5D(mesh, locations)

# Set smoothing parameter
lambda = 10^seq(-4,-2,by=0.25)

#### Test 2.1: Without GCV
output_CPP<-smooth.FEM(observations=data, 
                       locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, 
                       lambda=lambda[1])
plot(output_CPP$fit.FEM)
output_CPP$solution$beta

#### Test 2.2: grid with exact GCV
output_CPP<-smooth.FEM(observations=data, locations = projected_locations,
                       covariates = cov1,
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

#### Test 2.3: grid with stochastic GCV
output_CPP<-smooth.FEM(observations=data, locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

#### Test 2.4: Newton method with exact GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

#### Test 2.5: Newton_fd method with exact GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton_fd', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

#### Test 2.6: Newton_fd method with stochastic GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')

plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta
