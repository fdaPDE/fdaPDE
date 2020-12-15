##########################################
############## TEST SCRIPT ###############
##########################################

library(fdaPDE)

####### 3D ########

#### Test 1: square domain ####
#            locations != nodes 
#            no covariates
#            no BC
#            order FE = 1
library(fdaPDE)
rm(list=ls())
graphics.off()

data(cube3D)
mesh = cube3D
plot(mesh)

nnodes=dim(mesh$nodes)[1]

FEMbasis=create.FEM.basis(mesh)

# Test function
f = function(x, y, z){
  coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
  sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
}

# Exact solution (pointwise at nodes)
sol_exact=f(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3])
plot(FEM(sol_exact, FEMbasis))

# Set locations and add error to simulate data
x = seq(0,1, length.out = 11)
y = x
z = x
locations = expand.grid(x,y,z)
DatiEsatti = f(locations[,1], locations[,2], locations[,3])
set.seed(7893475)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda= 10^seq(-9,-7,by=1) # for simulation choose a finer grid 

#### Test 1.1: Without GCV
output_CPP<-smooth.FEM(locations = locations, 
                       observations=data, FEMbasis=FEMbasis, lambda=lambda[1])
plot(output_CPP$fit.FEM)

#### Test 1.2: grid with exact GCV
# it takes lot of time
output_CPP<-smooth.FEM(locations = locations, 
                       observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

#### Test 1.3: grid with stochastic GCV
lambda= 10^seq(-9,-5,by=0.5)
output_CPP<-smooth.FEM(locations = locations, 
                       observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

#No exact GCV computation for Newton's method, it takes too much time
#### Test 1.4: Newton_fd method with stochastic GCV, default initial lambda and tolerance
#monotone increasing GCV function
output_CPP<-smooth.FEM(locations = locations, 
                       observations=data, FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')

plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))


#### Test 2: c-shaped domain ####
#            locations = nodes
#            with covariates
#            no BC
#            order FE = 1
library(fdaPDE)
rm(list=ls())
graphics.off()

data(horseshoe3D)

mesh = horseshoe3D

FEMbasis=create.FEM.basis(mesh)

ndata = FEMbasis$nbasis

# Create covariates
set.seed(509875)
cov1=rnorm(ndata, mean = 1, sd = 2)
cov2=sin(mesh$nodes[,1])

# Exact solution (pointwise at nodes)
DatiEsatti=fs.test.3D(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]) + 2*cov1 - cov2
plot(FEM(DatiEsatti, FEMbasis))

# Add error to simulate data
set.seed(543663)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda= 10^seq(-3,3,by=0.25)

#### Test 2.1: Without GCV
output_CPP<-smooth.FEM(observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda[1])
# plot(output_CPP$fit.FEM)

#### Test 2.2: grid with exact GCV
output_CPP<-smooth.FEM(observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

#### Test 2.3: grid with stochastic GCV
output_CPP<-smooth.FEM(observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

#### Test 2.4: Newton method with exact GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis,
                       lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

#### Test 2.5: Newton_fd method with exact GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis,
                       lambda.selection.criterion='newton_fd', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

#### Test 2.6: Newton_fd method with stochastic GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis,
                       lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')

plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta


### Second order FEM in 3D ################################
# Compare the accuracy and robustness of first 
# and second order methods in 3D settings.
# Order 2 in 3D is a new functionality.

rm(list=ls())

# Function to generate random points in a sphere
rsphere <- function(n, r = 1.0, surface_only = FALSE) {
  phi       <- runif(n, 0.0, 2.0 * pi)
  cos_theta <- runif(n, -1.0, 1.0)
  sin_theta <- sqrt((1.0-cos_theta)*(1.0+cos_theta))
  radius <- r
  if (surface_only == FALSE) {
    radius <- r * runif(n, 0.0, 1.0)^(1.0/3.0)
  }
  
  x <- radius * sin_theta * cos(phi)
  y <- radius * sin_theta * sin(phi)
  z <- radius * cos_theta
  
  cbind(x, y, z)
}

# Percentage of noise standard deviation on data range
# Try also different values, e.g 0.01 or 1
noisepercent<-0.25

data("sphere3Ddata")

# Build spherical mesh of order 2 
mesh_sphere<-create.mesh.3D(sphere3Ddata$nodes,
                              sphere3Ddata$tetrahedrons,
                              order=2)
  
FEMbasis <- create.FEM.basis(mesh_sphere)
  
# Set smoothing parameters
lambda= 10^seq(-6,1,by=.25)
  
set.seed(5847947)
  
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)
  
nnodes<-nrow(mesh_sphere$nodes)
  
# Evaluate exact solution on mesh nodes
func_evaluation = a1* sin(2*pi*mesh_sphere$nodes[,1]) +  a2* sin(2*pi*mesh_sphere$nodes[,2]) +  a3*sin(2*pi*mesh_sphere$nodes[,3]) + 1
  
ran=range(func_evaluation)
  
exact_sol=func_evaluation
  
# Plot exact solution
plot(FEM(exact_sol, FEMbasis))
title3d(main="Exact solution, order 2", 
          col="black")
  
# Add noise to exact solution to generate data for the simulation
data= exact_sol + rnorm(nnodes,
                          mean=0,
                          sd=noisepercent * (ran[2]-ran[1]))
  
# Compute the solution for each lambda
output_CPP <- smooth.FEM(observations = data, 
                           FEMbasis = FEMbasis, 
                           lambda = lambda)

### PDE penalization in 3D ################################
# Include PDE parameters and use a penalization term with
# anisotropic diffusion and possibly advection and reaction
# terms. All this is a new functionality.

rm(list=ls())

# Function to generate random points in a sphere
rsphere <- function(n, r = 1.0, surface_only = FALSE) {
  phi       <- runif(n, 0.0, 2.0 * pi)
  cos_theta <- runif(n, -1.0, 1.0)
  sin_theta <- sqrt((1.0-cos_theta)*(1.0+cos_theta))
  radius <- r
  if (surface_only == FALSE) {
    radius <- r * runif(n, 0.0, 1.0)^(1.0/3.0)
  }
  
  x <- radius * sin_theta * cos(phi)
  y <- radius * sin_theta * sin(phi)
  z <- radius * cos_theta
  
  cbind(x, y, z)
}


# Build mesh: Sphere
data("sphere3Ddata")
mesh_sphere<-create.mesh.3D(sphere3Ddata$nodes,
                            sphere3Ddata$tetrahedrons,
                            order=1)

FEMbasis <- create.FEM.basis(mesh_sphere)

set.seed(5847947)

# Exact test function
nnodes = nrow(mesh_sphere$nodes)

# Set smoothing parameter
lambda=10^seq(-6, 1, by=.25)

# Set PDE parameters (in this case they are constant)
PDE_parameters_anys = list(K = diag(c(1,.5,1)), b = c(0,0,0), c = -4*pi^2)

# Evaluate exact solution on mesh nodes
exact_sol =  sin(2*pi*mesh_sphere$nodes[,1]) +  2 * sin(2*pi*mesh_sphere$nodes[,2]) +  sin(2*pi*mesh_sphere$nodes[,3])

# Plot exact solution
# plot(FEM(exact_sol,FEMbasis))
# title3d(main="Exact solution", 
# col="black")

# Add noise to generate data - 10% level of noise
data=exact_sol + rnorm(nrow(mesh_sphere$nodes), mean=0, sd=0.10*diff(range(exact_sol)))
  
# Compute the solution for each lambda
output_CPP <- smooth.FEM(observations = data, PDE_parameters = PDE_parameters_anys,
                           FEMbasis = FEMbasis, lambda = lambda)
  
