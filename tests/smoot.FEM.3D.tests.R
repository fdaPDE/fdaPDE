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
GCVFLAG=FALSE
output_CPP<-smooth.FEM(locations = locations, 
                       observations=data, FEMbasis=FEMbasis, lambda=lambda[1],
                       GCV=GCVFLAG)
plot(output_CPP$fit.FEM)

#### Test 1.2: With exact GCV
GCVFLAG=TRUE
GCVMETHODFLAG='Exact'
# it takes lot of time
output_CPP<-smooth.FEM(locations = locations, 
                       observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
plot(log10(lambda), output_CPP$GCV)
plot(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))

#### Test 1.3: With stochastic GCV
GCVFLAG=TRUE
GCVMETHODFLAG='Stochastic'
lambda= 10^seq(-9,-5,by=0.5)
output_CPP<-smooth.FEM(locations = locations, 
                       observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
plot(log10(lambda), output_CPP$GCV)
plot(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))


#### Test 2: c-shaped domain ####
#            locations = nodes
#            with covariates
#            no BC
#            order FE = 1
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
GCVFLAG=FALSE
output_CPP<-smooth.FEM(observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda[1],
                       GCV=GCVFLAG)
# plot(output_CPP$fit.FEM)

#### Test 2.2: With exact GCV
GCVFLAG=TRUE
GCVMETHODFLAG='Exact'
output_CPP<-smooth.FEM(observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
plot(log10(lambda), output_CPP$GCV)
plot(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))

output_CPP$beta[,which.min(output_CPP$GCV)]

#### Test 2.3: With stochastic GCV
GCVFLAG=TRUE
GCVMETHODFLAG='Stochastic'
output_CPP<-smooth.FEM(observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
plot(log10(lambda), output_CPP$GCV)
plot(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))

output_CPP$beta[,which.min(output_CPP$GCV)]


