##########################################
############## TEST SCRIPT ###############
##########################################

library(fdaPDE)

####### 2D ########

#### Test 1: square domain ####
#            locations = nodes 
#            laplacian
#            no covariates
#            no BC
#            order FE = 1
rm(list=ls())
graphics.off()

x = seq(0,1, length.out = 41)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)
plot(mesh)

nnodes=dim(mesh$nodes)[1]

FEMbasis=create.FEM.basis(mesh)

# Test function
f = function(x, y, z = 1){
  coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
  sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
}

# Exact solution (pointwise at nodes)
sol_exact=f(mesh$nodes[,1], mesh$nodes[,2])
image(FEM(sol_exact, FEMbasis))

# Add error to simulate data
set.seed(7893475)
ran=range(sol_exact)
data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda= 10^seq(-6,-3,by=0.25)

#### Test 1.1: Without GCV
GCVFLAG=FALSE
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       GCV=GCVFLAG)
image(output_CPP$fit.FEM)

#### Test 1.2: With exact GCV
GCVFLAG=TRUE
GCVMETHODFLAG='Exact'
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
plot(log10(lambda), output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))

#### Test 1.3: With stochastic GCV
GCVFLAG=TRUE
GCVMETHODFLAG='Stochastic'
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
plot(log10(lambda), output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))


#### Test 2: c-shaped domain ####
#            locations != nodes
#            laplacian
#            with covariates
#            no BC
#            order FE = 1
rm(list=ls())
graphics.off()

data(horseshoe2D)

mesh=create.mesh.2D(nodes=horseshoe2D$boundary_nodes, 
                    segments = horseshoe2D$boundary_segments)
locations=refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
mesh=refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
plot(mesh, pch = ".")
points(locations, pch = 16, col = "red")

FEMbasis=create.FEM.basis(mesh)

ndata = nrow(locations)

# Create covariates
set.seed(509875)
cov1=rnorm(ndata, mean = 1, sd = 2)
cov2=sin(locations[,1])

# Exact solution (pointwise at nodes)
DatiEsatti=fs.test(locations[,1], locations[,2]) + 2*cov1 - cov2

for(ind in 1:100){
  points(locations[which(round((DatiEsatti-min(DatiEsatti))/(max(DatiEsatti)-min(DatiEsatti))*100)==ind),1],locations[which(round((DatiEsatti-min(DatiEsatti))/(max(DatiEsatti)-min(DatiEsatti))*100)==ind),2],
         col=heat.colors(100)[ind], pch=16)
}

# Add error to simulate data
set.seed(543663)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda= 10^seq(-3,3,by=0.25)

#### Test 2.1: Without GCV
GCVFLAG=FALSE
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       GCV=GCVFLAG)
# image(output_CPP$fit.FEM)

#### Test 2.2: With exact GCV
GCVFLAG=TRUE
GCVMETHODFLAG='Exact'
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
plot(log10(lambda), output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))

output_CPP$beta[,which.min(output_CPP$GCV)]

#### Test 2.3: With stochastic GCV
GCVFLAG=TRUE
GCVMETHODFLAG='Stochastic'
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
plot(log10(lambda), output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))

output_CPP$beta[,which.min(output_CPP$GCV)]


#### Test 3: square domain ####
#            locations in nodes
#            PDE
#            no covariates
#            no BC
#            order FE = 2
rm(list=ls())
graphics.off()

x = seq(0,1, length.out = 11)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations, order = 2)
plot(mesh)

FEMbasis=create.FEM.basis(mesh)

# Test function
a1=1
a2=4
z<-function(p){  
  a1*sin(2*pi*p[,1])*cos(2*pi*p[,2])+a2*sin(3*pi*p[,1])}

# Exact solution (pointwise at nodes)
sol_exact=z(mesh$nodes)
image(FEM(sol_exact, FEMbasis))

DatiEsatti=z(locations)
ndati = length(DatiEsatti)

# Add error to simulate data
set.seed(7893475)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(ndati, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda= 10^seq(-6,-3,by=0.25)

# Set PDE parameters
PDE_parameters = list(K = matrix(c(1,0,0,4), nrow = 2), b = c(0,0), c = 0)

#### Test 3.1: Without GCV
GCVFLAG=FALSE
output_CPP<-smooth.FEM(observations=data, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda[7], 
                       PDE_parameters=PDE_parameters,
                       GCV=GCVFLAG)
image(output_CPP$fit.FEM)

#### Test 1.2: With exact GCV
GCVFLAG=TRUE
GCVMETHODFLAG='Exact'
output_CPP<-smooth.FEM(observations=data, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       PDE_parameters=PDE_parameters,
                       GCV=GCVFLAG,
                       GCVmethod = GCVMETHODFLAG)
plot(log10(lambda), output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))

#### Test 1.3: With stochastic GCV
GCVFLAG=TRUE
GCVMETHODFLAG='Stochastic'
output_CPP<-smooth.FEM(observations=data, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda, 
                       PDE_parameters=PDE_parameters,
                       GCV=GCVFLAG, 
                       GCVmethod = GCVMETHODFLAG)
plot(log10(lambda), output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))




#### Test 4: quasicircular domain ####
#            areal observations
#            PDE space varying
#            no covariates
#            with BC
#            order FE = 1

rm(list=ls())

data(quasicircle2Dareal)
mesh = quasicircle2Dareal$mesh
incidence_matrix = quasicircle2Dareal$incidence_matrix
DatiEsatti = quasicircle2Dareal$data

plot(mesh)
FEMbasis = create.FEM.basis(mesh)

# Add error to simulate data
set.seed(5839745)
data = DatiEsatti + rnorm(length(DatiEsatti), sd = 0.05*(max(DatiEsatti)-min(DatiEsatti)))

# Set smoothing parameter
lambda = 10^-3

# Set BC (equal to zero)
BC = NULL
BC$BC_indices = which(mesh$nodesmarkers == 1)
BC$BC_values = rep(0,length(BC$BC_indices))

# Set sv-PDE parameters
R = 2.8
K1 = 0.1
K2 = 0.2
beta = 0.5

K_func<-function(points)
{
  output = array(0, c(2, 2, nrow(points)))
  for (i in 1:nrow(points))
    output[,,i] = 10*rbind(c(points[i,2]^2 + K1*points[i,1]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2),
                             (K1-1)*points[i,1]*points[i,2]),
                           c((K1-1)*points[i,1]*points[i,2],
                             points[i,1]^2 + K1*points[i,2]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2)))
  output
}

b_func<-function(points)
{
  output = array(0, c(2, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = 10*beta*c(points[i,1],points[i,2])
  output
}

c_func<-function(points)
{
  rep(c(0), nrow(points))
}

u_func<-function(points)
{
  rep(c(0), nrow(points))
}
PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)

#### Test 4.1: Forcing term = 0
GCVFLAG=FALSE
output_CPP<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       GCV=GCVFLAG,
                       BC = BC, 
                       PDE_parameters = PDE_parameters)
plot(output_CPP$fit.FEM)

#### Test 4.2: Forcing term != 0
# forcing function != 0
u_func<-function(points)
{
  output = array(0, c(1, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = -ifelse((points[i,1]^2+points[i,2]^2)<1,100,0)
  output
}
# plot the forcing funcion
xgrid=seq(from=-3,to=3,by=0.1)
ygrid=seq(from=-3,to=3,by=0.1)
image(xgrid,ygrid,matrix(u_func(expand.grid(xgrid,ygrid)),nrow=length(xgrid),ncol=length(ygrid),byrow=FALSE),main='forcing function',asp=1)
lines(mesh$nodes[1:50,])

PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)

output_CPP<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       GCV=GCVFLAG,
                       BC = BC, 
                       PDE_parameters = PDE_parameters)
plot(output_CPP$fit.FEM)

#### Test 4.3: BC != 0

# Add a constat to the data to change true BC
data = data + 5

# Set new value for the BC
BC$BC_values = rep(5,length(BC$BC_indices))

output_CPP<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       GCV=GCVFLAG,
                       BC = BC, 
                       PDE_parameters = PDE_parameters)
plot(output_CPP$fit.FEM)
