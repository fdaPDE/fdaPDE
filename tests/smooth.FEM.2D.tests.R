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
library(fdaPDE)
rm(list=ls())
graphics.off()

x = seq(0,1, length.out = 60)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)
plot(mesh)

nnodes = dim(mesh$nodes)[1]

FEMbasis = create.FEM.basis(mesh)

# Test function
f = function(x, y, z = 1)
{
  coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
  sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
}

# Exact solution (pointwise at nodes)
sol_exact = f(mesh$nodes[,1], mesh$nodes[,2])
image(FEM(sol_exact, FEMbasis))

# Add error to simulate data
set.seed(7893475)
ran = range(sol_exact)
data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^seq(-6,-3,by=0.25)

#### Test 1.1: Without GCV
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda)
image(output_CPP$fit.FEM)

#### Test 1.2: grid with exact GCV
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda, lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

#### Test 1.3: grid with stochastic GCV
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda, lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV', DOF.stochastic.realizations = 1000)
plot(log10(lambda), output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

### Test 1.4: Newton exact method with exact GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

### Test 1.5: Newton_fd method with  exact GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda.selection.criterion='newton_fd', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

### Test 1.6: Newton_fd method with stochastic GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')

#### Test 2: c-shaped domain ####
#            locations != nodes
#            laplacian
#            with covariates
#            no BC
#            order FE = 1
library(fdaPDE)
rm(list=ls())
graphics.off()

data(horseshoe2D)

mesh = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
locations = refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
plot(mesh, pch = ".")
points(locations, pch = 16, col = "red")

FEMbasis = create.FEM.basis(mesh)

ndata = nrow(locations)

# Create covariates
set.seed(509875)
cov1 = rnorm(ndata, mean = 1, sd = 2)
cov2 = sin(locations[,1])

# Exact solution (pointwise at nodes)
DatiEsatti=fs.test(locations[,1], locations[,2]) + 2*cov1 - cov2

for(ind in 1:100)
{
  points(locations[which(round((DatiEsatti-min(DatiEsatti))/(max(DatiEsatti)-min(DatiEsatti))*100)==ind),1],
         locations[which(round((DatiEsatti-min(DatiEsatti))/(max(DatiEsatti)-min(DatiEsatti))*100)==ind),2],
         col=heat.colors(100)[ind], pch=16)
}

# Add error to simulate data
set.seed(543663)
ran = range(DatiEsatti)
data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^seq(-3,3,by=0.25)

#### Test 2.1: Without GCV
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda)

#### Test 2.2: grid with exact GCV
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

#### Test 2.3: grid with stochastic GCV
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

### Test 2.4: Newton exact method with exact  GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis,
                       lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

### Test 2.5: Newton_fd method with exact GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton_fd', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

### Test 2.6: Newton_fd method with stochastic GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')

image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

### Test 2.7: Inference on beta, hypothesis testing, Wald, Speckman, ESF, and enhanced ESF p_values
inf_obj<-inferenceDataObjectBuilder(test = c("sim", rep("oat",3)), dim = 2, n_cov = 2, type = c("w", "s", "esf", "enh-esf"), beta0 = c(2,-1))
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis,
                       lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', 
                       inference.data.object=inf_obj)

output_CPP$inference$beta$p_values

### Test 2.8: Inference on beta, hypothesis testing and confidence intervals of linear combinations, Wald and Speckman p_values  
inf_obj<-inferenceDataObjectBuilder(test = "oat", interval = "oat", dim = 2, n_cov = 2, type = c("w", "s"), coeff = matrix(data = c(1,1,1,-1), nrow = 2, byrow = T))
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', 
                       inference.data.object=inf_obj)

output_CPP$inference$beta$p_values

### Test 2.9: Inference on f, hypothesis testing, equality to f0, Wald, Sign-flip and ESF p_values
inf_obj<-inferenceDataObjectBuilder(test = "sim", dim = 2, n_cov = 2, type = c("w","sf","esf"), component = "nonparametric", f0 = fs.test)
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis,
                       lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', 
                       inference.data.object=inf_obj)

output_CPP$inference$f$p_values

### Test 2.10: Inference on f, hypothesis testing and confidence intervals, Wald with new locations  
mesh_loc = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
mesh_loc_ref<-refine.mesh.2D(mesh_loc, maximum_area = 0.05)

new_locs <- mesh_loc_ref$nodes[which(mesh_loc_ref$nodesmarkers!=1),]
new_locs[,1] <- new_locs[,1] + 0.2

inf_obj<-inferenceDataObjectBuilder(test = "sim", interval = "oat", dim = 2, n_cov = 2, type = "w", component = "nonparametric", locations = new_locs)
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis,
                       lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', 
                       inference.data.object=inf_obj)

output_CPP$inference$f$p_values
output_CPP$inference$f$CI

### Test 2.11: Inference on both beta and f, hypothesis testing: all implementations p_values 
inf_obj<-inferenceDataObjectBuilder(test = c("sim", "oat", "sim", "sim", "oat"), dim = 2, n_cov = 2, type = c("w","s","sf","esf","enh-esf"), 
                                    component = c("both", "parametric", "nonparametric", "both", "parametric"), f0 = fs.test, beta0 = c(2,-1))

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', 
                       inference.data.object=inf_obj)

output_CPP$inference$f$p_values
output_CPP$inference$beta$p_values


#### Test 3: square domain ####
#            locations in nodes
#            PDE
#            no covariates
#            no BC
#            order FE = 2
library(fdaPDE)
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
output_CPP<-smooth.FEM(observations=data, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda, 
                       PDE_parameters=PDE_parameters)

#### Test 3.2: grid with exact GCV
output_CPP<-smooth.FEM(observations=data, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       PDE_parameters=PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

#### Test 3.3: grid with stochastic GCV
output_CPP<-smooth.FEM(observations=data, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda, 
                       PDE_parameters=PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))


### Test 3.4: Newton exact method with exact GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, 
                       FEMbasis=FEMbasis, 
                       PDE_parameters=PDE_parameters,
                       lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

### Test 3.5: Newton_fd method with exact GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, 
                       FEMbasis=FEMbasis, 
                       PDE_parameters=PDE_parameters,
                       lambda.selection.criterion='newton_fd', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

### Test 3.6: Newton_fd method with stochastic GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, 
                       FEMbasis=FEMbasis, 
                       PDE_parameters=PDE_parameters,
                       lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')

image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

### Test 3.7: Test on inference: generating covariates in order to perform inference
set.seed(509875)
cov1 = rnorm(ndati, mean = 1, sd = 2)
cov2 = sin(locations[,1])

DatiEsatti = DatiEsatti + 2.5*cov1+1.5*cov2

set.seed(7893475)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(ndati, mean=0, sd=0.05*abs(ran[2]-ran[1]))

### Wald and Speckman CI and p_values for beta; Wald, Sign-flip and ESF p_values for f
inf_obj<-inferenceDataObjectBuilder(test = c(rep('oat',3),rep('sim',3)), interval = c('sim','oat','bonf',rep('none',3)),
                                    component=c(rep('parametric',3), rep('nonparametric',3)),
                                    dim = 2, n_cov = 2, type = c("w","w","s","w","sf","esf"),f0 = fs.test, beta0 = c(2.5,1.5))
output_CPP<-smooth.FEM(observations=data,
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, 
                       PDE_parameters=PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', lambda = lambda,
                       inference.data.object = inf_obj)

#### Test 4: quasicircular domain ####
#            areal observations
#            PDE space varying
#            no covariates
#            with BC
#            order FE = 1
library(fdaPDE)
rm(list=ls())
graphics.off()

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


#### Test 4.1.2: Forcing term = 0  grid without GCV 
output_CPP<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters)
plot(output_CPP$fit.FEM)


#### Test 4.1.2: Forcing term = 0  grid with exact  GCV 
output_CPP<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(output_CPP$fit.FEM)


#### Test 4.1.3: Forcing term = 0  grid with stochatic  GCV 
output_CPP<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
plot(output_CPP$fit.FEM)


#### Test 4.2.1: Forcing term != 0 without GCV 
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
                       BC = BC, 
                       PDE_parameters = PDE_parameters)
plot(output_CPP$fit.FEM)


#### Test 4.2.2: Forcing term != 0 grid with exact GCV 
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
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(output_CPP$fit.FEM)

#### Test 4.2.3: Forcing term != 0 grid with stochastic GCV 
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
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
plot(output_CPP$fit.FEM)


#### Test 4.2.4: Forcing term != 0 Newton finite differences with exact GCV 
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
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='newton_fd', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(output_CPP$fit.FEM)


#### Test 4.2.5: Forcing term != 0 Newton finite differences with stochastic GCV 
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
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
plot(output_CPP$fit.FEM)


#### Test 4.3.1: BC != 0      without GCV 
# Add a constat to the data to change true BC
data_backup=data #save a copy of original data
data = data + 5

# Set new value for the BC
BC$BC_values = rep(5,length(BC$BC_indices))

output_CPP<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters)
data=data_backup #restore original data for next tests
plot(output_CPP$fit.FEM)


#### Test 4.3.2: BC != 0      grid with exact GCV 
# Add a constat to the data to change true BC
data_backup=data #save a copy of original data
data = data + 5

# Set new value for the BC
BC$BC_values = rep(5,length(BC$BC_indices))

output_CPP<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
data=data_backup #restore original data for next tests
plot(output_CPP$fit.FEM)


#### Test 4.3.3: BC != 0      grid with stochastic GCV 
# Add a constat to the data to change true BC
data_backup=data #save a copy of original data
data = data + 5

# Set new value for the BC
BC$BC_values = rep(5,length(BC$BC_indices))

output_CPP<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
data=data_backup #restore original data for next tests
plot(output_CPP$fit.FEM)


#### Test 4.3.4: BC != 0     Newton Finite differences with exact GCV 
# Add a constat to the data to change true BC
data_backup=data #save a copy of original data
data = data + 5

# Set new value for the BC
BC$BC_values = rep(5,length(BC$BC_indices))

output_CPP<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix, 
                       FEMbasis=FEMbasis, 
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='newton_fd', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
data=data_backup #restore original data for next tests
plot(output_CPP$fit.FEM)

#### Test 4.3.5: BC != 0     Newton Finite differences with stochastic GCV 
# Add a constat to the data to change true BC
data_backup=data #save a copy of original data
data = data + 5

# Set new value for the BC
BC$BC_values = rep(5,length(BC$BC_indices))

output_CPP<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix, 
                       FEMbasis=FEMbasis, 
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
data=data_backup #restore original data for next tests
plot(output_CPP$fit.FEM)

#### Test on inference: creating covariates and special setting
# Create covariates
set.seed(509875)
cov1 = rnorm(length(DatiEsatti), mean = 2, sd = 1)

# Add error to simulate data
set.seed(5839745)
data = DatiEsatti + 3*cov1 
data = data + rnorm(length(DatiEsatti), sd = 0.1)
observations=data

# Inference obj: separate tests for beta and f. Wald, Speckman, ESF, Enhanced ESF p_values for beta; Wald, Sign-Flip and ESF p_values for f
inf_beta <- inferenceDataObjectBuilder (test='oat', interval='oat',  dim=2, n_cov=1, type=c('w', 's', 'esf', 'enh-esf'), beta0 = 3, component='parametric', n_flip=1000, f_var=T)
inf_f <- inferenceDataObjectBuilder(test = 'sim', type=c('w', 'sf', 'esf'), component = 'nonparametric', dim=2, n_cov=1)

#### Test 4.4: Forcing term = 0  grid with exact  GCV 
output_CPP<-smooth.FEM(observations=data, 
                       covariates = cov1,
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV',
                       inference.data.object = inf_beta)

output_CPP<-smooth.FEM(observations=data, 
                       covariates = cov1,
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV',
                       inference.data.object = inf_f)
