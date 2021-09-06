##########################################
############## TEST SCRIPT ###############
######## for space-time problems #########
##########################################

library(fdaPDE)

####### 2D ########

#### Test 1: square domain ####
#            locations = nodes 
#            laplacian
#            no covariates
#            no BC
#            order FE = 1
#            separate penalizations
rm(list=ls())
graphics.off()

x = seq(0,1, length.out = 21)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)
plot(mesh)

nnodes = dim(mesh$nodes)[1]

FEMbasis = create.FEM.basis(mesh)

# Test function
f = function(x, y, z)
{
  coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
  sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
}

NumTimePoints=11
TimePoints=seq(0,2,length.out = NumTimePoints)

SpacePoints=mesh$nodes
SpaceTimePoints=cbind(rep(mesh$nodes[,1],NumTimePoints),rep(mesh$nodes[,2],NumTimePoints),rep(TimePoints,each=nrow(mesh$nodes)))

# Exact solution (pointwise at nodes)
sol_exact=f(SpaceTimePoints[,1],SpaceTimePoints[,2],SpaceTimePoints[,3])
# plot(FEM.time(sol_exact, TimePoints, FEMbasis, FLAG_PARABOLIC = TRUE))

# Add error to simulate data
set.seed(7893475)
ran = range(sol_exact)
data = sol_exact + rnorm(length(sol_exact), mean=0, sd=0.05*abs(ran[2]-ran[1]))
observations=matrix(data,nrow(SpacePoints),NumTimePoints)

# Set smoothing parameter
lambdaS = 1e-2
lambdaT = 1e-2

#### Test 1.1: Without GCV
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                            FEMbasis=FEMbasis, lambdaS=lambdaS,
                            lambdaT=lambdaT)
plot(output_CPP$fit.FEM)

# RMSE evaluation on a fine grid
xeval=runif(1000,0,1)
yeval=runif(1000,0,1)
teval=runif(1000,0,2)
sol_eval=eval.FEM.time(output_CPP$fit.FEM.time,space.time.locations = cbind(teval,xeval,yeval))
sol_exact=f(xeval,yeval,teval)
RMSE<-function(f,g) sqrt(mean((f-g)^2))
RMSE(sol_eval,sol_exact)

#### Test 1.2: grid with exact GCV
# We test one value for speeding the tests, this is not what should be done in practice!
lambdaS = 1e-3
lambdaT = 1e-3
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                            FEMbasis=FEMbasis, lambdaS=lambdaS,lambdaT=lambdaT,
                            lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
bestlambdas = which(output_CPP$GCV == min(output_CPP$GCV), arr.ind = TRUE)
image(FEM.time(output_CPP$fit.FEM$coeff[,1,1],FEMbasis = FEMbasis, time_mesh = TimePoints), t=1)

#### Test 1.3: grid with stochastic GCV
lambdaS = 10^seq(-7, -5, length.out=3)
lambdaT = 10^seq(-7, -5, length.out=3)
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, 
                            time_locations = TimePoints,
                            observations=observations, 
                            FEMbasis=FEMbasis, 
                            lambdaS=lambdaS,lambdaT=lambdaT,
                            lambda.selection.criterion='grid', 
                            DOF.evaluation='stochastic', 
                            lambda.selection.lossfunction='GCV')
bestlambdas = which(output_CPP$GCV == min(output_CPP$GCV), arr.ind = TRUE)
image(FEM.time(output_CPP$fit.FEM$coeff[,1,1],FEMbasis = FEMbasis, time_mesh = TimePoints), t=1)

#### Test 2: c-shaped domain ####
#            locations != nodes
#            laplacian
#            with covariates
#            no BC
#            order FE = 1
#            separate penalizations 
library(fdaPDE)
rm(list=ls())
graphics.off()

data(horseshoe2D)

mesh = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
locations = refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)

# plot spatial locations
plot(mesh, pch = ".")
points(locations, pch = 16, col = "red")

FEMbasis = create.FEM.basis(mesh)

f<-function(x,y,t)
{
  K <- (y/0.1*as.double((abs(y)<=0.1 & x>-0.5))+as.double((abs(y)>0.1 | x<=-0.5)))^2
  res=numeric(length =length(x))
  for(i in 1:length(x))
  {
    if(x[i]>=0 && y[i]>0)
      res[i]=cos(t[i])*(0.25*pi+x[i])+(y[i]-0.5)^2
    if(x[i]>=0 && y[i]<=0)
      res[i]=cos(2*t[i])*(-0.25*pi-x[i])+(-y[i]-0.5)^2
    if(x[i]<0 && y[i]>0)
      res[i]=cos(t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
    if(x[i]<0 && y[i]<=0)
      res[i]=cos(2*t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
  }
  res
}

NumTimeInstants=5
TimePoints=seq(0,pi,length.out =NumTimeInstants)

space_time_locations = cbind(rep(TimePoints,each=nrow(locations)),rep(locations[,1],NumTimeInstants),rep(locations[,2],NumTimeInstants))
sol_exact = f(space_time_locations[,2],space_time_locations[,3],space_time_locations[,1])

ndata = length(sol_exact)

# Create covariates
set.seed(509875)
cov1 = rnorm(ndata, mean = 1, sd = 2)

# Add error to simulate data
set.seed(7893475)
data = sol_exact + 2*cov1 
data = data + rnorm(ndata, mean = 0, sd =  0.05*diff(range(sol_exact)))
observations = matrix(data,nrow(locations),NumTimeInstants)

# Set smoothing parameter
lambdaS = 10^-2
lambdaT = 10^-2

#### Test 2.1: Without GCV
output_CPP<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                            observations=observations, 
                            covariates = cov1,
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT)

plot(output_CPP$fit.FEM)
output_CPP$beta

#### Test 2.2: exact GCV
lambdaS = 10^(-1:1)
lambdaT = 10^(-6:-4)
output_CPP<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                            observations=observations, 
                            covariates = cov1,
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT, 
                            lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

bestlambdas = which(output_CPP$GCV == min(output_CPP$GCV), arr.ind = TRUE)
plot(FEM.time(output_CPP$fit.FEM$coeff[,bestlambdas[1],bestlambdas[2]],FEMbasis = FEMbasis, time_mesh = TimePoints), t=1)

output_CPP$beta

#### Test 2.3: stochastic GCV
lambdaS = 10^(-2:0)
lambdaT = 10^(-3:-1)
output_CPP<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                            observations=observations, 
                            covariates = cov1,
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT, 
                            lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')

bestlambdas = which(output_CPP$GCV == min(output_CPP$GCV), arr.ind = TRUE)
plot(FEM.time(output_CPP$fit.FEM$coeff[,bestlambdas[1],bestlambdas[2]],FEMbasis = FEMbasis, time_mesh = TimePoints), t=1)

output_CPP$beta



#### Test 3: square domain ####
#            locations in nodes
#            PDE
#            no covariates
#            no BC
#            order FE = 2
#            separate penalizations
library(fdaPDE)
rm(list=ls())
graphics.off()

x = seq(0,1, length.out = 21)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)
plot(mesh)

nnodes = dim(mesh$nodes)[1]

FEMbasis = create.FEM.basis(mesh)

# Test function
f = function(x, y, z)
{
  a1=1
  a2=4
  (a1*sin(2*pi*x)*cos(2*pi*y)+a2*sin(3*pi*x))*cos(pi*z)
}

NumTimePoints=11
TimePoints=seq(0,2,length.out = NumTimePoints)

SpacePoints=mesh$nodes
SpaceTimePoints=cbind(rep(mesh$nodes[,1],NumTimePoints),rep(mesh$nodes[,2],NumTimePoints),rep(TimePoints,each=nrow(mesh$nodes)))

# Exact solution (pointwise at nodes)
sol_exact=f(SpaceTimePoints[,1],SpaceTimePoints[,2],SpaceTimePoints[,3])

# Add error to simulate data
set.seed(7893475)
ran = range(sol_exact)
data = sol_exact + rnorm(length(sol_exact), mean=0, sd=0.05*abs(ran[2]-ran[1]))
observations=matrix(data,nrow(SpacePoints),NumTimePoints)

# Set PDE parameters
PDE_parameters = list(K = matrix(c(1,0,0,4), nrow = 2), b = c(0,0), c = 0)

# Set smoothing parameter
lambdaS = 1e-4
lambdaT = 1e-4

#### Test 3.1: Without GCV
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                            PDE_parameters=PDE_parameters)
image(output_CPP$fit.FEM, t = 1)
image(FEM.time(sol_exact, time_mesh = TimePoints, FEMbasis = FEMbasis, FLAG_PARABOLIC = TRUE), t=1)

plot(output_CPP$fit.FEM, locations = c(0.5,0.5))
plot(FEM.time(sol_exact, time_mesh = TimePoints, FEMbasis = FEMbasis, FLAG_PARABOLIC = TRUE), locations = c(0.5,0.5),add=TRUE,col=1)

#### Test 3.2: exact GCV
lambdaS = 10^(-5:-3)
lambdaT = 10^(-5:-3)
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                            PDE_parameters=PDE_parameters,
                            lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

#### Test 3.3: stochastic GCV
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                            PDE_parameters=PDE_parameters,
                            lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')


#### Test 4: quasicircular domain ####
#            areal observations
#            PDE space varying
#            no covariates
#            with BC
#            order FE = 1
#            parabolic smoothing
library(fdaPDE)
rm(list=ls())
graphics.off()

data(quasicircle2Dareal)
mesh = quasicircle2Dareal$mesh
incidence_matrix = quasicircle2Dareal$incidence_matrix
DatiEsatti = quasicircle2Dareal$data

plot(mesh)
FEMbasis = create.FEM.basis(mesh)
time_mesh=seq(0,4,length.out = 11)

DatiEsatti=rep(DatiEsatti,length(time_mesh))*rep(exp(time_mesh),each=length(DatiEsatti))

# Add error to simulate data
set.seed(5839745)
data = DatiEsatti + rnorm(length(DatiEsatti), sd = 0.1)
observations=matrix(data,nrow(incidence_matrix),length(time_mesh))

# Set smoothing parameter
lambdaS = 10^-6
lambdaT = 10^-6

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


#### Test 4.1.2: Forcing term = 0, estimated IC, without GCV
output_CPP<-smooth.FEM.time(observations=observations,
                            incidence_matrix = incidence_matrix,
                            FEMbasis = FEMbasis, time_mesh = time_mesh,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            BC = BC,
                            PDE_parameters = PDE_parameters,
                            FLAG_PARABOLIC = TRUE,
                            DOF.evaluation = NULL)
plot(output_CPP$fit.FEM.time, t = 1)

#### Test 4.1.2: Forcing term = 0 exact  GCV
output_CPP<-smooth.FEM.time(observations=observations,
                            incidence_matrix = incidence_matrix, areal.data.avg = TRUE,
                            FEMbasis = FEMbasis, time_mesh = time_mesh,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            BC = BC,
                            PDE_parameters = PDE_parameters,
                            FLAG_PARABOLIC = TRUE,
                            lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

#### Test 4.1.3: Forcing term = 0  stochastic  GCV
output_CPP<-smooth.FEM.time(observations=observations,
                            incidence_matrix = incidence_matrix, areal.data.avg = TRUE,
                            FEMbasis = FEMbasis, time_mesh = time_mesh,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            BC = BC,
                            PDE_parameters = PDE_parameters,
                            FLAG_PARABOLIC = TRUE,
                            lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')

#### Test 4.2.1: Forcing term != 0 without GCV
# forcing function != 0
u_func<-function(points)
{
  output = array(0, c(1, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = -ifelse((points[i,1]^2+points[i,2]^2)<1,100,0)
  output
}

PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)

output_CPP<-smooth.FEM.time(observations=observations,
                            incidence_matrix = incidence_matrix,
                            FEMbasis = FEMbasis, time_mesh = time_mesh,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            BC = BC,
                            PDE_parameters = PDE_parameters,
                            FLAG_PARABOLIC = TRUE)
plot(output_CPP$fit.FEM.time, t = 1)

#### Test 4.2.2: Forcing term != 0 exact GCV
output_CPP<-smooth.FEM.time(observations=observations,
                            incidence_matrix = incidence_matrix, areal.data.avg = TRUE,
                            FEMbasis = FEMbasis, time_mesh = time_mesh,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            BC = BC,
                            PDE_parameters = PDE_parameters,
                            FLAG_PARABOLIC = TRUE,
                            lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

#### Test 4.2.3: Forcing term != 0 stochastic GCV
# forcing function != 0
output_CPP<-smooth.FEM.time(observations=observations,
                            incidence_matrix = incidence_matrix, areal.data.avg = TRUE,
                            FEMbasis = FEMbasis, time_mesh = time_mesh,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            BC = BC,
                            PDE_parameters = PDE_parameters,
                            FLAG_PARABOLIC = TRUE,
                            lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')


#### Test 4.3.1: BC != 0  without GCV
u_func<-function(points)
{
  rep(c(0), nrow(points))
}
PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)

# Add a constat to the data to change true BC
obs_backup=observations #save a copy of original data
observations = observations + 5

# Set new value for the BC
BC$BC_values = rep(5,length(BC$BC_indices))

output_CPP<-smooth.FEM.time(observations=observations,
                            incidence_matrix = incidence_matrix, areal.data.avg = TRUE,
                            FEMbasis = FEMbasis, time_mesh = time_mesh,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            BC = BC,
                            PDE_parameters = PDE_parameters,
                            FLAG_PARABOLIC = TRUE)
plot(output_CPP$fit.FEM, t = 1)


#### Test 4.3.2: BC != 0 exact GCV
output_CPP<-smooth.FEM.time(observations=observations,
                            incidence_matrix = incidence_matrix, areal.data.avg = TRUE,
                            FEMbasis = FEMbasis, time_mesh = time_mesh,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            BC = BC,
                            PDE_parameters = PDE_parameters,
                            FLAG_PARABOLIC = TRUE,
                            lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

#### Test 4.3.3: BC != 0 stochastic GCV
output_CPP<-smooth.FEM.time(observations=observations,
                            incidence_matrix = incidence_matrix, areal.data.avg = TRUE,
                            FEMbasis = FEMbasis, time_mesh = time_mesh,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            BC = BC,
                            PDE_parameters = PDE_parameters,
                            FLAG_PARABOLIC = TRUE,
                            lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')


####### 2.5D ########

#### hub pointwise (examples with and without covariates) ####
library(fdaPDE)
rm(list=ls())

data(hub2.5D)
mesh <- create.mesh.2.5D(nodes = hub2.5D$hub2.5D.nodes,triangles = hub2.5D$hub2.5D.triangles)
FEMbasis <- create.FEM.basis(mesh)

# Locations at nodes
nodesLocations=mesh$nodes

# Exact data - Locations at nodes
nnodes = nrow(mesh$nodes)
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)
TimeNodes = 0:4

locations = cbind(rep(TimeNodes,each=nnodes),rep(nodesLocations[,1],length(TimeNodes)),rep(nodesLocations[,2],length(TimeNodes)),rep(nodesLocations[,3],length(TimeNodes)))

func = function(x)
{
  (a1*sin(2*pi*x[,2])+a2*sin(2*pi*x[,3])+a3*sin(2*pi*x[,4])+1)*cos(x[,1])
}

func_evaluation = func(locations)
# Plot the exact solution
plot(FEM.time(coeff=array(func_evaluation,dim=c(length(func_evaluation),1,1)), FEMbasis = FEMbasis,time_mesh=TimeNodes,FLAG_PARABOLIC=T),TimeNodes)


lambdaS=10^seq(-9, -7, 0.5)
lambdaT=10^seq(-6, -4, 0.5)

lambdaS_par=10^seq(-4, -3, 0.25)
lambdaT_par=10^seq(1, 1.8, 0.2)

cov1=4*sin(2*pi*locations[,2])*cos(2*pi*locations[,3])
cov2=rnorm(nnodes*length(TimeNodes), mean=3, sd=0.1)*rep(exp(-TimeNodes/length(TimeNodes)),each=nnodes)
W=cbind(cov1,cov2)

# plot(FEM(coeff = cov1[1:nnodes], FEMbasis = FEMbasis))
# plot(FEM(coeff = cov2[1:nnodes], FEMbasis = FEMbasis))

# Fix betas
beta_exact=c(0.45,0.3)

ran=range(W%*%beta_exact + func_evaluation)
ran=range(func_evaluation)

# Plot exact solution
plot(FEM.time(coeff=array(W%*%beta_exact + func_evaluation, dim = c(nnodes*length(TimeNodes),1,1)),FEMbasis=FEMbasis,time_mesh = TimeNodes,FLAG_PARABOLIC = T),TimeNodes)

ran = range(func_evaluation)
data = func_evaluation +rnorm(nnodes*length(TimeNodes),mean=0,sd=0.05*(ran[2]-ran[1]))

ran = range(func_evaluation+ W%*%beta_exact)
datacov=func_evaluation+ W%*%beta_exact +rnorm(nnodes*length(TimeNodes),mean=0,sd=0.05*(ran[2]-ran[1]))

data = matrix(data,nrow(mesh$nodes),length(TimeNodes))
datacov = matrix(datacov,nrow(mesh$nodes),length(TimeNodes))

#########################################SEPARABLE####################################################
solSep = smooth.FEM.time(observations=data,
                         FEMbasis = FEMbasis, time_mesh = TimeNodes, time_locations = TimeNodes,
                         lambdaS = lambdaS, lambdaT = lambdaT,
                         FLAG_PARABOLIC = FALSE)

solSepCov = smooth.FEM.time(observations=datacov, covariates = W,
                         FEMbasis = FEMbasis, time_mesh = TimeNodes,
                         lambdaS = lambdaS, lambdaT = lambdaT,
                         FLAG_PARABOLIC = FALSE)


##########################################PARABOLIC####################################################
solPar = smooth.FEM.time(observations=data,
                         FEMbasis = FEMbasis, time_mesh = TimeNodes, time_locations = TimeNodes,
                         lambdaS = lambdaS, lambdaT = lambdaT,
                         FLAG_PARABOLIC = TRUE)

solParCov = smooth.FEM.time(observations=datacov[,2:length(TimeNodes)], covariates = W[(1+nrow(mesh$nodes)):(length(TimeNodes)*nrow(mesh$nodes)),],
                            FEMbasis = FEMbasis, time_mesh = TimeNodes,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            IC=func_evaluation[1:nrow(mesh$nodes)],
                            FLAG_PARABOLIC = TRUE)
  
# Example of RMSE computation
TimeNodesEval=seq(0,4,length.out = 9)
eval_locations = cbind(rep(TimeNodesEval,each=nnodes),rep(nodesLocations[,1],length(TimeNodesEval)),rep(nodesLocations[,2],length(TimeNodesEval)),rep(nodesLocations[,3],length(TimeNodesEval)))
sol_eval=eval.FEM.time(solSep$fit.FEM.time,locations = nodesLocations,time.instants = TimeNodesEval, lambdaS = solSep$bestlambda[1],lambdaT = solSep$bestlambda[2])
sol_exact = func(eval_locations)
RMSE<-function(f,g) sqrt(mean((f-g)^2))
RMSE(sol_eval,sol_exact)


#### hub areal (examples with and without covariates) ####
library(fdaPDE)
rm(list=ls())

data(hub2.5Dareal)

nodesLocations=mesh$nodes

nnodes = nrow(mesh$nodes)
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)
TimeNodes = 0:4
TimeNodesRMSE = seq(0,4,length.out = 15)

locations = cbind(rep(TimeNodesRMSE,each=nnodes),rep(nodesLocations[,1],length(TimeNodesRMSE)),rep(nodesLocations[,2],length(TimeNodesRMSE)),rep(nodesLocations[,3],length(TimeNodesRMSE)))

func = function(x)
{
  (a1*sin(2*pi*x[,2])+a2*sin(2*pi*x[,3])+a3*sin(2*pi*x[,4])+1)*cos(x[,1])
}

func_evaluation = func(locations)
FEMbasis=create.FEM.basis(mesh)

plot(FEM.time(coeff=array(func_evaluation,dim=c(length(func_evaluation),1,1)),time_mesh = TimeNodesRMSE,FEMbasis,FLAG_PARABOLIC = T),3)

sol_exact=func_evaluation

W_areal=cbind(rep(cov_areal,length(TimeNodes))*rep(exp(-TimeNodes/length(TimeNodes)),each=RDD_groups))

beta_exact=c(1)

lambdaS=10^seq(-9, -7, 0.5)
lambdaT=10^seq(-6, -4, 0.5)

lambdaS_par=10^seq(-5.2, -4.8, 0.1)
lambdaT_par=10^seq(1, 1.8, 0.2)

obs_areal = rep(obs_areal,length(TimeNodes))*rep(cos(TimeNodes),each=RDD_groups)

ran = range(obs_areal)
data = obs_areal +rnorm(RDD_groups*length(TimeNodes),mean=0,sd=0.02*(ran[2]-ran[1]))

ran = range(obs_areal + W_areal%*%beta_exact)
datacov=obs_areal + W_areal%*%beta_exact + rnorm(RDD_groups*length(TimeNodes),mean=0,sd=0.02*(ran[2]-ran[1]))

data = matrix(data,RDD_groups,length(TimeNodes))
datacov = matrix(datacov,RDD_groups,length(TimeNodes))

###########################SEPARABLE###########################################
solSep = smooth.FEM.time(observations = data,incidence_matrix = incidence_matrix, time_locations=TimeNodes,
                         time_mesh = TimeNodes,FEMbasis = FEMbasis, 
                         lambdaS = lambdaS, lambdaT = lambdaT)

solSep = smooth.FEM.time(observations = datacov,time_mesh = TimeNodes, covariates = W_areal,incidence_matrix = incidence_matrix,
                         FEMbasis = FEMbasis, lambdaS = lambdaS, lambdaT = lambdaT)

##########################################PARABOLIC####################################################
solPar = smooth.FEM.time(observations = data[,2:length(TimeNodes)],time_mesh = TimeNodes, incidence_matrix = incidence_matrix,
                         FEMbasis = FEMbasis, lambdaS = lambdaS_par, lambdaT = lambdaT_par, FLAG_PARABOLIC = TRUE,
                         IC=func_evaluation[1:nrow(mesh$nodes)])

solPar = smooth.FEM.time(observations = datacov[,2:length(TimeNodes)],time_mesh = TimeNodes, incidence_matrix = incidence_matrix,covariates = W_areal[(1+RDD_groups):(length(TimeNodes)*RDD_groups),],
                         FEMbasis = FEMbasis, lambdaS = lambdaS_par, lambdaT = lambdaT_par, FLAG_PARABOLIC = TRUE, 
                         IC=func_evaluation[1:nrow(mesh$nodes)])


######### 3D (These tests are slow!) #########

#### sphere 3D pointwise (with or without covariates + locations at nodes or not + stochastic GCV) ####
library(fdaPDE)
rm(list=ls())

# Build mesh: Sphere
data(sphere3Ddata)
sphere3D<-create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)
# plot(sphere3D)
FEMbasis <- create.FEM.basis(sphere3D)
nodesLocations=sphere3D$nodes
nnodes = nrow(sphere3D$nodes)
TimeLocations = seq(0,1,length.out = 5)
Locations = cbind(rep(TimeLocations,each=nnodes),rep(nodesLocations[,1],length(TimeLocations)),rep(nodesLocations[,2],length(TimeLocations)),rep(nodesLocations[,3],length(TimeLocations)))

# Exact test function
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)
a4 = rnorm(1,mean = 1, sd = 1)

func = function(x)
{
  a1*sin(2*pi*(x[,1]*x[,2]))+a2*cos(2*pi*x[,2])+a3*cos(2*pi*x[,3])+a4*sin(2*pi*x[,4])
}

func_evaluation = func(Locations)
ran=range(func_evaluation)
#
# plot(FEM(func_evaluation[1:nnodes],FEMbasis))
# Set smoothing parameter

# Generate locations
nloc = 1000
loc=matrix(data=runif(3*nloc, min=-1,max=1),nrow=nloc,ncol=3,byrow=T)

ind=NULL
for(row in 1:nloc){
  normvec = (loc[row,1]^2+loc[row,2]^2+loc[row,3]^2)
  if(normvec>0.975)   # check points outside the sphere and remove them
    ind = c(ind,row)
}

loc=loc[-ind,]
nloc=dim(loc)[1]
timeloc = seq(0,1,length.out=5)
loc = cbind(rep(timeloc,each=nloc),rep(loc[,1],length(timeloc)),rep(loc[,2],length(timeloc)),rep(loc[,3],length(timeloc)))


# Exact test function - locations different from nodes
func_evaluation2=func(loc)


cov1=(4*sin(2*pi*Locations[,2])+6*sin((2*pi*Locations[,3])^2))*(1-exp(-Locations[,1]))/3
cov2=cos(-2*pi*Locations[,4])+2*Locations[,1]*sin(2*pi*Locations[,2])/6

cov1_nonod=(4*sin(2*pi*loc[,2])+6*sin((2*pi*loc[,3])^2))*(1-exp(-loc[,1]))/3
cov2_nonod=cos(-2*pi*loc[,4])+2*loc[,1]*sin(2*pi*loc[,2])/6

W=cbind(cov1,cov2)
W2=cbind(cov1_nonod,cov2_nonod)



lambdaS=10^seq(-5.0, -4.0, 0.25)
lambdaT=10^seq(-1.5, -0.5, 0.25)

lambdaS2=10^seq(-5.5, -4.5, 0.25)
lambdaT2=10^seq(-1.5, -0.5, 0.25)

lambdaS_par=10^seq(-4.8, -4.4, 0.1)
lambdaT_par=10^seq(1.4, 1.8, 0.1)

lambdaS_par2=10^seq(-4.4, -4.0, 0.1)
lambdaT_par2=10^seq(1.4, 1.8, 0.1)

beta_exact= c(0.7,2.0)
ran = range(func_evaluation)
data = func_evaluation +rnorm(length(func_evaluation),mean=0,sd=0.05*(ran[2]-ran[1]))

ran = range(func_evaluation2)
data_noloc = func_evaluation2 +rnorm(length(func_evaluation2),mean=0,sd=0.05*(ran[2]-ran[1]))

ran = range(func_evaluation+ W%*%beta_exact)
datacov=func_evaluation+ W%*%beta_exact +rnorm(length(func_evaluation),mean=0,sd=0.05*(ran[2]-ran[1]))

data = matrix(data,nnodes,length(TimeLocations))
data_noloc = matrix(data_noloc,nloc,length(timeloc))
datacov = matrix(datacov,nnodes,length(TimeLocations))
###########################SEPARABLE###########################################

solSep = smooth.FEM.time(observations = data,time_mesh = TimeLocations,
                         FEMbasis = FEMbasis, lambdaS = lambdaS, lambdaT = lambdaT)

solSepNoNodes = smooth.FEM.time(locations=loc[1:nloc,2:4],time_locations = timeloc,
                                observations = data_noloc,
                                time_mesh = timeloc,FEMbasis = FEMbasis, 
                                lambdaS = lambdaS2, lambdaT = lambdaT2)

solSepCov = smooth.FEM.time(observations = datacov,time_mesh = TimeLocations, covariates = W,
                            FEMbasis = FEMbasis, lambdaS = lambdaS, lambdaT = lambdaT)

##########################################PARABOLIC####################################################
solPar = smooth.FEM.time(observations = data,time_mesh = TimeLocations,
                         FEMbasis = FEMbasis, lambdaS = lambdaS_par, lambdaT = lambdaT_par, FLAG_PARABOLIC = TRUE)

solParNoNodes = smooth.FEM.time(locations=loc[1:nloc,2:4],observations = data_noloc,time_mesh = timeloc,
                                FEMbasis = FEMbasis, lambdaS = lambdaS_par2, lambdaT = lambdaT_par2, FLAG_PARABOLIC = TRUE)

solParCov = smooth.FEM.time(observations = datacov[,2:length(TimeLocations)],time_mesh = TimeLocations, covariates = W[(1+nnodes):(length(TimeLocations)*nnodes),],
                            FEMbasis = FEMbasis, lambdaS = lambdaS_par, lambdaT = lambdaT_par, FLAG_PARABOLIC = TRUE,
                            IC=func_evaluation[1:nnodes])

#### Example of RMSE computation

sol_eval=eval.FEM.time(solParNoNodes$fit.FEM.time,locations = loc[1:nloc,2:4], time.instants = timeloc, lambdaS = solParNoNodes$bestlambda[1],lambdaT = solParNoNodes$bestlambda[2])
sol_exact = func_evaluation2
RMSE<-function(f,g) sqrt(mean((f-g)^2))
RMSE(sol_eval,sol_exact)
