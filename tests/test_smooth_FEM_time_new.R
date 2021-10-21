##########################################
############## TEST SCRIPT ###############
######## for space-time problems #########
########## using Evaluator_new ###########

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
fdaPDE::plot.mesh.2D(mesh)

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

xeval=runif(1000,0,1)
yeval=runif(1000,0,1)
teval=runif(1000,0,2)

sol_eval=eval.FEM.time(output_CPP$fit.FEM.time,space.time.locations = cbind(teval,xeval,yeval))
sol_eval.new= eval.FEM.time.new(output_CPP$fit.FEM.time,space.time.locations = cbind(teval,xeval,yeval)) 

diff.max = max( abs(sol_eval-sol_eval.new) )

sol_exact=f(xeval,yeval,teval)

RMSE<-function(f,g) sqrt(mean((f-g)^2))
RMSE(sol_eval,sol_exact)
RMSE(sol_eval.new,sol_exact)

#### Test 1.2: grid with exact GCV (salto -> no evaluation)
#### Test 1.3: grid with stochastic GCV (salto -> no evaluation)

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

plot(output_CPP$fit.FEM.time)

xeval=runif(1000,-0.5,0.5)
yeval=runif(1000,-0.5,0.5)
teval=runif(1000,0,2)

sol_eval=eval.FEM.time(output_CPP$fit.FEM.time,space.time.locations = cbind(teval,xeval,yeval))
sol_eval.new= eval.FEM.time.new(output_CPP$fit.FEM.time,space.time.locations = cbind(teval,xeval,yeval)) 


sol_eval_ = sol_eval[ is.finite(sol_eval)]
sol_eval.new_ = sol_eval.new[is.finite(sol_eval.new)]
diff.max = max( abs(sol_eval_ - sol_eval.new_))
diff.max

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

xeval=runif(1000,0,1)
yeval=runif(1000,0,1)
teval=runif(1000,0,2)

sol_eval=eval.FEM.time(output_CPP$fit.FEM.time,space.time.locations = cbind(teval,xeval,yeval))
sol_eval.new= eval.FEM.time.new(output_CPP$fit.FEM.time,space.time.locations = cbind(teval,xeval,yeval)) 

diff.max = max( abs(sol_eval-sol_eval.new) )
diff.max

#### Test 3.2: exact GCV (salto)
#### Test 3.2: exact GCV (salto)
#### Test 3.3: stochastic GCV (salto)

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

xeval=runif(1000,0,1)
yeval=runif(1000,0,1)
teval=runif(1000,0,2)

sol_eval=eval.FEM.time(output_CPP$fit.FEM.time,time.instants = teval,incidence_matrix = incidence_matrix)
sol_eval.new= eval.FEM.time.new(output_CPP$fit.FEM.time,time.instants = teval, incidence_matrix = incidence_matrix) 

sol_eval_ = sol_eval[ is.finite(sol_eval)]
sol_eval.new_ = sol_eval.new[is.finite(sol_eval.new)]
diff.max = max( abs(sol_eval_ - sol_eval.new_))
diff.max



                            