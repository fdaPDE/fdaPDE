########################################
############## TEST SCRIPT #############
#######  space-time GLR problems #######
########################################

####### 2D ####### 

### Test 1: C-shaped domain 
#           locations != nodes
#           with covariates
#           no BC
#           order FE = 1

library(fdaPDE)
FAMILY = "gamma"
inv.link <-  function(mu){ return(-1 / mu)} 

data("horseshoe2D")
mesh = create.mesh.2D(nodes = horseshoe2D$boundary_nodes,
                      segments = horseshoe2D$boundary_segments)
mesh = refine.mesh.2D(mesh, maximum_area =  0.025, minimum_angle = 30 )
FEMbasis = create.FEM.basis(mesh)
plot(mesh)

f <- function(x, y, t) {
  
  return((-1.5 * sin(2 * pi * x) * cos(2 * pi * y) +
              2 / 5 * sin(3 * pi * x * t) - 2) * (t + 1))
}

is.in.horseshoe <- function(x, y) {
  r <- .5
  r0 <- .1
  
  if ((x - 3)^2 + (y - r)^2 < (r - r0)^2 & x > 3) {
    return(TRUE)
  }
  if (r0^2 <= x^2 + y^2 & x^2 + y^2 <= (2 * r - r0)^2 & x <= 0) {
    return(TRUE)
  }
  if ((x - 3)^2 + (y + r)^2 < (r - r0)^2 & x > 3) {
    return(TRUE)
  }
  if (abs(y) > r0 & abs(y) < 2 * r - r0 & 0 < x & x <= 3) {
    return(TRUE)
  }
  return(FALSE)
}

is.p.in.horseshoe <- function(p) {
  return(is.in.horseshoe(p[1], p[2]))
}

### generating data ###
M = 4
m = M
n = 400

time_mesh <- seq(0, 1, length.out = M)
time_locations = time_mesh

set.seed(32)

loc = cbind(runif(2 * n, min = -1, max = 4),
            runif(2 * n, min = -1, max = 1))
ww <- apply(loc, 1, is.p.in.horseshoe) # ! is.na(fs.test(loc[,1],
                                                # loc[ ,2], exclude = T))
loc <- loc[ww, ]  

space_time_locations <-cbind(rep(time_locations, each = nrow(loc)),
                             rep(loc[, 1], length(time_locations)),
                             rep(loc[, 2], length(time_locations)) )
### covariates ###
betas = c(-.2, .3)
desmat <- matrix(0, nrow = nrow(space_time_locations), ncol = 2)
desmat[, 1] <- rbeta(n = nrow(desmat), shape1 = 1.5, shape2 = 2)
desmat[, 2] <- rbeta(n = nrow(desmat), shape1 = 3, shape2 = 2)

true.field <- f(space_time_locations[, 2], space_time_locations[, 3],space_time_locations[, 1])  
field = true.field + desmat %*% betas 

scale = 1
data <- rgamma(n = nrow(space_time_locations),
               shape = -1 / field / scale, scale = scale)
data <- matrix(data, nrow = nrow(loc), ncol = m)

### smoothing parameters ###
lambdaS <- 10^seq(-3, -1, by=0.5)
lambdaT <- 10^seq(-3, -2, by=0.5)

########### PARABOLIC ###########

#### Test 1.1: without GCV
output_CPP <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                           locations = loc, observations = data,  FEMbasis = FEMbasis, 
                           covariates = desmat,  
                           lambdaS = lambdaS[2], lambdaT = lambdaT[2],
                           family = FAMILY, FLAG_PARABOLIC = T)

#### Test 1.2: with stochastic GCV
output_CPP <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                           locations = loc, observations = data,  FEMbasis = FEMbasis, 
                           covariates = desmat,  
                           DOF.evaluation = "stochastic", 
                           lambda.selection.lossfunction = "GCV",
                           lambdaS = lambdaS, lambdaT = lambdaT,
                           family = FAMILY, FLAG_PARABOLIC = T)

########### SEPARABLE ############

#### Test 1.3: without GCV
output_CPP <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                          locations = loc, observations = data, FEMbasis = FEMbasis, 
                          covariates = desmat, 
                          lambdaS = lambdaS[2], lambdaT = lambdaT[2],
                          family = FAMILY, FLAG_PARABOLIC = F)

#### Test 1.4: with stochastic GCV
output_CPP <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                          locations = loc, observations = data, FEMbasis = FEMbasis, 
                          covariates = desmat, 
                          DOF.evaluation = "stochastic", lambda.selection.lossfunction = "GCV",
                          lambdaS = lambdaS, lambdaT = lambdaT,
                          family = FAMILY, FLAG_PARABOLIC = F)

####### 2.5D ####### 

#### Test 2: Hub domain 
#            locations != nodes
#            with covariates
#            no BC
#            order FE = 1

library(fdaPDE)
FAMILY = "poisson"

l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

data("hub2.5D")
mesh <- create.mesh.2.5D(nodes = hub2.5D$hub2.5D.nodes,
                         triangles = hub2.5D$hub2.5D.triangles)
FEMbasis <- create.FEM.basis(mesh)

### generating data ###
nodesLocations=mesh$nodes
set.seed(32)
x = runif(400, min = -0.6, max = 0.6)
y = runif(400, min = -0.6, max = 0.6)
z = runif(400, min =  0.0, max = 1.0)
points_ = cbind(x,y,z)
loc = projection.points.2.5D(mesh, points_)

nloc = nrow(loc) 
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)
time_mesh = as.numeric(0:4)

time_locations = time_mesh
space_time_locations = cbind(rep(time_locations,each=nloc),
                  rep(loc[,1],length(time_locations)),
                  rep(loc[,2],length(time_locations)),
                  rep(loc[,3],length(time_locations)))

func = function(x)
{
  (a1*sin(2*pi*x[,2])+a2*sin(2*pi*x[,3])+a3*sin(2*pi*x[,4])+1)*cos(x[,1])
}

func_evaluation = func(space_time_locations)

lambdaS=10^seq(from=-3, to=0, by=1)
lambdaT=10^seq(from=-2, to=1, by=1)

cov1=4*sin(2*pi*space_time_locations[,2])*cos(2*pi*space_time_locations[,3])
cov2=rnorm(nloc*length(time_locations), mean=3, sd=0.1)*rep(exp(-time_locations/length(time_locations)),each=nloc)
W=cbind(cov1,cov2)

# Fix betas
beta_exact=c(0.45,0.3)

ran=range(W%*%beta_exact + func_evaluation)
ran=range(func_evaluation)

ran = range(func_evaluation+ W%*%beta_exact)
datacov=func_evaluation+ W%*%beta_exact +rnorm(nloc*length(time_locations),mean=0,sd=0.05*(ran[2]-ran[1]))

mu<-inv.link(datacov)
range(mu)
response <- rpois(length(mu), lambda = mu)

response <- matrix(response, nloc, length(time_locations))
storage.mode(response) <- "numeric"

########### PARABOLIC ###########

#### Test 2.1: without GCV 
output_CPP <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                              locations = loc, observations = response, FEMbasis = FEMbasis, 
                              covariates = W, 
                              lambdaS = lambdaS[1], lambdaT =  lambdaT[1],
                              family = FAMILY, FLAG_PARABOLIC = T)

#### Test 2.1: with exact GCV 
output_CPP <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                              locations = loc, observations = response, FEMbasis = FEMbasis, 
                              covariates = W,
                              lambda.selection.lossfunction = "GCV",
                              DOF.evaluation = "exact",
                              lambdaS = lambdaS, lambdaT =  lambdaT,
                              family = FAMILY, FLAG_PARABOLIC = T)

########### SEPARABLE ###########

### Test 2.3: without GCV 
output_CPP <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                              locations = loc, observations = response, FEMbasis = FEMbasis, 
                              covariates = W, 
                              lambdaS = lambdaS[1], lambdaT =  lambdaT[1],
                              family = FAMILY, FLAG_PARABOLIC = F)

### Test 2.4; with exact GCV
output_CPP <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                              locations = loc, observations = response, FEMbasis = FEMbasis, 
                              covariates = W, 
                              lambda.selection.lossfunction = "GCV",
                              DOF.evaluation = "exact",
                              lambdaS = lambdaS, lambdaT =  lambdaT,
                              family = FAMILY, FLAG_PARABOLIC = F)

####### 3D ####### 

#### Test 3: Sphere domain 
#            locations at nodes
#            no covariates
#            no BC
#            order FE = 1

library(fdaPDE)
FAMILY = "poisson"
l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

data(sphere3Ddata)
sphere3D<-create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)

FEMbasis <- create.FEM.basis(sphere3D)

nodesLocations=sphere3D$nodes
nnodes = nrow(sphere3D$nodes)
time_locations = seq(0,1,length.out = 3)
space_time_locations = cbind(rep(time_locations,each=nnodes),rep(nodesLocations[,1],length(time_locations)),rep(nodesLocations[,2],length(time_locations)),rep(nodesLocations[,3],length(time_locations)))

# Exact test function
a1 = -1
a2 = -2
a3 = -3

func = function(x){
  (a1* sin(x[,2]) +  a2* sin(x[,3]) +  a3*sin(x[,4])) * (sin(x[,1]))
}

func_evaluation = func(space_time_locations)
ran=range(func_evaluation)

# generating data
data = func_evaluation +rnorm(length(func_evaluation),mean=0,sd=0.05*(ran[2]-ran[1]))
mu = inv.link(data)
response  <- rpois(length(mu), lambda = mu)
range(response)

response = matrix(response, nnodes, length(time_locations))
storage.mode(response) <- "numeric"

# Set smoothing parameter
lambdaS=10^seq(0, 1, 0.5)
lambdaT=10^seq(0, 1, 0.5)

########### SEPARABLE ###########

#### Test 3.1: without GCV 
output_CPP <- smooth.FEM.time(locations = NULL, observations = response, FEMbasis =FEMbasis, 
                              covariates = NULL, time_mesh = time_locations, time_locations = time_locations,
                              max.steps=15, family=FAMILY,
                              FLAG_PARABOLIC = F,
                              lambdaS = lambdaS[1], lambdaT = lambdaT[1])

#### Test 3.2: with stochastic GCV
output_CPP <- smooth.FEM.time(locations = NULL, observations = response, FEMbasis =FEMbasis, 
                              covariates = NULL, time_mesh = time_locations, time_locations = time_locations,
                              max.steps=15, family=FAMILY,
                              lambda.selection.lossfunction = "GCV",
                              DOF.evaluation = "stochastic",
                              FLAG_PARABOLIC = F,
                              lambdaS = lambdaS, lambdaT = lambdaT)


#### Test 3.3: with exact GCV
output_CPP <- smooth.FEM.time(locations = NULL, observations = response, FEMbasis =FEMbasis, 
                              covariates = NULL, time_mesh = time_locations, time_locations = time_locations,
                              max.steps=15, family=FAMILY,
                              lambda.selection.lossfunction = "GCV",
                              DOF.evaluation = "exact",
                              FLAG_PARABOLIC = F,
                              lambdaS = lambdaS, lambdaT = lambdaT)


####### 1.5D ####### 

#### Test 4: C-shaped domain 
#            locations at nodes
#            no covariates
#            no BC
#            order FE = 1

library(fdaPDE)
FAMILY = "poisson"
l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

eps = 1 / 2
x = c(0., 1)
y = c(0.,eps)
vertices = expand.grid(x,y)
vertices = cbind(vertices[,1], vertices[,2])
edges = matrix(c(1,2,1,3,3,4), nrow=3,ncol=2, byrow=T)

mesh = create.mesh.1.5D(vertices, edges)
mesh = refine.mesh.1.5D(mesh,delta=0.0125)

nnodes=dim(mesh$nodes)[1]
FEMbasis=create.FEM.basis(mesh)

# Locations at nodes
nodesLocations=mesh$nodes

nnodes = nrow(mesh$nodes)
time_locations = seq(0,pi,length.out = 6)
space_time_locations = cbind(rep(time_locations,each=nnodes),rep(nodesLocations[,1],length(time_locations)),rep(nodesLocations[,2],length(time_locations)) )

aux.4 = function(x,y){
  h = 1
  source = 4 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,1] - mesh$nodes[source,1])
    if(delta < h ){
      coef[i] = 3*exp(1/(delta^2 / h^2 -1)+1 ) - 1
      
    }
  }
  
  return(coef)
}
aux.3 = function(x,y){
  
  h = eps
  source = 3 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,1]-mesh$nodes[source,1]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,2] - mesh$nodes[source,2])
    if(delta < h ){
      coef[i] = -1 - 1/h*delta
    }
    
    
  }
  return(coef)
}
aux.1 = function(x,y){
  
  h = 1
  source = 1 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,1] - mesh$nodes[source,1])
    if(delta <= h ){
      coef[i] = -2 - 1/h*delta
    }
  }
  return(coef)
  
}  
func = function(x,y,t){
  
  res = (aux.1(x,y) + aux.3(x,y) + aux.4(x,y)) * cos(t)
  return(res)
}

func_evaluation = func(space_time_locations[,2], 
                       space_time_locations[,3], 
                       space_time_locations[,1])

ran = range(func_evaluation)
data = func_evaluation +rnorm(nnodes*length(time_locations),mean=0,sd=0.05*(ran[2]-ran[1]))
mu = inv.link(data)
response = rpois(nnodes*length(time_locations), lambda = mu)

data = matrix(data,nrow(mesh$nodes),length(time_locations))

storage.mode(response)<-"numeric"
response <- matrix(response, nrow(mesh$nodes), length(time_locations))

# Set smoothing parameters
lambdaS = 10^seq(1, 2, by=0.25)
lambdaT = 10^seq(0, 1, by=0.25)

########### PARABOLIC ###########

#### Test 4.1: with stochastic GCV
output_CPP <- smooth.FEM.time(observations=response,
                         FEMbasis = FEMbasis, time_mesh = time_locations, 
                         time_locations = time_locations,
                         lambdaS = lambdaS, lambdaT = lambdaT, 
                         family = FAMILY, 
                         lambda.selection.criterion = "grid",
                         DOF.evaluation = "stochastic", 
                         lambda.selection.lossfunction = "GCV",
                         FLAG_PARABOLIC = TRUE)

#### Test 4.2: with exact GCV
output_CPP <- smooth.FEM.time(observations=response,
                           FEMbasis = FEMbasis, time_mesh = time_locations, 
                           time_locations = time_locations,
                           lambdaS = lambdaS, lambdaT = lambdaT, 
                           family = FAMILY, 
                           lambda.selection.criterion = "grid",
                           DOF.evaluation = "exact", 
                           lambda.selection.lossfunction = "GCV",
                           FLAG_PARABOLIC = TRUE)

########### SEPARABLE ###########

#### Test 4.3: with stochastic GCV
output_CPP <- smooth.FEM.time(observations=response,
                           FEMbasis = FEMbasis, time_mesh = time_locations, 
                           time_locations = time_locations,
                           lambdaS = lambdaS, lambdaT = lambdaT, 
                           family = FAMILY, 
                           lambda.selection.criterion = "grid",
                           DOF.evaluation = "stochastic", 
                           lambda.selection.lossfunction = "GCV",
                           FLAG_PARABOLIC = FALSE)

#### Test 4.4: with exact GCV
output_CPP <- smooth.FEM.time(observations=response,
                           FEMbasis = FEMbasis, time_mesh = time_locations, 
                           time_locations = time_locations,
                           lambdaS = lambdaS, lambdaT = lambdaT, 
                           family = FAMILY, 
                           lambda.selection.criterion = "grid",
                           DOF.evaluation = "exact", 
                           lambda.selection.lossfunction = "GCV",
                           FLAG_PARABOLIC = FALSE)
