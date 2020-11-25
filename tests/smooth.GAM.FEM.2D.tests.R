###############################################
############## TEST SCRIPT GAM  ###############
###############################################
rm(list=ls())
graphics.off()

library(fdaPDE)

####### 2D ########

#### Test 1: square domain BINOMIAL family ####
#            locations = nodes 
#            laplacian
#            no covariates
#            no BC
#            order FE = 1
#
library(fdaPDE)
FAMILY = "binomial"
library(purrr)

logit <- function(x){qlogis(x)}
inv.logit <- function(x){plogis(x)}
link = logit
inv.link = inv.logit

# mesh
x = seq(0,1, length.out = 40)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)
plot(mesh)

nnodes = dim(mesh$nodes)[1]

FEMbasis = create.FEM.basis(mesh)

set.seed(42)
# locations
nloc = 800
xobs=runif(min=0,max=1,n=nloc)
yobs=runif(min=0,max=1,n=nloc)
loc=cbind(xobs,yobs)

# 2D random field (function f) 
a1=-2.5
a2=0.8

z<-function(p)
{
  a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1])
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

range(sol_exact) 
param = sol_exact
mu<-inv.link(param)
range(mu)
# sampling response:
response <- rbernoulli(length(loc[,1]),p = mu)

# Set smoothing parameter
lambda = 10^seq(-5,0,length.out = 20)

#### Test 1.1: Without GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis = FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda)
image(FEM(sol_nodes, FEMbasis))
image(output_CPP$fit.FEM)


#### Test 1.2: With exact GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))
image(FEM(sol_nodes, FEMbasis))

#### Test 1.3: With stochastic GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))
image(FEM(sol_nodes, FEMbasis))


#### Test 2: square domain EXPONENTIAL family ####
#            locations = nodes 
#            laplacian
#            no covariates
#            no BC
#            order FE = 1
#
library(fdaPDE)
rm(list=ls())
graphics.off()

FAMILY = "exponential"

link<-function(x){-1/x}
inv.link<-link 

# mesh
x = seq(0,1, length.out = 40)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)
plot(mesh)

nnodes = dim(mesh$nodes)[1]

FEMbasis = create.FEM.basis(mesh)

# locations -----------------------------------------------
set.seed(42)
# locations
nloc = 800
xobs=runif(min=0,max=1,n=nloc)
yobs=runif(min=0,max=1,n=nloc)
loc=cbind(xobs,yobs)

# 2D random field (function f) ------------------------------
a1=-1.5
a2=0.4

z<-function(p)
{
  a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1]) - 2
  
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

range(sol_exact) 
param = sol_exact
mu<-inv.link(param)
range(mu)
# sampling response:
set.seed(95)
response <- response <- rexp(nloc, rate = 1/mu)

# Set smoothing parameter
lambda = 10^seq(-5,0,length.out = 20)

#### Test 2.1: Without GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda=lambda)
image(output_CPP$fit.FEM)
image(FEM(sol_nodes, FEMbasis))

#### Test 2.2: With exact GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))
image(FEM(sol_nodes, FEMbasis))

#### Test 2.3: With stochastic GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))
image(FEM(sol_nodes, FEMbasis))


#### Test 3: square domain GAMMA family ####
#            locations = nodes 
#            laplacian
#            no covariates
#            no BC
#            order FE = 1
#
library(fdaPDE)
rm(list=ls())
graphics.off()

FAMILY = "gamma"

link<-function(x){-1/x}
inv.link<-link

# mesh
x = seq(0,1, length.out = 40)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)
plot(mesh)

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
  a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1]) - 2
  
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

range(sol_exact) 
param = sol_exact
mu<-inv.link(param)
range(mu)
# sampling response:
set.seed(95)
scale.param=1

response <- rgamma(length(loc[,1]), shape=mu/scale.param, scale=scale.param)

# Set smoothing parameter
lambda = 10^seq(-5,0,length.out = 20)

#### Test 3.1: Without GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda)
image(output_CPP$fit.FEM)

#### Test 3.2: With exact GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))
image(FEM(sol_nodes, FEMbasis))

#### Test 3.3: With stochastic GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))
image(FEM(sol_nodes, FEMbasis))


#### Test 4: square domain POISSON family ####
#            locations = nodes 
#            laplacian
#            no covariates
#            no BC
#            order FE = 1
#
library(fdaPDE)
rm(list=ls())
graphics.off()

FAMILY = "poisson"

l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

# mesh
x = seq(0,1, length.out = 40)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)
plot(mesh)

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


range(sol_exact) 
param = sol_exact
mu<-inv.link(param)
range(mu)
# sampling response:
set.seed(95)

response <- rpois(length(loc[,1]), lambda = mu)

# Set smoothing parameter
lambda = 10^seq(-5,0,length.out = 20)

#### Test 4.1: Without GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda=lambda)
image(output_CPP$fit.FEM)

#### Test 4.2: With exact GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))
image(FEM(sol_nodes, FEMbasis))

#### Test 4.3: With stochastic GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV' )
plot(log10(lambda),output_CPP$GCV)
image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))
image(FEM(sol_nodes, FEMbasis))

