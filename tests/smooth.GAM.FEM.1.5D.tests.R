###############################################
############ TEST SCRIPT GAM 1.5D #############
###############################################

####### 1.5D ########
#### Test 1: C-shaped domain (1.5D) POISSON family ####
#            locations = nodes
#            with covariates
#            no BC
library(fdaPDE)
library(purrr)
rm(list=ls())
graphics.off()

# family
FAMILY = "poisson"

l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

# beta
beta1 = 0.5
beta2 = -0.2
betas_truth = c(beta1,beta2)

# lambda 
lambda = 10^seq(-2,0,length.out = 10)

# scale param
scale.param = 1

# mesh
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
plot(mesh, asp=1)

# locations 
loc = mesh$nodes
nloc = dim(loc)[1]

# 1.5D spatial field (function f)
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
f = function(x,y){
  
  res = aux.1(x,y) + aux.3(x,y) + aux.4(x,y)
  return(res)
}

sol_exact = f(loc[,1] , loc[,2])
range(sol_exact) 

# covariates
set.seed(42)

desmat=matrix(0,nrow=nloc, ncol=2)
desmat[,1]= rnorm(nloc, mean=2, sd=0.5)
desmat[,2]= rnorm(nloc, mean=0,sd=1)

# samopling response
ran=range(desmat%*%betas_truth + sol_exact) 
param=sol_exact+beta1*desmat[,1]+beta2*desmat[,2]

mu<-inv.link(param)
range(mu)
response <- rpois(nloc, lambda = mu)

# Plot the True Field
plot(FEM(sol_exact, FEMbasis))

#### Test 1.1: Without GCV
output_CPP<- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                        max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                        lambda = lambda)

#### Test 1.2: grid with exact GCV
output_CPP<- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                        max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                        lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$optimization$GCV_vector)

beta = output_CPP$solution$beta
func_estimation = output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position]
plot(FEM(func_estimation,FEMbasis))

#### Test 1.3: grid with stochastic GCV
output_CPP <- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                         max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                         lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$optimization$GCV_vector)

beta = output_CPP$solution$beta
func_estimation = output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position]
plot(FEM(func_estimation,FEMbasis))
