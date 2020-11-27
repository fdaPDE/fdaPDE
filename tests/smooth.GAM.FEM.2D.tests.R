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
plot(log10(lambda),output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position],FEMbasis))
image(FEM(sol_nodes, FEMbasis))

#### Test 1.3: With stochastic GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position],FEMbasis))
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
plot(log10(lambda),output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position],FEMbasis))
image(FEM(sol_nodes, FEMbasis))

#### Test 2.3: With stochastic GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position],FEMbasis))
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
plot(log10(lambda),output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position],FEMbasis))
image(FEM(sol_nodes, FEMbasis))

#### Test 3.3: With stochastic GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position],FEMbasis))
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
plot(log10(lambda),output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position],FEMbasis))
image(FEM(sol_nodes, FEMbasis))

#### Test 4.3: With stochastic GCV
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV' )
plot(log10(lambda),output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position],FEMbasis))
image(FEM(sol_nodes, FEMbasis))

#### Test 5: horseshoe domain POISSON family ####
#            laplacian
#            no covariates
#            no BC
#            areal data
#
library(fdaPDE)
library(purrr)
library(pracma)
library(ggplot2)
rm(list=ls())
graphics.off()

# family
FAMILY = "poisson"

l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

# beta
beta1 = 5

# lambda
lambda = 10^seq(-8,3,length.out = 20)

# scale param
scale.param = 1

# mesh data
load("../data/horseshoe2D_areal.RData")
nodes = mesh$nodes

plot(mesh, lwd=3, cex = 1.9)

nnodes = dim(mesh$nodes)[1]
FEMbasis = fdaPDE::create.FEM.basis(mesh)

# 2D random field (function f)
a_fun <- function(p){
  
  if(p[1]>= 0 && p[2] > 0){
    pi/4 + p[1]
  }else{
    if(p[1]>= 0 && p[2] <= 0){
      -pi/4 -p[1]
    }else{
      if(p[1] < 0){
        -0.5*atan(p[2]/p[1])
      }
    }
  }
}

d_fun <- function(p){
  
  if(p[1]>= 0 && p[2] > 0){
    -0.5 + p[2]
  }else{
    if(p[1]>= 0 && p[2] <= 0){
      -0.5 - p[2]
    }else{
      if(p[1] < 0){
        sqrt(p[1]^2 + p[2]^2) - 0.5
      }
    }
  }
}

z <- function(p){
  a_fun(p) + d_fun(p)^2
}

tri = mesh$triangles
dim(incidence_matrix)
dim(tri)
integration_nodes = data.frame(matrix(nrow = dim(incidence_matrix)[2],ncol = 11))
names(integration_nodes) = c("T1","T2","N1","N2","N3","N4", "xmin", "xmax", "ymin", "ymax", "label")
nodi_plot = NULL

for(i in 1:dim(incidence_matrix)[2])
{
  tri_used = which(incidence_matrix[,i] == 1)
  integration_nodes$T1[i] = tri_used[1]
  integration_nodes$T2[i] = tri_used[2]
  nodes_used = unique(c(tri[tri_used[1],],tri[tri_used[2],]))
  integration_nodes$N1[i] = nodes_used[1]
  integration_nodes$N2[i] = nodes_used[2]
  integration_nodes$N3[i] = nodes_used[3]
  integration_nodes$N4[i] = nodes_used[4]
  integration_nodes$label[i] = i
  xvec = c(nodes[integration_nodes$N1[i],1],nodes[integration_nodes$N2[i],1],nodes[integration_nodes$N3[i],1],nodes[integration_nodes$N4[i],1])
  yvec = c(nodes[integration_nodes$N1[i],2],nodes[integration_nodes$N2[i],2],nodes[integration_nodes$N3[i],2],nodes[integration_nodes$N4[i],2])
  integration_nodes$xmin[i] = min(xvec)
  integration_nodes$xmax[i] = max(xvec)
  integration_nodes$ymin[i] = min(yvec)
  integration_nodes$ymax[i] = max(yvec)
  
  nodi_tmp = rbind(mesh$nodes[integration_nodes$N1[i],],mesh$nodes[integration_nodes$N2[i],],mesh$nodes[integration_nodes$N3[i],],mesh$nodes[integration_nodes$N4[i],])
  nodi_tmp = cbind(nodi_tmp,rep(i,4))
  nodi_plot = rbind(nodi_plot,nodi_tmp)
}

nodi_plot <- as.data.frame(nodi_plot)
names(nodi_plot) = c("X","Y","Tri")
nodi_plot$Tri = as.factor(nodi_plot$Tri)

# plot areal mesh
nnodes = dim(mesh$nodes)[1]
plot_values = numeric(nnodes)

for(i in 1:dim(incidence_matrix)[2])
{
  plot_values[integration_nodes$N1[i]] = 15
  plot_values[integration_nodes$N2[i]] = 15
  plot_values[integration_nodes$N3[i]] = 15
  plot_values[integration_nodes$N4[i]] = 15
}

image(FEM(plot_values,FEMbasis))

nnodes = dim(mesh$nodes)[1]
sol_exact <- numeric(nnodes)
for(i in 1:nnodes){
  sol_exact[i] <- z(mesh$nodes[i,])
}

z_xy <- function(x,y){z(c(x,y))}
sol_integrand = matrix(nrow = 20 , ncol = 1)
for(i in 1:20){
  res = integral2(fun = z_xy, xmin = integration_nodes$xmin[i], xmax = integration_nodes$xmax[i] , ymin = integration_nodes$ymin[i], ymax = integration_nodes$ymax[i])
  sol_integrand[i] = res$Q
}

range(sol_integrand)

set.seed(2020)

# covariates
desmat=matrix(0,nrow=length(sol_integrand),ncol=1)
desmat[,1]=rbeta(length(sol_integrand),shape1=2,shape2=2)  # sampling covariates from beta distr.

theta = sol_integrand + desmat%*%beta1
range(theta)

mu = inv.link(theta)
range(mu)

# sampling response
response <- rpois(length(mu), lambda = mu)

#### Test 5.1: Without GCV
output_CPP <- smooth.FEM(observations = as.numeric(response), FEMbasis = FEMbasis, covariates = desmat,
                         max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                         incidence_matrix = t(incidence_matrix), areal.data.avg = FALSE,
                         lambda = lambda)

#### Test 5.2: grid with exact GCV
output_CPP <- smooth.FEM(observations = as.numeric(response), FEMbasis = FEMbasis, covariates = desmat,
                         max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                         incidence_matrix = t(incidence_matrix), areal.data.avg = FALSE,
                         lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$optimization$GCV_vector)

image(FEM(output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position], FEMbasis))

#### Test 5.3: grid with stochastic GCV
output_CPP <- smooth.FEM(observations = as.numeric(response), FEMbasis = FEMbasis, covariates = desmat,
                         max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                         incidence_matrix = t(incidence_matrix), areal.data.avg = FALSE,
                         lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position], FEMbasis))