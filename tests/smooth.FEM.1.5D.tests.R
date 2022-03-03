library(fdaPDE)
rm(list=ls())
graphics.off()

# Create mesh
delta = 1 / 2
x = c(0., 1)
y = c(0., delta)
vertices = expand.grid(x,y)
vertices = cbind(vertices[,1], vertices[,2])
edges = matrix(c(1,2,1,3,3,4), nrow=3,ncol=2, byrow=T)

mesh = create.mesh.1.5D(vertices, edges)
mesh = refine.mesh.1.5D(mesh,delta=0.0125)

nnodes=dim(mesh$nodes)[1]
FEMbasis=create.FEM.basis(mesh)
plot(mesh, asp=1, pch=16)

# True spatial field
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
  
  h = delta
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
f= function(x,y){
  
  res = aux.1(x,y) + aux.3(x,y) + aux.4(x,y)
  return(res)
}

sol_exact=f(mesh$nodes[,1],mesh$nodes[,2])
plot(FEM(sol_exact, FEMbasis))

# Generate covariate and data
locations = refine.by.splitting.mesh.1.5D(mesh)$nodes
nloc = nrow(locations)

cov1 = rnorm(nloc, mean=2,sd=1)
DatiEsatti = f(locations[,1],locations[,2]) + cov1

# Add error to simulate data
set.seed(32)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(nloc, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^seq(-4,-1,by=0.25)

#### Test 1: Without GCV
output_CPP<-smooth.FEM(observations=data, 
                       locations = locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, 
                       lambda=lambda)
plot(output_CPP$fit.FEM)
output_CPP$solution$beta

#### Test 2: grid with exact GCV
output_CPP<-smooth.FEM(observations=data, locations = locations,
                       covariates = cov1,
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector, ylab="CV", main="Exact GCV")
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))
output_CPP$solution$beta

#### Test 3: grid with stochastic GCV
output_CPP<-smooth.FEM(observations=data, locations = locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector, ylab="CV", main="Stochastic GCV")
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))
output_CPP$solution$beta
