### Confronto Mattina - fdaPDE ###

vertices_x = as.matrix(simplenet$vertices$x)
vertices_y = as.matrix(simplenet$vertices$y)
vertices = as.matrix(cbind(vertices_x, vertices_y))
start_coord = cbind(simplenet$lines$ends$x0, simplenet$lines$ends$y0)
end_coord = cbind(simplenet$lines$ends$x1, simplenet$lines$ends$y1)
start = simplenet$from 
end = simplenet$to
edges = cbind(start, end)
M = simplenet$m
find.boundary(M)
x11()
plot(simplenet)
L = as.linnet(simplenet)
# delta > 0.0001 altrimenti PHI...
delta = 0.01

mesh = create.mesh.1D.vertices(vertices, edges, delta)
mesh$nnodes
nodes = mesh$nodes
nodes_x = mesh$nodes[,1]
nodes_y = mesh$nodes[,2]
points = nodes 
points_x = nodes_x
points_y = nodes_y

nobs = mesh$nnodes

FEMbasis = create.FEM.basis.1D(mesh) 


### fdaPDE ### 
mesh.fdaPDE = fdaPDE::create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)

## controllo se ho stessi nodi
dim(mesh.fdaPDE$nodes)[1]
sum( abs(mesh.fdaPDE$nodes - mesh$nodes) > .Machine$double.eps )
###


# CAMPO f 

aux = function(x, y, seg, tp) { cos((2*x+y)/4) +  ((x+y))^2 }

my_dens <- linfun(aux, L)

x11()
plot(my_dens)

W1 = rnorm(nobs, 1, sd=1.5)
W2 = rnorm(nobs, 2, sd=2)

W1 = rnorm(nobs,0,1)
W2 = rnorm(nobs,1,1.5)*cos(nodes_x)*sin(nodes_y)
W = NULL
W = cbind(W1,W2)

beta1 = -0.5
beta2 = 0.2

temp = range(beta1*W1 + beta2*W2 + my_dens(points_x, points_y))
my_range = temp[2]-temp[1]
EPS = rnorm(nobs, 0, sd=sqrt(0.05*my_range) )

Z = beta1*W1 + beta2*W2 + my_dens(points_x, points_y) + EPS

PHI=diag(nobs)

epsilon = 0.0001

perm = sample(1:nobs, nobs)
points_perm = points[perm,]
PHI_perm = PHI[perm,]
Z_perm = Z[perm]
W_perm = W[perm,]

lambda_vector = c(10^-10,10^-9,10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,10,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10)
output.fdaPDE = smooth.FEM( observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W, 
                            lambda=lambda_vector, lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV' )

lambdas = as.matrix( output.fdaPDE$optimization$lambda_solution * c(1:5) ) 
beta_ex = rbind(beta1,beta2)

# || beta.fdaPDE - beta || 
beta.err = matrix(0,nrow=length(lambdas),ncol=1)
# || beta.* - beta_ex ||
err_ = matrix(0, nrow=length(lambdas), ncol=2) 
times = matrix(0,nrow=length(lambdas), ncol=2)

options(digits.sec=6)

for(i in 1:length(lambdas)){
  start.time <- Sys.time()
  output = covariate.FEM.1D(  W, Z, FEMbasis, PHI, lambdas[i])
  end.time  <-Sys.time()
  
  
  output.fdaPDE = smooth.FEM( observations=Z, 
                              FEMbasis=FEMbasis.fdaPDE, 
                              covariates=W, lambda=lambdas[i] )
  
  
  beta.err[i]   = norm(output$beta_estimate - output.fdaPDE$solution$beta, type = "2")
  err_[i,1]     = norm(beta_ex - output$beta_estimate, type="2")
  err_[i,2]     = norm(beta_ex - output.fdaPDE$solution$beta,type="2")
  
  times[i,1] = end.time - start.time
  times[i,2] = output.fdaPDE$time 
  
  }

precision   = sum( abs(err_[,2] -  err_[,1]) < 10 * .Machine$double.eps) /length(lambdas)  
performance = sum( times[,2] < times[,1]) /length(lambdas)

precision
performance

