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
FEMbasis = create.FEM.basis.1D(mesh)

nodes_x = mesh$nodes[,1]
nodes_y = mesh$nodes[,2]

#locations == nodes
nobs = mesh$nnodes
points_x = nodes_x
points_y = nodes_y

### fdaPDE - MESH ### 
mesh.fdaPDE = fdaPDE::create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)
# check
dim(mesh.fdaPDE$nodes)[1]
sum( abs(mesh.fdaPDE$nodes - mesh$nodes) > .Machine$double.eps )
###


### CAMPO f ### 
aux = function(x, y, seg, tp) { cos((2*x+y)/4) +  ((x+y))^2 }
my_dens <- linfun(aux, L)

x11()
plot(my_dens)

### observations ###
W1 = rnorm(nobs, 1, sd=1.5)
W2 = rnorm(nobs, 2, sd=2)
W = NULL
W = cbind(W1,W2)

beta1 = -0.5
beta2 = 0.2

temp = range(beta1*W1 + beta2*W2 + my_dens(points_x, points_y))
my_range = temp[2]-temp[1]
EPS = rnorm(nobs, 0, sd=sqrt(0.05*my_range) )
Z = beta1*W1 + beta2*W2 + my_dens(points_x, points_y) + EPS
PHI=diag(nobs)

### cross - validation ###
lambda_vector = c(10^-10,10^-9,10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,10,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10)
output.fdaPDE = smooth.FEM( observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W, 
                            lambda=lambda_vector, lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV' )

lambdas = as.matrix( output.fdaPDE$optimization$lambda_solution * c(1:5) ) 
beta_ex = rbind(beta1,beta2)

# || beta.fdaPDE - beta || 
beta.err = matrix(0,nrow=length(lambdas),ncol=1)
# || beta.* - beta_ex ||
err_ = matrix(0, nrow=length(lambdas), ncol=2) 

times = matrix(0,nrow=length(lambdas), ncol=3)
options(digits.sec=6)

for(i in 1:length(lambdas)){
  start.time <- Sys.time()
  output = covariate.FEM.1D(  W, Z, FEMbasis, PHI, lambdas[i])
  end.time  <-Sys.time()
  
  start.time.fdaPDE <- Sys.time()
  output.fdaPDE = smooth.FEM( observations=Z, 
                              FEMbasis=FEMbasis.fdaPDE, 
                              covariates=W, lambda=lambdas[i] )
  end.time.fdaPDE  <-Sys.time()
  
  beta.err[i]   = norm(output$beta_estimate - output.fdaPDE$solution$beta, type = "2")
  err_[i,1]     = norm(beta_ex - output$beta_estimate, type="2")
  err_[i,2]     = norm(beta_ex - output.fdaPDE$solution$beta,type="2")
  
  times[i,1] =  end.time - start.time
  times[i,2] =  end.time.fdaPDE - start.time.fdaPDE
  times[i,3] =  output.fdaPDE$time
  }



x11()
plot(1:length(lambdas), times[,1],type="b",pch = 15,col="red", xlab='',ylab='times', 
        ylim=c(min(times),max(times)) )
points(1:length(lambdas), times[,2],type="b",pch = 16, col="blue")
points(1:length(lambdas), times[,3],type="b",pch = 17, col="magenta")

legend("topleft",legend=c("time Mattina", "time fdaPDE R","time fdaPDE C++"), 
      col=c("red","blue","magenta"),pch=c(15,16,17),cex = 0.8,text.font=4 )

x11()
plot(1:length(lambdas), err_[,1],type="b", pch=5 ,col="red",xlab='',ylab='|| beta_ex - beta_hat ||',
        ylim=c(min(err_),max(err_)))
points(1:length(lambdas), err_[,2],type="b",pch=20, col="blue")
legend("topleft", legend("Mattina","fdaPDE"), col=c("red","blue"), pch=c(5,20))
as.matrix(abs(err_[,1] - err_[,2]))

########################
### decreasing delta ###
## fixed lambda ##
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
L = as.linnet(simplenet)


#num subdivsion
N = 5

# || beta.fdaPDE - beta || 
beta.err = matrix(0,nrow=N,ncol=1)
# || beta.* - beta_ex ||
err_ = matrix(0, nrow=N, ncol=2) 

times = matrix(0,nrow=N, ncol=3)
options(digits.sec=6)
lambda = 1e-4


deltas= 0.0625 * 0.5^(1:N)
segments_= as.matrix(0,nrow=N,ncol=1)

for( i in 1:N){
  delta = deltas[i]
  mesh = create.mesh.1D.vertices(vertices, edges, delta)
  FEMbasis = create.FEM.basis.1D(mesh)
  
  nodes_x = mesh$nodes[,1]
  nodes_y = mesh$nodes[,2]
  
  #locations == nodes
  nobs = mesh$nnodes
  points_x = nodes_x
  points_y = nodes_y
  
  ### fdaPDE - MESH ### 
  mesh.fdaPDE = fdaPDE::create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
  FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)
  
  
  ### CAMPO f ### 
  aux = function(x, y, seg, tp) { cos((2*x+y)/4) +  ((x+y))^2 }
  my_dens <- linfun(aux, L)
  
  ### observations ###
  W1 = rnorm(nobs, 1, sd=1.5)
  W2 = rnorm(nobs, 2, sd=2)
  W = NULL
  W = cbind(W1,W2)
  
  beta1 = -0.5
  beta2 = 0.2
  
  temp = range(beta1*W1 + beta2*W2 + my_dens(points_x, points_y))
  my_range = temp[2]-temp[1]
  EPS = rnorm(nobs, 0, sd=sqrt(0.05*my_range) )
  Z = beta1*W1 + beta2*W2 + my_dens(points_x, points_y) + EPS
  PHI=diag(nobs)
  beta_ex = rbind(beta1,beta2)
  
  start.time <- Sys.time()
  output = covariate.FEM.1D(  W, Z, FEMbasis, PHI, lambda)
  end.time  <-Sys.time()
  
  start.time.fdaPDE <- Sys.time()
  output.fdaPDE = smooth.FEM( observations=Z, 
                              FEMbasis=FEMbasis.fdaPDE, 
                              covariates=W, lambda=lambda )
  end.time.fdaPDE  <-Sys.time()
  
  beta.err[i]   = norm(output$beta_estimate - output.fdaPDE$solution$beta, type = "2")
  err_[i,1]     = norm(beta_ex - output$beta_estimate, type="2")
  err_[i,2]     = norm(beta_ex - output.fdaPDE$solution$beta,type="2")
  
  times[i,1] =  end.time - start.time
  times[i,2] =  end.time.fdaPDE - start.time.fdaPDE
  times[i,3] =  output.fdaPDE$time
  
  segments_[i]=dim(mesh.fdaPDE$edges)[1]
}

x11()
plot(1:N, times[,1],type="b",pch = 15,col="red", xlab='',ylab='times', 
     ylim=c(min(times),max(times)),asp=2)
points(1:N, times[,2],type="b",pch = 16, col="blue")
points(1:N, times[,3],type="b",pch = 17, col="magenta")
legend("topleft",legend=c("time Mattina", "time fdaPDE R","time fdaPDE C++"), 
       col=c("red","blue","magenta"),pch=c(15,16,17),cex = 0.8,text.font=4 )

x11()
plot(1:length(1:N), err_[,1],type="b", pch=5 ,col="red",xlab='',ylab='|| beta_ex - beta_hat ||',
     ylim=c(min(err_),max(err_)))
points(1:N, err_[,2],type="b",pch=20, col="blue")
legend("topleft", legend=c("Mattina","fdaPDE"), col=c("red","blue"), pch=c(5,20))
as.matrix(abs(err_[,1] - err_[,2]))

# order (?)
p.fdaPDE = log2(err_[1:N-1,2]/err_[2:N,2])
p.fdaPDE

### solo fda.PDE ### 

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
L = as.linnet(simplenet)


#num subdivsion
N = 8

# || beta.fdaPDE - beta || 
beta.err = matrix(0,nrow=N,ncol=1)
# || beta.* - beta_ex ||
err_ = matrix(0, nrow=N, ncol=1) 
mise = matrix(0,nrow=N,ncol=1)
times = matrix(0,nrow=N, ncol=2)
options(digits.sec=6)
lambda = 1e-4


deltas= 0.0625 * 0.5^(1:N)
segments_= as.matrix(0,nrow=N,ncol=1)

for( i in 1:N){
  delta = deltas[i]
  mesh = create.mesh.1D.vertices(vertices, edges, delta)
  #locations == nodes
  nobs = mesh$nnodes
  points_x = mesh$nodes[,1]
  points_y = mesh$nodes[,2]
  
  ### fdaPDE - MESH ### 
  mesh.fdaPDE = fdaPDE::create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
  FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)
  
  
  ### CAMPO f ### 
  aux = function(x, y, seg, tp) { cos((2*x+y)/4) +  ((x+y))^2 }
  my_dens <- linfun(aux, L)
  
  ### observations ###
  W1 = rnorm(nobs, 1, sd=1.5)
  W2 = rnorm(nobs, 2, sd=2)
  W = NULL
  W = cbind(W1,W2)
  
  beta1 = -0.5
  beta2 = 0.2
  
  temp = range(beta1*W1 + beta2*W2 + my_dens(points_x, points_y))
  my_range = temp[2]-temp[1]
  EPS = rnorm(nobs, 0, sd=sqrt(0.05*my_range) )
  Z = beta1*W1 + beta2*W2 + my_dens(points_x, points_y) + EPS
  beta_ex = rbind(beta1,beta2)
  
  start.time.fdaPDE <- Sys.time()
  output.fdaPDE = smooth.FEM( observations=Z, 
                              FEMbasis=FEMbasis.fdaPDE, 
                              covariates=W, lambda=lambda )
  end.time.fdaPDE  <-Sys.time()
  

  err_[i,1]     = norm(beta_ex - output.fdaPDE$solution$beta,type="2")
  mise[i,1]     = norm( my_dens(points_x,points_y) - output.fdaPDE$fit.FEM$coeff,type="2")
  times[i,1] =  end.time.fdaPDE - start.time.fdaPDE
  times[i,2] =  output.fdaPDE$time
  
  segments_[i]=dim(mesh.fdaPDE$edges)[1]
}

x11()
plot(1:N, times[,1],type="b",pch = 15,col="red", xlab='',ylab='times', 
     ylim=c(min(times),max(times)),asp=1)
points(1:N, times[,2],type="b",pch = 16, col="blue")
legend("topleft",legend=c( "time fdaPDE R","time fdaPDE C++"), 
       col=c("red","blue","magenta"),pch=c(15,16),cex = 0.8,text.font=4 )

x11()
plot(1:length(1:N), err_[,1],type="b", pch=5 ,col="red",xlab='',ylab='|| beta_ex - beta_hat ||',
     ylim=c(min(err_),max(err_)))
legend("topleft", legend="fdaPDE", col="red", pch=5)
# order (?)
p.fdaPDE = log2(err_[1:N-1,1]/err_[2:N,1])
p.fdaPDE

