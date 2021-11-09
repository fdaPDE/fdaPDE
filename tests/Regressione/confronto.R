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
delta = 0.03
mesh = create.mesh.1D.vertices(vertices, edges, delta)
FEMbasis = create.FEM.basis.1D(mesh)

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

### cross - validation ###
lambda_vector = c(10^-10,10^-9,10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,10,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10)

beta1=-0.5
beta2=0.2
beta_ex = rbind(beta1,beta2)

nobs=c(50,75,100,150)
N=length(nobs)
M=30
# || beta.fdaPDE - beta || 
beta.err = matrix(0,nrow=N,ncol=1)
# || beta.* - beta_ex ||
err_ = matrix(0, nrow=N, ncol=2) 

times = matrix(0,nrow=N, ncol=3)
options(digits.sec=6)
epsilon=0.0001
for(i in 1:N){
  PPP = rlpp(nobs[i], my_dens)
  points_x = PPP$data$x
  points_y = PPP$data$y
  points = NULL 
  points = cbind(points_x, points_y)
  points = as.matrix(points)
  colnames(points) = NULL
  
  ### observations ###
  W1 = rnorm(nobs[i], 0, sd=1.)
  W2 = rnorm(nobs[i], 0, sd=0.5)
  W = NULL
  W = cbind(W1,W2)
  temp = range(beta1*W1 + beta2*W2 + my_dens(points_x, points_y))
  my_range = temp[2]-temp[1]
  EPS = rnorm(nobs[i], 0, sd=sqrt(0.05*my_range) )
  Z = beta1*W1 + beta2*W2 + my_dens(points_x, points_y) + EPS
 
  start.time.fdaPDE <- Sys.time()
  output.fdaPDE = smooth.FEM( observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W,
                              locations=points,
                              lambda=lambda_vector, 
                              lambda.selection.criterion='grid', 
                              DOF.evaluation='exact', 
                              lambda.selection.lossfunction='GCV' )
  end.time.fdaPDE  <-Sys.time()
  
  lambda = output.fdaPDE$optimization$lambda_solution 
  result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
  PHI = result$PHI
  
  start.time <- Sys.time()
  output = covariate.FEM.1D(  W, Z, FEMbasis, PHI, lambda)
  end.time  <-Sys.time()
  
  beta.err[i]   = norm(output$beta_estimate - output.fdaPDE$solution$beta, type = "2")
  err_[i,1]     = norm(beta_ex - output$beta_estimate, type="2")
  err_[i,2]     = norm(beta_ex - output.fdaPDE$solution$beta,type="2")
  
  times[i,1] =  end.time - start.time
  times[i,2] =  end.time.fdaPDE - start.time.fdaPDE
  times[i,3] =  output.fdaPDE$time
  }



x11()
plot(nobs, times[,1],type="b",pch = 15,col="red", xlab='num obs',ylab='times', 
        ylim=c(min(times),max(times)) )
points(nobs, times[,2],type="b",pch = 16, col="blue")
points(nobs, times[,3],type="b",pch = 17, col="magenta")
legend("topleft",legend=c("time Mattina", "time fdaPDE R","time fdaPDE C++"), 
      col=c("red","blue","magenta"),pch=c(15,16,17),cex = 0.8,text.font=4 )

x11()
plot(nobs, err_[,1],type="b", pch=5 ,col="red",xlab='num obs',ylab='|| beta_ex - beta_hat ||',
        ylim=c(min(err_),max(err_)))
points(nobs, err_[,2],type="b",pch=20, col="blue")
legend("topleft", legend=c("Mattina","fdaPDE"), col=c("red","blue"), pch=c(5,20))
as.matrix(abs(err_[,1] - err_[,2]))

########################

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

delta = 0.03
mesh = create.mesh.1D.vertices(vertices, edges, delta)
FEMbasis = create.FEM.basis.1D(mesh)

### fdaPDE - MESH ### 
mesh.fdaPDE = fdaPDE::create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)


### CAMPO f ### 
aux = function(x, y, seg, tp) { cos((2*x+y)/4) +  ((x+y))^2 }
my_dens <- linfun(aux, L)

#num subdivsion
#nobs = c(200,300,400,500,600,700,800)
nobs  = c(50,100,150,200,250,300,350)
N=length(nobs)
M=30

beta1 = -0.5
beta2 = 0.2
beta_ex = rbind(beta1,beta2)
# || beta.fdaPDE - beta || 
rmse.beta.1 = matrix(0,nrow=M,ncol=N)
rmse.beta.2 = matrix(0,nrow=M,ncol=N)
rmse.fdaPDE.beta.1 = matrix(0,nrow=M,ncol=N)
rmse.fdaPDE.beta.2 = matrix(0,nrow=M,ncol=N)

times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)
options(digits.sec=6)
lambda_vector = c(10^-10,10^-9,10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,10,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10)
epsilon=0.0001
for( i in 1:N){
  for( j in 1:M){
    PPP = rlpp(nobs[i], my_dens)
    points_x = PPP$data$x
    points_y = PPP$data$y
    points = NULL 
    points = cbind(points_x, points_y)
    points = as.matrix(points)
    colnames(points) = NULL
  
   ### observations ###
    W1 = rnorm(nobs[i], 0, sd=1.)
    W2 = rnorm(nobs[i], 0, sd=0.5)
    W = NULL
    W = cbind(W1,W2)
  
    temp = range(beta1*W1 + beta2*W2 + my_dens(points_x, points_y))
    my_range = temp[2]-temp[1]
    EPS = rnorm(nobs[i], 0, sd=sqrt(0.05*my_range) )
    Z = beta1*W1 + beta2*W2 + my_dens(points_x, points_y) + EPS
  
    output.fdaPDE = smooth.FEM( observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W,
                              locations=points,
                              lambda=lambda_vector, 
                              lambda.selection.criterion='grid', 
                              DOF.evaluation='exact', 
                              lambda.selection.lossfunction='GCV' )
   
    lambda = output.fdaPDE$optimization$lambda_solution
    
    output.fdaPDE = smooth.FEM(observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W,
                               locations=points,
                               lambda=lambda)
    
    #result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    #PHI = result$PHI
  
    start.time <- Sys.time()
    result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    PHI = result$PHI
    output = covariate.FEM.1D(  W, Z, FEMbasis, PHI, lambda)
    end.time  <-Sys.time()
    
    rmse.beta.1[j,i] = norm(beta_ex[1]- output$beta_estimate[1],type="2")
    rmse.beta.2[j,i] = norm(beta_ex[2]- output$beta_estimate[2],type="2")
    rmse.fdaPDE.beta.1[j,i] = norm(beta_ex[1]- output.fdaPDE$solution$beta[1],type="2")
    rmse.fdaPDE.beta.2[j,i] = norm(beta_ex[2]- output.fdaPDE$solution$beta[2],type="2")
  
    times[j,i] = end.time - start.time
    times.fdaPDE[j,i] = output.fdaPDE$time
  }
}

x11()
boxplot(rmse.beta.1[,1],rmse.fdaPDE.beta.1[,1],
        rmse.beta.1[,2],rmse.fdaPDE.beta.1[,2],
        rmse.beta.1[,3],rmse.fdaPDE.beta.1[,3],
        rmse.beta.1[,4],rmse.fdaPDE.beta.1[,4],
        #main="beta1",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11),
        ##rmse.beta.1[,5],rmse.fdaPDE.beta.1[,5],
        ##rmse.beta.1[,6],rmse.fdaPDE.beta.1[,6],
        ##rmse.beta.1[,7],rmse.fdaPDE.beta.1[,7],
        ##main="beta2",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11,13,15,17,18,20,21),
        col=c("red","blue"))
legend("topright", legend=c("Mattina","fdaPDE"), col=c("red","blue"),pch=c(15,15) )

x11()
boxplot(rmse.beta.2[,1],rmse.fdaPDE.beta.2[,1],
        rmse.beta.2[,2],rmse.fdaPDE.beta.2[,2],
        rmse.beta.2[,3],rmse.fdaPDE.beta.2[,3],
        rmse.beta.2[,4],rmse.fdaPDE.beta.2[,4],
        main="beta2",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11),
         # rmse.beta.2[,5],rmse.fdaPDE.beta.2[,5],
         # rmse.beta.2[,6],rmse.fdaPDE.beta.2[,6],
         # rmse.beta.2[,7],rmse.fdaPDE.beta.2[,7],
       # main="beta2",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11,13,15,17,18,20,21),
        col=c("red","blue"))
legend("topright", legend=c("Mattina","fdaPDE"), col=c("red","blue"),pch=c(15,15) )

time.plot = matrix(0,nrow=1,ncol=N)
time.fdaPDE.plot = matrix(0,nrow=1,ncol=N)

for(i in 1:N){
  time.plot[i] = mean(times[,i])
  time.fdaPDE.plot[i] = mean(times.fdaPDE[,i])
}

x11()
plot(nobs, time.plot,type="b",col="red",pch=15 ,xlab = "nobs",ylab = "times", main="mean time." ,ylim=c(min(time.plot),max(time.plot)))
points(nobs, time.fdaPDE.plot,type="b",pch = 16, col="blue")
legend("topleft",legend=c("eval+algo Mattina", "fdaPDE"), 
       col=c("red","blue"),pch=c(15,16),cex = 0.8,text.font=4 )


x11()
boxplot(times[,1],times.fdaPDE[,1],
        times[,2],times.fdaPDE[,2],
        times[,3],times.fdaPDE[,3],
        times[,4],times.fdaPDE[,4],
        #main="time",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11),
        times[,5],times.fdaPDE[,5],
        times[,6],times.fdaPDE[,6],
        times[,7],times.fdaPDE[,7],
        main="times",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20),
        col=c("red","blue"))
legend("topright", legend=c("eval+algo Mattina","fdaPDE"), col=c("red","blue"),pch=c(15,15) )