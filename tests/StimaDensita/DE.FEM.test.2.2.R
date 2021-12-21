# RightCV vs SimplifiedCV #
# RightCV Gradient descent with fixed step vs BFGS with fixed step #

# DE test 2 #

source("~/Scrivania/fdaPDE/tests/integrate_f.R")
library(spatstat)
data("dendrite")

vertices = cbind(dendrite$domain$vertices$x, dendrite$domain$vertices$y) 
edges = cbind(dendrite$domain$from, dendrite$domain$to)
L = as.linnet(dendrite$domain)
M = dendrite$domain$m

mesh = create.mesh.1.5D(nodes=vertices,edges=edges)
FEMbasis = create.FEM.basis(mesh=mesh)

plot(mesh)
points(mesh$nodes[600,1],mesh$nodes[600,2],pch=16)
points(mesh$nodes[70,1],mesh$nodes[70,2],pch=16)
points(mesh$nodes[10,1],mesh$nodes[10,2],pch=16)
# density
mu1 = mesh$nodes[70,]
mu2 = mesh$nodes[600,]
mu3 = mesh$nodes[10,]
sigmax = 20
sigmay = 20
aux = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], x3=mu3[1], y3=mu3[2], sigma_x=sigmax, sigma_y=sigmay) 
{   ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))+ 
    ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x3)^2 + (1/(2*sigma_y^2))*(y-y3)^2 ))}


integral_ex = integrate_f( FEM(aux(mesh$nodes[,1],mesh$nodes[,2]),FEMbasis) )

aux = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], x3=mu3[1], y3=mu3[2], sigma_x=sigmax, sigma_y=sigmay) 
{   1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))+
    1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x3)^2 + (1/(2*sigma_y^2))*(y-y3)^2 ))}
my_dens <- linfun(aux, L) # FUNZIONE SU UN NETWORK

coef.ex = aux(mesh$nodes[,1],mesh$nodes[,2])
plot(my_dens)

niter = 25
nfolds=5
 
nnodes = dim(mesh$nodes)[1]
###
nobs = c(100,150,200,250,500)
#nobs = c(150)
N = length(nobs)
#cvec=matrix(0,nrow=mesh$nnodes,ncol=1)
# 
lambda = c(1e-4,1e-3,1e-2,1e-1,1) 

M=30
mise.RightCV.Grad = matrix(0,nrow=M,ncol=N)
mise.RightCV.BFGS = matrix(0,nrow=M,ncol=N)
mise.SimplifiedCV.Grad = matrix(0,nrow=M,ncol=N) # RightCV vs Smplified con Gradient (?)
mise.SimplifiedCV.BFGS = matrix(0,nrow=M,ncol=N)

err.L2.RightCV.Grad.vs.SimplifiedCV.Grad= matrix(0,nrow=M,ncol=N)
err.L2.RightCV.BFGS.vs.SimplifiedCV.BFGS = matrix(0,nrow=M,ncol=N)

times.RightCV.Grad = matrix(0,nrow=M, ncol=N)
times.RightCV.BFGS = matrix(0,nrow=M,ncol=N)
times.SimplifiedCV.Grad = matrix(0,nrow=M,ncol=N)
times.SimplifiedCV.BFGS = matrix(0,nrow=M,ncol=N)

sols.RightCV.Grad = array(0,dim=c(M,N,nnodes) ) #contiene vettore dei coeff
sols.RightCV.BFGS = array(0,dim=c(M,N,nnodes) )
sols.SimplifiedCV.Grad = array(0,dim=c(M,N,nnodes) )
sols.SimplifiedCV.BFGS = array(0,dim=c(M,N,nnodes) )

lambda.RightCV.Grad = matrix(0,nrow=M,ncol=N)
lambda.RightCV.BFGS = matrix(0,nrow=M,ncol=N)
lambda.SimplifiedCV.Grad = matrix(0,nrow=M,ncol=N)
lambda.SimplifiedCV.BFGS = matrix(0,nrow=M,ncol=N)

tot.start = Sys.time()
for(i in 1:N){
  for(j in 1:M){
    print(paste("##### nobs = ",nobs[i]," #####",sep=""))
    print(paste("### ",j,"/",M," ###",sep=""))
    PPP = rlpp(nobs[i], my_dens)
    points_x = PPP$data$x
    points_y = PPP$data$y
    points = NULL 
    points = cbind(points_x, points_y)
    points = as.matrix(points)
    colnames(points) = NULL
    
    # RightCV Gradient
    print("# RightCV - Gradient #")
    start <- Sys.time()
    sol = DE.FEM(data=points,FEMbasis=FEMbasis,lambda=lambda, 
                        nsimulations=niter, 
                        #fvec = exp(cvec),
                        tol1=1e-5,
                        nfolds=nfolds, 
                        step_method="Fixed_Step", 
                        direction_method="Gradient",
                        preprocess_method = "RightCV")
    end <- Sys.time()
    coef = exp(sol$g)
    sols.RightCV.Grad[j,i,] = coef 
    mise.RightCV.Grad[j,i]   = integrate_f( FEM( (coef - coef.ex)^2, FEMbasis)) 
    times.RightCV.Grad[j,i] = difftime(end, start, units="secs")
    lambda.RightCV.Grad[j,i] = sol$lambda
    
    # RightCV BFGS
    print("# RightCV - BFGS #")
    start <- Sys.time()
    sol = DE.FEM(data=points,FEMbasis=FEMbasis,lambda=lambda, 
                 nsimulations=niter, 
                 #fvec = exp(cvec),
                 tol1 = 1e-5,
                 nfolds=nfolds, 
                 step_method="Fixed_Step", 
                 direction_method="BFGS",
                 preprocess_method = "RightCV")
    end <- Sys.time()
    coef = exp(sol$g)
    sols.RightCV.BFGS[j,i,] = coef 
    mise.RightCV.BFGS[j,i]   = integrate_f( FEM( (coef - coef.ex)^2, FEMbasis)) 
    times.RightCV.BFGS[j,i] = difftime(end, start, units="secs")
    lambda.RightCV.BFGS[j,i] = sol$lambda
    
    # SimplifiedCV Gradient
    print("# SimplifiedCV - Gradient #")
    start <- Sys.time()
    sol = DE.FEM(data=points,FEMbasis=FEMbasis,lambda=lambda, 
                 nsimulations=niter, 
                 #fvec = exp(cvec),
                 tol1 = 1e-5,
                 nfolds=nfolds, 
                 step_method="Fixed_Step", 
                 direction_method="Gradient",
                 preprocess_method = "SimplifiedCV")
    end <- Sys.time()
    coef = exp(sol$g)
    sols.SimplifiedCV.Grad[j,i,] = coef 
    mise.SimplifiedCV.Grad[j,i]   = integrate_f( FEM( (coef - coef.ex)^2, FEMbasis)) 
    times.SimplifiedCV.Grad[j,i] = difftime(end, start, units="secs")
    lambda.SimplifiedCV.Grad[j,i] = sol$lambda
    
    # SimplifiedCV BFGS
    print("# SimplifiedCV - BFGS #")
    start <- Sys.time()
    sol = DE.FEM(data=points,FEMbasis=FEMbasis,lambda=lambda, 
                 nsimulations=niter, 
                 #fvec = exp(cvec),
                 tol1 = 1e-5,
                 nfolds=nfolds, 
                 step_method="Fixed_Step", 
                 direction_method="BFGS",
                 preprocess_method = "SimplifiedCV")
    end <- Sys.time()
    coef = exp(sol$g)
    sols.SimplifiedCV.BFGS[j,i,] = coef 
    mise.SimplifiedCV.BFGS[j,i]   = integrate_f( FEM( (coef - coef.ex)^2, FEMbasis)) 
    times.SimplifiedCV.BFGS[j,i] = difftime(end, start, units="secs")
    lambda.SimplifiedCV.BFGS[j,i] = sol$lambda
    
    
    err.L2.RightCV.Grad.vs.SimplifiedCV.Grad[j,i] = 
            integrate_f( FEM( (sols.RightCV.Grad[j,i,] - sols.SimplifiedCV.Grad[j,i,])^2, FEMbasis))
                    
    err.L2.RightCV.BFGS.vs.SimplifiedCV.BFGS[j,i] = 
            integrate_f( FEM( (sols.RightCV.BFGS[j,i,] - sols.SimplifiedCV.BFGS[j,i,])^2, FEMbasis))
  
  }
}
tot.end = Sys.time()
tot.time = difftime(tot.end,tot.start, units="secs")


today_ = Sys.Date()
init_ = "-fdaPDE" # "-null_vector"
ntest_ = "-test.2.3"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")


save(nobs, N, M, coef.ex,niter,
     mise.RightCV.Grad, mise.RightCV.BFGS,
     mise.SimplifiedCV.Grad, mise.SimplifiedCV.BFGS,
     err.L2.RightCV.Grad.vs.SimplifiedCV.Grad,
     err.L2.RightCV.BFGS.vs.SimplifiedCV.BFGS,
     times.RightCV.Grad, times.RightCV.BFGS,
     times.SimplifiedCV.Grad, times.SimplifiedCV.BFGS,
     sols.RightCV.Grad, sols.RightCV.BFGS,
     sols.SimplifiedCV.Grad, sols.SimplifiedCV.BFGS,
     lambda.RightCV.Grad, lambda.RightCV.BFGS,
     lambda.SimplifiedCV.Grad, lambda.SimplifiedCV.BFGS,
     file = save.file)
