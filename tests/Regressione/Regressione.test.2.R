source("~/Scrivania/fdaPDE/tests/Regressione/Regressione_Mattina.R")

data("spiders")
vertices = cbind(spiders$domain$vertices$x, spiders$domain$vertices$y)
edges = cbind(spiders$domain$from, spiders$domain$to)

M = spiders$domai$m
find.boundary(M)
L = as.linnet(spiders$domain)

#412 nodes -> delta = 40
delta =40
mesh = create.mesh.1D.vertices(vertices, edges, delta)
mesh$nnodes
FEMbasis = create.FEM.basis.1D(mesh)

### fdaPDE - MESH ### 
mesh.fdaPDE = create.mesh.1.5D(mesh$nodes,mesh$segments)
FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)



### CAMPO f ### 
source("~/Scrivania/fdaPDE/tests/integrate_f.R")

plot(mesh.fdaPDE)
points(mesh$nodes[350,1],mesh$nodes[350,2],pch=16)
points(mesh$nodes[35,1],mesh$nodes[35,2],pch=16)
points(mesh$nodes[196,1],mesh$nodes[196,2],pch=16)
# density
mu1 = mesh$nodes[350,]
mu2 = mesh$nodes[35,]
mu3 = mesh$nodes[196,]
sigmax = 100
sigmay = 100
aux = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], x3=mu3[1], y3=mu3[2], sigma_x=sigmax, sigma_y=sigmay) 
{   exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))+
    exp(- ( (1/(2*sigma_x^2))*(x-x3)^2 + (1/(2*sigma_y^2))*(y-y3)^2 ))}

my_dens <- linfun(aux, L) # FUNZIONE SU UN NETWORK
plot(my_dens)

coef.ex = aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])

nobs = c(300)
N=length(nobs)
M=30

beta_ex = 0.8 #test-2

# || beta.fdaPDE - beta || 
rmse.beta = matrix(0,nrow=M,ncol=N)
rmse.fdaPDE.beta = matrix(0,nrow=M,ncol=N)

mise = matrix(0,nrow=M,ncol=N)
mise.fdaPDE = matrix(0,nrow=M,ncol=N)
# (||f_{fdaPDE} - f{R} ||_{L2} )^2
err.L2 = matrix(0,nrow=M,ncol=N)

norms.beta.inf = matrix(0,nrow=M,ncol=N)
norms.beta.2 =matrix(0,nrow=M,ncol=N)

sols = array(0,dim=c(M,N,mesh$nnodes) ) #contiene vettore dei coeff
sols.fdaPDE = array(0,dim=c(M,N,mesh$nnodes) )

times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)
times.fdaPDE.2 = matrix(0,nrow=M,ncol=N)

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
    
    W = matrix( rnorm(nobs[i], 0, sd=0.5), nrow=nobs[i], ncol=1)
    temp = range(beta_ex*W +  my_dens(points_x, points_y))
    my_range = temp[2]-temp[1]
    EPS = rnorm(nobs[i], 0, sd=sqrt(0.05*my_range) )
    Z = beta_ex*W + my_dens(points_x, points_y) + EPS
    
    output.fdaPDE = smooth.FEM( observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W,
                                locations=points,
                                lambda=lambda_vector, 
                                lambda.selection.criterion='grid', 
                                DOF.evaluation='exact', 
                                lambda.selection.lossfunction='GCV' )
    
    lambda = output.fdaPDE$optimization$lambda_solution
    
    start.fdaPDE = Sys.time()
    output.fdaPDE = smooth.FEM(observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W,
                               locations=points,
                               lambda=lambda)
    end.fdaPDE = Sys.time()
    
    start.time <- Sys.time()
    result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    PHI = result$PHI
    output = covariate.FEM.1D(  W, Z, FEMbasis, PHI, lambda)
    end.time  <-Sys.time()
    
    rmse.beta[j,i] = norm(beta_ex[1]- output$beta_estimate[1],type="2")
    rmse.fdaPDE.beta[j,i] = norm(beta_ex[1]- output.fdaPDE$solution$beta[1],type="2")
    
    norms.beta.inf[j,i] = abs( output$beta_estimate[1] - output.fdaPDE$solution$beta[1])
    norms.beta.2[j,i]   = norm(output$beta_estimate[1] - output.fdaPDE$solution$beta[1],type="2")
    
    times[j,i] = difftime(end.time, start.time, units="secs")
    
    times.fdaPDE[j,i] = output.fdaPDE$time
    times.fdaPDE.2[j,i] = difftime(end.fdaPDE, start.fdaPDE, units="secs")
    coef = output$f_coeff
    coef.fdaPDE     = output.fdaPDE$fit.FEM$coeff
    
    sols[j,i,] =  coef
    sols.fdaPDE[j,i,] = coef.fdaPDE
    
    mise[j,i]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
    mise.fdaPDE[j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
    err.L2[j,i] = integrate_f( FEM((coef.fdaPDE-coef)^2,FEMbasis.fdaPDE) )
  }
}

today_ = Sys.Date()
ntest_ = "-test-2-sample-2"  

file.name = paste("Regression-",today_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

save(nobs,M,N,
     rmse.beta,
     rmse.fdaPDE.beta,
     mise, mise.fdaPDE, err.L2,
     sols,sols.fdaPDE,
     norms.beta.inf,norms.beta.2,
     times,times.fdaPDE,times.fdaPDE.2, file =save.file)

