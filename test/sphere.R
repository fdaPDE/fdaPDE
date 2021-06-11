# Function to generate random points in a sphere
rsphere <- function(n, r = 0.9, surface_only = FALSE) {
  phi       <- runif(n, 0.0, 2.0 * pi)
  cos_theta <- runif(n, -1.0, 1.0)
  sin_theta <- sqrt((1.0-cos_theta)*(1.0+cos_theta))
  radius <- r
  if (surface_only == FALSE) {
    radius <- r * runif(n, 0.0, 1.0)^(1.0/3.0)
  }
  
  x <- radius * sin_theta * cos(phi)
  y <- radius * sin_theta * sin(phi)
  z <- radius * cos_theta
  
  cbind(x, y, z)
}

# Build mesh: Sphere
data("sphere3Ddata")
mesh<-create.mesh.3D(sphere3Ddata$nodes,
                            sphere3Ddata$tetrahedrons,
                            order=1)

FEMbasis <- create.FEM.basis(mesh)

set.seed(5847947)

a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

# Exact test function
nodes=mesh$nodes
nnodes = nrow(mesh$nodes)

# Set smoothing parameter
lambda=10^(-12:3)

##### no covariate case #####

# Evaluate exact solution on mesh nodes
coeffCov =  sin(2*pi*mesh$nodes[,1]) +  2 * sin(2*pi*mesh$nodes[,2]) +  sin(2*pi*mesh$nodes[,3])

# Plot exact solution
plot(FEM(coeffCov,FEMbasis))

dataCov = coeffCov
dataCov1 = dataCov + rnorm(nrow(mesh$nodes), mean=0, sd=0.001*diff(range(coeffCov)))


func = function(x)
{
  a1*sin(2*pi*(x[,1]*x[,2]))+a2*cos(2*pi*x[,2])+a3*cos(2*pi*x[,3])
}

solnoCov = smooth.FEM(observations = dataCov, FEMbasis = FEMbasis, 
                    lambda = lambda)
solnoCov1 = smooth.FEM(observations = dataCov1, FEMbasis = FEMbasis, 
                     lambda = lambda)

sol_approxnoCov = solnoCov$fit.FEM
femsol_approxnoCov<-t(eval.FEM(sol_approxnoCov,nodes))
sol_approx1noCov = solnoCov1$fit.FEM
femsol_approx1noCov<-t(eval.FEM(sol_approx1noCov,nodes))

errornoCov <- matrix(nrow=length(dataCov), ncol= length(lambda))
for (i in 1:length(lambda))
  errornoCov[,i] <- abs(coeffCov-t(femsol_approxnoCov)[,i])
errnoCov2 <- colMeans(errornoCov^2)
plot(log10(lambda), errnoCov1, pch=4, xlab='log_10(lambda)', ylab='error', col=1)
grid()
abline(v=log10(best_lambdanoCov1), lty=2, col='red')
lines(-12:3, errnoCov1, col=1)
points(log10(best_lambdanoCov1), errnoCov1[11], pch=19, col='red')
points(-12:3, errnoCov2, pch=4, col=2)
lines(-12:3, errnoCov2, col=2)
points(log10(best_lambdanoCov2), errnoCov2[11], pch=3, col=2)
lines(-12:3, errnoCov3, col=3)
points(-12:3, errnoCov4, pch=10, col=4)
points(log10(best_lambdanoCov3), errnoCov3[11], pch=1, col=3)
legend('topleft',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

residualnoCov <- matrix(nrow=dim(nodes)[1], ncol= length(lambda))
residualnoCov <- abs(t(femsol_approx1noCov)-t(femsol_approxnoCov))
resnoCov5 <- colMeans(residualnoCov)
plot(-12:3, resnoCov1, pch=4, xlab='log_10(lambda)', ylab='residual', col=1, ylim=c(min(resnoCov1),0.01))
grid()
abline(v=log10(best_lambdanoCov1), lty=2, col='red')
lines(-12:3, resnoCov1, col=1)
points(log10(best_lambdanoCov1), resnoCov1[10], pch=19, col=1)
points(-12:3, resnoCov5, pch=4, col=2)
lines(-12:3, resnoCov2, col=2)
points(log10(best_lambdanoCov2), resnoCov2[10], pch=19, col=2)
points(-12:3, resnoCov3, pch=4, col=3)
lines(-12:3, resnoCov3, col=3)
points(log10(best_lambdanoCov3), resnoCov3[10], pch=1, col=3)
legend('topright',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

solnoCov <- smooth.FEM(observations = dataCov,
                     FEMbasis = FEMbasis,
                     lambda.selection.criterion = 'newton', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
sol_approxnoCov <- solnoCov$fit.FEM
best_lambdanoCov5 <- solnoCov$optimization$lambda_solution
RMSE<-function(f,g) sqrt(mean((f-g)^2))
loc_eval = rsphere(5000)
rmsenoCov =NULL
for (i in 1:20)
{
  solnoCov <- smooth.FEM(observations = dataCov+rnorm(nrow(mesh$nodes), mean=0, sd=0.01*diff(range(coeffCov))),
                         FEMbasis = FEMbasis,
                         lambda.selection.criterion = 'newton', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
  sol_approxnoCov <- solnoCov$fit.FEM
  best_lambdanoCov5 <- c(best_lambdanoCov5,solnoCov$optimization$lambda_solution)
  rmsenoCov <- c(rmsenoCov, RMSE(func(loc_eval),
                             eval.FEM(solnoCov$fit.FEM,loc_eval)))
}
boxplot(rmsenoCov)
rmsenoCovE = rmsenoCov

boxplot(rmsenoCovA,rmsenoCovB,rmsenoCovC,rmsenoCovD,names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4))


##### covariate case ######

cov1=(4*sin(2*pi*nodes[,2])+6*sin((2*pi*nodes[,3])^2))*(1-exp(-nodes[,1]))/3
cov2=1+2*nodes[,1]*sin(2*pi*nodes[,2])/6

W=cbind(cov1,cov2)

beta=c(0.7,0.2)

solCov = smooth.FEM(observations = dataCov+beta[1]*cov1+beta[2]*cov2, covariates = W, FEMbasis = FEMbasis, 
                                lambda = lambda)
solCov1 = smooth.FEM(observations = dataCov1+beta[1]*cov1+beta[2]*cov2, covariates = W, FEMbasis = FEMbasis, 
                      lambda = lambda)

sol_approxCov = solCov$fit.FEM
femsol_approxCov<-t(eval.FEM(sol_approxCov,nodes))
sol_approx1Cov = solCov1$fit.FEM
femsol_approx1Cov<-t(eval.FEM(sol_approx1Cov,nodes))

errorCov <- matrix(nrow=length(dataCov), ncol= length(lambda))
for (i in 1:length(lambda))
  errorCov[,i] <- abs(coeffCov-t(femsol_approxCov)[,i])
errCov2 <- colMeans(errorCov^2)
plot(log10(lambda), errCov1, pch=4, xlab='log_10(lambda)', ylab='error', col=1)
grid()
abline(v=log10(best_lambdaCov1), lty=2, col='red')
lines(-12:3, errCov1, col=1)
points(log10(best_lambdaCov1), errCov1[11], pch=19, col='red')
points(-12:3, errCov2, pch=4, col=2)
lines(-12:3, errCov2, col=2)
points(log10(best_lambdaCov2), errCov2[11], pch=3, col=2)
lines(-12:3, errCov3, col=3)
points(-12:3, errCov3, pch=4, col=3)
points(log10(best_lambdaCov3), errCov3[11], pch=1, col=3)
legend('topleft',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

residualCov <- matrix(nrow=dim(nodes)[1], ncol= length(lambda))
residualCov <- abs(t(femsol_approx1Cov)-t(femsol_approxCov))
resCov2 <- colMeans(residualCov)
plot(-12:3, resCov1, pch=4, xlab='log_10(lambda)', ylab='residual', col=1, ylim=c(0,0.0085))
grid()
abline(v=log10(best_lambdaCov1), lty=2, col='red')
lines(-12:3, resCov1, col=1)
points(log10(best_lambdaCov1), resCov1[10], pch=19, col=1)
points(-12:3, resCov2, pch=4, col=2)
lines(-12:3, resCov2, col=2)
points(log10(best_lambdaCov2), resCov2[10], pch=19, col=2)
points(-12:3, resCov3, pch=4, col=4)
lines(-12:3, resCov3, col=4)
points(log10(best_lambdaCov1), resCov3[10], pch=1, col=4)
legend('topright',legend=c('no preconditioner', 'mass lumping',' block preconditioner'), cex=0.62,col=c(1,2,4), pch=19)

solCov <- smooth.FEM(observations = dataCov+beta[1]*cov1+beta[2]*cov2, covariates = W,
                  FEMbasis = FEMbasis,
                  lambda.selection.criterion = 'newton', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
sol_approxCov <- solCov$fit.FEM
best_lambdaCov2 <- solCov$optimization$lambda_solution
RMSE<-function(f,g) sqrt(mean((f-g)^2))
rmseCov =NULL
for (i in 1:20)
{
  solCov <- smooth.FEM(observations = dataCov+beta[1]*cov1+beta[2]*cov2+rnorm(nrow(mesh$nodes), mean=0, sd=0.01*diff(range(coeffCov))), covariates = W,
                       FEMbasis = FEMbasis,
                       lambda.selection.criterion = 'newton', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
  sol_approxCov <- solCov$fit.FEM
  best_lambdaCov2 <- c(best_lambdaCov2,solCov$optimization$lambda_solution)
  rmseCov <- c(rmseCov, RMSE(func(loc_eval),
                       eval.FEM(solCov$fit.FEM,loc_eval)))
}
boxplot(rmseCov)
rmseCovB = rmseCov

boxplot(rmseCovA,rmseCovB,rmseCovC,rmseCovD, names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4),ylab='RMSE')


##### PDE penalization #####

# Set PDE parameters (in this case they are constant)
PDE_parameters_anys = list(K = diag(c(1,.5,1)), b = c(0,0,0), c = -4*pi^2)

# Evaluate exact solution on mesh nodes
coeff =  sin(2*pi*mesh$nodes[,1]) +  2 * sin(2*pi*mesh$nodes[,2]) +  sin(2*pi*mesh$nodes[,3])

# Plot exact solution
plot(FEM(coeff,FEMbasis))

data = coeff
data1 = data + rnorm(nrow(mesh$nodes), mean=0, sd=0.001*diff(range(coeff)))


# Compute the solution for each lambda
output_CPP <- smooth.FEM(observations = data, PDE_parameters = PDE_parameters_anys,
                         FEMbasis = FEMbasis, lambda = lambda)
output_CPP1 <- smooth.FEM(observations = data1, PDE_parameters = PDE_parameters_anys,
                         FEMbasis = FEMbasis, lambda = lambda)

sol_approx = output_CPP$fit.FEM
femsol_approx<-t(eval.FEM(sol_approx,nodes))
sol_approx1 = output_CPP1$fit.FEM
femsol_approx1<-t(eval.FEM(sol_approx1,nodes))

error <- matrix(nrow=length(data), ncol= length(lambda))
for (i in 1:length(lambda))
  error[,i] <- abs(coeff-t(femsol_approx)[,i])
err2 <- colMeans(error^2)
plot(log10(lambda), err1, pch=4, xlab='log_10(lambda)', ylab='error', col=1)
grid()
abline(v=log10(best_lambda1), lty=2, col='red')
lines(-12:3, err1, col=1)
points(log10(best_lambda1), err1[11], pch=19, col='red')
points(-12:3, err2, pch=4, col=2)
lines(-12:3, err2, col=2)
points(log10(best_lambda2), err2[11], pch=3, col=2)
lines(-12:3, err3, col=3)
points(-12:3, err3, pch=4, col=3)
points(log10(best_lambda3), err3[11], pch=1, col=3)
legend('topleft',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

residual <- matrix(nrow=dim(nodes)[1], ncol= length(lambda))
residual <- abs(t(femsol_approx1)-t(femsol_approx))
res2 <- colMeans(residual)
plot(-12:3, res1, pch=4, xlab='log_10(lambda)', ylab='residual', col=1)
grid()
abline(v=log10(best_lambda1), lty=2, col='red')
lines(-12:3, res1, col=1)
points(log10(best_lambda1), res1[11]/2+res1[10]/2, pch=19, col=1)
points(-12:3, res2, pch=4, col=2)
lines(-12:3, res2, col=2)
points(log10(best_lambda2), res2[11]/2+res2[10]/2, pch=1, col=2)
points(-12:3, res3, pch=4, col=3)
lines(-12:3, res3, col=3)
points(log10(best_lambda1), res2[11], pch=19, col=3)
legend('bottomleft',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

sol <- smooth.FEM(observations = data, PDE_parameters = PDE_parameters_anys,
                  FEMbasis = FEMbasis,
                  lambda.selection.criterion = 'newton', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
sol_approx <- sol$fit.FEM
best_lambda2 <- sol$optimization$lambda_solution
RMSE<-function(f,g) sqrt(mean((f-g)^2))
locations_eval = rsphere(5000)
rmse =NULL
for (i in 1:20)
{
  index_eval = sample(dim(locations_eval)[1])[1:500]
  loc_eval = locations_eval[index_eval,]
  # plot(mesh)
  # points3d(loc_eval,pch=19,col='red')
  rmse <- c(rmse, RMSE(sin(2*pi*loc_eval[,1])+sin(2*pi*loc_eval[,2])+sin(2*pi*loc_eval[,3]),
                       eval.FEM(sol$fit.FEM,loc_eval)))
}
boxplot(rmse)
rmseB = rmse

boxplot(rmseA,rmseB,rmseC,rmseD,names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4))

###### execution times #####
tD=microbenchmark(smooth.FEM(observations = data, FEMbasis = FEMbasis, lambda=best_lambda1))
tcovD=microbenchmark(smooth.FEM(observations = dataCov, covariates = W, FEMbasis = FEMbasis, lambda=best_lambdaCov1))
tnoncovD=microbenchmark(smooth.FEM(observations = dataCov,FEMbasis = FEMbasis, lambda=best_lambdanoCov1))

boxplot(log(tA$time),log(tB$time),log(tC$time),log(tD$time), names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')
boxplot(log(tcovA$time),log(tcovB$time),log(tcovC$time),log(tcovD$time), names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')
boxplot(log(tnoncovA$time),log(tnoncovB$time),log(tnoncovC$time),log(tnoncovD$time), names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')

boxplot(tA$time,tB$time,tC$time,tD$time, names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='time')
boxplot(tcovA$time,tcovB$time,tcovC$time,tcovD$time, names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='time')
boxplot(tnoncovA$time,tnoncovB$time,tnoncovC$time,tnoncovD$time, names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='time')

###############################################
############# TEST SCRIPT GAM 3D ##############
###############################################

####### 3D ########
#### Test 1: Sphere 3D GAMMA family ####
#            locations not at nodes
#            with covariates
#            no BC

library(fdaPDE)
library(purrr)
# family
FAMILY = "gamma"

link<-function(x){-1/x}
inv.link<-link

# beta
beta1=0.8
beta2=1.2
beta_exact=c(beta1,beta2)

# lambda
lambdaGLM= 10^seq(-1,3,length.out = 10)

# scale param
scale.param = 1

data(sphere3Ddata)

sphere3D<-create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)
FEMbasis <- fdaPDE::create.FEM.basis(sphere3D)
nodesLocations=sphere3D$nodes

set.seed(5847947)
set.seed(42)

# add evaluation points
new_grid3d = cbind(runif(1000,min = -1 , max = 1),runif(1000,min = -1 , max = 1),runif(1000,min = -1 , max = 1))
for(i in 1:length(new_grid3d[,1])){
  if(new_grid3d[i,1]^2 + new_grid3d[i,2]^2 + new_grid3d[i,3]^2 >= 0.98){
    new_grid3d[i,] <- c(NA,NA,NA)
  }
}

new_loc <- matrix(ncol = 3, nrow = sum(!is.na(new_grid3d[,1])))
new_loc[,1] <- new_grid3d[!is.na(new_grid3d[,1]),1] 
new_loc[,2] <- new_grid3d[!is.na(new_grid3d[,1]),2] 
new_loc[,3] <- new_grid3d[!is.na(new_grid3d[,1]),3] 

loc <- rbind(nodesLocations, new_loc)
dim(loc)

# Exact test function - locations at nodes
nnodes = sphere3D$nnodes
nloc = dim(loc)[1]
a1 = rnorm(1,mean = -2, sd = 1/8)
a2 = rnorm(1,mean = -2, sd = 1/8)
a3 = rnorm(1,mean = -2, sd = 1/8)

a1 = -1
a2 = -2
a3 = -3
func_evaluation = numeric(nloc)

for (i in 0:(nloc-1)){
  func_evaluation[i+1] = a1* sin(loc[i+1,1]) +  a2* sin(loc[i+1,2]) +  a3*sin(loc[i+1,3]) - 7
}
ran=range(func_evaluation) 
ran

# covariates
cov1_nonod=sin(2*pi*loc[,1])+sin((2*pi*loc[,2])^2)
cov2_nonod=cos(-2*pi*loc[,3])

plot(FEM(cov1_nonod[1:nnodes], FEMbasis))
plot(FEM(cov2_nonod[1:nnodes], FEMbasis))

W2=cbind(cov1_nonod,cov2_nonod)

# # Evaluation on a planar cut [?? purpose]
# {
#   grid_planar_ = cbind(runif(1000,min = -1 , max = 1),runif(1000,min = -1 , max = 1))
#   radius = 1
#   for(i in 1:length(grid_planar_[,1])){
#     
#     if(grid_planar_[i,1]^2+grid_planar_[i,2]^2 >= radius){
#       grid_planar_[i,] <- c(NA,NA)
#     }
#     
#   }    
#   
#   grid_planar <- matrix(ncol=2, nrow = sum(!is.na(grid_planar_[,1])))
#   grid_planar[,1] <- grid_planar_[!is.na(grid_planar_[,1]),1]
#   grid_planar[,2] <- grid_planar_[!is.na(grid_planar_[,2]),2]
#   
#   cov1_planar=sin(2*pi*grid_planar[,1])+sin((2*pi*grid_planar[,2])^2)
#   meshgrid2d <- create.mesh.2D(nodes=grid_planar)
#   plot(meshgrid2d)
#   FEM_basis_planar = fdaPDE::create.FEM.basis(meshgrid2d)
#   
#   image(FEM(cov1_planar, FEM_basis_planar))
# }

range(W2%*%beta_exact + func_evaluation) 

# data generation
theta = func_evaluation + W2%*%beta_exact
plot(FEM(theta[1:nnodes], FEMbasis))

mu = inv.link(theta)
plot(FEM(mu[1:nnodes], FEMbasis))

range(mu)

response <- rgamma(length(loc[,1]), shape=mu/scale.param, scale=scale.param)

# Plot the True Field
plot(FEM(response[1:nnodes],FEMbasis))

#### Test 1.2: grid with exact GCV
output_CPP <- smooth.FEM(locations = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = W2,
                         max.steps=15, family =FAMILY, mu0=NULL, scale.param=NULL,
                         lambda = lambdaGLM, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV', GCV.inflation.factor = 1)

best_lambdaGLM4<- output_CPP$bestlambda

RMSE<-function(f,g) sqrt(mean((f-g)^2))
rmseGLM =NULL
for (i in 1:20)
{
  output_CPP<- smooth.FEM(locations = loc, observations = as.numeric(rgamma(length(loc[,1]), shape=mu/scale.param, scale=scale.param)), FEMbasis =FEMbasis, covariates = W2,
                          max.steps=15, family =FAMILY, mu0=NULL, scale.param=NULL,
                          lambda = lambdaGLM, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV', GCV.inflation.factor = 1)
  best_lambdaGLM4 <- c(best_lambdaGLM4,output_CPP$bestlambda)
  rmseGLM<- c(rmseGLM, RMSE(a1* sin(loc_eval[,1]) +  a2* sin(loc_eval[,2]) +  a3*sin(loc_eval[,3]) - 7,eval.FEM(output_CPP$fit.FEM,loc_eval)))
}
boxplot(rmseGLM)

rmseGLMD <- rmseGLM

boxplot(rmseGLMA,rmseGLMB,rmseGLMC,rmseGLMD, names=c('fdaPDE','mass lumping','lambda','block'), col=c('grey',2,3,4),ylab='RMSE')

tGLMD=microbenchmark(smooth.FEM(locations = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = W2,
                                max.steps=15, family =FAMILY, mu0=NULL, scale.param=NULL,
                                lambda = best_lambdaGLM1[1]),times = 10)

boxplot(log(tGLMA$time),log(tGLMB$time),log(tGLMC$time),log(tGLMD$time), names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')
