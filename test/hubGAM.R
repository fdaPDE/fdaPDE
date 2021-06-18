#### hub pointwise (examples with and without covariates) ####
library(fdaPDE)
rm(list=ls())

data(hub2.5D)
mesh <- create.mesh.2.5D(nodes = hub2.5D$hub2.5D.nodes,triangles = hub2.5D$hub2.5D.triangles)
FEMbasis <- create.FEM.basis(mesh)

# Locations at nodes
nodes=mesh$nodes

# Exact data - Locations at nodes
nnodes = nrow(mesh$nodes)
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

func = function(x)
{
  a1*sin(2*pi*x[,2])+a2*sin(2*pi*x[,3])
}

coeff = func(nodes)
plot(FEM(coeff,FEMbasis))

lambdaS=10^(-12:3)
ran=range(coeff)

data=coeff
data1=data+rnorm(length(data)[1],sd=0.001*(ran[2]-ran[1]))

sol = smooth.FEM(observations=data,
                         FEMbasis = FEMbasis, lambda=lambdaS)
sol1 = smooth.FEM(observations=data1,
                      FEMbasis = FEMbasis, lambda=lambdaS)

sol_approx = sol$fit.FEM
femsol_approx<-t(eval.FEM(sol_approx,nodes))
sol_approx1 = sol1$fit.FEM
femsol_approx1<-t(eval.FEM(sol_approx1,nodes))

error <- matrix(nrow=length(data), ncol= length(lambdaS))
for (i in 1:length(lambdaS))
  error[,i] <- abs(coeff-t(femsol_approx)[,i])
err8 <- colMeans(error^2)
plot(log10(lambdaS), err, pch=4, xlab='log_10(lambda)', ylab='error', col=1)
grid()
abline(v=log10(best_lambda1), lty=2, col='red')
lines(-12:3, err1, col=1)
points(log10(best_lambda1), err1[14], pch=19, col='red')
points(-12:3, err2, pch=4, col=2)
lines(-12:3, err2, col=2)
points(log10(best_lambda1), err2[14], pch=3, col=2)
lines(-12:3, err3, col=3)
points(-12:3, err3, pch=3, col=3)
points(log10(best_lambda3), err3[11], pch=1, col=3)
legend('topleft',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

residual <- matrix(nrow=dim(nodes)[1], ncol= length(lambdaS))
residual <- abs(t(femsol_approx1)-t(femsol_approx))
res8 <- colMeans(residual)
plot(-12:3, res5, pch=4, xlab='log_10(lambda)', ylab='residual', col=1)
grid()
abline(v=log10(best_lambda8[1]), lty=2, col='red')
lines(-12:3, res1, col=1)
points(log10(best_lambda1), res1[14], pch=19, col='red')
points(-12:3, res6, pch=4, col=2)
lines(-12:3, res6, col=2)
points(log10(best_lambda2), res2[14], pch=19, col=2)
points(-12:3, res8, pch=1, col=4)
lines(-12:3, res4, col=3)
points(log10(best_lambda3), res3[14], pch=1, col=3)
legend('bottomleft',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

sol <- smooth.FEM(observations = data,FEMbasis = FEMbasis,
                    lambda.selection.criterion = 'newton', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
sol_approx <- sol$fit.FEM
best_lambda8 <- sol$optimization$lambda_solution
RMSE<-function(f,g) sqrt(mean((f-g)^2))
loc_eval = projection.points.2.5D(mesh, cbind(runif(10*dim(nodes)[1],range(nodes[,1])[1],range(nodes[,1][2])),
                                                    runif(10*dim(nodes)[1],range(nodes[,2])[1],range(nodes[,2][2])),
                                                    runif(10*dim(nodes)[1],range(nodes[,3])[1],range(nodes[,3][2]))))
rmse =NULL
for (i in 1:20)
{
  sol <- smooth.FEM(observations = data+rnorm(length(data)[1],sd=0.01*(ran[2]-ran[1])),FEMbasis = FEMbasis,
                    lambda.selection.criterion = 'newton', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
  sol_approx <- sol$fit.FEM
  best_lambda8 <- c(best_lambda8,sol$optimization$lambda_solution)
  rmse <- c(rmse, RMSE(func(loc_eval),eval.FEM(sol$fit.FEM,loc_eval)))
}

boxplot(rmse)
rmseH = rmse

boxplot(rmseA,rmseB,rmseC,rmseD, names=c('fdaPDE', 'mass lumping', 'lambda','block'), col=c('grey',2,3,4))

##### covariate case ####

cov1 = 4*sin(2*pi*nodes[,2])*cos(2*pi*nodes[,3])
cov2 = rnorm(nodes[,1],mean=3,sd=0.1)
plot(FEM(cov1,FEMbasis))
plot(FEM(cov2,FEMbasis))
beta=c(0.45,0.3)

dataCov <- coeff +beta[1]*cov1+beta[2]*cov2
data1Cov <- dataCov + rnorm (nnodes, sd=0.001*(ran[2]-ran[1]))

solCov <- smooth.FEM(observations = dataCov,
                     FEMbasis = FEMbasis, covariates=cbind(cov1,cov2), lambda=lambdaS)
sol_approxCov <- solCov$fit.FEM

solCov1 <- smooth.FEM(observations = data1Cov,
                     FEMbasis = FEMbasis, covariates=cbind(cov1,cov2), lambda=lambdaS)
sol_approxCov1 <- solCov1$fit.FEM

femsol_approxCov <- t(eval.FEM(sol_approxCov,nodes))
femsol_approxCov1 <- t(eval.FEM(sol_approxCov1,nodes))

errorCov <- matrix(nrow=dim(nodes)[1], ncol= length(lambdaS))
for (i in 1:length(lambdaS))
  errorCov[,i] <- abs(coeff-t(femsol_approxCov)[,i])
errCov8<- sqrt(colMeans(errorCov^2))
plot(-12:3, errCov1, pch=4, xlab='log_10(lambda)', ylab='error', col=1)
grid()
abline(v=log10(best_lambdaCov1), lty=2, col='red')
lines(-12:3, errCov1, col=1)
points(log10(best_lambdaCov1), errCov1[8], pch=19, col=1)
points(-12:3, errCov2, pch=4, col=2)
lines(-12:3, errCov2, col=2)
points(log10(best_lambdaCov1), errCov2[8], pch=19, col=2)
lines(-12:3, errCov4, col=4)
points(-12:3, errCov4, pch=4, col=4)
points(log10(best_lambdaCov3), errCov3[8], pch=1, col=3)
legend('topleft',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

residualCov <- matrix(nrow=dim(nodes)[1], ncol= length(lambdaS))
residualCov <- abs(t(femsol_approxCov1)-t(femsol_approxCov))
resCov8 <- colMeans(residualCov)
plot(-12:3, resCov5, pch=4, xlab='log_10(lambda)', ylab='residual', col=1, ylim=c(min(resCov4),range(resCov2)[2]))
grid()
abline(v=log10(best_lambdaCov1[1]), lty=2, col='red')
lines(-12:3, resCov5, col=1)
points(log10(best_lambdaCov1), resCov1[8], pch=19, col=1)
points(-12:3, resCov6, pch=4, col=2)
lines(-12:3, resCov6, col=2)
points(log10(best_lambdaCov2), resCov2[8], pch=19, col=2)
points(-12:3, resCov7, pch=4, col=3)
lines(-12:3, resCov4, col=4)
points(-12:3, resCov8, pch=1,col=4)
points(log10(best_lambdaCov3), resCov4[8], pch=19, col=4)
legend('topright',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner','block preconditioner'), cex=0.62,col=1:4, pch=19)

solCov <- smooth.FEM(observations = dataCov,
                     FEMbasis = FEMbasis, covariates=cbind(cov1,cov2),
                     lambda.selection.criterion = 'newton', DOF.evaluation = 'exact')
best_lambdaCov8 <- solCov$optimization$lambda_solution

RMSE<-function(f,g) sqrt(mean((f-g)^2))
rmseCov =NULL
for (i in 1:20)
{
  solCov <- smooth.FEM(observations = dataCov+ rnorm (nnodes, sd=0.01*(ran[2]-ran[1])),
                       FEMbasis = FEMbasis, covariates=cbind(cov1,cov2),
                       lambda.selection.criterion = 'newton', DOF.evaluation = 'exact')
  best_lambdaCov8 <- c(best_lambdaCov8,solCov$optimization$lambda_solution)
  rmseCov <- c(rmseCov, RMSE(func(loc_eval),eval.FEM(solCov$fit.FEM,loc_eval)))
}
boxplot(rmseCov)

rmseCovH <- rmseCov

boxplot(rmseCovE,rmseCovF,rmseCovG,rmseCovH, names=c('fdaPDE','mass lumping','lambda','block'), col=c('grey',2,3,4))

###### execution times #####
tG=microbenchmark(smooth.FEM(observations = data, FEMbasis = FEMbasis,lambda=best_lambda1[1]))
tcovG=microbenchmark(smooth.FEM(observations = dataCov, covariates=cbind(cov1,cov2), FEMbasis = FEMbasis, lambda=best_lambdaCov1[1]))

boxplot(log(tE$time),log(tF$time),log(tG$time),log(tH$time), names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')
boxplot(log(tcovE$time),log(tcovF$time),log(tcovG$time),log(tcovH$time), names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')

boxplot(tA$time,tB$time,tC$time, names=c('fdaPDE','mass lumping','lambda'),col=c('grey',2,3), ylab='time')
boxplot(tcovA$time,tcovB$time,tcovC$time, names=c('fdaPDE','mass lumping','lambda'),col=c('grey',2,3), ylab='time')


###############################################
############ TEST SCRIPT GAM 2.5D #############
###############################################

####### 2.5D ########
#### Test 1: HUB (2.5D) POISSON family ####
#            locations = nodes
#            with covariates
#            no BC
library(fdaPDE)
library(purrr)

# family
FAMILY = "poisson"

l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

# beta
beta1 = 2
beta2 = -3
betas_truth = c(beta1,beta2)

# lambda 
lambdaGLM = c(0.0001,0.001,0.005,0.006,0.007,0.008,0.009,0.0095,
           0.00975,0.01,0.0125,0.015,0.02,0.03,0.04,0.05,0.1,
           1,10,100,1000)

# scale param
scale.param = 1

# mesh
data(hub2.5D)
mesh <- create.mesh.2.5D(nodes = hub2.5D$hub2.5D.nodes,triangles = hub2.5D$hub2.5D.triangles)
FEMbasis <- create.FEM.basis(mesh)

# locations 
loc = mesh$nodes
nloc = dim(loc)[1]

# 2.5D random field (function f)
a1 = rnorm(1,mean = 0, sd = 1)
a2 = rnorm(1,mean = 0, sd = 1)
a3 = rnorm(1,mean = 0, sd = 1)

sol_exactGLM = numeric(nloc)
for (i in 0:(nloc-1)){
  sol_exactGLM[i+1] = a1* sin(2*pi*loc[i+1,1]) +  a2* sin(2*pi*loc[i+1,2]) +  a3*sin(2*pi*loc[i+1,3]) + 7
}
range(sol_exactGLM) 

# covariates
set.seed(42)

desmat=matrix(0,nrow=nloc, ncol=2)
desmat[,1]= sin(2*pi*loc[,1])*cos(2*pi*loc[,2])
desmat[,2]= rnorm(nloc, mean=2, sd=0.1)

# samopling response
ranGLM=range(desmat%*%betas_truth + sol_exactGLM) 
param=sol_exactGLM+beta1*desmat[,1]+beta2*desmat[,2]

mu<-inv.link(param)
range(mu)
response <- rpois(nloc, lambda = mu)

# Plot the True Field
plot(FEM(sol_exactGLM, FEMbasis))

# #### Test 1.1: Without GCV
# output_CPP<- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
#                         max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
#                         lambda = lambda)

#### Test 1.2: grid with exact GCV
output_CPP<- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                        max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                        lambda = lambdaGLM, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')

best_lambdaGLM4<- output_CPP$bestlambda

RMSE<-function(f,g) sqrt(mean((f-g)^2))
rmseGLM =NULL
for (i in 1:20)
{
  output_CPP<- smooth.FEM(location = NULL, observations = as.numeric(rpois(nloc, lambda = mu)), FEMbasis =FEMbasis, covariates = desmat,
                          max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                          lambda = lambdaGLM, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
  best_lambdaGLM4 <- c(best_lambdaGLM4,output_CPP$bestlambda)
  rmseGLM<- c(rmseGLM, RMSE(a1* sin(2*pi*loc_eval[i+1,1]) +  a2* sin(2*pi*loc_eval[i+1,2]) +  a3*sin(2*pi*loc_eval[i+1,3]) + 7,eval.FEM(output_CPP$fit.FEM,loc_eval)))
}
boxplot(rmseGLM)

rmseGLMD <- rmseGLM

boxplot(rmseGLMA,rmseGLMB,rmseGLMC,rmseGLMD, names=c('fdaPDE','mass lumping','lambda','block'), col=c('grey',2,3,4),ylab='rmse')

tGLMD=microbenchmark(smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                                max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                                lambda = best_lambdaGLM1[1]))

boxplot(log(tGLME$time),log(tGLMF$time),log(tGLMG$time),log(tGLMD$time), names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')
