data(horseshoe2D)
boundary_nodes = horseshoe2D$boundary_nodes
boundary_segments = horseshoe2D$boundary_segments
locations = horseshoe2D$locations

mesh <- create.mesh.2D(nodes = rbind (boundary_nodes, locations), segments = boundary_segments)
#mesh <- create.mesh.2D(nodes=boundary_nodes, segments=boundary_segments)
nodes <- mesh$nodes
locations = refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
fembasis <- create.FEM.basis(mesh)

coeff <- fs.test (nodes[,1],nodes[,2])

femfun = FEM(coeff, fembasis)
femsol <- eval.FEM(femfun, locations)
plot(femfun)
image(femfun)

#### no covariates case ####

lambda = 10^(-12:3)
sol_exact <- fs.test(locations[,1],locations[,2])
#data <- coeff + rnorm (nrow(nodes), sd=0.01)
#data <- coeff
data <- sol_exact
ran = range(sol_exact)
data1 <- data + rnorm (length(data), sd=0.01*(ran[2]-ran[1]))

sol <- smooth.FEM(observations = data, FEMbasis = fembasis, locations=locations,lambda=lambda)#,
                  #lambda.selection.criterion = 'newton', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
sol_approx <- sol$fit.FEM
best_lambda8 <- sol$optimization$lambda_solution
#plot(sol_approx)

sol1 <- smooth.FEM(observations = data1, FEMbasis = fembasis, location=locations,
                   lambda = lambda,lambda.selection.lossfunction = NULL)

sol_approx1 <- sol1$fit.FEM
#plot(sol_approx1)

femsol_approx <- t(eval.FEM(sol_approx,locations))
femsol_approx1 <- t(eval.FEM(sol_approx1,locations))

error <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
for (i in 1:length(lambda))
  error[,i] <- abs(sol_exact-t(femsol_approx)[,i])
err8 <- sqrt(colMeans(error^2))
plot(-12:3, err1, pch=4, xlab='log_10(lambda)', ylab='error', col=1)
grid()
abline(v=log10(best_lambda1), lty=2, col='red')
lines(-12:3, err1, col=1)
points(log10(best_lambda1), err1[10], pch=19, col='red')
points(-12:3, err2, pch=4, col=2)
lines(-12:3, err2, col=2)
points(log10(best_lambda2), err2[10], pch=19, col='red')
lines(-12:3, err3, col=3)
points(-12:3, err6, pch=4, col=6)
points(log10(best_lambda3), err3[10], pch=19, col='red')
legend('topleft',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

residual <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
residual <- abs(t(femsol_approx1)-t(femsol_approx))
res8 <- colMeans(residual)
plot(-12:3, res5, pch=4, xlab='log_10(lambda)', ylab='residual', col=1, ylim=c(0,0.1))
grid()
abline(v=log10(best_lambda1[1]), lty=2, col='red')
lines(-12:3, res5, col=1)
points(log10(best_lambda1[1]), res1[10], pch=19, col=1)
points(-12:3, res6, pch=4, col=2)
lines(-12:3, res6, col=2)
points(log10(best_lambda2[1]), res2[10], pch=19, col=2)
lines(-12:3, res7, col=3)
lines(-12:3, res8, pch=4, col=4)
points(log10(best_lambda3), res3[10], pch=1, col=3)
legend('topright',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

sol <- smooth.FEM(observations = data, FEMbasis = fembasis, locations=locations,
                lambda.selection.criterion = 'newton', DOF.evaluation = 'exact',
                lambda.selection.lossfunction = 'GCV')
sol_approx <- sol$fit.FEM
best_lambda5 <- sol$optimization$lambda_solution
RMSE<-function(f,g) sqrt(mean((f-g)^2))
loc_eval = refine.mesh.2D(mesh, minimum_angle = 30, maximum_area = 0.005)$nodes
rmse1 =NULL
for (i in 1:20)
{
  sol <- smooth.FEM(observations = data+rnorm(length(data), sd=0.01*(ran[2]-ran[1])),
                    FEMbasis = fembasis, locations=locations,
                    lambda.selection.criterion = 'newton', DOF.evaluation = 'exact',
                    lambda.selection.lossfunction = 'GCV')
  sol_approx <- sol$fit.FEM
  best_lambda5 <- c(best_lambda5,sol$optimization$lambda_solution)
  rmse1 <- c(rmse1, RMSE(fs.test(loc_eval[,1],loc_eval[,2]),eval.FEM(sol$fit.FEM,loc_eval)))
}
boxplot(rmse1)

rmse1E <- rmse1

boxplot(rmse1E,rmse1F,rmse1G,rmse1EH, names=c('fdaPDE', 'mass lumping', 'lambda','block'), col=c('grey', 2,3,4),ylab='RMSE')

locations_eval = refine.mesh.2D(mesh, minimum_angle = 30, maximum_area = 0.005)$nodes
sol_eval=eval.FEM(sol_approx,locations_eval)
sol_exact=fs.test(locations_eval[,1],locations_eval[,2])
rmse1=NULL
for (i in 1:length(lambda))
  rmse1 <- c(rmse1,RMSE(sol_exact, sol_eval[,i]))
plot(log10(lambda), rmse1, pch=4, xlab='log_10(lambda)', ylab='rmse', col=1)
grid()
abline(v=log10(best_lambda1), lty=2, col='red')
lines(log10(lambda), rmse1, col=1)
points(log10(best_lambda1), rmse1[10], pch=19, col='red')
points(log10(lambda), rmse2, pch=4, col=2)
lines(log10(lambda), rmse2, col=2)
points(log10(best_lambda1), rmse3[10], pch=3, col=2)
lines(log10(lambda), rmse3, col=4)
points(log10(lambda), rmse3, pch=4, col=4)
points(log10(best_lambda3), rmse3[10], pch=1, col=3)


##### covariates case ####

#load covariates as covs.test

# covariate = covs.test(nodes[,1], nodes[,2])
# coeffCov = fs.test(nodes[,1], nodes[,2]) + 2*covariate
ndata=length(coeff)
cov1 = rnorm(ndata, mean = 1, sd = 2)
cov2 = sin(nodes[,1])
coeffCov = fs.test(nodes[,1], nodes[,2])

cov1_exactCov = rnorm(ndata, mean = 1, sd = 2)
cov2_exactCov = sin(locations[,1])
sol_exactCov = fs.test(locations[,1], locations[,2])

femfunCov = FEM(coeffCov, fembasis)
plot(femfunCov)
image(femfunCov)

lambda = 10^(-12:3)

#sol_exactCov <- fs.test(locations[,1],locations[,2]) + 2*covs.test(locations[,1],locations[,2])
#data <- coeff + rnorm (nrow(nodes), sd=0.01)
dataCov <- sol_exactCov +2*cov1_exactCov -cov2_exactCov
ranCov = range(sol_exactCov)
data1Cov <- dataCov + rnorm (ndata, sd=0.001*(ranCov[2]-ranCov[1]))

solCov <- smooth.FEM(locations=locations,observations = dataCov,
                     FEMbasis = fembasis, covariates=cbind(cov1_exactCov,cov2_exactCov),lambda=lambda)
                     #lambda.selection.criterion = 'newton', DOF.evaluation = 'exact')
sol_approxCov <- solCov$fit.FEM
best_lambdaCov8 <- solCov$optimization$lambda_solution
#image(solCov$fit.FEM)


sol1Cov <- smooth.FEM(locations=locations, observations = data1Cov,
                      FEMbasis = fembasis, covariates = cbind(cov1_exactCov,cov2_exactCov), lambda = lambda)
sol_approx1Cov <- sol1Cov$fit.FEM
#plot(sol_approx1Cov)

femsol_approxCov <- t(eval.FEM(sol_approxCov,locations))
femsol_approx1Cov <- t(eval.FEM(sol_approx1Cov,locations))

errorCov <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
for (i in 1:length(lambda))
  errorCov[,i] <- abs(sol_exactCov-t(femsol_approxCov)[,i])
errCov8<- sqrt(colMeans(errorCov^2))
plot(-12:3, errCov1, pch=4, xlab='log_10(lambda)', ylab='error', col=1, ylim=c(min(errCov4),1.2))
grid()
abline(v=log10(best_lambdaCov1), lty=2, col='red')
lines(-12:3, errCov1, col=1)
points(log10(best_lambdaCov1), errCov1[10], pch=19, col=1)
points(-12:3, errCov2, pch=4, col=2)
lines(-12:3, errCov2, col=2)
points(log10(best_lambdaCov1), errCov2[10], pch=19, col=2)
lines(-12:3, errCov4, col=4)
points(-12:3, errCov3, pch=4, col=3)
points(log10(best_lambdaCov3), errCov3[10], pch=1, col=3)
points(-12:3, errCov4, pch=4, col=4)
#lines(-12:3, errCov3, col=3)
points(log10(best_lambdaCov4), errCov4[10], pch=19, col=4)
legend('topleft',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner','block preconditioner'), cex=0.62,col=1:4, pch=19)

residualCov <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
residualCov <- abs(t(femsol_approx1Cov)-t(femsol_approxCov))
resCov8 <- colMeans(residualCov)
plot(-12:3, resCov7, pch=4, xlab='log_10(lambda)', ylab='residual', col=1)
grid()
abline(v=log10(best_lambdaCov1[1]), lty=2, col='red')
lines(-12:3, resCov7, col=1)
points(log10(best_lambdaCov1[1]), resCov1[10], pch=19, col=1)
points(-12:3, resCov2, pch=4, col=2)
lines(-12:3, resCov2, col=2)
points(log10(best_lambdaCov2), resCov2[10], pch=19, col=2)
points(-12:3, resCov7, pch=4, col=3)
#lines(-12:3, resCov3, col=3)
points(log10(best_lambdaCov3), resCov3[10], pch=19, col=3)
points(-12:3, resCov8, pch=1, col=4)
lines(-12:3, resCov8, col=4)
points(log10(best_lambdaCov4), resCov4[10], pch=19, col=4)
legend('topright',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner','block preconditioner'), cex=0.62,col=1:4, pch=19)

solCov <- smooth.FEM(locations=locations,observations = dataCov,
                     FEMbasis = fembasis, covariates=cbind(cov1_exactCov,cov2_exactCov),
                    lambda.selection.criterion = 'newton', DOF.evaluation = 'exact')
sol_approxCov <- solCov$fit.FEM
best_lambdaCov5 <- solCov$optimization$lambda_solution
RMSE<-function(f,g) sqrt(mean((f-g)^2))
rmseCov =NULL
for (i in 1:20)
{
  solCov <- smooth.FEM(locations=locations,observations = dataCov+rnorm(ndata, sd=0.01*(ranCov[2]-ranCov[1])),
                       FEMbasis = fembasis, covariates=cbind(cov1_exactCov,cov2_exactCov),
                      lambda.selection.criterion = 'newton', DOF.evaluation = 'exact')
  sol_approxCov <- solCov$fit.FEM
  best_lambdaCov5 <- c(best_lambdaCov5,solCov$optimization$lambda_solution)
  rmseCov <- c(rmseCov, RMSE(fs.test(loc_eval[,1],loc_eval[,2]), eval.FEM(solCov$fit.FEM,loc_eval)))
}
boxplot(rmseCov)

rmseCovE <- rmseCov

boxplot(rmseCovA,rmseCovB,rmseCovC,rmseCovD, names=c('fdaPDE','mass lumping','lambda','block'), col=c('grey',2,3,4),ylab='RMSE')

locations_eval = refine.mesh.2D(mesh, minimum_angle = 30, maximum_area = 0.005)$nodes
sol_evalCov=eval.FEM(sol_approx,locations_eval)
sol_exactCov=fs.test(locations_eval[,1],locations_eval[,2]) + 2*covs.test(locations_eval[,1],locations_eval[,2])
rmseCov1=NULL
for (i in 1:length(lambda))
  rmseCov1 <- c(rmseCov1,RMSE(sol_exact, sol_eval[,i]))
plot(log10(lambda), rmseCov1, pch=4, xlab='log_10(lambda)', ylab='residual', col=1)
grid()
abline(v=log10(best_lambdaCov1), lty=2, col='red')
lines(log10(lambda), rmseCov1, col=1)
points(log10(best_lambdaCov1), rmseCov1[10], pch=19, col='red')
points(log10(lambda), rmseCov2, pch=4, col=2)
lines(log10(lambda), rmseCov2, col=2)
points(log10(best_lambdaCov1), rmseCov3[10], pch=3, col=2)
lines(log10(lambda), rmseCov3, col=4)
points(log10(lambda), rmseCov3, pch=4, col=4)
points(log10(best_lambdaCov3), rmseCov3[10], pch=1, col=3)

#### 3D case ####
data(horseshoe2.5D)
mesh3d = horseshoe2.5D

nodes3d <- mesh3d$nodes
fembasis3d <- create.FEM.basis(mesh3d)
coeff3d <- fs.test.3D (nodes3d[,1],nodes3d[,2],nodes3d[,3])


femfun3d = FEM(coeff3d, fembasis3d)
plot(femfun3d)

#### no covariates case ####

lambda = 10^(-12:3)
sol_exact3d <- coeff3d
#data <- coeff + rnorm (nrow(nodes), sd=0.01)
data3d <- coeff3d
ran3d = range(sol_exact3d)
data13d <- data3d + rnorm (nrow(nodes3d), sd=0.001*(ran3d[2]-ran3d[1]))

sol3d <- smooth.FEM(observations = data3d,FEMbasis = fembasis3d, lambda=lambda)
                    #lambda.selection.criterion = 'newton', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
sol_approx3d <- sol3d$fit.FEM
#plot(sol_approx3d)

sol13d <- smooth.FEM(observations = data13d, FEMbasis = fembasis3d, lambda = lambda,lambda.selection.lossfunction = NULL)
sol_approx13d <- sol13d$fit.FEM
#plot(sol_approx13d)

femsol_approx3d <- t(eval.FEM(sol_approx3d,nodes3d))
femsol_approx13d <- t(eval.FEM(sol_approx13d,nodes3d))

error3d <- matrix(nrow=length(sol_exact3d), ncol= length(lambda))
for (i in 1:length(lambda))
  error3d[,i] <- abs(sol_exact3d-t(femsol_approx3d)[,i])
err3d4 <- colMeans(error3d^2)
plot(-12:3, err3d1, pch=4, xlab='log_10(lambda)', ylab='error', col=1)
grid()
abline(v=log10(best_lambda3d), lty=2, col='red')
lines(-12:3, err3d1, col=1)
points(log10(best_lambda3d), err3d1[10], pch=19, col='red')
points(-12:3, err3d2, pch=4, col=2)
lines(-12:3, err3d2, col=2)
points(log10(best_lambda3d), err3d2[10], pch=3, col=2)
lines(-12:3, err3d3, col=3)
points(-12:3, err3d4, pch=4, col=4)
points(log10(best_lambda3d), err3d3[10], pch=1, col=3)

residual3d <- matrix(nrow=dim(nodes3d)[1], ncol= length(lambda))
residual3d <- abs(t(femsol_approx13d)-t(femsol_approx3d))
res3d4 <- colMeans(residual3d)
plot(-12:3, res3d1, pch=4, xlab='log_10(lambda)', ylab='residual', col=1,  ylim=c(0,0.018))
grid()
abline(v=log10(best_lambda3d2), lty=2, col='red')
lines(-12:3, res3d1, col=1)
points(log10(best_lambda3d2), res3d1[10]/2+res3d1[9]/2, pch=19, col=1)
points(-12:3, res3d2, pch=4, col=2)
lines(-12:3, res3d2, col=2)
points(log10(best_lambda3d2), res3d2[10]/2+res3d2[9]/2, pch=19, col=2)
points(-12:3, res3d3, pch=4, col=3)
lines(-12:3, res3d3, col=3)
points(log10(best_lambda3d3), res3d3[10]/2+res3d3[9]/2, pch=1, col=3)
legend('topright',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

sol3d <- smooth.FEM(observations = data3d+ rnorm (nrow(nodes3d), sd=0.001*(ran3d[2]-ran3d[1])),FEMbasis = fembasis3d, #lambda=lambda)
                    lambda.selection.criterion = 'newton', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
sol_approx3d <- sol3d$fit.FEM
best_lambda3d1 <- sol3d$optimization$lambda_solution
RMSE<-function(f,g) sqrt(mean((f-g)^2))
loc_eval3D = projection.points.2.5D(mesh3d, cbind(runif(dim(nodes3d)[1],-2.25,6),runif(dim(nodes3d)[1],-1,1),runif(dim(nodes3d)[1],-2.25,2.25)))#cbind(locations_eval,runif(dim(locations_eval)[1]),-2,2))#rep(0,dim(locations_eval)[1])))
rmse3D =NULL
for (i in 1:20)
{
  sol3d <- smooth.FEM(observations = data3d,FEMbasis = fembasis3d,
                      lambda.selection.criterion = 'newton', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
  sol_approx3d <- sol3d$fit.FEM
  best_lambda3d1 <- c(best_lambda3d1,sol3d$optimization$lambda_solution)
  rmse3D <- c(rmse3D, RMSE(fs.test.3D(loc_eval3d[,1],loc_eval3d[,2],loc_eval3d[,3]),eval.FEM(sol3d$fit.FEM,loc_eval3d)))
}

boxplot(rmse3D)
rmse3DA = rmse3D

boxplot(rmse3DA,rmse3DB,rmse3DC,rmse3DD, names=c('fdaPDE','mass lumping','lambda','block'), col=c('grey',2,3,4))

locations_eval3D = refine.mesh.2D(mesh, minimum_angle = 30, maximum_area = 0.005)$nodes
sol_eval3D=eval.FEM(sol_approx3D,locations_eval)
sol_exact3D=fs.test(locations_eval3D[,1],locations_eval3D[,2])
rmse3D1=NULL
for (i in 1:length(lambda))
  rmse3D1 <- c(rmse3D1,RMSE(sol_exact3D, sol_eval3D[,i]))
plot(log10(lambda), rmse3D1, pch=4, xlab='log_10(lambda)', ylab='residual', col=1)
grid()
abline(v=log10(best_lambda3D1), lty=2, col='red')
lines(log10(lambda), rmse3D1, col=1)
points(log10(best_lambda3D1), rmse3D1[10], pch=19, col='red')
points(log10(lambda), rmse3D2, pch=4, col=2)
lines(log10(lambda), rmse3D2, col=2)
points(log10(best_lambda3D1), rmse3D3[10], pch=3, col=2)
lines(log10(lambda), rmse3D3, col=4)
points(log10(lambda), rmse3D3, pch=4, col=4)
points(log10(best_lambda3D3), rmse3D3[10], pch=1, col=3)

###### execution times #####
tG=microbenchmark(smooth.FEM(observations = data, FEMbasis = fembasis, locations=locations,lambda=best_lambda1[1]))
tcovG=microbenchmark(smooth.FEM(observations = dataCov, FEMbasis = fembasis, covariates = cbind(cov1_exactCov,cov2_exactCov), locations=locations,lambda=best_lambdaCov1[1]),times=30)
t3dH=microbenchmark(smooth.FEM(observations = data3d,FEMbasis = fembasis3d, lambda=best_lambda3d),times=10)

boxplot(log(tE$time),log(tF$time),log(tG$time),log(tH$time), names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')
boxplot(log(tcovE$time),log(tcovF$time),log(tcovG$time),log(tcovH$time),names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')
boxplot(log(t3dA$time),log(t3dB$time),log(t3dC$time),log(t3dD$time), names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')

boxplot(tA$time,tB$time,tC$time,tD$time, names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='time')
boxplot(tcovA$time,tcovB$time,tcovC$time,tcovD$time, names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='time')
boxplot(t3dA$time,t3dB$time,t3dC$time, names=c('fdaPDE','mass lumping','lambda'),col=c('grey',2,3), ylab='time')