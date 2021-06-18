x = seq(0,1, length.out = 30)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)
plot(mesh)
nodes=mesh$nodes
nnodes = dim(mesh$nodes)[1]

FEMbasis = create.FEM.basis(mesh)

# Test function
f = function(x, y, z = 1)
{
  coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
  sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
}

# Exact solution (pointwise at nodes)
sol_exact = f(mesh$nodes[,1], mesh$nodes[,2])
image(FEM(sol_exact, FEMbasis))

# Add error to simulate data
set.seed(7893475)
ran = range(sol_exact)
data = sol_exact
data1 = sol_exact + rnorm(nnodes, mean=0, sd=0.01*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^(-12:3)

#### Test 1.1: Without GCV
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda)
#image(output_CPP$fit.FEM)
sol_approx = output_CPP$fit.FEM
femsol_approx<-t(eval.FEM(sol_approx,nodes))

output_CPP1<-smooth.FEM(observations=data1, FEMbasis=FEMbasis, lambda=lambda)
sol_approx1 = output_CPP1$fit.FEM
femsol_approx1<-t(eval.FEM(sol_approx1,nodes))

# #### Test 1.2: grid with exact GCV
# output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda, lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
# plot(log10(lambda), output_CPP$optimization$GCV_vector)
# image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

# #### Test 1.3: grid with stochastic GCV
# output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda, lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV', DOF.stochastic.realizations = 1000)
# plot(log10(lambda), output_CPP$optimization$GCV_vector)
# image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

### Test 1.4: Newton exact method with exact GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
best_lambda5<-output_CPP$optimization$lambda_solution
image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

# ### Test 1.5: Newton_fd method with  exact GCV, default initial lambda and tolerance
# output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda.selection.criterion='newton_fd', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
# 
# ### Test 1.6: Newton_fd method with stochastic GCV, default initial lambda and tolerance
# output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')

error <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
for (i in 1:length(lambda))
  error[,i] <- abs(sol_exact-t(femsol_approx)[,i])
err5 <- colMeans(error^2)
plot(log10(lambda), err1, pch=4, xlab='log_10(lambda)', ylab='error', col=1, ylim=c(0,0.45))
grid()
abline(v=log10(best_lambda1), lty=2, col='red')
lines(-12:3, err1, col=1)
points(log10(best_lambda1), err1[12]/2+err1[13]/2, pch=19, col='red')
points(-12:3, err2, pch=4, col=2)
lines(-12:3, err2, col=2)
points(log10(best_lambda1), err2[12]/2+err2[13]/2, pch=3, col='red')
lines(-12:3, err3, col=3)
points(-12:3, err3, pch=4, col=3)
points(log10(best_lambda3), err3[12]/2+err3[13]/2, pch=19, col='red')
legend('topright',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner'), cex=0.62,col=1:3, pch=19)

residual <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
residual <- abs(t(femsol_approx1)-t(femsol_approx))
res5 <- colMeans(residual)
plot(log10(lambda), res5, pch=4, xlab='log_10(lambda)', ylab='residual', col=1)#, ylim=c(0,0.026))
grid()
abline(v=log10(best_lambda5[1]), lty=2, col='red')
lines(-12:3, res5, col=1)
points(log10(best_lambda5[1]), res5[12], pch=19, col='red')
points(-12:3, res6, pch=4, col=2)
lines(-12:3, res6, col=2)
points(log10(best_lambda1), res2[12], pch=19, col='red')
lines(-12:3, res8, col=4)
points(-12:3, res7, pch=4, col=3)
points(log10(best_lambda3), res3[12], pch=19, col='red')
legend('topright',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner'), cex=0.62,col=1:3, pch=19)

par(mfrow=c(1,2))

plot(log10(lambda), err1, pch=4, xlab='log_10(lambda)', ylab='error', col=1)
grid()
abline(v=log10(best_lambda1), lty=2, col='red')
lines(-12:3, err1, col=1)
points(log10(best_lambda1), err1[12]/2+err1[13]/2, pch=19, col='red')

plot(log10(lambda), res1, pch=4, xlab='log_10(lambda)', ylab='residual', col=1)
grid()
abline(v=log10(best_lambda1), lty=2, col='red')
lines(-12:3, res1, col=1)
points(log10(best_lambda1), res1[12], pch=19, col='red')


RMSE<-function(f,g) sqrt(mean((f-g)^2))
xeval=runif(10000,0,1)
yeval=runif(10000,0,1)
rmse <- NULL
for (i in 1:20)
{
  output_CPP<-smooth.FEM(observations=data+rnorm(nnodes, mean=0, sd=0.01*abs(ran[2]-ran[1])), FEMbasis=FEMbasis, lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
  best_lambda5<-c(best_lambda5,output_CPP$optimization$lambda_solution)
  sol_approx=output_CPP$fit.FEM
  rmse <- c(rmse, RMSE(f(xeval,yeval),eval.FEM(sol_approx,locations=cbind(xeval,yeval))))
}

boxplot(rmse)
rmseA <- rmse

boxplot(rmseA,rmseB,rmseC,rmseD, names=c('fdaPDE','mass lumping','lambda','block'), col=c('grey',2,3,4))

xeval=runif(10000,0,1)
yeval=runif(10000,0,1)
sol_eval=eval.FEM(sol_approx,locations = cbind(xeval,yeval))
sol_exact=f(xeval,yeval)
rmse = NULL
for (i in 1:length(lambda))
  rmse <- c(rmse, RMSE(sol_exact,sol_eval[,i]))
rmse3 <- rmse
plot(log10(lambda), rmse1, pch=4, xlab='log_10(lambda)', ylab='rmse', col=1, ylim=c(0,0.8))
grid()
abline(v=log10(best_lambda1), lty=2, col='red')
lines(-12:3, rmse1, col=1)
points(log10(best_lambda1), rmse1[12]/2+rmse1[13]/2, pch=19, col='red')
points(-12:3, rmse2, pch=4, col=2)
lines(-12:3, rmse2, col=2)
points(log10(best_lambda2), rmse2[12]/2+rmse2[13]/2, pch=3, col='red')
lines(-12:3, rmse3, col=3)
points(-12:3, rmse3, pch=4, col=3)
points(log10(best_lambda3), rmse3[12]/2+rmse3[13]/2, pch=19, col='red')
legend('topright',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner'), cex=0.62,col=1:3, pch=19)


# Order=2
mesh2 = create.mesh.2D(locations)
plot(mesh2)
nodes2=mesh2$nodes
nnodes2 = dim(mesh2$nodes)[1]

FEMbasis2 = create.FEM.basis(mesh2)

# Test function
a1=1
a2=4
z<-function(p){  
  a1*sin(2*pi*p[,1])*cos(2*pi*p[,2])+a2*sin(3*pi*p[,1])}

# Exact solution (pointwise at nodes)
sol_exact2 = z(mesh2$nodes)
image(FEM(sol_exact2, FEMbasis2))

# Add error to simulate data
set.seed(7893475)
ran2 = range(sol_exact2)
data2 = sol_exact2
data12 = sol_exact2 + rnorm(nnodes2, mean=0, sd=0.01*abs(ran2[2]-ran2[1]))

# Set smoothing parameter
lambda = 10^(-12:3)

#### Test 1.1: Without GCV
output_CPP2<-smooth.FEM(observations=data2, FEMbasis=FEMbasis2, lambda=lambda)
#image(output_CPP$fit.FEM)
sol_approx2 = output_CPP2$fit.FEM
femsol_approx2<-t(eval.FEM(sol_approx2,nodes2))

output_CPP12<-smooth.FEM(observations=data12, FEMbasis=FEMbasis2, lambda=lambda)
sol_approx12 = output_CPP12$fit.FEM
femsol_approx12<-t(eval.FEM(sol_approx12,nodes2))

# #### Test 1.2: grid with exact GCV
# output_CPP2<-smooth.FEM(observations=data2, FEMbasis=FEMbasis, lambda=lambda, lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
# plot(log10(lambda), output_CPP2$optimization$GCV_vector)
# image(FEM(output_CPP2$fit.FEM$coeff,FEMbasis))

# #### Test 1.3: grid with stochastic GCV
# output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda, lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV', DOF.stochastic.realizations = 1000)
# plot(log10(lambda), output_CPP$optimization$GCV_vector)
# image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

### Test 1.4: Newton exact method with exact GCV, default initial lambda and tolerance
output_CPP2<-smooth.FEM(observations=data2, FEMbasis=FEMbasis2, lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
best_lambda21<-output_CPP2$optimization$lambda_solution
image(FEM(output_CPP2$fit.FEM$coeff,FEMbasis))

# ### Test 1.5: Newton_fd method with  exact GCV, default initial lambda and tolerance
# output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda.selection.criterion='newton_fd', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
# 
# ### Test 1.6: Newton_fd method with stochastic GCV, default initial lambda and tolerance
# output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')


error2 <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
for (i in 1:length(lambda))
  error2[,i] <- abs(sol_exact2-t(femsol_approx2)[,i])
err24 <- colMeans(error2^2)
plot(log10(lambda), err21, pch=4, xlab='log_10(lambda)', ylab='error', col=1, ylim=c(0,12))
grid()
abline(v=log10(best_lambda21), lty=2, col='red')
lines(-12:3, err21, col=1)
points(log10(best_lambda21), err21[11]/2, pch=19, col='red')
points(-12:3, err22, pch=4, col=2)
lines(-12:3, err22, col=2)
points(log10(best_lambda22), err22[11]/2, pch=3, col='red')
lines(-12:3, err23, col=3)
points(-12:3, err23, pch=4, col=3)
points(log10(best_lambda23), err23[11]/2, pch=19, col='red')
legend('topright',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner'), cex=0.62,col=1:3, pch=19)

residual2 <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
residual2 <- abs(t(femsol_approx12)-t(femsol_approx2))
res24 <- colMeans(residual2)
plot(log10(lambda), res5, pch=4, xlab='log_10(lambda)', ylab='residual', col=1)
grid()
abline(v=log10(best_lambda21), lty=2, col='red')
lines(-12:3, res5, col=1)
points(log10(best_lambda21), res21[10]/2+res21[11]/2, pch=19, col='red')
points(-12:3, res6, pch=4, col=2)
lines(-12:3, res6, col=2)
points(log10(best_lambda21), res22[10]/2+res22[11]/2, pch=19, col='red')
lines(-12:3, res7, col=3)
points(-12:3, res7, pch=4, col=3)
points(log10(best_lambda23), res23[10]/2+res23[11]/2, pch=19, col='red')
legend('topright',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner'), cex=0.62,col=1:3, pch=19)

xeval=runif(10000,0,1)
yeval=runif(10000,0,1)
sol_eval2=eval.FEM(sol_approx2,locations = cbind(xeval,yeval))
sol_exact2=z(xeval,yeval)
RMSE<-function(f,g) sqrt(mean((f-g)^2))
rmse2 <- NULL
for (i in 1:20)
{
  xeval=runif(10000,0,1)
  yeval=runif(10000,0,1)
  rmse2 <- c(rmse2, RMSE(z(cbind(xeval,yeval)),eval.FEM(sol_approx2,locations=cbind(xeval,yeval))))
}
boxplot(rmse2)
rmse2D <- rmse

rmse <- NULL
for (i in 1:length(lambda))
  rmse <- c(rmse, RMSE(sol_exact,sol_eval[,i]))
rmse23 <- rmse
plot(log10(lambda), rmse21, pch=4, xlab='log_10(lambda)', ylab='rmse', col=1, ylim=c(0,0.8))
grid()
abline(v=log10(best_lambda21), lty=2, col='red')
lines(-12:3, rmse21, col=1)
points(log10(best_lambda21), rmse21[10]/2+rmse21[11]/2, pch=19, col='red')
points(-12:3, rmse22, pch=4, col=2)
lines(-12:3, rmse22, col=2)
points(log10(best_lambda21), rmse22[10]/2+rmse22[11]/2, pch=3, col='red')
lines(-12:3, rmse23, col=3)
points(-12:3, rmse23, pch=4, col=3)
points(log10(best_lambda23), rmse23[10]/2+rmse23[11]/2, pch=19, col='red')
legend('topright',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner'), cex=0.62,col=1:3, pch=19)

##### covariate case #####

a1=runif(1,min=-1.5,max=1.5)
a2=runif(1,min=-1.5,max=1.5)
zCov<-function(p){  
  a1*sin(2*pi*p[,1])*cos(2*pi*p[,2])+a2*sin(3*pi*p[,1])}

coeffCov = zCov(mesh$nodes)
cov=cos(3*pi*nodes[,2])
beta=1.2

image(FEM(coeffCov,FEMbasis))
image(FEM(cov,FEMbasis))

dataCov = coeffCov + beta*cov
ranCov = range(coeffCov)
data1Cov <- dataCov + rnorm (nnodes, sd=0.001*(ranCov[2]-ranCov[1]))

solCov <- smooth.FEM(observations = dataCov,
                     FEMbasis = FEMbasis, covariates=cov,lambda=lambda)
                     #lambda.selection.criterion = 'newton', DOF.evaluation = 'exact')
sol_approxCov <- solCov$fit.FEM

#image(solCov$fit.FEM)


sol1Cov <- smooth.FEM(observations = data1Cov,
                      FEMbasis = FEMbasis, covariates = cov, lambda = lambda)
sol_approx1Cov <- sol1Cov$fit.FEM
#plot(sol_approx1Cov)

femsol_approxCov <- t(eval.FEM(sol_approxCov,nodes))
femsol_approx1Cov <- t(eval.FEM(sol_approx1Cov,nodes))

errorCov <- matrix(nrow=dim(nodes)[1], ncol= length(lambda))
for (i in 1:length(lambda))
  errorCov[,i] <- abs(coeffCov-t(femsol_approxCov)[,i])
errCov5<- sqrt(colMeans(errorCov^2))
plot(-12:3, errCov1, pch=4, xlab='log_10(lambda)', ylab='error', col=1)
grid()
abline(v=log10(best_lambdaCov1), lty=2, col='red')
lines(-12:3, errCov1, col=1)
points(log10(best_lambdaCov1), errCov1[11], pch=19, col='red')
points(-12:3, errCov2, pch=4, col=2)
lines(-12:3, errCov2, col=2)
points(log10(best_lambdaCov1), errCov2[11], pch=19, col=2)
lines(-12:3, errCov4, col=4)
points(-12:3, errCov3, pch=4, col=3)
points(log10(best_lambdaCov3), errCov3[11], pch=1, col=3)
points(-12:3, errCov5, pch=1, col=4)
#lines(-12:3, errCov3, col=3)
points(log10(best_lambdaCov4), errCov4[11], pch=19, col=4)
legend('topleft',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner','block preconditioner'), cex=0.62,col=1:4, pch=19)

residualCov <- matrix(nrow=dim(nodes)[1], ncol= length(lambda))
residualCov <- abs(t(femsol_approx1Cov)-t(femsol_approxCov))
resCov5 <- colMeans(residualCov)
plot(-12:3, resCov5, pch=4, xlab='log_10(lambda)', ylab='residual', col=1)
grid()
abline(v=log10(best_lambdaCov1[1]), lty=2, col='red')
lines(-12:3, resCov5, col=1)
points(log10(best_lambdaCov1[1]), resCov1[13], pch=19, col='red')
points(-12:3, resCov6, pch=4, col=2)
lines(-12:3, resCov6, col=2)
points(log10(best_lambdaCov2), resCov2[10], pch=19, col=2)
points(-12:3, resCov3, pch=4, col=3)
#lines(-12:3, resCov3, col=3)
points(log10(best_lambdaCov3), resCov3[10], pch=19, col=3)
points(-12:3, resCov4, pch=1, col=4)
lines(-12:3, resCov4, col=4)
points(log10(best_lambdaCov4), resCov4[10], pch=19, col=4)
legend('topright',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner','block preconditioner','lumping + block preconditioner'), cex=0.62,col=1:5, pch=19)

solCov <- smooth.FEM(observations = dataCov, FEMbasis = FEMbasis, covariates=cov,
                     lambda.selection.criterion = 'newton', DOF.evaluation = 'exact')
sol_approxCov <- solCov$fit.FEM
best_lambdaCov5 <- solCov$optimization$lambda_solution

RMSE<-function(f,g) sqrt(mean((f-g)^2))
xeval=runif(10000,0,1)
yeval=runif(10000,0,1)
rmseCov <- NULL
for (i in 1:20)
{
  output_CPP<-smooth.FEM(observations=dataCov+rnorm(nnodes, mean=0, sd=0.01*abs(ranCov[2]-ranCov[1])),
                         FEMbasis=FEMbasis, covariates=cov, lambda.selection.criterion='newton',
                         DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
  best_lambdaCov5<-c(best_lambdaCov5,output_CPP$optimization$lambda_solution)
  sol_approxCov=output_CPP$fit.FEM
  rmseCov <- c(rmseCov, RMSE(zCov(cbind(xeval,yeval)),eval.FEM(sol_approxCov,locations=cbind(xeval,yeval))))
}
boxplot(rmseCov)

rmseCovE <- rmseCov

boxplot(rmseCovA,rmseCovB,rmseCovC,rmseCovD, names=c('fdaPDE','mass lumping','lambda','block'), col=c('grey',2,3,4),ylab='RMSE')


###### execution times #####
tA=microbenchmark(smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=best_lambda7[1]),times=30)
t2D=microbenchmark(smooth.FEM(observations=data2, FEMbasis=FEMbasis2, lambda=best_lambda21))
tcovA=microbenchmark(smooth.FEM(observations=dataCov, FEMbasis=FEMbasis, covariates=cov, lambda=best_lambdaCov1[1]),times=30)

boxplot(log(tA$time),log(tB$time),log(tC$time),log(tD$time), names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')
boxplot(log(t2A$time),log(t2B$time),log(t2C$time), names=c('fdaPDE','mass lumping','lambda'),col=c('grey',2,3), ylab='log(time)')
boxplot(log(tcovA$time),log(tcovF$time),log(tcovC$time), log(tcovD$time),names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')

boxplot(tA$time,tB$time,tC$time,tD$time, names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='time')
boxplot(t2A$time,t2B$time,t2C$time,t2D$time, names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='time')
boxplot(tcovA$time,tcovB$time,tcovC$time, names=c('fdaPDE','mass lumping','lambda'),col=c('grey',2,3), ylab='time')
