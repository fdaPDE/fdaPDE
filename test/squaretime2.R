x = seq(0,1, length.out = 21)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)
plot(mesh)

nnodes = dim(mesh$nodes)[1]

FEMbasis = create.FEM.basis(mesh)

# Test function
f = function(x, y, z)
{
  a1=1
  a2=4
  (a1*sin(2*pi*x)*cos(2*pi*y)+a2*sin(3*pi*x))*cos(pi*z)
}

NumTimePoints=11
TimePoints=seq(0,2,length.out = NumTimePoints)

SpacePoints=mesh$nodes
SpaceTimePoints=cbind(rep(mesh$nodes[,1],NumTimePoints),rep(mesh$nodes[,2],NumTimePoints),rep(TimePoints,each=nrow(mesh$nodes)))

# Exact solution (pointwise at nodes)
sol_exact=f(SpaceTimePoints[,1],SpaceTimePoints[,2],SpaceTimePoints[,3])

# Add error to simulate data
set.seed(7893475)
ran = range(sol_exact)
data =sol_exact
data1 = data + rnorm(nnodes, mean=0, sd=0.01*abs(ran[2]-ran[1]))
observations=matrix(data,nrow(SpacePoints),NumTimePoints)
observations1=matrix(data1,nrow(SpacePoints),NumTimePoints)

# Set PDE parameters
PDE_parameters = list(K = matrix(c(1,0,0,4), nrow = 2), b = c(0,0), c = 0)

# Set smoothing parameter
lambdaS = 1e-4
lambdaT = 1e-4

#### Test 3.1: Without GCV
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                            PDE_parameters=PDE_parameters)
sol_approx = output_CPP$fit.FEM.time
femsol_approx=t(eval.FEM.time(sol_approx,space.time.locations = cbind(TimePoints,x,y)))

output_CPP1<-smooth.FEM.time(time_mesh = TimePoints, observations=observations1, 
                             FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                             PDE_parameters=PDE_parameters)
sol_approx1 = output_CPP1$fit.FEM.time
femsol_approx1<-t(eval.FEM.time(sol_approx1,space.time.locations = cbind(TimePoints,x,y)))

output_CPP<-smooth.FEM(observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

lambdaS = 10^(-5:-3)
lambdaT = 10^(-5:-3)
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                            PDE_parameters=PDE_parameters,
                            lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
best_lambda = which(output_CPP$GCV == min(output_CPP$GCV), arr.ind = TRUE)

error <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
for (i in 1:length(lambda))
  error[,i] <- abs(sol_exact-t(femsol_approx)[,i])
err3 <- colMeans(error^2)
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
res3 <- colMeans(residual)
plot(log10(lambda), res1, pch=4, xlab='log_10(lambda)', ylab='residual', col=1, ylim=c(0,0.026))
grid()
abline(v=log10(best_lambda1), lty=2, col='red')
lines(-12:3, res1, col=1)
points(log10(best_lambda1), res1[12], pch=19, col='red')
points(-12:3, res2, pch=4, col=2)
lines(-12:3, res2, col=2)
points(log10(best_lambda1), res2[12], pch=19, col='red')
lines(-12:3, res3, col=3)
points(-12:3, res3, pch=4, col=3)
points(log10(best_lambda3), res3[12], pch=19, col='red')
legend('topright',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner'), cex=0.62,col=1:3, pch=19)
