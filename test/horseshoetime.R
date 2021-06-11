data(horseshoe2D)

mesh = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
locations = refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)

# plot spatial locations
plot(mesh, pch = ".")
points(locations, pch = 16, col = "red")

FEMbasis = create.FEM.basis(mesh)

f<-function(x,y,t)
{
  K <- (y/0.1*as.double((abs(y)<=0.1 & x>-0.5))+as.double((abs(y)>0.1 | x<=-0.5)))^2
  res=numeric(length =length(x))
  for(i in 1:length(x))
  {
    if(x[i]>=0 && y[i]>0)
      res[i]=cos(t[i])*(0.25*pi+x[i])+(y[i]-0.5)^2
    if(x[i]>=0 && y[i]<=0)
      res[i]=cos(2*t[i])*(-0.25*pi-x[i])+(-y[i]-0.5)^2
    if(x[i]<0 && y[i]>0)
      res[i]=cos(t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
    if(x[i]<0 && y[i]<=0)
      res[i]=cos(2*t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
  }
  res
}

NumTimeInstants=5
TimePoints=seq(0,pi,length.out =NumTimeInstants)

space_time_locations = cbind(rep(TimePoints,each=nrow(locations)),rep(locations[,1],NumTimeInstants),rep(locations[,2],NumTimeInstants))
sol_exact = f(space_time_locations[,2],space_time_locations[,3],space_time_locations[,1])

ndata = length(sol_exact)

# Create covariates
set.seed(509875)
cov1 = rnorm(ndata, mean = 1, sd = 2)

# Add error to simulate data
set.seed(7893475)
data = sol_exact + 2*cov1 
data1 = data + rnorm(length(sol_exact), mean = 0, sd =  0.01*diff(range(sol_exact)))
observations = matrix(data,nrow(locations),NumTimeInstants)
observations1 = matrix(data1,nrow(locations),NumTimeInstants)

# Set smoothing parameter
lambdaS = 10^-2
lambdaT = 10^-2

#### Test 2.1: Without GCV
output_CPP<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                            observations=observations, 
                            covariates = cov1,
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT)

#plot(output_CPP$fit.FEM)
sol_approx = output_CPP$fit.FEM.time
femsol_approx<-t(eval.FEM.time(sol_approx,space.time.locations = cbind(TimePoints,locations)))

output_CPP1<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                        observations=observations1, 
                        covariates = cov1,
                        FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT)
sol_approx1 = output_CPP1$fit.FEM.time
femsol_approx1<-t(eval.FEM.time(sol_approx1,space.time.locations = cbind(TimePoints,locations)))


#### Test 2.2: exact GCV
lambdaS = 10^(-1:1)
lambdaT = 10^(-6:-4)
output_CPP<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                            observations=observations, 
                            covariates = cov1,
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT, 
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
