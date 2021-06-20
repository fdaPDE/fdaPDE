## SQUARE ##

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
  coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
  sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
}

NumTimePoints=11
TimePoints=seq(0,2,length.out = NumTimePoints)

SpacePoints=mesh$nodes
SpaceTimePoints=cbind(rep(mesh$nodes[,1],NumTimePoints),rep(mesh$nodes[,2],NumTimePoints),rep(TimePoints,each=nrow(mesh$nodes)))

# Exact solution (pointwise at nodes)
coeff=f(SpaceTimePoints[,1],SpaceTimePoints[,2],SpaceTimePoints[,3])
# plot(FEM.time(coeff, TimePoints, FEMbasis, FLAG_PARABOLIC = TRUE))

sol_exact=NULL
for (t in TimePoints)
  sol_exact <- c(sol_exact, f(x,y,t))

# Add error to simulate data
set.seed(7893475)
ran = range(coeff)
data = coeff
data1 = data + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))
observations=matrix(data,nrow(SpacePoints),NumTimePoints)
observations1=matrix(data1,nrow(SpacePoints),NumTimePoints)

# Set smoothing parameter
lambdaS = 1e-2
lambdaT = 1e-2

#### Test 1.1: Without GCV
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                            FEMbasis=FEMbasis, lambdaS=lambdaS,
                            lambdaT=lambdaT)
sol_approx = output_CPP$fit.FEM.time
femsol_approx=t(eval.FEM.time(sol_approx,locations = cbind(x,y), time.instants=TimePoints))

output_CPP1<-smooth.FEM.time(time_mesh = TimePoints, observations=observations1, 
                        FEMbasis=FEMbasis, lambdaS=lambdaS,
                        lambdaT=lambdaT)
sol_approx1 = output_CPP1$fit.FEM.time
femsol_approx1<-t(eval.FEM.time(sol_approx1,locations = cbind(x,y), time.instants=TimePoints))


#### Test 1.3: grid with stochastic GCV
lambdaS = 10^seq(-7, -5, length.out=3)
lambdaT = 10^seq(-7, -5, length.out=3)
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, 
                            time_locations = TimePoints,
                            observations=observations, 
                            FEMbasis=FEMbasis, 
                            lambdaS=lambdaS,lambdaT=lambdaT,
                            lambda.selection.criterion='grid', 
                            DOF.evaluation='exact', 
                            lambda.selection.lossfunction='GCV')
best_lambda= c(lambdaS[which(output_CPP$GCV == min(output_CPP$GCV), arr.ind = TRUE)[1]],
               lambdaT[which(output_CPP$GCV == min(output_CPP$GCV), arr.ind = TRUE)[2]])

error <- matrix(nrow=length(sol_exact), ncol= length(lambdaS))
for (i in 1:length(lambdaS))
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
