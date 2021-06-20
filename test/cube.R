data(cube3D)
mesh = cube3D
plot(mesh)
nodes=mesh$nodes
nnodes=dim(mesh$nodes)[1]

FEMbasis=create.FEM.basis(mesh)

# Test function
f = function(x, y, z){
  coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
  sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
}

# Exact solution (pointwise at nodes)
coeff=f(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3])
plot(FEM(coeff, FEMbasis))

# Set locations and add error to simulate data
x = seq(0,1, length.out = 11)
y = x
z = x
locations = expand.grid(x,y,z)
sol_exact = f(locations[,1], locations[,2], locations[,3])
set.seed(7893475)
ran=range(sol_exact)
data = sol_exact
data1 = data + rnorm(length(sol_exact), mean=0, sd=0.01*abs(ran[2]-ran[1]))

lambda = 10^(-12:3)

#### Test 1.1: Without GCV
output_CPP<-smooth.FEM(locations=locations,observations=data, FEMbasis=FEMbasis, lambda=lambda[1])
#image(output_CPP$fit.FEM)
sol_approx = output_CPP$fit.FEM
femsol_approx<-t(eval.FEM(sol_approx,nodes))

output_CPP1<-smooth.FEM(observations=data1, FEMbasis=FEMbasis, lambda=lambda)
sol_approx1 = output_CPP1$fit.FEM
femsol_approx1<-t(eval.FEM(sol_approx1,nodes))


### Test 1.4: Newton exact method with exact GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
best_lambda3<-output_CPP$optimization$lambda_solution


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
