best_lambdaO =NULL
err=NULL
res=NULL
rmseOOO =NULL
tD = NULL


for (k in seq(1,7,by=2))
{
x = seq(0,1, length.out =k*10)
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

# Add error to simulate data
set.seed(7893475)
ran = range(sol_exact)
data = sol_exact
data1 = sol_exact + rnorm(nnodes, mean=0, sd=0.01*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^(-12:3)

#### Test 1.1: Without GCV
# output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda)
# #image(output_CPP$fit.FEM)
# sol_approx = output_CPP$fit.FEM
# femsol_approx<-t(eval.FEM(sol_approx,nodes))
# 
# output_CPP1<-smooth.FEM(observations=data1, FEMbasis=FEMbasis, lambda=lambda)
# sol_approx1 = output_CPP1$fit.FEM
# femsol_approx1<-t(eval.FEM(sol_approx1,nodes))

### Test 1.4: Newton exact method with exact GCV, default initial lambda and tolerance
# output_CPP<-smooth.FEM(observations=data1, FEMbasis=FEMbasis, lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
# best_lambda8<-c(best_lambda8,output_CPP$optimization$lambda_solution)

# error <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
# for (i in 1:length(lambda))
#   error[,i] <- abs(sol_exact-t(femsol_approx)[,i])
# err1 <- cbind(err1,colMeans(error^2))

# residual <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
# residual <- abs(t(femsol_approx1)-t(femsol_approx))
# res1 <- cbind(res1,colMeans(residual))

xeval=runif(10000,0,1)
yeval=runif(10000,0,1)
RMSE<-function(f,g) sqrt(mean((f-g)^2))
rmse <- NULL
#rmse <- c(rmse, RMSE(f(xeval,yeval),eval.FEM(output_CPP$fit.FEM,locations=cbind(xeval,yeval))))

# for (i in 1)
# {
#   output_CPP<-smooth.FEM(observations=data+ rnorm(nnodes, mean=0, sd=0.01*abs(ran[2]-ran[1])),
#                          FEMbasis=FEMbasis, lambda.selection.criterion='newton',
#                          DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
#   rmse <- c(rmse, RMSE(f(xeval,yeval),eval.FEM(output_CPP$fit.FEM,locations=cbind(xeval,yeval))))
# }
# 
# boxplot(rmse)
# rmseOOO <- cbind(rmseOOO,rmse)

#tD <- cbind(tD,microbenchmark(smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=best_lambda7[length(best_lambda7)]),times=10)$time)

}

for (i in 1:4)
{
  plot(log10(lambda), err[,i], pch=4, xlab='log_10(lambda)', ylab='error', col='grey')
  grid()
  lines(-12:3, err[,i], lty=2, col='grey')
  points(-12:3, err1[,i], pch=4, col=1)
  abline(v=log10(best_lambda1[i]), lty=2, col='red')
  lines(-12:3, err1[,i], col=1)
  points(log10(best_lambda1[i]), err1[12,i]/2+err1[13,i]/2, pch=19, col='red')
  points(-12:3, err2[,i], pch=4, col=2)
  lines(-12:3, err2[,i], col=2)
  points(log10(best_lambda1[i]), err2[12,i]/2+err2[13,i]/2, pch=3, col='red')
  #lines(-12:3, err3[,i], col=3)
  points(-12:3, err3[,i], pch=4, col=3)
  points(log10(best_lambda3[i]), err3[12,i]/2+err3[13,i]/2, pch=19, col='red')
  points(-12:3, err4[,i], pch=10, col=4)
  points(log10(best_lambda4[i]), err4[12,i]/2+err4[13,i]/2, pch=19, col='red')
  points(-12:3, err5[,i], pch=19, col=5)
  points(log10(best_lambda5[i]), err5[12,i]/2+err5[13,i]/2, pch=19, col='red')
  legend('topleft',legend=c('original','no preconditioner', 'mass lumping', 'lambda preconditioner','block preconditioner'), cex=0.62,col=c('grey',1,2,3,4), pch=19)
  
  plot(log10(lambda),res[,i], pch=4, xlab='log_10(lambda)', ylab='residual', col='grey')
  grid()
  lines(-12:3, res[,i], lty=2, col='grey')
  points(-12:3, res[,i], pch=4, col=1)
  abline(v=log10(best_lambda1[i]), lty=2, col='red')
  lines(-12:3, res1[,i], col=1)
  points(log10(best_lambda1[i]), res1[12,i]/2+res1[13,i]/2, pch=19, col='red')
  points(-12:3, res2[,i], pch=4, col=2)
  lines(-12:3, res2[,i], col=2)
  points(log10(best_lambda1[i]), res2[12,i]/2+res2[13,i]/2, pch=3, col='red')
  #lines(-12:3, res3[,i], col=3)
  points(-12:3, res3[,i], pch=4, col=3)
  points(log10(best_lambda3[i]), res3[12,i]/2+res3[13,i]/2, pch=19, col='red')
  points(-12:3, res4[,i], pch=10, col=4)
  points(log10(best_lambda4[i]), res4[12,i]/2+res4[13,i]/2, pch=19, col='red')
  points(-12:3, res5[,i], pch=19, col=5)
  points(log10(best_lambda5[i]), res5[12,i]/2+res5[13,i]/2, pch=19, col='red')
  legend('topright',legend=c('original','no preconditioner', 'mass lumping', 'lambda preconditioner','block preconditioner'), cex=0.62,col=c('grey',1,2,3,4), pch=19)
}


plot((seq(1,7,by=2)*10)^2, log(colMeans(rmseOO[6:11,])), xlab='number of nodes', ylab='log(RMSE)',pch=19, col='red')
grid()
points((seq(1,7,by=2)*10)^2, log(colMeans(rmseA)), xlab='number of nodes', ylab='log(mse)',pch=1, col=1)
points((seq(1,7,by=2)*10)^2, log(colMeans(rmseB)), xlab='number of nodes', ylab='log(mse)',pch=19, col=2)
points((seq(1,7,by=2)*10)^2, log(colMeans(rmseC)), xlab='number of nodes', ylab='log(mse)',pch=19, col=3)
points((seq(1,7,by=2)*10)^2, log(colMeans(rmseD)), xlab='number of nodes', ylab='log(mse)',pch=4, col=4)
points((c((1:4)*10)^2), log(colMeans(rmseE)), xlab='number of nodes', ylab='log(mse)',pch=4, col=5)
lines((seq(1,7,by=2)*10)^2, log(colMeans(rmseO)), xlab='number of nodes', ylab='log(mse)', col='grey')
lines((seq(1,7,by=2)*10)^2, log(colMeans(rmseA)), xlab='number of nodes', ylab='log(mse)', col=1)
lines((seq(1,7,by=2)*10)^2, log(colMeans(rmseB)), xlab='number of nodes', ylab='log(mse)', col=2)
lines((seq(1,7,by=2)*10)^2, log(colMeans(rmseC)), xlab='number of nodes', ylab='log(mse)', col=3)
lines((seq(1,7,by=2)*10)^2, log(colMeans(rmseD)), xlab='number of nodes', ylab='log(mse)', col=4)
lines((seq(1,7,by=2)*10)^2, log(colMeans(rmseE)), xlab='number of nodes', ylab='log(mse)', col=5)
legend('topright',legend=c('original','no preconditioner', 'mass lumping', 'lambda preconditioner','block preconditioner'), cex=0.62,col=c('grey',1,2,3,4), pch=19)

plot(seq(1,11,by=2)^2, log(colMeans(tA[,1:4])), xlab='number of nodes', ylab='log(time)', pch=19,col='grey', ylim=c(min(log(colMeans(tB))),max(log(colMeans(tD)))+5))
grid()
plot(seq(1,7,by=2)^2, log(colMeans(tA[,1:4])), xlab='number of nodes', ylab='log(time)',pch=19,col=5)
points(seq(1,7,by=2)^2, log(colMeans(tB[,1:4])), xlab='number of nodes', ylab='log(time)',pch=19,col=2)
points(seq(1,7,by=2)^2, log(colMeans(tC[,1:4])), xlab='number of nodes', ylab='log(time)',pch=19,col=3)
points(seq(1,7,by=2)^2, log(tD[1:4]), xlab='number of nodes', ylab='log(time)',pch=19,col=4)
points(seq(1,7,by=2)^2, log(tF), xlab='number of nodes', ylab='log(time)',pch=19,col=5)
lines(c(((1:4)*10)^2), log(colMeans(t)), xlab='number of nodes', ylab='log(time)',pch=19,col='grey')
lines(c(((1:4)*10)^2), log(colMeans(tA)), xlab='number of nodes', ylab='log(time)',pch=19,col=1)
lines(c(((1:4)*10)^2), log(colMeans(tB)), xlab='number of nodes', ylab='log(time)',pch=19,col=2)
lines(c(((1:4)*10)^2), log(colMeans(tC)), xlab='number of nodes', ylab='log(time)',pch=19,col=3)
lines(c(((1:4)*10)^2), log(colMeans(tD)), xlab='number of nodes', ylab='log(time)',pch=19,col=4)
lines(c(((1:4)*10)^2), log(colMeans(tE)), xlab='number of nodes', ylab='log(time)',pch=19,col=5)
legend('topleft',legend=c('original','no preconditioner', 'mass lumping', 'lambda preconditioner', 'block preconditioner','block no class'), cex=0.62,col=c('grey',1,2,3,4,5), pch=19)

