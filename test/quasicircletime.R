data(quasicircle2Dareal)
mesh = quasicircle2Dareal$mesh
incidence_matrix = quasicircle2Dareal$incidence_matrix
data = quasicircle2Dareal$data

plot(mesh)
FEMbasis = create.FEM.basis(mesh)
time_mesh=seq(0,4,length.out = 11)

data=rep(data,length(time_mesh))*rep(exp(time_mesh),each=length(data))

# Add error to simulate data
set.seed(5839745)
data1 = data + rnorm(length(data), sd = 0.01)
observations=matrix(data,nrow(incidence_matrix),length(time_mesh))
observations1=matrix(data1,nrow(incidence_matrix),length(time_mesh))

# Set smoothing parameter
lambdaS = 10^-6
lambdaT = 10^-6

# Set BC (equal to zero)
BC = NULL
BC$BC_indices = which(mesh$nodesmarkers == 1)
BC$BC_values = rep(0,length(BC$BC_indices))

# Set sv-PDE parameters
R = 2.8
K1 = 0.1
K2 = 0.2
beta = 0.5

K_func<-function(points)
{
  output = array(0, c(2, 2, nrow(points)))
  for (i in 1:nrow(points))
    output[,,i] = 10*rbind(c(points[i,2]^2 + K1*points[i,1]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2),
                             (K1-1)*points[i,1]*points[i,2]),
                           c((K1-1)*points[i,1]*points[i,2],
                             points[i,1]^2 + K1*points[i,2]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2)))
  output
}

b_func<-function(points)
{
  output = array(0, c(2, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = 10*beta*c(points[i,1],points[i,2])
  output
}

c_func<-function(points)
{
  rep(c(0), nrow(points))
}

u_func<-function(points)
{
  rep(c(0), nrow(points))
}
PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)


#### Test 4.1.2: Forcing term = 0, estimated IC, without GCV
output_CPP<-smooth.FEM.time(observations=observations,
                            incidence_matrix = incidence_matrix,
                            FEMbasis = FEMbasis, time_mesh = time_mesh,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            BC = BC,
                            PDE_parameters = PDE_parameters,
                            FLAG_PARABOLIC = TRUE,
                            DOF.evaluation = NULL)
sol_approx = output_CPP$fit.FEM.time
femsol_approx<-t(eval.FEM.time(sol_approx,space.time.locations = cbind(TimePoints,mesh$nodes)))

output_CPP1<-smooth.FEM.time(observations=observations1,
                            incidence_matrix = incidence_matrix,
                            FEMbasis = FEMbasis, time_mesh = time_mesh,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            BC = BC,
                            PDE_parameters = PDE_parameters,
                            FLAG_PARABOLIC = TRUE,
                            DOF.evaluation = NULL)
sol_approx1 = output_CPP$fit.FEM.time
femsol_approx1<-t(eval.FEM.time(sol_approx1,space.time.locations = cbind(TimePoints,mesh$nodes)))

output_CPP<-smooth.FEM.time(observations=observations,
                            incidence_matrix = incidence_matrix, areal.data.avg = TRUE,
                            FEMbasis = FEMbasis, time_mesh = time_mesh,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            BC = BC,
                            PDE_parameters = PDE_parameters,
                            FLAG_PARABOLIC = TRUE,
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
