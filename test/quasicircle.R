# QUASICIRCLE2D

data(quasicircle2D)
boundary_nodes = quasicircle2D$boundary_nodes
boundary_segments = quasicircle2D$boundary_segments
locations = quasicircle2D$locations
data = quasicircle2D$data
ran = range(data)
data1 = data + rnorm (length(data), sd=0.01*(ran[2]-ran[1]))
mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
FEMbasis = create.FEM.basis(mesh)

mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
fembasis = create.FEM.basis(mesh)
lambda = 10^(-12:3)

R = 2.8 
K1 = 0.1
K2 = 0.2
beta = 0.5

K_func<-function(points)
{
  output = array(0, c(2, 2, nrow(points)))
  for (i in 1:nrow(points))
    output[,,i]=10*rbind(c(points[i,2]^2+K1*points[i,1]^2+K2*(R^2-points[i,1]^2-points[i,2]^2),
                           (K1-1)*points[i,1]*points[i,2]),
                         c((K1-1)*points[i,1]*points[i,2],
                           points[i,1]^2+K1*points[i,2]^2+K2*(R^2-points[i,1]^2-points[i,2]^2)))
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

BC = NULL
BC$BC_indices = which(mesh$nodesmarkers == 1) # b.c. on the complete boundary
BC$BC_values = rep(0,length(BC$BC_indices)) # homogeneus b.c.
# 
# dataNA = rep(NA, fembasisC$nbasis)
# dataNA[mesh$nodesmarkers == 0] = dataC

sol = smooth.FEM(observations = data,  FEMbasis = FEMbasis, locations = locations, PDE_parameters = PDE_parameters, BC = BC,
                  #lambda.selection.criterion = 'newton_fd',DOF.evaluation = 'exact')
                  lambda=lambda)#,lambda.selection.criterion = 'grid', DOF.evaluation = 'exact')
sol_approx <- sol$fit.FEM
best_lambda8 <- sol$optimization$lambda_solution

sol1 = smooth.FEM(observations = data1,  FEMbasis = FEMbasis, lambda=lambda,locations = locations,
                   PDE_parameters = PDE_parameters, BC = BC)
sol_approx1 <- sol1$fit.FEM

femsol_approx <- t(eval.FEM(sol_approx,locations))
femsol_approx1 <- t(eval.FEM(sol_approx1,locations))

error <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
for (i in 1:length(lambda))
  error[,i] <- abs(data-t(femsol_approx)[,i])
err8 <- colMeans(error^2)
plot(log10(lambda), err5, pch=4, xlab='log_10(lambda)', ylab='error', col=1,ylim=c(0,0.031))
grid()
abline(v=log10(best_lambda1), lty=2, col=1)
lines(-12:3, err5, col=1)
points(log10(best_lambda1), err1[10], pch=19, col=1)
points(-12:3, err6, pch=4, col=2)
abline(v=log10(best_lambda2), lty=2, col=2)
lines(-12:3, err6, col=2)
points(log10(best_lambda2), err2[10], pch=19, col=2)
lines(-12:3, err7, col=3)
points(-12:3, err3, pch=3, col=3)
abline(v=log10(best_lambda3), lty=2, col=3)
points(log10(best_lambda3), err3[10], pch=1, col=3)
lines(-12:3, err8, col=4)
points(-12:3, err4, pch=1, col=4)
abline(v=log10(best_lambda4), lty=2, col=4)
points(log10(best_lambda4), err4[7], pch=19, col=4)
legend('topleft',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner','block preconditioner'), cex=0.62,col=1:4, pch=19)

residual <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
residual <- abs(t(femsol_approx1)-t(femsol_approx))
res8 <- colMeans(residual)
plot(-12:3, res1, pch=4, xlab='log_10(lambda)', ylab='residual', col=1)
grid()
abline(v=log10(best_lambda1), lty=2, col=1)
lines(-12:3, res1, col=1)
points(log10(best_lambda1), res1[8], pch=19, col=1)
points(-12:3, res2, pch=4, col=2)
abline(v=log10(best_lambda2), lty=2, col=2)
lines(-12:3, res2, col=2)
points(log10(best_lambda2), 2*res2[9]/3+res2[8]/3, pch=19, col=2)
lines(-12:3, res8, col=4)
points(-12:3, res3, pch=4, col=3)
abline(v=log10(best_lambda4), lty=2, col=4)
points(log10(best_lambda4), res4[8], pch=19, col=4)
legend('topright',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner','block preconditioner'), cex=0.62,col=1:4, pch=19)



# AREAL DATA

data(quasicircle2Dareal)
incidence_matrix = quasicircle2Dareal$incidence_matrix
data_areal = quasicircle2Dareal$data
ran_areal = range(data_areal)
data1_areal <- data_areal + rnorm (length(data_areal), sd=0.05*(ran_areal[2]-ran_areal[1]))

mesh_areal = quasicircle2Dareal$mesh
nodes_areal = mesh_areal$nodes

fembasis_areal = create.FEM.basis(mesh_areal)
lambda = 10^(-12:3)

R_areal = 2.8 
K1_areal = 0.1
K2_areal = 0.1
beta_areal = 0.5

K_func_areal<-function(points)
{
  output = array(0, c(2, 2, nrow(points)))
  for (i in 1:nrow(points))
    output[,,i]=10*rbind(c(points[i,2]^2+K1_areal*points[i,1]^2+K2_areal*(R_areal^2-points[i,1]^2-points[i,2]^2),
                           (K1_areal-1)*points[i,1]*points[i,2]),
                         c((K1_areal-1)*points[i,1]*points[i,2],
                           points[i,1]^2+K1_areal*points[i,2]^2+K2_areal*(R_areal^2-points[i,1]^2-points[i,2]^2)))
  output
}

b_func_areal<-function(points)
{
  output = array(0, c(2, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = 10*beta_areal*c(points[i,1],points[i,2])
  output
}

c_func_areal<-function(points)
{
  rep(c(0), nrow(points))
}

u_func_areal<-function(points)
{
  rep(c(0), nrow(points))
}
PDE_parameters_areal = list(K = K_func_areal, b = b_func_areal, c = c_func_areal, u = u_func_areal)

BC_areal = NULL
BC_areal$BC_indices = which(mesh_areal$nodesmarkers == 1) # b.c. on the complete boundary
BC_areal$BC_values = rep(0,length(BC_areal$BC_indices)) # homogeneus b.c.

sol_areal = smooth.FEM(observations = data_areal, incidence_matrix = incidence_matrix, FEMbasis = fembasis_areal,
                   #lambda.selection.criterion = 'newton_fd', DOF.evaluation = 'exact', PDE_parameters = PDE_parameters_areal, BC = BC_areal)
                   lambda=lambda, PDE_parameters = PDE_parameters_areal, BC = BC_areal)
sol_approx_areal <- sol_areal$fit.FEM
best_lambda_areal8 <- sol_areal$optimization$lambda_solution

sol1_areal = smooth.FEM(observations = data1_areal, incidence_matrix = incidence_matrix, FEMbasis = fembasis_areal, lambda = lambda,
                    PDE_parameters = PDE_parameters_areal, BC = BC_areal)
sol_approx1_areal <- sol1_areal$fit.FEM

femsol_approx_areal <- t(eval.FEM(sol_approx_areal,mesh_areal$nodes))
femsol_approx1_areal <- t(eval.FEM(sol_approx1_areal,mesh_areal$nodes))

residual_areal <- matrix(nrow=dim(nodes_areal)[1], ncol= length(lambda))
residual_areal <- abs(t(femsol_approx1_areal)-t(femsol_approx_areal))
res_areal8 <- colMeans(residual_areal)
plot(-12:3,res_areal, pch=4, xlab='log_10(lambda)', ylab='reisdual', col=1)
grid()
abline(v=log10(best_lambda_areal1), lty=2, col=1)
lines(-12:3, res_areal1, col=1)
points(log10(best_lambda_areal1), res_areal1[8], pch=19, col=1)
points(-12:3, res_areal2, pch=4, col=2)
abline(v=log10(best_lambda_areal2), lty=2, col=2)
lines(-12:3, res_areal2, col=2)
points(log10(best_lambda_areal2), res_areal2[11]/2+res_areal2[10]/2, pch=19, col=2)
lines(-12:3, res_areal3, col=3)
points(-12:3, res_areal3, pch=4, col=3)
points(log10(best_lambda_areal3), res_areal3[8], pch=1, col=3)
lines(-12:3, res_areal4, col=4)
points(-12:3, res_areal4, pch=4, col=4)
points(log10(best_lambda_areal4), res_areal4[8], pch=19, col=4)
legend('bottomleft',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner','block preconditioner'), cex=0.62,col=1:4, pch=19)

###### execution times #####
tD=microbenchmark(smooth.FEM(observations = data,  FEMbasis = FEMbasis, locations = locations,
                  lambda=best_lambda2, PDE_parameters = PDE_parameters, BC = BC))
t2D=microbenchmark(smooth.FEM(observations = data_areal, incidence_matrix = incidence_matrix, FEMbasis = fembasis_areal,
                      lambda=best_lambda_areal2, PDE_parameters = PDE_parameters_areal, BC = BC_areal))

boxplot(log(tE$time),log(tF$time),log(tG$time), log(tD$time),names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')
boxplot(log(t2E$time),log(t2F$time),log(t2G$time),log(t2D$time), names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')

boxplot(tA$time,tB$time,tC$time,tD$time, names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='time')
boxplot(t2A$time,t2B$time,t2C$time,t2D$time, names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='time')
