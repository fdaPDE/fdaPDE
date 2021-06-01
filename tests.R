#### No prior information about anysotropy/non-stationarity (laplacian smoothing) ####
data(horseshoe2D)
boundary_nodes = horseshoe2D$boundary_nodes
boundary_segments = horseshoe2D$boundary_segments
locations = horseshoe2D$locations

mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
FEMbasis = create.FEM.basis(mesh)
lambda = 10^-1
# no covariate
data = fs.test(mesh$nodes[,1], mesh$nodes[,2]) + rnorm(nrow(mesh$nodes), sd = 0.5)

solution = smooth.FEM(observations = data, FEMbasis = FEMbasis, lambda = lambda)
plot(solution$fit.FEM)

# with covariates
covariate = covs.test(mesh$nodes[,1], mesh$nodes[,2])
data = fs.test(mesh$nodes[,1], mesh$nodes[,2]) + 2*covariate + rnorm(nrow(mesh$nodes), sd = 0.5)

solution = smooth.FEM(observations = data, covariates = covariate, 
                      FEMbasis = FEMbasis, lambda = lambda)

# beta estimate:
solution$solution$beta
# non-parametric estimate:
plot(solution$fit.FEM)

# Choose lambda with GCV - stochastic grid evaluation:
lambda = 10^(-12:3)
solution = smooth.FEM(observations = data,
                      covariates = covariate,
                      FEMbasis = FEMbasis,
                      lambda = lambda)
bestLambda = solution$optimization$lambda_solution
# Choose lambda with GCV - Newton finite differences stochastic evaluation -:
solution = smooth.FEM(observations = data,
                      covariates = covariate,
                      FEMbasis = FEMbasis,
                      DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV')
bestLambda = solution$optimization$lambda_solution

data1 <- data + rnorm (nrow(mesh$nodes), sd=0.01)

solution1 = smooth.FEM(observations = data1,
                      covariates = covariate,
                      FEMbasis = FEMbasis, lambda=lambda)

sol_approx <- solution$fit.FEM
sol_approx1 <- solution1$fit.FEM

sol_exact <- fs.test(locations[,1],locations[,2])

femsol_approx <- t(eval.FEM(sol_approx,locations))
error <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
for (i in 1:length(lambda))
  error[,i] <- abs(sol_exact-t(femsol_approx)[,i])
err <- colMeans(error^2)
plot(err)

femsol_approx1 <- t(eval.FEM(sol_approx1,locations))
residual <- matrix(nrow=dim(locations)[1], ncol= length(lambda))
residual <- abs(t(femsol_approx1)-t(femsol_approx))
res <- colMeans(residual^2)
plot(res)




#### Smoothing with prior information about anysotropy/non-stationarity and boundary conditions ####
# See Azzimonti et al. for reference to the current exemple
data(quasicircle2D)
boundary_nodes = quasicircle2D$boundary_nodes
boundary_segments = quasicircle2D$boundary_segments
locations = quasicircle2D$locations
data = quasicircle2D$data

mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
FEMbasis = create.FEM.basis(mesh)
lambda = 10^-2

# Set the PDE parameters
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

# Set the boundary conditions
BC = NULL
BC$BC_indices = which(mesh$nodesmarkers == 1) # b.c. on the complete boundary
BC$BC_values = rep(0,length(BC$BC_indices)) # homogeneus b.c.

# Since the data locations are a subset of the mesh nodes for a faster solution use:
dataNA = rep(NA, FEMbasis$nbasis)
dataNA[mesh$nodesmarkers == 0] = data
#grid evaluation
solution = smooth.FEM(observations = dataNA,
                      FEMbasis = FEMbasis,
                      lambda = lambda,
                      PDE_parameters = PDE_parameters,
                      BC = BC)
plot(solution$fit.FEM)
image(solution$fit.FEM)
# Newton's method
solution = smooth.FEM(observations = dataNA,
                      FEMbasis = FEMbasis,
                      PDE_parameters = PDE_parameters,
                      BC = BC)
plot(solution$fit.FEM)
image(solution$fit.FEM)

#### Smoothing with areal data ####
# See Azzimonti et al. for reference to the current exemple
data(quasicircle2Dareal)
incidence_matrix = quasicircle2Dareal$incidence_matrix
data = quasicircle2Dareal$data
mesh = quasicircle2Dareal$mesh

FEMbasis = create.FEM.basis(mesh)
lambda = 10^-4

# Set the PDE parameters
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

# Set the boundary conditions
BC = NULL
BC$BC_indices = which(mesh$nodesmarkers == 1) # b.c. on the complete boundary
BC$BC_values = rep(0,length(BC$BC_indices)) # homogeneus b.c.
#grid evaluation
solution = smooth.FEM(observations = data,
                      incidence_matrix = incidence_matrix,
                      FEMbasis = FEMbasis,
                      lambda = lambda,
                      PDE_parameters = PDE_parameters,
                      BC = BC)
plot(solution$fit.FEM)
image(solution$fit.FEM)
#Newton's method
solution = smooth.FEM(observations = data,
                      incidence_matrix = incidence_matrix,
                      FEMbasis = FEMbasis,
                      PDE_parameters = PDE_parameters,
                      BC = BC)
plot(solution$fit.FEM)
image(solution$fit.FEM)


