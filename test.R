library( fdaPDE)

# mesh generation
x = seq (0, 1, length.out = 30)
y = x
locations = expand.grid(x, y)
mesh = create.mesh.2D( locations )
plot(mesh )
nodes=mesh$nodes
nnodes = dim(mesh$nodes) [ 1 ]
FEMbasis = create.FEM.basis(mesh)

# Test function
a1 = runif (1 ,min=1.5, max=1.5)
a2 = runif (1 ,min=1.5, max=1.5)
f = function (p)
{
  a1* sin (2*pi *p[ , 1]) * cos(2*pi *p[ , 2])+ a2* sin (3* pi *p[ , 1])
}

# Exact solution ( pointwise at nodes )
sol_exact = f(mesh$nodes)
image (FEM(sol_exact, FEMbasis))
# Add error to simulate data
set.seed (7893475)
ran = range (sol_exact)
data = sol_exact
data1 = sol_exact + rnorm(nnodes, mean=0,
                                 sd=0.001* abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^( -12:3)

### Computation of the residual
# no GCV

# no preconditioner
output_CPP <- smooth.FEM(observations=data,
                         FEMbasis=FEMbasis,
                         lambda=lambda)
image (output_CPP$fit.FEM)
sol_approx = output_CPP$fit.FEM
femsol_approx <- eval.FEM(sol_approx, nodes)

output_CPP1 <- smooth.FEM(observations=data1,
                          FEMbasis=FEMbasis,
                          lambda=lambda)
sol_approx1 = output_CPP1$fit.FEM
femsol_approx1 <- eval.FEM(sol_approx1, nodes)

residual <- matrix(nrow=dim(locations)[1], ncol=length(lambda))
residual <- abs(femsol_approx1 - femsol_approx)
res1 <- colMeans(residual)

# mass lumping
output_CPP <- smooth.FEM(observations=data,
                         FEMbasis=FEMbasis,
                         lambda=lambda, solver.options="mass_lumping")
image (output_CPP$fit.FEM)
sol_approx = output_CPP$fit.FEM
femsol_approx <- eval.FEM(sol_approx, nodes)

output_CPP1 <- smooth.FEM(observations=data1,
                          FEMbasis=FEMbasis,
                          lambda=lambda, solver.options="mass_lumping")
sol_approx1 = output_CPP1$fit.FEM
femsol_approx1 <- eval.FEM(sol_approx1, nodes)

residual <- matrix(nrow=dim(locations)[1], ncol=length(lambda))
residual <- abs(femsol_approx1 - femsol_approx)
res2 <- colMeans(residual)

# lambda preconditoner
output_CPP <- smooth.FEM(observations=data,
                         FEMbasis=FEMbasis,
                         lambda=lambda, solver.options="lambda_preconditioner")
image (output_CPP$fit.FEM)
sol_approx = output_CPP$fit.FEM
femsol_approx <- eval.FEM(sol_approx, nodes)

output_CPP1 <- smooth.FEM(observations=data1,
                          FEMbasis=FEMbasis,
                          lambda=lambda, solver.options="lambda_preconditioner")
sol_approx1 = output_CPP1$fit.FEM
femsol_approx1 <- eval.FEM(sol_approx1, nodes)

residual <- matrix(nrow=dim(locations)[1], ncol=length(lambda))
residual <- abs(femsol_approx1 - femsol_approx)
res3 <- colMeans(residual)

# block preconditoner
output_CPP <- smooth.FEM(observations=data,
                         FEMbasis=FEMbasis,
                         lambda=lambda, solver.options="block_preconditioner")
image (output_CPP$fit.FEM)
sol_approx = output_CPP$fit.FEM
femsol_approx <- eval.FEM(sol_approx, nodes)

output_CPP1 <- smooth.FEM(observations=data1,
                          FEMbasis=FEMbasis,
                          lambda=lambda, solver.options="block_preconditioner")
sol_approx1 = output_CPP1$fit.FEM
femsol_approx1 <- eval.FEM(sol_approx1, nodes)

residual <- matrix(nrow=dim(locations)[1], ncol=length(lambda))
residual <- abs(femsol_approx1 - femsol_approx)
res4 <- colMeans(residual)

### Best lambda corresponding to exact data
### Newton exact method with exact GCV
### default initial lambda and tolerance
# no preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV')
best_lambda1 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# mass lumping
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "mass_lumping")
best_lambda2 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# lambda preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "lambda_preconditoner")
best_lambda3 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# block preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "block_preconditoner")
best_lambda4 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )


### RMSE corresponding to the best lambda according to exact Newton
RMSE <- function (f,g) sqrt (mean((f-g)^2))
xeval = runif(1000,0,1)
yeval = runif(1000,0,1)

# no preconditioner
rmse <- NULL
for ( i in 1:30 )
{
  output_CPP <- smooth.FEM(observations=data+rnorm(nnodes, mean=0, sd=0.1*abs(ran[2]-ran[1])),
                           FEMbasis=FEMbasis, lambda.selection.criterion='newton',
                           DOF.evaluation='exact',
                           lambda.selection.lossfunction='GCV')
  best_lambda1 <- c (best_lambda1, output_CPP$optimization$lambda_solution)
  sol_approx=output_CPP$fit.FEM
  rmse <- c (rmse, RMSE( f (cbind(xeval, yeval)),
                           eval.FEM( sol_approx, locations=cbind(xeval,yeval))))
}
boxplot (rmse)
rmseA <- rmse

# mass lumping
rmse <- NULL
for ( i in 1:30 )
{
  output_CPP <- smooth.FEM(observations=data+rnorm(nnodes, mean=0, sd=0.1*abs(ran[2]-ran[1])),
                           FEMbasis=FEMbasis, lambda.selection.criterion='newton',
                           DOF.evaluation='exact',
                           lambda.selection.lossfunction='GCV',
                           solver.options = "mass_lumping")
  best_lambda2 <- c (best_lambda2, output_CPP$optimization$lambda_solution)
  sol_approx=output_CPP$fit.FEM
  rmse <- c (rmse, RMSE( f (cbind(xeval, yeval)),
                         eval.FEM( sol_approx, locations=cbind(xeval,yeval))))
}
boxplot (rmse)
rmseB <- rmse

# lambda preconditoner
rmse <- NULL
for ( i in 1:30 )
{
  output_CPP <- smooth.FEM(observations=data+rnorm(nnodes, mean=0, sd=0.1*abs(ran[2]-ran[1])),
                           FEMbasis=FEMbasis, lambda.selection.criterion='newton',
                           DOF.evaluation='exact',
                           lambda.selection.lossfunction='GCV',
                           solver.options = "lambda_preconditioner")
  best_lambda3 <- c (best_lambda3, output_CPP$optimization$lambda_solution)
  sol_approx=output_CPP$fit.FEM
  rmse <- c (rmse, RMSE( f (cbind(xeval, yeval)),
                         eval.FEM( sol_approx, locations=cbind(xeval,yeval))))
}
boxplot (rmse)
rmseC <- rmse

# block preconditioner
rmse <- NULL
for ( i in 1:30 )
{
  output_CPP <- smooth.FEM(observations=data+rnorm(nnodes, mean=0, sd=0.1*abs(ran[2]-ran[1])),
                           FEMbasis=FEMbasis, lambda.selection.criterion='newton',
                           DOF.evaluation='exact',
                           lambda.selection.lossfunction='GCV',
                           solver.options = "block_preconditioner")
  best_lambda4 <- c (best_lambda4, output_CPP$optimization$lambda_solution)
  sol_approx=output_CPP$fit.FEM
  rmse <- c (rmse, RMSE( f (cbind(xeval, yeval)),
                         eval.FEM( sol_approx, locations=cbind(xeval,yeval))))
}
rmseD <- rmse

# comparison
boxplot(rmseA,rmseB,rmseC,rmseD, names=c('fdaPDE','mass lumping','lambda','block'), col=c('grey',2,3,4), ylab="RMSE")


### RMSE on a grid of lambdas

rmse_gridA <- rmse_gridB <- rmse_gridC <- rmse_gridD <- NULL

# no preconditioner
for ( l in lambda )
{
  set.seed (1527)
  rmse <- NULL
  for(i in 1:30)
  {
    output_CPP <- smooth.FEM(observations = data+rnorm (length(nodes[1,]),
                                            sd=0.1*(ran[2]=ran[1])),
                                            FEMbasis=FEMbasis,
                                            lambda=l )
    rmse <- c(rmse, RMSE( f(cbind(xeval, yeval)) ,
                             eval.FEM( output_CPP$fit.FEM, locations=cbind(xeval,yeval))))
  }
  rmse_gridA <- cbind( rmse_gridA, rmse )
}

# mass lumping
for ( l in lambda )
{
  set.seed (1527)
  rmse <- NULL
  for(i in 1:30)
  {
    output_CPP <- smooth.FEM(observations = data+rnorm (length(nodes[1,]),
                                                        sd=0.1*(ran[2]=ran[1])),
                             FEMbasis=FEMbasis,
                             lambda=l, solver.options = "mass_lumping" )
    rmse <- c(rmse, RMSE( f(cbind(xeval, yeval)) ,
                          eval.FEM( output_CPP$fit.FEM, locations=cbind(xeval,yeval))))
  }
  rmse_gridB <- cbind( rmse_gridB, rmse )
}

# lambda preconditioner
for ( l in lambda )
{
  set.seed (1527)
  rmse <- NULL
  for(i in 1:30)
  {
    output_CPP <- smooth.FEM(observations = data+rnorm (length(nodes[1,]),
                                                        sd=0.1*(ran[2]=ran[1])),
                             FEMbasis=FEMbasis,
                             lambda=l, solver.options = "lambda_preconditioner" )
    rmse <- c(rmse, RMSE( f(cbind(xeval, yeval)) ,
                          eval.FEM( output_CPP$fit.FEM, locations=cbind(xeval,yeval))))
  }
  rmse_gridC <- cbind( rmse_gridC, rmse )
}

# block preconditioner
for ( l in lambda )
{
  set.seed (1527)
  rmse <- NULL
  for(i in 1:30)
  {
    output_CPP <- smooth.FEM(observations = data+rnorm (length(nodes[1,]),
                                                        sd=0.1*(ran[2]=ran[1])),
                             FEMbasis=FEMbasis,
                             lambda=l, solver.options = "block_preconditioner" )
    rmse <- c(rmse, RMSE( f(cbind(xeval, yeval)) ,
                          eval.FEM( output_CPP$fit.FEM, locations=cbind(xeval,yeval))))
  }
  rmse_gridD <- cbind( rmse_gridD, rmse )
}

# comparison
boxplot(rmse_gridA,rmse_gridB,rmse_gridC,rmse_gridD, names=c('fdaPDE','mass lumping','lambda','block'), col=c('grey',2,3,4), ylab="RMSE")


### Computational times
# smooth.fem given exact data, lambda fixed at its optimal value
# for exact data according to exact Newton
library(microbenchmark)

# no preconditioner
tA=microbenchmark(smooth.FEM(observations=data, FEMbasis=FEMbasis,
    lambda=best_lambda1[1]), times=30)

# mass lumping
tB=microbenchmark(smooth.FEM(observations=data, FEMbasis=FEMbasis,
                             lambda=best_lambda2[1], solver.options = "mass_lumping"), times=30)

# lambda preconditioner
tC=microbenchmark(smooth.FEM(observations=data, FEMbasis=FEMbasis,
                             lambda=best_lambda3[1], solver.options = "lambda_preconditioner"), times=30)

# block preconditioner
tD=microbenchmark(smooth.FEM(observations=data, FEMbasis=FEMbasis,
                             lambda=best_lambda4[1], solver.options = "block_preconditioner"), times=30)

# comparison
boxplot(log(tA$time),log(tB$time),log(tC$time),log(tD$time), names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')