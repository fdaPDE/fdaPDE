##########################################
############## TEST SCRIPT ###############
##########################################

##### Examples on how to perform the computation of
##### residual, RMSE and computational time
##### are reported only for the square 2D domain
##### both with and without covariates.
##### The extension to all the other cases is trivial

library(fdaPDE)

####### 2D ########

#### Test 1: square domain ####
#            locations = nodes 
#            laplacian
#            no covariates
#            no BC
#            order FE = 1
library(fdaPDE)
rm(list=ls())
graphics.off()

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

#### Test 2: square domain ####
#            locations = nodes 
#            laplacian
#            covariates
#            no BC
#            order FE = 1

# keep the same mesh and basis functions as before

# covariate definition
cov=cos(3*pi*nodes[,2])
beta=1.2
image(FEM(cov,FEMbasis))

data = sol_exact + beta*cov
ran = range(sol_exact)
data1 <- data + rnorm (nnodes, sd=0.001*(ran[2]-ran[1]))

### Computation of the residual
# no GCV

# no preconditioner
output_CPP <- smooth.FEM(observations=data,
                         FEMbasis=FEMbasis,
                         lambda=lambda, covariates = cov)
image (output_CPP$fit.FEM)
sol_approx = output_CPP$fit.FEM
femsol_approx <- eval.FEM(sol_approx, nodes)

output_CPP1 <- smooth.FEM(observations=data1,
                          FEMbasis=FEMbasis,
                          lambda=lambda, covariates = cov)
sol_approx1 = output_CPP1$fit.FEM
femsol_approx1 <- eval.FEM(sol_approx1, nodes)

residual <- matrix(nrow=dim(locations)[1], ncol=length(lambda))
residual <- abs(femsol_approx1 - femsol_approx)
res1 <- colMeans(residual)

# mass lumping
output_CPP <- smooth.FEM(observations=data,
                         FEMbasis=FEMbasis,
                         lambda=lambda, covariates = cov,
                         solver.options="mass_lumping")
image (output_CPP$fit.FEM)
sol_approx = output_CPP$fit.FEM
femsol_approx <- eval.FEM(sol_approx, nodes)

output_CPP1 <- smooth.FEM(observations=data1,
                          FEMbasis=FEMbasis,
                          lambda=lambda, covariates = cov,
                          solver.options="mass_lumping")
sol_approx1 = output_CPP1$fit.FEM
femsol_approx1 <- eval.FEM(sol_approx1, nodes)

residual <- matrix(nrow=dim(locations)[1], ncol=length(lambda))
residual <- abs(femsol_approx1 - femsol_approx)
res2 <- colMeans(residual)

# lambda preconditoner
output_CPP <- smooth.FEM(observations=data,
                         FEMbasis=FEMbasis,
                         lambda=lambda, covariates = cov,
                         solver.options="lambda_preconditioner")
image (output_CPP$fit.FEM)
sol_approx = output_CPP$fit.FEM
femsol_approx <- eval.FEM(sol_approx, nodes)

output_CPP1 <- smooth.FEM(observations=data1,
                          FEMbasis=FEMbasis,
                          lambda=lambda, covariates = cov,
                          solver.options="lambda_preconditioner")
sol_approx1 = output_CPP1$fit.FEM
femsol_approx1 <- eval.FEM(sol_approx1, nodes)

residual <- matrix(nrow=dim(locations)[1], ncol=length(lambda))
residual <- abs(femsol_approx1 - femsol_approx)
res3 <- colMeans(residual)

# block preconditoner
output_CPP <- smooth.FEM(observations=data,
                         FEMbasis=FEMbasis,
                         lambda=lambda, covariates = cov,
                         solver.options="block_preconditioner")
image (output_CPP$fit.FEM)
sol_approx = output_CPP$fit.FEM
femsol_approx <- eval.FEM(sol_approx, nodes)

output_CPP1 <- smooth.FEM(observations=data1,
                          FEMbasis=FEMbasis,
                          lambda=lambda, covariates = cov,
                          solver.options="block_preconditioner")
sol_approx1 = output_CPP1$fit.FEM
femsol_approx1 <- eval.FEM(sol_approx1, nodes)

residual <- matrix(nrow=dim(locations)[1], ncol=length(lambda))
residual <- abs(femsol_approx1 - femsol_approx)
res4 <- colMeans(residual)

### Best lambda corresponding to exact data
### Newton exact method with exact GCV
### default initial lambda and tolerance
# no preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV')
best_lambda1 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# mass lumping
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "mass_lumping")
best_lambda2 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# lambda preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "lambda_preconditoner")
best_lambda3 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# block preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates = cov,
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
                           FEMbasis=FEMbasis,  covariates = cov,
                           lambda.selection.criterion='newton',
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
                           FEMbasis=FEMbasis,  covariates = cov,
                           lambda.selection.criterion='newton',
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
                           FEMbasis=FEMbasis,  covariates = cov,
                           lambda.selection.criterion='newton',
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
                           FEMbasis=FEMbasis,  covariates = cov,
                           lambda.selection.criterion='newton',
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
    output_CPP <- smooth.FEM(observations = data+rnorm(length(nodes[1,]),
                                                        sd=0.1*(ran[2]=ran[1])),
                             FEMbasis=FEMbasis,
                             lambda=l, covariates = cov)
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
    output_CPP <- smooth.FEM(observations = data+rnorm(length(nodes[1,]),
                                                        sd=0.1*(ran[2]=ran[1])),
                             FEMbasis=FEMbasis,
                             lambda=l, covariates = cov,
                             solver.options = "mass_lumping" )
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
                             lambda=l, covariates = cov,
                             solver.options = "lambda_preconditioner" )
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
                             lambda=l, covariates = cov,
                             solver.options = "block_preconditioner" )
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
tA=microbenchmark(smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates=cov,
                             lambda=best_lambda1[1]), times=30)

# mass lumping
tB=microbenchmark(smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates=cov,
                             lambda=best_lambda2[1], solver.options = "mass_lumping"), times=30)

# lambda preconditioner
tC=microbenchmark(smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates=cov,
                             lambda=best_lambda3[1], solver.options = "lambda_preconditioner"), times=30)

# block preconditioner
tD=microbenchmark(smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates=cov,
                             lambda=best_lambda4[1], solver.options = "block_preconditioner"), times=30)

# comparison
boxplot(log(tA$time),log(tB$time),log(tC$time),log(tD$time), names=c('fdaPDE','mass lumping','lambda','block'),col=c('grey',2,3,4), ylab='log(time)')


#### Test 3: c-shaped domain ####
#            locations != nodes
#            laplacian
#            with covariates
#            no BC
#            order FE = 1
library(fdaPDE)
rm(list=ls())
graphics.off()

data(horseshoe2D)
boundary_nodes = horseshoe2D$boundary_nodes
boundary_segments = horseshoe2D$boundary_segments
locations = horseshoe2D$locations

mesh <- create.mesh.2D(nodes = rbind (boundary_nodes, locations), segments = boundary_segments)
nodes <- mesh$nodes
nnodes=length(nodes[,1])
locations = refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
FEMbasis <- create.FEM.basis(mesh)

coeff <- fs.test (nodes[,1],nodes[,2])

femfun = FEM(coeff, FEMbasis)
plot(femfun)
image(femfun)

lambda = 10^(-12:3)
sol_exact <- fs.test(locations[,1],locations[,2])
ran=range(sol_exact)

data = sol_exact
data1 = sol_exact + rnorm(nnodes, mean=0,
                          sd=0.001* abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^( -12:3)

### Best lambda corresponding to exact data
### Newton exact method with exact GCV
### default initial lambda and tolerance
# no preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, locations=locations,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV')
best_lambda1 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# mass lumping
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, locations=locations,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "mass_lumping")
best_lambda2 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# lambda preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, locations=locations,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "lambda_preconditoner")
best_lambda3 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# block preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, locations=locations,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "block_preconditoner")
best_lambda4 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

#### Test 4: c-shaped domain ####
#            locations != nodes
#            laplacian
#            with covariates
#            no BC
#            order FE = 1

# keep the same mesh and basis functions as before

ndata=length(locations[,1])
cov1 = rnorm(ndata, mean = 1, sd = 2)
cov2 = sin(nodes[,1])
coeff = fs.test(nodes[,1], nodes[,2])

cov1_loc = rnorm(ndata, mean = 1, sd = 2)
cov2_loc = sin(locations[,1])
cov = cbind(cov1_loc,cov2_loc)
sol_exact = fs.test(locations[,1], locations[,2])

data <- sol_exact +2*cov1_loc -cov2_loc
ran = range(sol_exact)
data1 <- data + rnorm (ndata, sd=0.001*(ran[2]-ran[1]))

### Best lambda corresponding to exact data
### Newton exact method with exact GCV
### default initial lambda and tolerance
# no preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, locations=locations,
                         covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV')
best_lambda1 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# mass lumping
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, locations=locations,
                         covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "mass_lumping")
best_lambda2 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# lambda preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, locations=locations,
                         covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "lambda_preconditoner")
best_lambda3 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# block preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, locations=locations,
                         covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "block_preconditoner")
best_lambda4 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

#### Test 5: quasicircular domain ####
#            areal observations
#            PDE space varying
#            no covariates
#            with BC
#            order FE = 1
library(fdaPDE)
rm(list=ls())
graphics.off()

data(quasicircle2Dareal)
incidence_matrix = quasicircle2Dareal$incidence_matrix
data_areal = quasicircle2Dareal$data
ran_areal = range(data_areal)
data1_areal <- data_areal + rnorm (length(data_areal), sd=0.05*(ran_areal[2]-ran_areal[1]))

mesh_areal = quasicircle2Dareal$mesh
nodes_areal = mesh_areal$nodes

fembasis_areal = create.FEM.basis(mesh_areal)
lambda = 10^(-12:3)

# PDE parameters
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

### Best lambda corresponding to exact data
### Newton exact method with exact GCV
### default initial lambda and tolerance
# no preconditioner
sol_areal = smooth.FEM(observations = data_areal, incidence_matrix = incidence_matrix, FEMbasis = fembasis_areal,
                       lambda.selection.criterion = 'newton_fd', DOF.evaluation = 'exact',
                       PDE_parameters = PDE_parameters_areal, BC = BC_areal)
sol_approx_areal <- sol_areal$fit.FEM
best_lambda_areal1 <- sol_areal$optimization$lambda_solution

# mass lumping
sol_areal = smooth.FEM(observations = data_areal, incidence_matrix = incidence_matrix, FEMbasis = fembasis_areal,
                       lambda.selection.criterion = 'newton_fd', DOF.evaluation = 'exact',
                       solver.options="mass_lumping",
                       PDE_parameters = PDE_parameters_areal, BC = BC_areal)
sol_approx_areal <- sol_areal$fit.FEM
best_lambda_areal2 <- sol_areal$optimization$lambda_solution

# lambda preconditioner
sol_areal = smooth.FEM(observations = data_areal, incidence_matrix = incidence_matrix, FEMbasis = fembasis_areal,
                       lambda.selection.criterion = 'newton_fd', DOF.evaluation = 'exact',
                       solver.options="lambda_preconditioner",
                       PDE_parameters = PDE_parameters_areal, BC = BC_areal)
sol_approx_areal <- sol_areal$fit.FEM
best_lambda_areal3 <- sol_areal$optimization$lambda_solution

# block preconditioner
sol_areal = smooth.FEM(observations = data_areal, incidence_matrix = incidence_matrix, FEMbasis = fembasis_areal,
                       lambda.selection.criterion = 'newton_fd', DOF.evaluation = 'exact',
                       solver.options="block_preconditioner",
                       PDE_parameters = PDE_parameters_areal, BC = BC_areal)
sol_approx_areal <- sol_areal$fit.FEM
best_lambda_areal4 <- sol_areal$optimization$lambda_solution

####### 2.5D ########

#### Test 1: hub domain ####
#            locations = nodes 
#            no covariates
#            no BC
#            order FE = 1

library(fdaPDE)
rm(list=ls())

data(hub2.5D)
mesh <- create.mesh.2.5D(nodes = hub2.5D$hub2.5D.nodes,triangles = hub2.5D$hub2.5D.triangles)
FEMbasis <- create.FEM.basis(mesh)

# Locations at nodes
nodes=mesh$nodes

# Exact data - Locations at nodes
nnodes = nrow(mesh$nodes)
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

func = function(x)
{
  a1*sin(2*pi*x[,2])+a2*sin(2*pi*x[,3])
}

coeff = func(nodes)
plot(FEM(coeff,FEMbasis))

lambda=10^(-12:3)
ran=range(coeff)

data=coeff
data1=data+rnorm(length(data)[1],sd=0.001*(ran[2]-ran[1]))

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

#### Test 2: hub domain ####
#            locations = nodes 
#            covariates
#            no BC
#            order FE = 1

# keep the same mesh and basis functions as before

cov1 = 4*sin(2*pi*nodes[,2])*cos(2*pi*nodes[,3])
cov2 = rnorm(nodes[,1],mean=3,sd=0.1)
cov = cbind(cov1,cov2)
plot(FEM(cov1,FEMbasis))
plot(FEM(cov2,FEMbasis))
beta=c(0.45,0.3)

data <- coeff +beta[1]*cov1+beta[2]*cov2
data1 <- data + rnorm (nnodes, sd=0.001*(ran[2]-ran[1]))

### Best lambda corresponding to exact data
### Newton exact method with exact GCV
### default initial lambda and tolerance
# no preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV')
best_lambda1 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# mass lumping
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "mass_lumping")
best_lambda2 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# lambda preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "lambda_preconditoner")
best_lambda3 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# block preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "block_preconditoner")
best_lambda4 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )


####### 3D ########

#### Test 1: sphere domain ####
#            locations = nodes 
#            covariates
#            no BC
#            order FE = 1

library(fdaPDE)
rm(list=ls())
graphics.off()

# # Function to generate random points in a sphere
# rsphere <- function(n, r = 0.9, surface_only = FALSE) {
#   phi       <- runif(n, 0.0, 2.0 * pi)
#   cos_theta <- runif(n, -1.0, 1.0)
#   sin_theta <- sqrt((1.0-cos_theta)*(1.0+cos_theta))
#   radius <- r
#   if (surface_only == FALSE) {
#     radius <- r * runif(n, 0.0, 1.0)^(1.0/3.0)
#   }
#   
#   x <- radius * sin_theta * cos(phi)
#   y <- radius * sin_theta * sin(phi)
#   z <- radius * cos_theta
#   
#   cbind(x, y, z)
# }

# Build mesh: Sphere
data("sphere3Ddata")
mesh<-create.mesh.3D(sphere3Ddata$nodes,
                     sphere3Ddata$tetrahedrons,
                     order=1)

FEMbasis <- create.FEM.basis(mesh)

set.seed(5847947)

a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)


nodes=mesh$nodes
nnodes = nrow(mesh$nodes)

# Set smoothing parameter
lambda=10^(-12:3)

# Evaluate exact solution on mesh nodes
coeff =  sin(2*pi*mesh$nodes[,1]) +  2 * sin(2*pi*mesh$nodes[,2]) +  sin(2*pi*mesh$nodes[,3])

# Plot exact solution
plot(FEM(coeff,FEMbasis))

data = coeff
data1 = data + rnorm(nrow(mesh$nodes), mean=0, sd=0.001*diff(range(coeff)))


func = function(x)
{
  a1*sin(2*pi*(x[,1]*x[,2]))+a2*cos(2*pi*x[,2])+a3*cos(2*pi*x[,3])
}

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

#### Test 2: sphere domain ####
#            locations = nodes 
#            covariates
#            no BC
#            order FE = 1

# keep the same mesh and basis functions as before

cov1=(4*sin(2*pi*nodes[,2])+6*sin((2*pi*nodes[,3])^2))*(1-exp(-nodes[,1]))/3
cov2=1+2*nodes[,1]*sin(2*pi*nodes[,2])/6
cov=cbind(cov1,cov2)
beta=c(0.7,0.2)

### Best lambda corresponding to exact data
### Newton exact method with exact GCV
### default initial lambda and tolerance
# no preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV')
best_lambda1 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# mass lumping
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "mass_lumping")
best_lambda2 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# lambda preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "lambda_preconditoner")
best_lambda3 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )

# block preconditioner
output_CPP <- smooth.FEM(observations=data, FEMbasis=FEMbasis, covariates = cov,
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV', solver.options = "block_preconditoner")
best_lambda4 <- output_CPP$optimization$lambda_solution
image (FEM( output_CPP$fit.FEM$coeff , FEMbasis ) )


####### GLM ########

library(fdaPDE)
library(purrr)

#### Test 1: square domain BINOMIAL family ####
#            locations = nodes 
#            laplacian
#            no covariates
#            no BC
#            order FE = 1
#

rm(list=ls())
graphics.off()

library(fdaPDE)
library(purrr)

FAMILY1 = "binomial"

logit <- function(x){qlogis(x)}
inv.logit <- function(x){plogis(x)}
link = logit
inv.link = inv.logit

# mesh
x = seq(0,1, length.out = 40)
y = x
locations = expand.grid(x,y)

mesh1 = create.mesh.2D(locations)
plot(mesh1)

nnodes1 = dim(mesh1$nodes)[1]

FEMbasis = create.FEM.basis(mesh1)

set.seed(42)
# locations
nloc = 800
xobs=runif(min=0,max=1,n=nloc)
yobs=runif(min=0,max=1,n=nloc)
loc=cbind(xobs,yobs)

# 2D random field (function f) 
a1=-2.5
a2=0.8

z<-function(p)
{
  a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1])
}

# exact solution
sol_exact1=rep(0,length(loc[,1]))
for(i in 1:length(loc[,1])){
  sol_exact1[i] <- z(loc[i,])
}

nnodes1 = dim(mesh1$nodes)[1]
sol_nodes1 = numeric(nnodes1)
for(i in 1:nnodes1){
  sol_nodes1[i] = z(mesh1$nodes[i,])
}

range(sol_exact1) 
param1 = sol_exact1
mu1<-inv.link(param1)
range(mu1)
# sampling response:
response1 <- rbernoulli(length(loc[,1]),p = mu1)

# Set smoothing parameter
lambda = 10^seq(-5,0,length.out = 20)

### Best lambda corresponding to exact data
### Newton exact method with exact GCV
### default initial lambda and tolerance
# no preconditioner
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response1), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY1, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda1 <- output_CPP$bestlambda

# mass lumping
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response1), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY1, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, solver.options="mass_lumping",
                                 lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda2 <- output_CPP$bestlambda

# lambda preconditioner
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response1), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY1, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, solver.options="lambda_preconditioner",
                                 lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda3 <- output_CPP$bestlambda

# block preconditioner
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response1), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY1, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, solver.options="block_preconditioner",
                                 lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda4 <- output_CPP$bestlambda

#### Test 2: square domain EXPONENTIAL family ####
#            locations = nodes 
#            laplacian
#            no covariates
#            no BC
#            order FE = 1
#

rm(list=ls())
graphics.off()

library(fdaPDE)
library(purrr)

FAMILY2 = "exponential"

link<-function(x){-1/x}
inv.link<-link 

# 2D random field (function f)
a1=-1.5
a2=0.4

z<-function(p)
{
  a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1]) - 2
  
}

# mesh
x = seq(0,1, length.out = 40)
y = x
locations = expand.grid(x,y)

mesh1 = create.mesh.2D(locations)
plot(mesh1)

nnodes1 = dim(mesh1$nodes)[1]

FEMbasis = create.FEM.basis(mesh1)

set.seed(42)
# locations
nloc = 800
xobs=runif(min=0,max=1,n=nloc)
yobs=runif(min=0,max=1,n=nloc)
loc=cbind(xobs,yobs)

# exact solution
sol_exact2=rep(0,length(loc[,1]))
for(i in 1:length(loc[,1])){
  sol_exact2[i] <- z(loc[i,])
}

sol_nodes2 = numeric(nnodes1)
for(i in 1:nnodes1){
  sol_nodes2[i] = z(mesh1$nodes[i,])
}

range(sol_exact2) 
param2 = sol_exact2
mu2<-inv.link(param2)
range(mu2)
# sampling response:
set.seed(95)
response2 <- response2 <- rexp(nloc, rate = 1/mu2)

# Set smoothing parameter
lambda = 10^seq(-5,0,length.out = 20)

### Best lambda corresponding to exact data
### Newton exact method with exact GCV
### default initial lambda and tolerance
# no preconditioner
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response2), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY2, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda1 <- output_CPP$bestlambda

# mass lumping
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response2), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY2, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, solver.options="mass_lumping",
                                 lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda2 <- output_CPP$bestlambda

# lambda preconditioner
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response2), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY2, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, solver.options="lambda_preconditioner",
                                 lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda3 <- output_CPP$bestlambda

# block preconditioner
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response2), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY2, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, solver.options="block_preconditioner",
                                 lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda4 <- output_CPP$bestlambda

#### Test 3: square domain GAMMA family ####
#            locations = nodes 
#            laplacian
#            no covariates
#            no BC
#            order FE = 1
#

rm(list=ls())
graphics.off()

library(fdaPDE)
library(purrr)

FAMILY3 = "gamma"

link<-function(x){-1/x}
inv.link<-link


# 2D random field (function f) 
a1=-1.5
a2=0.4

z<-function(p)
{
  a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1]) - 2
  
}

# mesh
x = seq(0,1, length.out = 40)
y = x
locations = expand.grid(x,y)

mesh1 = create.mesh.2D(locations)
plot(mesh1)

nnodes1 = dim(mesh1$nodes)[1]

FEMbasis = create.FEM.basis(mesh1)

set.seed(42)
# locations
nloc = 800
xobs=runif(min=0,max=1,n=nloc)
yobs=runif(min=0,max=1,n=nloc)
loc=cbind(xobs,yobs)


# exact solution
sol_exact3=rep(0,length(loc[,1]))
for(i in 1:length(loc[,1])){
  sol_exact3[i] <- z(loc[i,])
}

sol_nodes3 = numeric(nnodes1)
for(i in 1:nnodes1){
  sol_nodes3[i] = z(mesh1$nodes[i,])
}

range(sol_exact3) 
param3 = sol_exact3
mu3<-inv.link(param3)
range(mu3)
# sampling response:
set.seed(95)
scale.param3=1

response3 <- rgamma(length(loc[,1]), shape=mu3/scale.param3, scale=scale.param3)

# Set smoothing parameter
lambda = 10^seq(-5,0,length.out = 20)

### Best lambda corresponding to exact data
### Newton exact method with exact GCV
### default initial lambda and tolerance
# no preconditioner
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response3), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY3, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda1 <- output_CPP$bestlambda

# mass lumping
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response3), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY3, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, solver.options="mass_lumping",
                                 lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda2 <- output_CPP$bestlambda

# lambda preconditioner
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response3), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY3, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, solver.options="lambda_preconditioner",
                                 lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda3 <- output_CPP$bestlambda

# block preconditioner
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response3), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY3, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, solver.options="block_preconditioner",
                                 lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda4 <- output_CPP$bestlambda

#### Test 4: square domain POISSON family ####
#            locations = nodes 
#            laplacian
#            no covariates
#            no BC
#            order FE = 1
#
rm(list=ls())
graphics.off()

library(fdaPDE)
library(purrr)

FAMILY4 = "poisson"

l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv


# 2D random field (function f) 
a1=-1.5
a2=0.4

z<-function(p)
{
  a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1]) + 2
  
}

# mesh
x = seq(0,1, length.out = 40)
y = x
locations = expand.grid(x,y)

mesh1 = create.mesh.2D(locations)
plot(mesh1)

nnodes1 = dim(mesh1$nodes)[1]

FEMbasis = create.FEM.basis(mesh1)

set.seed(42)
# locations
nloc = 800
xobs=runif(min=0,max=1,n=nloc)
yobs=runif(min=0,max=1,n=nloc)
loc=cbind(xobs,yobs)


# exact solution
sol_exact4=rep(0,length(loc[,1]))
for(i in 1:length(loc[,1])){
  sol_exact4[i] <- z(loc[i,])
}


sol_nodes4 = numeric(nnodes1)
for(i in 1:nnodes1){
  sol_nodes4[i] = z(mesh1$nodes[i,])
}


range(sol_exact4) 
param4 = sol_exact4
mu4<-inv.link(param4)
range(mu4)
# sampling response:
set.seed(95)

response4 <- rpois(length(loc[,1]), lambda = mu4)

# Set smoothing parameter
lambda = 10^seq(-5,0,length.out = 20)

### Best lambda corresponding to exact data
### Newton exact method with exact GCV
### default initial lambda and tolerance
# no preconditioner
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response4), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY4, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda1 <- output_CPP$bestlambda

# mass lumping
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response4), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY4, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, solver.options="mass_lumping",
                                 lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda2 <- output_CPP$bestlambda

# lambda preconditioner
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response4), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY4, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, solver.options="lambda_preconditioner",
                                 lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda3 <- output_CPP$bestlambda

# block preconditioner
output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response4), FEMbasis =FEMbasis, covariates = NULL,
                                 max.steps=15, fam=FAMILY4, mu0=NULL, scale.param=NULL,
                                 lambda = lambda, solver.options="block_preconditioner",
                                 lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
best_lambda4 <- output_CPP$bestlambda