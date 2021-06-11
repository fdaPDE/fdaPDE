####### 2D ########

#### Test 1: square domain ####
#            locations = nodes 
#            laplacian
#            no covariates
#            no BC
#            order FE = 1
#            separate penalizations

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
sol_exact=f(SpaceTimePoints[,1],SpaceTimePoints[,2],SpaceTimePoints[,3])
# plot(FEM.time(sol_exact, TimePoints, FEMbasis, FLAG_PARABOLIC = TRUE))

# Add error to simulate data
set.seed(7893475)
ran = range(sol_exact)
data = sol_exact
data1 = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))
observations=matrix(data,nrow(SpacePoints),NumTimePoints)
observations1=matrix(data1,nrow(SpacePoints),NumTimePoints)

# Set smoothing parameter
lambdaS = 10^-3
lambdaT = 10^(-12:3)

# #### Test 1.1: Without GCV
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations,
                            FEMbasis=FEMbasis, lambdaS=lambdaS,
                            lambdaT=lambdaT, FLAG_PARABOLIC = FALSE)
output_CPP1<-smooth.FEM.time(time_mesh = TimePoints, observations=observations1,
                            FEMbasis=FEMbasis, lambdaS=lambdaS,
                            lambdaT=lambdaT, FLAG_PARABOLIC = FALSE)


residual <- matrix(nrow=length(TimePoints)*dim(SpacePoints)[1], ncol= length(lambdaT))
for (i in 1:length(lambdaT))
{
  sol_approx = output_CPP$fit.FEM
  femsol_approx<-eval.FEM.time(sol_approx,locations = SpacePoints,time.instants = TimePoints,lambdaT=i)
  sol_approx1 = output_CPP1$fit.FEM
  femsol_approx1<-eval.FEM.time(sol_approx1,locations = SpacePoints,time.instants = TimePoints,lambdaT=i)
  residual[,i] <- abs((femsol_approx1)-(femsol_approx))
}
res5 <- colMeans(residual)
plot(-12:3, res1, pch=4, xlab='log_10(lambda)', ylab='residual', col=1,ylim=c(min(res3),max(res2)+0.003))
grid()
abline(v=-13+bestlambdas1[2], lty=3, col=1)
abline(v=-13+bestlambdas2[2], lty=2, col=2)
abline(v=-13+bestlambdas3[2], lty=2, col=3)
abline(v=-13+bestlambdas4[2], lty=4, col=4)
lines(-12:3, res1, col=1)
points(-12:3, res2, pch=4, col=2)
lines(-12:3, res2, col=2)
lines(-12:3, res3, col=3)
points(-12:3, res3, pch=4, col=3)
lines(-12:3, res4, col=4)
points(-12:3, res4, pch=4, col=4)
lines(-12:3, res5, col=5)
points(-12:3, res5, pch=4, col=5)
legend('topright',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner','block preconditioner','lambdaT'), cex=0.62,col=1:5, pch=19)


# #### Test 1.2: grid with exact GCV
# # We test one value for speeding the tests, this is not what should be done in practice!
lambdaS = 1e-3
lambdaT = 10^(-6:0)
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations,
                            FEMbasis=FEMbasis, lambdaS=lambdaS,lambdaT=lambdaT,
                            lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV',
                            FLAG_PARABOLIC = FALSE)
bestlambdas5 = which(output_CPP$GCV == min(output_CPP$GCV), arr.ind = TRUE)
# image(FEM.time(output_CPP$fit.FEM$coeff[,bestlambdas[1],bestlambdas[2]],FEMbasis = FEMbasis, time_mesh = TimePoints), t=1)

# RMSE evaluation on a fine grid
xeval=runif(1000,0,1)
yeval=runif(1000,0,1)
teval=runif(1000,0,2)
sol_exact=f(xeval,yeval,teval)
RMSE<-function(f,g) sqrt(mean((f-g)^2))
rmse =NULL
for (i in 1:20)
{
  output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=matrix(data+ rnorm(nnodes, mean=0, sd=0.01*abs(ran[2]-ran[1])),nrow(SpacePoints),NumTimePoints), 
                              FEMbasis=FEMbasis, lambdaS=lambdaS,
                              lambdaT=bestlambdas3[2], lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV',
                              FLAG_PARABOLIC = FALSE)
  bestlambdas5 <- c(bestlambdas5, which(output_CPP$GCV == min(output_CPP$GCV), arr.ind = TRUE))
  sol_eval=eval.FEM.time(output_CPP$fit.FEM.time,space.time.locations = cbind(teval,xeval,yeval))
  rmse <- c(rmse,RMSE(sol_eval,sol_exact))
}

rmseF <- rmse

#### Test 2: c-shaped domain ####
#            locations != nodes
#            laplacian
#            with covariates
#            no BC
#            order FE = 1
#            separate penalizations 
library(fdaPDE)

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
data1 = data + rnorm(length(sol_exact), mean = 0, sd =  0.05*diff(range(sol_exact)))
observations = matrix(data,nrow(locations),NumTimeInstants)

#### Test 2.2: exact GCV
lambdaS = 10^(-1:1)
lambdaT = 10^(-6:0)
output_CPP<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                            observations=observations, 
                            covariates = cov1,
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT, 
                            lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV',
                            FLAG_PARABOLIC = FALSE)

best_lambda1 <- which(output_CPP$GCV == min(output_CPP$GCV), arr.ind = FALSE)

teval=runif(1000,0,pi)
sol_exact=f(xeval,yeval,teval)
RMSE<-function(f,g) sqrt(mean((f-g)^2))
loc_eval = refine.mesh.2D(mesh, minimum_angle = 30, maximum_area = 0.005)$nodes
rmse1 =NULL
for (i in 1:20)
{
  output_CPP <- smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                    observations=matrix(data + rnorm(length(sol_exact), mean = 0, sd =  0.05*diff(range(sol_exact))),nrow(locations),NumTimeInstants), 
                    covariates = cov1,
                    FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT, 
                    lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV',
                    FLAG_PARABOLIC = FALSE)
  best_lambda1 <- c(best_lambda1, which(output_CPP$GCV == min(output_CPP$GCV), arr.ind = TRUE))
  sol_eval=eval.FEM.time(output_CPP$fit.FEM.time,space.time.locations = cbind(teval,loc_eval))
  rmse1 <- c(rmse1,RMSE(sol_eval,sol_exact))
}
boxplot(rmse1)

rmse1A <- rmse1

boxplot(rmse1A,rmse1B,rmse1C,rmse1D,rmse1E, names=c('fdaPDE','mass lumping','lambda','block','lambdaT'),
        col=c('grey',2,3,4,5),ylab='RMSE')

# Set smoothing parameter
lambdaT = best_lambda2[43]
lambdaS = 10^(-12:3)

# #### Test 1.1: Without GCV
output_CPP<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                            observations=matrix(data,nrow(locations),NumTimeInstants), 
                            covariates = cov1,
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                            FLAG_PARABOLIC = FALSE)
output_CPP1<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                             observations=matrix(data + rnorm(length(data), mean = 0, sd =  0.05*diff(range(sol_exact))),nrow(locations),NumTimeInstants), 
                             covariates = cov1,
                             FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                             FLAG_PARABOLIC = FALSE)

residual <- matrix(nrow=dim(space_time_locations)[1], ncol= length(lambdaS))
for (i in 1:length(lambdaS))
{
  sol_approx = output_CPP$fit.FEM
  femsol_approx<-eval.FEM.time(sol_approx,space.time.locations=space_time_locations,lambdaS=i)
  sol_approx1 = output_CPP1$fit.FEM
  femsol_approx1<-eval.FEM.time(sol_approx1,space.time.locations=space_time_locations,lambdaS=i)
  residual[,i] <- abs((femsol_approx1)-(femsol_approx))
}
reshs5S <- colMeans(residual)
plot(-12:3, reshs1, pch=4, xlab='log_10(lambda)', ylab='residual', col=1,ylim=c(min(reshs4),max(reshs1)+0.025))
grid()
abline(v=-13+best_lambda1[2], lty=3, col=1)
abline(v=-13+best_lambda2[2], lty=2, col=2)
abline(v=-13+best_lambda3[2], lty=4, col=3)
abline(v=-13+best_lambda1[2], lty=2, col=4)
lines(-12:3, reshs1, col=1)
points(-12:3, reshs2, pch=4, col=2)
lines(-12:3, reshs2, col=2)
lines(-12:3, reshs3, col=3)
points(-12:3, reshs3, pch=4, col=3)
lines(-12:3, reshs4, col=4)
points(-12:3, reshs4, pch=4, col=4)
lines(-12:3, reshs5, col=5)
points(-12:3, reshs5, pch=4, col=5)
legend('topright',legend=c('no preconditioner', 'mass lumping', 'lambda preconditioner','block preconditioner','lambdaT'), cex=0.62,col=1:5, pch=19)


# ####### 2.5D ########
# 
# #### hub pointwise (examples with and without covariates) ####
# library(fdaPDE)
# rm(list=ls())
# 
# data(hub2.5D)
# mesh <- create.mesh.2.5D(nodes = hub2.5D$hub2.5D.nodes,triangles = hub2.5D$hub2.5D.triangles)
# FEMbasis <- create.FEM.basis(mesh)
# 
# # Locations at nodes
# nodesLocations=mesh$nodes
# 
# # Exact data - Locations at nodes
# nnodes = mesh$nnodes
# a1 = rnorm(1,mean = 1, sd = 1)
# a2 = rnorm(1,mean = 1, sd = 1)
# a3 = rnorm(1,mean = 1, sd = 1)
# TimeNodes = 0:4
# 
# locations = cbind(rep(TimeNodes,each=nnodes),rep(nodesLocations[,1],length(TimeNodes)),rep(nodesLocations[,2],length(TimeNodes)),rep(nodesLocations[,3],length(TimeNodes)))
# 
# func = function(x)
# {
#   (a1*sin(2*pi*x[,2])+a2*sin(2*pi*x[,3])+a3*sin(2*pi*x[,4])+1)*cos(x[,1])
# }
# 
# func_evaluation = func(locations)
# # Plot the exact solution
# plot(FEM.time(coeff=array(func_evaluation,dim=c(length(func_evaluation),1,1)), FEMbasis = FEMbasis,time_mesh=TimeNodes,FLAG_PARABOLIC=T),TimeNodes)
# 
# 
# lambdaS=10^seq(-9, -7, 0.5)
# lambdaT=10^seq(-6, -4, 0.5)
# 
# lambdaS_par=10^seq(-4, -3, 0.25)
# lambdaT_par=10^seq(1, 1.8, 0.2)
# 
# cov1=4*sin(2*pi*locations[,2])*cos(2*pi*locations[,3])
# cov2=rnorm(nnodes*length(TimeNodes), mean=3, sd=0.1)*rep(exp(-TimeNodes/length(TimeNodes)),each=nnodes)
# W=cbind(cov1,cov2)
# 
# # plot(FEM(coeff = cov1[1:nnodes], FEMbasis = FEMbasis))
# # plot(FEM(coeff = cov2[1:nnodes], FEMbasis = FEMbasis))
# 
# # Fix betas
# beta_exact=c(0.45,0.3)
# 
# ran=range(W%*%beta_exact + func_evaluation)
# ran=range(func_evaluation)
# 
# # Plot exact solution
# plot(FEM.time(coeff=array(W%*%beta_exact + func_evaluation, dim = c(nnodes*length(TimeNodes),1,1)),FEMbasis=FEMbasis,time_mesh = TimeNodes,FLAG_PARABOLIC = T),TimeNodes)
# 
# ran = range(func_evaluation)
# data = func_evaluation +rnorm(nnodes,mean=0,sd=0.05*(ran[2]-ran[1]))
# 
# ran = range(func_evaluation+ W%*%beta_exact)
# datacov=func_evaluation+ W%*%beta_exact +rnorm(nnodes,mean=0,sd=0.05*(ran[2]-ran[1]))
# 
# data = matrix(data,mesh$nnodes,length(TimeNodes))
# datacov = matrix(datacov,mesh$nnodes,length(TimeNodes))
# 
# #########################################SEPARABLE####################################################
# solSep = smooth.FEM.time(observations=data,
#                          FEMbasis = FEMbasis, time_mesh = TimeNodes, time_locations = TimeNodes,
#                          lambdaS = lambdaS, lambdaT = lambdaT,
#                          FLAG_PARABOLIC = FALSE)
# 
# solSepCov = smooth.FEM.time(observations=datacov, covariates = W,
#                             FEMbasis = FEMbasis, time_mesh = TimeNodes,
#                             lambdaS = lambdaS, lambdaT = lambdaT,
#                             FLAG_PARABOLIC = FALSE)
# 
# 
# ##########################################PARABOLIC####################################################
# solPar = smooth.FEM.time(observations=data,
#                          FEMbasis = FEMbasis, time_mesh = TimeNodes, time_locations = TimeNodes,
#                          lambdaS = lambdaS, lambdaT = lambdaT,
#                          FLAG_PARABOLIC = TRUE)
# 
# solParCov = smooth.FEM.time(observations=datacov[,2:length(TimeNodes)], covariates = W[(1+mesh$nnodes):(length(TimeNodes)*mesh$nnodes),],
#                             FEMbasis = FEMbasis, time_mesh = TimeNodes,
#                             lambdaS = lambdaS, lambdaT = lambdaT,
#                             IC=func_evaluation[1:mesh$nnodes],
#                             FLAG_PARABOLIC = TRUE)
# 
# # Example of RMSE computation
# TimeNodesEval=seq(0,4,length.out = 9)
# eval_locations = cbind(rep(TimeNodesEval,each=nnodes),rep(nodesLocations[,1],length(TimeNodesEval)),rep(nodesLocations[,2],length(TimeNodesEval)),rep(nodesLocations[,3],length(TimeNodesEval)))
# sol_eval=eval.FEM.time(solSep$fit.FEM.time,locations = nodesLocations,time.instants = TimeNodesEval, lambdaS = solSep$bestlambda[1],lambdaT = solSep$bestlambda[2])
# sol_exact = func(eval_locations)
# RMSE<-function(f,g) sqrt(mean((f-g)^2))
# RMSE(sol_eval,sol_exact)
# 
# 
# #### hub areal (examples with and without covariates) ####
# library(fdaPDE)
# rm(list=ls())
# 
# data(hub2.5Dareal)
# 
# nodesLocations=mesh$nodes
# 
# nnodes = mesh$nnodes
# a1 = rnorm(1,mean = 1, sd = 1)
# a2 = rnorm(1,mean = 1, sd = 1)
# a3 = rnorm(1,mean = 1, sd = 1)
# TimeNodes = 0:4
# TimeNodesRMSE = seq(0,4,length.out = 15)
# 
# locations = cbind(rep(TimeNodesRMSE,each=nnodes),rep(nodesLocations[,1],length(TimeNodesRMSE)),rep(nodesLocations[,2],length(TimeNodesRMSE)),rep(nodesLocations[,3],length(TimeNodesRMSE)))
# 
# func = function(x)
# {
#   (a1*sin(2*pi*x[,2])+a2*sin(2*pi*x[,3])+a3*sin(2*pi*x[,4])+1)*cos(x[,1])
# }
# 
# func_evaluation = func(locations)
# FEMbasis=create.FEM.basis(mesh)
# 
# plot(FEM.time(coeff=array(func_evaluation,dim=c(length(func_evaluation),1,1)),time_mesh = TimeNodesRMSE,FEMbasis,FLAG_PARABOLIC = T),3)
# 
# sol_exact=func_evaluation
# 
# W_areal=cbind(rep(cov_areal,length(TimeNodes))*rep(exp(-TimeNodes/length(TimeNodes)),each=RDD_groups))
# 
# beta_exact=c(1)
# 
# lambdaS=10^seq(-9, -7, 0.5)
# lambdaT=10^seq(-6, -4, 0.5)
# 
# lambdaS_par=10^seq(-5.2, -4.8, 0.1)
# lambdaT_par=10^seq(1, 1.8, 0.2)
# 
# obs_areal = rep(obs_areal,length(TimeNodes))*rep(cos(TimeNodes),each=RDD_groups)
# 
# ran = range(obs_areal)
# data = obs_areal +rnorm(RDD_groups*length(TimeNodes),mean=0,sd=0.02*(ran[2]-ran[1]))
# 
# ran = range(obs_areal + W_areal%*%beta_exact)
# datacov=obs_areal + W_areal%*%beta_exact + rnorm(RDD_groups*length(TimeNodes),mean=0,sd=0.02*(ran[2]-ran[1]))
# 
# data = matrix(data,RDD_groups,length(TimeNodes))
# datacov = matrix(datacov,RDD_groups,length(TimeNodes))
# 
# ###########################SEPARABLE###########################################
# solSep = smooth.FEM.time(observations = data,incidence_matrix = incidence_matrix, time_locations=TimeNodes,
#                          time_mesh = TimeNodes,FEMbasis = FEMbasis, 
#                          lambdaS = lambdaS, lambdaT = lambdaT)
# 
# solSep = smooth.FEM.time(observations = datacov,time_mesh = TimeNodes, covariates = W_areal,incidence_matrix = incidence_matrix,
#                          FEMbasis = FEMbasis, lambdaS = lambdaS, lambdaT = lambdaT)
# 
# ##########################################PARABOLIC####################################################
# solPar = smooth.FEM.time(observations = data[,2:length(TimeNodes)],time_mesh = TimeNodes, incidence_matrix = incidence_matrix,
#                          FEMbasis = FEMbasis, lambdaS = lambdaS_par, lambdaT = lambdaT_par, FLAG_PARABOLIC = TRUE,
#                          IC=func_evaluation[1:mesh$nnodes])
# 
# solPar = smooth.FEM.time(observations = datacov[,2:length(TimeNodes)],time_mesh = TimeNodes, incidence_matrix = incidence_matrix,covariates = W_areal[(1+RDD_groups):(length(TimeNodes)*RDD_groups),],
#                          FEMbasis = FEMbasis, lambdaS = lambdaS_par, lambdaT = lambdaT_par, FLAG_PARABOLIC = TRUE, 
#                          IC=func_evaluation[1:mesh$nnodes])
# 
