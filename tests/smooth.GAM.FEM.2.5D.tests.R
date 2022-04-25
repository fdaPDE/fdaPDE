###############################################
############ TEST SCRIPT GAM 2.5D #############
###############################################

####### 2.5D ########
#### Test 1: HUB (2.5D) POISSON family ####
#            locations = nodes
#            with covariates
#            no BC
library(fdaPDE)
library(purrr)
rm(list=ls())
graphics.off()

# family
FAMILY = "poisson"

l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

# beta
beta1 = 2
beta2 = -3
betas_truth = c(beta1,beta2)

# lambda 
lambda = c(0.0001,0.001,0.005,0.006,0.007,0.008,0.009,0.0095,
           0.00975,0.01,0.0125,0.015,0.02,0.03,0.04,0.05,0.1,
           1,10,100,1000)

# scale param
scale.param = 1

# mesh
data(hub2.5D)
mesh <- create.mesh.2.5D(nodes = hub2.5D$hub2.5D.nodes,triangles = hub2.5D$hub2.5D.triangles)
FEMbasis <- create.FEM.basis(mesh)

# locations 
loc = mesh$nodes
nloc = dim(loc)[1]

# 2.5D random field (function f)
a1 = rnorm(1,mean = 0, sd = 1)
a2 = rnorm(1,mean = 0, sd = 1)
a3 = rnorm(1,mean = 0, sd = 1)

sol_exact = numeric(nloc)
for (i in 0:(nloc-1)){
  sol_exact[i+1] = a1* sin(2*pi*loc[i+1,1]) +  a2* sin(2*pi*loc[i+1,2]) +  a3*sin(2*pi*loc[i+1,3]) + 7
}
range(sol_exact) 

# covariates
set.seed(42)

desmat=matrix(0,nrow=nloc, ncol=2)
desmat[,1]= sin(2*pi*loc[,1])*cos(2*pi*loc[,2])
desmat[,2]= rnorm(nloc, mean=2, sd=0.1)

# samopling response
ran=range(desmat%*%betas_truth + sol_exact) 
param=sol_exact+beta1*desmat[,1]+beta2*desmat[,2]

mu<-inv.link(param)
range(mu)
response <- rpois(nloc, lambda = mu)

# Plot the True Field
plot(FEM(sol_exact, FEMbasis))

#### Test 1.1: Without GCV
output_CPP<- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                        max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                        lambda = lambda)

# Test 1.1.1: mass lumping
output_CPP<- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                        max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                        lambda = lambda, preconditioner='mass_lumping')

# Test 1.1.2: diagonal preconditoner
output_CPP<- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                        max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                        lambda = lambda, preconditioner='lambda_preconditioner')

# Test 1.1.3: block preconditoner
output_CPP<- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                        max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                        lambda = lambda, preconditioner='block_preconditioner')

#### Test 1.2: grid with exact GCV
output_CPP<- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                         max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                         lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$optimization$GCV_vector)

beta = output_CPP$solution$beta
func_estimation = output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position]
plot(FEM(func_estimation,FEMbasis))

# Test 1.2.1: mass lumping
output_CPP<- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                        max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                        lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV',
                        preconditioner='mass_lumping')
plot(log10(lambda),output_CPP$optimization$GCV_vector)

beta = output_CPP$solution$beta
func_estimation = output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position]
plot(FEM(func_estimation,FEMbasis))

# Test 1.2.2: diagonal preconditoner
output_CPP<- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                        max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                        lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV',
                        preconditioner='lambda_preconditioner')
plot(log10(lambda),output_CPP$optimization$GCV_vector)

beta = output_CPP$solution$beta
func_estimation = output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position]
plot(FEM(func_estimation,FEMbasis))


# Test 1.2.3: block preconditoner
output_CPP<- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                        max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                        lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV',
                        preconditioner='block_preconditioner')
plot(log10(lambda),output_CPP$optimization$GCV_vector)

beta = output_CPP$solution$beta
func_estimation = output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position]
plot(FEM(func_estimation,FEMbasis))


#### Test 1.3: grid with stochastic GCV
output_CPP <- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                         max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                         lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda),output_CPP$optimization$GCV_vector)

beta = output_CPP$solution$beta
func_estimation = output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position]
plot(FEM(func_estimation,FEMbasis))


# Test 1.3.1: mass lumping
output_CPP <- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                         max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                         lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV',
                         preconditioner = 'mass_lumping')
plot(log10(lambda),output_CPP$optimization$GCV_vector)

beta = output_CPP$solution$beta
func_estimation = output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position]
plot(FEM(func_estimation,FEMbasis))

# Test 1.3.2: diagonal preconditoner
output_CPP <- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                         max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                         lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV',
                         preconditioner = 'lambda_preconditioner')
plot(log10(lambda),output_CPP$optimization$GCV_vector)

beta = output_CPP$solution$beta
func_estimation = output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position]
plot(FEM(func_estimation,FEMbasis))


# Test 1.3.3: block preconditoner

output_CPP <- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                         max.steps=15, family=FAMILY, mu0=NULL, scale.param=NULL,
                         lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV',
                         preconditioner = 'block_preconditioner')
plot(log10(lambda),output_CPP$optimization$GCV_vector)

beta = output_CPP$solution$beta
func_estimation = output_CPP$fit.FEM$coeff[,output_CPP$optimization$lambda_position]
plot(FEM(func_estimation,FEMbasis))
