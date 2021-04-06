##########################################
############## TEST SCRIPT ###############
##########################################

library(fdaPDE)

####### 2.5D ########

#### Test 2: sphere domain ####
#            locations != nodes 
#            with covariates
#            no BC
#            order FE = 1
rm(list=ls())
graphics.off()

data(sphere2.5D)
mesh = sphere2.5D

FEMbasis=create.FEM.basis(mesh)

# Test function 
f = function(x, y, z){
  phi = (asin(y/sqrt(x^2+y^2)))
  theta = acos(z/sqrt(x^2+y^2+z^2))
  # rho = 1
  
  sin(4*(1/2*sin(theta)*exp(-sin(theta)^2)+1)*theta)*cos(2*(1/2*cos(phi)*exp(-cos(phi)^2)+1)*phi)
}

# Exact solution (pointwise at nodes)
sol_exact=f(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3])
plot(FEM(sol_exact, FEMbasis))

# Generate data locations on the sphere
set.seed(598944)
ndata = 500
locations = matrix(rnorm(ndata*3), ncol = 3)
locations = locations/sqrt(locations[,1]^2 + locations[,2]^2 + locations[,3]^2)

# Generate covariate and data
cov1 = runif(ndata, min = -1, max = 1)
DatiEsatti = f(locations[,1],locations[,2], locations[,3]) + cov1

# Add error to simulate data
set.seed(7893475)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(ndata, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Project locations on the mesh
projected_locations = projection.points.2.5D(mesh, locations)

# Set smoothing parameter
lambda = 10^seq(-4,-2,by=0.25)

#### Test 2.1.1: Without GCV, only first lambda, one-at-the-time tests, Wald-type, exact computation
R_Inference_Object <- inferenceDataObjectBuilder(test = "one-at-the-time", exact="True", dim=1)

output_CPP<-smooth.FEM(observations=data, 
                       locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, 
                       lambda=lambda[1], 
                       R_Inference_Data_Object = R_Inference_Object)
plot(output_CPP$fit.FEM)
output_CPP$solution$beta
output_CPP$inference$p_vals

#### Test 2.1.1: Without GCV, only first lambda, one-at-the-time test, Wald-type, exact computation
R_Inference_Object <- inferenceDataObjectBuilder(test = "one-at-the-time", exact="True", dim=1)

output_CPP<-smooth.FEM(observations=data, 
                       locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, 
                       lambda=lambda[1], 
                       R_Inference_Data_Object = R_Inference_Object)
plot(output_CPP$fit.FEM)
output_CPP$solution$beta
output_CPP$inference$p_vals


#### Test 2.1.2: Without GCV, only first lambda, one-at-the-time test, Speckman-type, exact computation
R_Inference_Object <- inferenceDataObjectBuilder(test = "one-at-the-time", type = "speckman", exact="True", dim=1)

output_CPP<-smooth.FEM(observations=data, 
                       locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, 
                       lambda=lambda[1], 
                       R_Inference_Data_Object = R_Inference_Object)
plot(output_CPP$fit.FEM)
output_CPP$solution$beta
output_CPP$inference$p_vals

#### Test 2.1.3: Without GCV, only first lambda, one-at-the-time test, Permutational-type, exact computation
R_Inference_Object <- inferenceDataObjectBuilder(test = "one-at-the-time", type = "permutational", exact="True", dim=1)

output_CPP<-smooth.FEM(observations=data, 
                       locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, 
                       lambda=lambda[1], 
                       R_Inference_Data_Object = R_Inference_Object)
plot(output_CPP$fit.FEM)
output_CPP$solution$beta
output_CPP$inference$p_vals

#### Test 2.2.1: grid with exact GCV, one-at-the-time interval, Wald-type, exact computation
R_Inference_Object <- inferenceDataObjectBuilder(interval = "one-at-the-time", exact="True", dim=1)

output_CPP<-smooth.FEM(observations=data, locations = projected_locations,
                       covariates = cov1,
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', 
                       R_Inference_Data_Object = R_Inference_Object)
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta
output_CPP$inference$CI

#### Test 2.2.2: grid with exact GCV, one-at-the-time interval, Speckman-type, exact computation
R_Inference_Object <- inferenceDataObjectBuilder(interval = "one-at-the-time", type = "speckman", exact="True", dim=1)

output_CPP<-smooth.FEM(observations=data, locations = projected_locations,
                       covariates = cov1,
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', 
                       R_Inference_Data_Object = R_Inference_Object)
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta
output_CPP$inference$CI

#### Test 2.3.1: grid with stochastic GCV, one-at-the-time test, Wald-type, beta0, exact computation
R_Inference_Object <- inferenceDataObjectBuilder(test = "one-at-the-time", beta0 = 0.99, exact="True", dim=1)

output_CPP<-smooth.FEM(observations=data, locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV',
                       R_Inference_Data_Object = R_Inference_Object)
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta
output_CPP$inference$p_vals

#### Test 2.3.2: grid with stochastic GCV, one-at-the-time test, Speckman-type, beta0, exact computation
R_Inference_Object <- inferenceDataObjectBuilder(test = "one-at-the-time", type = "speckman", beta0 = 0.99, exact="True", dim=1)

output_CPP<-smooth.FEM(observations=data, locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV',
                       R_Inference_Data_Object = R_Inference_Object)
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta
output_CPP$inference$p_vals

#### Test 2.3.3: grid with stochastic GCV, one-at-the-time test, Permutational-type, beta0, exact computation
R_Inference_Object <- inferenceDataObjectBuilder(test = "one-at-the-time", type = "permutational", beta0 = 0.99, exact="True", dim=1)

output_CPP<-smooth.FEM(observations=data, locations = projected_locations, 
                       covariates = cov1,
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV',
                       R_Inference_Data_Object = R_Inference_Object)
plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta
output_CPP$inference$p_vals
