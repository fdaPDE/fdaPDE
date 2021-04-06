##########################################
############## TEST SCRIPT ###############
##########################################

library(fdaPDE)

####### 2D ########
#### Test 2: c-shaped domain ####
#            locations != nodes
#            laplacian
#            with covariates
#            no BC
#            order FE = 1
rm(list=ls())
graphics.off()

data(horseshoe2D)

mesh = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
locations = refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
plot(mesh, pch = ".")
points(locations, pch = 16, col = "red")

FEMbasis = create.FEM.basis(mesh)

ndata = nrow(locations)

# Create covariates
set.seed(509875)
cov1 = rnorm(ndata, mean = 1, sd = 2)
cov2 = sin(locations[,1])

# Exact solution (pointwise at nodes)
DatiEsatti=fs.test(locations[,1], locations[,2]) + 2*cov1 - cov2

for(ind in 1:100)
{
  points(locations[which(round((DatiEsatti-min(DatiEsatti))/(max(DatiEsatti)-min(DatiEsatti))*100)==ind),1],
         locations[which(round((DatiEsatti-min(DatiEsatti))/(max(DatiEsatti)-min(DatiEsatti))*100)==ind),2],
         col=heat.colors(100)[ind], pch=16)
}

# Add error to simulate data
set.seed(543663)
ran = range(DatiEsatti)
data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^seq(-3,3,by=0.25)


#### Test 2.1.1: grid with stochastic GCV, one-at-the-time tests, Wald-type, exact computation
R_inference_object = inferenceDataObjectBuilder(test = "one-at-the-time", exact = "True", dim = 2)

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda, lambda.selection.lossfunction = "GCV",
                       R_Inference_Data_Object = R_inference_object
)
output_CPP$inference$p_vals

#### Test 2.1.2: grid with stochastic GCV, one-at-the-time tests, Speckman-type, exact computation
R_inference_object = inferenceDataObjectBuilder(test = "one-at-the-time", exact = "True", type = "speckman", dim = 2)

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2), lambda.selection.criterion = "grid", lambda.selection.lossfunction = "GCV",
                       FEMbasis=FEMbasis, lambda = lambda, R_Inference_Data_Object = R_inference_object
)
output_CPP$inference$p_vals

#### Test 2.1.3: grid with stochastic GCV, one-at-the-time tests, Permutational-type, exact computation
R_inference_object = inferenceDataObjectBuilder(test = "one-at-the-time", exact = "True", type = "permutational", dim = 2)

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2), lambda.selection.lossfunction = "GCV",
                       FEMbasis=FEMbasis, lambda=lambda, R_Inference_Data_Object = R_inference_object
)
output_CPP$inference$p_vals

#### Test 2.2.1: grid with stochastic GCV, simultaneous test, Wald-type, exact computation
R_inference_object = inferenceDataObjectBuilder(test = "simultaneous", exact = "True", dim = 2)

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2), lambda.selection.lossfunction = "GCV",
                       FEMbasis=FEMbasis, lambda=lambda, R_Inference_Data_Object = R_inference_object
)

output_CPP$inference$p_vals

#### Test 2.2.2: grid with stochastic GCV, simultaneous test, Speckman-type, exact computation
R_inference_object = inferenceDataObjectBuilder(test = "simultaneous", exact = "True", type = "speckman", dim = 2)

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2), lambda.selection.lossfunction = "GCV",
                       FEMbasis=FEMbasis, lambda=lambda, R_Inference_Data_Object = R_inference_object
)

output_CPP$inference$p_vals

#### Test 2.2.3: grid with stochastic GCV, simultaneous test, Permutational-type, exact computation
R_inference_object = inferenceDataObjectBuilder(test = "simultaneous", exact = "True", type = "permutational", dim = 2)

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2), lambda.selection.lossfunction = "GCV",
                       FEMbasis=FEMbasis, lambda=lambda, R_Inference_Data_Object = R_inference_object
)

output_CPP$inference$p_vals

#### Test 2.3.1: grid with stochastic GCV, one-at-the-time intervals, Wald-type, exact computation, default level
R_inference_object = inferenceDataObjectBuilder(interval = "one-at-the-time", exact = "True", dim = 2)

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2), lambda.selection.lossfunction = "GCV",
                       FEMbasis=FEMbasis, lambda=lambda, R_Inference_Data_Object = R_inference_object
)

output_CPP$inference$CI

#### Test 2.3.2: grid with stochastic GCV, one-at-the-time intervals, Speckman-type, exact computation, default level
R_inference_object = inferenceDataObjectBuilder(interval = "one-at-the-time", exact = "True", type = "speckman", dim = 2)

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2), lambda.selection.lossfunction = "GCV",
                       FEMbasis=FEMbasis, lambda=lambda, R_Inference_Data_Object = R_inference_object
)

output_CPP$inference$CI


#### Test 2.4.1: grid with stochastic GCV, simultaneous intervals, Wald-type, exact computation, custom level
R_inference_object = inferenceDataObjectBuilder(interval = "simultaneous", exact = "True", dim = 2, level = 0.1)

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2), lambda.selection.lossfunction = "GCV",
                       FEMbasis=FEMbasis, lambda=lambda, R_Inference_Data_Object = R_inference_object
)

output_CPP$inference$CI

#### Test 2.4.2: grid with stochastic GCV, simultaneous intervals, Speckman-type, exact computation, custom level
R_inference_object = inferenceDataObjectBuilder(interval = "simultaneous", exact = "True", type = "speckman" , dim = 2, level = 0.1)

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2), lambda.selection.lossfunction = "GCV",
                       FEMbasis=FEMbasis, lambda=lambda, R_Inference_Data_Object = R_inference_object
)

output_CPP$inference$CI

#### Test 2.5.1: grid with stochastic GCV, one-at-the-time tests, one-at-the-time intervals, Wald-type, exact computation, default level
R_inference_object = inferenceDataObjectBuilder(test = "one-at-the-time", interval = "one-at-the-time", exact = "True", dim = 2)

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda, R_Inference_Data_Object = R_inference_object
)

output_CPP$inference$p_vals
output_CPP$inference$CI

#### Test 2.5.2: grid with stochastic GCV, one-at-the-time tests, one-at-the-time intervals, Speckman-type, exact computation, default level
R_inference_object = inferenceDataObjectBuilder(test = "one-at-the-time", interval = "one-at-the-time", type = "speckman", exact = "True", dim = 2)

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda, R_Inference_Data_Object = R_inference_object
)

output_CPP$inference$p_vals
output_CPP$inference$CI

output_CPP$solution$beta
output_CPP$optimization$lambda_position
#### Test 2.6: Without GCV, one-at-the-time tests, one-at-the-time intervals, beta0, Wald-type, exact computation, default level
R_inference_object = inferenceDataObjectBuilder(test = "one-at-the-time", interval = "one-at-the-time", exact = "True", dim = 2, beta0 = c(1.9,-3.55))

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda, R_Inference_Data_Object = R_inference_object
)

output_CPP$inference$p_vals #we expect to accept H0 in this case
output_CPP$inference$CI

#### Test 2.7: Without GCV, simultaneous test, one-at-the-time intervals, linear combination, Wald-type, exact computation, default level
R_inference_object = inferenceDataObjectBuilder(test = "simultaneous", interval = "one-at-the-time", exact = "True", dim = 2, coeff = matrix(data = c(1,1,-1,1), nrow = 2, ncol = 2, byrow = T))

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda, R_Inference_Data_Object = R_inference_object
)

output_CPP$inference$p_vals 
output_CPP$inference$CI

#### Test 2.8: Without GCV, one-at-the-time tests, bonferroni intervals, linear combination, beta0, Wald-type, exact computation, default level
R_inference_object = inferenceDataObjectBuilder(test = "one-at-the-time", interval = "bonferroni", exact = "True", dim = 2, coeff = matrix(data = c(1,0,0,1), nrow = 2, ncol = 2, byrow = T), beta0 = c(1.9, -3.55))

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda, R_Inference_Data_Object = R_inference_object
)

output_CPP$inference$p_vals 
output_CPP$inference$CI
