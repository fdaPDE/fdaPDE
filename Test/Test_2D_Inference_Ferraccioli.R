library(MASS)
library(flip)
library(plyr)
library(fdaPDE)
library(Matrix)
library(ggplot2)
library(RandomFields)
library(reshape2)
library(viridis)
library(soobench)

# Mesh generation
bounds <- cbind(c(0:9, rep(10, 10), 10:1, rep(0, 10)), c(rep(0, 10), 0:9, rep(10, 10), 10:1))
mesh <- create.mesh.2D(nodes = bounds, order = 1)
plot(mesh)

# Mesh refinement
mesh_ref <- refine.mesh.2D(mesh, maximum_area = 0.3)
plot(mesh_ref)

FEMbasis <- create.FEM.basis(mesh_ref)

locations=mesh_ref$nodes
n <- nrow(mesh_ref$nodes) # Both n_observations and n_nodes since they coincide here

# Generation of the baseline function for data
spat_field <- function(x, y) {
  cos((2*x + y)/4) + ((x + y)/15)^2
}

plot(FEM(coeff = spat_field(mesh_ref$nodes[, 1], mesh_ref$nodes[, 2]), FEMbasis = FEMbasis))

c_true <- spat_field(mesh_ref$nodes[, 1], mesh_ref$nodes[, 2])
image(FEM(coeff = c_true, FEMbasis = FEMbasis))

# Covariates(random fields) definition
set.seed(1)
model <- RMgauss(scale = 0.3)
S <- RandomFields::RFsimulate(model, x = SpatialPoints(mesh_ref$nodes), n = 1)@data
S <- unlist(S)
plot(S)
image(FEM(coeff = S, FEMbasis = FEMbasis))

set.seed(4)
model <- RMexp(var = 1, scale = 2)
S3 <- RandomFields::RFsimulate(model, x = SpatialPoints(mesh_ref$nodes), n = 1)@data
S3 <- unlist(S3)
#plot(FEM(coeff = S3, FEMbasis = FEMbasis))
image(FEM(coeff = S3, FEMbasis = FEMbasis))
#rgl.snapshot("Sim_RDD/Simulations_mesh_centered/S3.png", fmt = "png")

set.seed(3)
model <- RMmatern(nu = 1, var = 8, scale = 0.5)
S2 <- RandomFields::RFsimulate(model, x = SpatialPoints(mesh_ref$nodes), n = 1)@data
S2 <- unlist(S2)
image(FEM(coeff = S2, FEMbasis = FEMbasis))

set.seed(2)
model <- RMmatern(nu = 5, var = 2, scale = 1)
S1 <- RandomFields::RFsimulate(model, x = SpatialPoints(mesh_ref$nodes), n = 1)@data
S1 <- unlist(S1)
image(FEM(coeff = S1, FEMbasis = FEMbasis))

#Setting Parameters for smoothing.R
lambda<-0.1
R_Inference_Object <- inferenceDataObjectBuilder(test='one-at-the-time', exact='True', dim=1)
sd <- 0.25

# Real simulations: 1
beta_H1_list <- seq(from = 0, by = 0.02, length.out = 11)

covariates=S
res1_S <- list()

for (i in 1:length(beta_H1_list)) {
  
  set.seed(35895389)
  rand= rnorm(n, 0, sd = sd)
  observations <- covariates * beta_H1_list[i] + diag(n) %*% c_true + rand # Build observations in H1
  
  out_CPP <- smooth.FEM(locations = locations, observations=observations, 
                        covariates = covariates,
                        FEMbasis=FEMbasis, lambda=lambda,
                        lambda.selection.criterion='grid',
                        R_Inference_Data_Object = R_Inference_Object)
  
  res1_S[[i]] <- out_CPP$inference$p_vals
  print(i)
}

# Real simulations: 2
beta_H1_list <- seq(from = 0, by = 0.02, length.out = 11)

covariates=S1
res1_S1 <- list()

for (i in 1:length(beta_H1_list)) {
  
  set.seed(35895389)
  rand= rnorm(n, 0, sd = sd)
  observations <- covariates * beta_H1_list[i] + diag(n) %*% c_true + rand # Build observations in H1
  
  out_CPP <- smooth.FEM(locations = locations, observations=observations, 
                        covariates = covariates,
                        FEMbasis=FEMbasis, lambda=lambda,
                        lambda.selection.criterion='grid',
                        R_Inference_Data_Object = R_Inference_Object)
  
  res1_S1[[i]] <- out_CPP$inference$p_vals
  print(i)
}

# Real simulations: 3
beta_H1_list <- seq(from = 0, by = 0.02, length.out = 11)

covariates=S2
res1_S2 <- list()

for (i in 1:length(beta_H1_list)) {
  
  set.seed(35895389)
  rand= rnorm(n, 0, sd = sd)
  observations <- covariates * beta_H1_list[i] + diag(n) %*% c_true + rand # Build observations in H1
  
  out_CPP <- smooth.FEM(locations = locations, observations=observations, 
                        covariates = covariates,
                        FEMbasis=FEMbasis, lambda=lambda,
                        lambda.selection.criterion='grid',
                        R_Inference_Data_Object = R_Inference_Object)
  
  res1_S2[[i]] <- out_CPP$inference$p_vals
  print(i)
}


# Real simulations: 4
beta_H1_list <- seq(from = 0, by = 0.02, length.out = 11)

covariates=S3
res1_S3 <- list()

for (i in 1:length(beta_H1_list)) {
  set.seed(35895389)
  rand= rnorm(n, 0, sd = sd)
  observations <- covariates * beta_H1_list[i] + diag(n) %*% c_true + rand # Build observations in H1
  
  out_CPP <- smooth.FEM(locations = locations, observations=observations, 
                        covariates = covariates,
                        FEMbasis=FEMbasis, lambda=lambda,
                        lambda.selection.criterion='grid',
                        R_Inference_Data_Object = R_Inference_Object)
  
  res1_S3[[i]] <- out_CPP$inference$p_vals
  print(i)
}

