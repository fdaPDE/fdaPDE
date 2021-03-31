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
sd <- 0.25
#inference object for Wald test
R_Inference_Object_Wald <- inferenceDataObjectBuilder(test='one-at-the-time', exact='True', dim=1)
#inference object for Speckman test
R_Inference_Object_Speckman <- inferenceDataObjectBuilder(test='one-at-the-time', type = 'speckman', exact='True', dim=1)


# Real simulations: 1
beta_H1_list <- seq(from = 0, by = 0.02, length.out = 11)

covariates=S
covariates = scale(covariates)

res1_S <- list()

SEED=15061997
set.seed(SEED)
rand= rnorm(n, 0, sd = sd)

init = Sys.time()
for (i in 1:length(beta_H1_list)) {
  
  #set.seed(SEED)
  #rand= rnorm(n, 0, sd = sd)
  observations <- covariates * beta_H1_list[i] + diag(n) %*% c_true + rand # Build observations in H1
  
  out_CPP_Wald <- smooth.FEM(locations = locations, observations=observations, 
                        covariates = covariates,
                        FEMbasis=FEMbasis, lambda=lambda,
                        lambda.selection.criterion='grid',
                        R_Inference_Data_Object = R_Inference_Object_Wald)
  
  out_CPP_Speckman <- smooth.FEM(locations = locations, observations=observations, 
                             covariates = covariates,
                             FEMbasis=FEMbasis, lambda=lambda,
                             lambda.selection.criterion='grid',
                             R_Inference_Data_Object = R_Inference_Object_Speckman)
  
  res1_S[[i]] <- c(out_CPP_Wald$inference$p_vals, out_CPP_Speckman$p_vals)
  print(i)
}

# Real simulations: 2
beta_H1_list <- seq(from = 0, by = 0.02, length.out = 11)

covariates=S1
covariates = scale(covariates)

res1_S1 <- list()

for (i in 1:length(beta_H1_list)) {
  
  #set.seed(SEED)
  #rand= rnorm(n, 0, sd = sd)
  observations <- covariates * beta_H1_list[i] + diag(n) %*% c_true + rand # Build observations in H1
  
  oout_CPP_Wald <- smooth.FEM(locations = locations, observations=observations, 
                              covariates = covariates,
                              FEMbasis=FEMbasis, lambda=lambda,
                              lambda.selection.criterion='grid',
                              R_Inference_Data_Object = R_Inference_Object_Wald)
  
  out_CPP_Speckman <- smooth.FEM(locations = locations, observations=observations, 
                                 covariates = covariates,
                                 FEMbasis=FEMbasis, lambda=lambda,
                                 lambda.selection.criterion='grid',
                                 R_Inference_Data_Object = R_Inference_Object_Speckman)
  
  res1_S1[[i]] <- c(out_CPP_Wald$inference$p_vals, out_CPP_Speckman$p_vals)
  print(i)
}

# Real simulations: 3
beta_H1_list <- seq(from = 0, by = 0.02, length.out = 11)

covariates=S2
covariates = scale(covariates)

res1_S2 <- list()

for (i in 1:length(beta_H1_list)) {
  
  #set.seed(SEED)
  #rand= rnorm(n, 0, sd = sd)
  observations <- covariates * beta_H1_list[i] + diag(n) %*% c_true + rand # Build observations in H1
  
  out_CPP_Wald <- smooth.FEM(locations = locations, observations=observations, 
                             covariates = covariates,
                             FEMbasis=FEMbasis, lambda=lambda,
                             lambda.selection.criterion='grid',
                             R_Inference_Data_Object = R_Inference_Object_Wald)
  
  out_CPP_Speckman <- smooth.FEM(locations = locations, observations=observations, 
                                 covariates = covariates,
                                 FEMbasis=FEMbasis, lambda=lambda,
                                 lambda.selection.criterion='grid',
                                 R_Inference_Data_Object = R_Inference_Object_Speckman)
  
  res1_S2[[i]] <- c(out_CPP_Wald$inference$p_vals, out_CPP_Speckman$p_vals)
  print(i)
}


# Real simulations: 4
beta_H1_list <- seq(from = 0, by = 0.02, length.out = 11)

covariates=S3
covariates = scale(covariates)

res1_S3 <- list()

for (i in 1:length(beta_H1_list)) {
  #set.seed(SEED)
  #rand= rnorm(n, 0, sd = sd)
  observations <- covariates * beta_H1_list[i] + diag(n) %*% c_true + rand # Build observations in H1
  
  out_CPP_Wald <- smooth.FEM(locations = locations, observations=observations, 
                             covariates = covariates,
                             FEMbasis=FEMbasis, lambda=lambda,
                             lambda.selection.criterion='grid',
                             R_Inference_Data_Object = R_Inference_Object_Wald)
  
  out_CPP_Speckman <- smooth.FEM(locations = locations, observations=observations, 
                                 covariates = covariates,
                                 FEMbasis=FEMbasis, lambda=lambda,
                                 lambda.selection.criterion='grid',
                                 R_Inference_Data_Object = R_Inference_Object_Speckman)
  
  res1_S3[[i]] <- c(out_CPP_Wald$inference$p_vals, out_CPP_Speckman$p_vals)
  print(i)
}
end = Sys.time()
timing = end - init

plot_res1_S <-  matrix(data = NA, ncol=2, nrow=11)
plot_res1_S1 <- matrix(data = NA, ncol=2, nrow=11)
plot_res1_S2 <- matrix(data = NA, ncol=2, nrow=11)
plot_res1_S3 <- matrix(data = NA, ncol=2, nrow=11)

for (i in 1:11){
  plot_res1_S[i,]=c(1-res1_S[[i]][1], 1-res1_S[[i]][2])
  plot_res1_S1[i,]=c(1-res1_S1[[i]][1], 1-res1_S1[[i]][2])
  plot_res1_S2[i,]=c(1-res1_S2[[i]][1], 1-res1_S2[[i]][2])
  plot_res1_S3[i,]=c(1-res1_S3[[i]][1], 1-res1_S3[[i]][2])
}

x11()
par(mfrow=c(2,2))
plot(beta_H1_list, res1_S[,1], main="a. Gaussian random field", type='b', pch=15, col='red',xlab=expression(beta),ylab="Power",ylim=c(0,1))
lines(beta_H1_list, res1_S[,2], type = 'b', pch = 15, col = 'g')
grid()
abline(h=0.05, type='l', lty=2)


plot(beta_H1_list, res1_S1, main="b. Matern random field", type='b', pch=15, col='red',xlab=expression(beta),ylab="Power",ylim=c(0,1))
lines(beta_H1_list, res1_S1[,2], type = 'b', pch = 15, col = 'g')
grid()
abline(h=0.05, type='l', lty=2)

plot(beta_H1_list, res1_S2, main="c. Exponential random field", type='b', pch=15, col='red',xlab=expression(beta),ylab="Power",ylim=c(0,1))
lines(beta_H1_list, res1_S2[,2], type = 'b', pch = 15, col = 'g')
grid()
abline(h=0.05, type='l', lty=2)

plot(beta_H1_list, res1_S3, main="b. Matern random field", type='b', pch=15, col='red',xlab=expression(beta),ylab="Power",ylim=c(0,1))
lines(beta_H1_list, res1_S3[,2], type = 'b', pch = 15, col = 'g')
grid()
abline(h=0.05, type='l', lty=2)

graphics.off()




