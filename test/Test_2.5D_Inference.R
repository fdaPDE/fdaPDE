##########################################
############## TEST SCRIPT ###############
##########################################

### Installation of fdaPDE 
### Refer to Report_PACS_Cavazzutti_Galiberti for further details about installation
### This is just a quick recap

### If you are on Windows you need to install the right version of R tools before following the steps below
### Launch the commands by typing ctrl + l once selected the right line
### install.packages("rgl")
### install.packages("plot3D")
### install.packages("plot3Drgl")
### install.packages("geometry")
### install.packages("RcppEigen")
### download fdaPDE-dev or fdaPDE-no-fspai-dev from GitHub
### install.packages("path/to/fdaPDE", type='source', repos=NULL)

library(fdaPDE)

####### 2.5D ########
#### Test 1: sphere domain ####
#            locations != nodes 
#            with covariates
#            no BC
#            order FE = 1
rm(list=ls())
graphics.off()

# set the working directory to the data fold 
#setwd("~/Scrivania/project/fdaPDE/data")

load("sphere2.5D.RData")
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
set.seed(30101997)
ndata = 500
locations = matrix(rnorm(ndata*3), ncol = 3)
locations = locations/sqrt(locations[,1]^2 + locations[,2]^2 + locations[,3]^2)

# Generate covariate and data
cov1 = runif(ndata, min = -1, max = 1)
DatiEsatti = f(locations[,1],locations[,2], locations[,3]) + cov1

# Add error to simulate data
set.seed(9121997)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(ndata, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Project locations on the mesh
projected_locations = projection.points.2.5D(mesh, locations)

# Set smoothing parameter
lambda = 10^seq(-4,-2,by=0.25)

# Create inferenceDataObjects for the tests

inf_Obj_1 <- inferenceDataObjectBuilder(test = c('one-at-the-time', 'simultaneous', 'one-at-the-time'), interval = c('one-at-the-time','bonferroni','none'),
                                        type=c('wald', 'speckman','eigen-sign-flip'), exact='True', dim=1, beta0 = 1)

### The following example should give a warning in no-fspai-dev library, performing exact inference in turn of fspai. 
inf_Obj_2 <- inferenceDataObjectBuilder(test = c('one-at-the-time', 'one-at-the-time'), interval = c('simultaneous','one-at-the-time'),
                                        type=c('wald', 'speckman'), exact='False', dim=1, beta0 = 1, tol_fspai = 0.005)

inf_Obj_3 <- inferenceDataObjectBuilder(test = c('one-at-the-time', 'simultaneous', 'one-at-the-time'), interval = c('one-at-the-time','bonferroni','none'),
                                        type=c('wald', 'speckman','eigen-sign-flip'), exact='True', dim=1, beta0 = 0)

#### Test 1.1: Exact inference, accepting H0
output_CPP_1<-smooth.FEM(locations = projected_locations, observations=data, 
                         covariates = cov1,
                         FEMbasis=FEMbasis, lambda=lambda,
                         lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', R_Inference_Data_Object = inf_Obj_1)

# retrieving p_values, H0 is expected to be accepted, hence p_values are expected to be high (e.g. > 0.05) 
output_CPP_1$inference$p_values$wald
output_CPP_1$inference$p_values$speckman
output_CPP_1$inference$p_values$eigen_sign_flip

# retrieving confidence intervals, they are expected to contain the real values of beta, i.e. 1
output_CPP_1$inference$CI$wald
output_CPP_1$inference$CI$speckman

#### Test 1.2: Non-Exact inference, accepting H0 (it should give a warning in no-fspai-dev library, performing exact inference in turn of fspai)
output_CPP_2<-smooth.FEM(locations = projected_locations, observations=data, 
                         covariates = cov1,
                         FEMbasis=FEMbasis, lambda=lambda,
                         lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', R_Inference_Data_Object = inf_Obj_2)

# retrieving p_values, H0 is expected to be accepted, hence p_values are expected to be high (e.g. > 0.05) 
output_CPP_2$inference$p_values$wald
output_CPP_2$inference$p_values$speckman

# retrieving confidence intervals, they are expected to contain the real values of combinations of beta, i.e. 1
output_CPP_2$inference$CI$wald
output_CPP_2$inference$CI$speckman


#### Test 1.3: Exact inference, rejecting H0 
output_CPP_3<-smooth.FEM(locations = projected_locations, observations=data, 
                         covariates = cov1,
                         FEMbasis=FEMbasis, lambda=lambda,
                         lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', R_Inference_Data_Object = inf_Obj_3)

# retrieving p_values, H0 is expected to be rejected, hence p_values are expected to be low (e.g. < 0.05) 
output_CPP_3$inference$p_values$wald
output_CPP_3$inference$p_values$speckman

# retrieving confidence intervals, they are expected to contain the real values of beta, i.e. 1
output_CPP_3$inference$CI$wald
output_CPP_3$inference$CI$speckman
