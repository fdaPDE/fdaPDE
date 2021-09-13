##########################################
############## TEST SCRIPT ###############
##########################################

library(fdaPDE)

####### 2D ########
#### Test 1: c-shaped domain ####
#            locations != nodes
#            laplacian
#            with covariates
#            no BC
#            order FE = 1
rm(list=ls())
graphics.off()

# set the working directory to the data fold 
#setwd("~/Scrivania/project/fdaPDE/data")

data(horseshoe2D)

mesh = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
locations = refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
plot(mesh, pch = ".")
points(locations, pch = 16, col = "red")

FEMbasis = create.FEM.basis(mesh)

ndata = nrow(locations)

c_true <- fs.test(mesh$nodes[, 1], mesh$nodes[, 2])
image(FEM(coeff = c_true, FEMbasis = FEMbasis))

# Create covariates
set.seed(30101997)
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
set.seed(9121997)
ran = range(DatiEsatti)
data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^seq(-3,3,by=0.25)

# Create inferenceDataObjects for the tests

inf_Obj_1 <- inferenceDataObjectBuilder(test = c('one-at-the-time', 'simultaneous', 'one-at-the-time'), interval = c('one-at-the-time','bonferroni','none'),
                                  type=c('wald', 'speckman','eigen-sign-flip'), exact='True', dim=2, beta0 = c(2,-1))

### The following example should give a warning in no-fspai-dev library, performing exact inference in turn of fspai. 
inf_Obj_2 <- inferenceDataObjectBuilder(test = c('one-at-the-time', 'one-at-the-time'), interval = c('simultaneous','one-at-the-time'),
                                    type=c('wald', 'speckman'), exact='False', coeff=rbind(c(1,1),c(1,-1)), dim=2, beta0 = c(1,3), tol_fspai = 0.005)

inf_Obj_3 <- inferenceDataObjectBuilder(test = c('one-at-the-time', 'simultaneous', 'one-at-the-time'), interval = c('one-at-the-time','bonferroni','none'),
                                    type=c('wald', 'speckman','eigen-sign-flip'), exact='True', dim=2, beta0 = c(0,0))

#### Test 1.1: Exact inference, accepting H0
output_CPP_1<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', R_Inference_Data_Object = inf_Obj_1)

# retrieving p_values, H0 is expected to be accepted, hence p_values are expected to be high (e.g. > 0.05) 
output_CPP_1$inference$p_values$wald
output_CPP_1$inference$p_values$speckman
output_CPP_1$inference$p_values$eigen_sign_flip

# retrieving confidence intervals, they are expected to contain the real values of beta, i.e. c(2,-1)
output_CPP_1$inference$CI$wald
output_CPP_1$inference$CI$speckman

#### Test 1.2: Non-Exact inference, accepting H0 for linear combinations (it should give a warning in no-fspai-dev library, performing exact inference in turn of fspai)
output_CPP_2<-smooth.FEM(locations = locations, observations=data, 
                         covariates = cbind(cov1, cov2),
                         FEMbasis=FEMbasis, lambda=lambda,
                         lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', R_Inference_Data_Object = inf_Obj_2)

# retrieving p_values, H0 is expected to be accepted, hence p_values are expected to be high (e.g. > 0.05) 
output_CPP_2$inference$p_values$wald
output_CPP_2$inference$p_values$speckman

# retrieving confidence intervals, they are expected to contain the real values of combinations of beta, i.e. c(1,3)
output_CPP_2$inference$CI$wald
output_CPP_2$inference$CI$speckman


#### Test 1.3: Exact inference, rejecting H0 for linear combinations 
output_CPP_3<-smooth.FEM(locations = locations, observations=data, 
                         covariates = cbind(cov1, cov2),
                         FEMbasis=FEMbasis, lambda=lambda,
                         lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', R_Inference_Data_Object = inf_Obj_3)

# retrieving p_values, H0 is expected to be rejected, hence p_values are expected to be low (e.g. < 0.05) 
output_CPP_3$inference$p_values$wald
output_CPP_3$inference$p_values$speckman

# retrieving confidence intervals, they are expected to contain the real values of beta, i.e. c(2,-1)
output_CPP_3$inference$CI$wald
output_CPP_3$inference$CI$speckman
