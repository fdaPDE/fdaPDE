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

####### 3D ########
#### Test 1: c-shaped domain ####
#            locations = nodes
#            with covariates
#            no BC
#            order FE = 1
rm(list=ls())
graphics.off()

# set the working directory to the data fold 
#setwd("~/Scrivania/project/fdaPDE/data")

load("horseshoe3D.RData")

mesh = horseshoe3D

FEMbasis=create.FEM.basis(mesh)

ndata = FEMbasis$nbasis

# Create covariates
set.seed(30101997)
cov1=rnorm(ndata, mean = 1, sd = 2)
cov2=sin(mesh$nodes[,1])

# Exact solution (pointwise at nodes)
DatiEsatti=fs.test.3D(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]) + 2*cov1 - cov2
plot(FEM(DatiEsatti, FEMbasis))

# Add error to simulate data
set.seed(9121997)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda= 10^seq(-3,3,by=0.25)

# Create inferenceDataObjects for the tests

inf_Obj_1 <- inferenceDataObjectBuilder(test = c('one-at-the-time', 'simultaneous', 'one-at-the-time'), interval = c('one-at-the-time','bonferroni','none'),
                                        type=c('wald', 'speckman','eigen-sign-flip'), exact='True', dim=2, beta0 = c(2,-1))

### The following example should give a warning in no-fspai-dev library, performing exact inference in turn of fspai. 
inf_Obj_2 <- inferenceDataObjectBuilder(test = c('one-at-the-time', 'one-at-the-time'), interval = c('simultaneous','one-at-the-time'),
                                        type=c('wald', 'speckman'), exact='False', coeff=rbind(c(1,1),c(1,-1)), dim=2, beta0 = c(1,3), tol_fspai = 0.005)

inf_Obj_3 <- inferenceDataObjectBuilder(test = c('one-at-the-time', 'simultaneous', 'one-at-the-time'), interval = c('one-at-the-time','bonferroni','none'),
                                        type=c('wald', 'speckman','eigen-sign-flip'), exact='True', dim=2, beta0 = c(0,0))

#### Test 1.1: Exact inference, accepting H0
output_CPP_1<-smooth.FEM(observations=data, 
                         covariates = cbind(cov1, cov2),
                         FEMbasis=FEMbasis, lambda=lambda,
                         lambda.selection.criterion='grid', lambda.selection.lossfunction = "GCV", DOF.evaluation = "exact",
                         R_Inference_Data_Object = inf_Obj_1)

# retrieving p_values, H0 is expected to be accepted, hence p_values are expected to be high (e.g. > 0.05) 
output_CPP_1$inference$p_values$wald
output_CPP_1$inference$p_values$speckman
output_CPP_1$inference$p_values$eigen_sign_flip

# retrieving confidence intervals, they are expected to contain the real values of beta, i.e. c(2,-1)
output_CPP_1$inference$CI$wald
output_CPP_1$inference$CI$speckman

#### Test 1.2: Non-Exact inference, accepting H0 for linear combinations (it should give a warning in no-fspai-dev library, performing exact inference in turn of fspai)
output_CPP_2<-smooth.FEM(observations=data, 
                         covariates = cbind(cov1, cov2),
                         FEMbasis=FEMbasis, lambda=lambda,
                         lambda.selection.criterion='grid', lambda.selection.lossfunction = "GCV", DOF.evaluation = "exact",
                         R_Inference_Data_Object = inf_Obj_2)

# retrieving p_values, H0 is expected to be accepted, hence p_values are expected to be high (e.g. > 0.05) 
output_CPP_2$inference$p_values$wald
output_CPP_2$inference$p_values$speckman

# retrieving confidence intervals, they are expected to contain the real values of combinations of beta, i.e. c(1,3)
output_CPP_2$inference$CI$wald
output_CPP_2$inference$CI$speckman


#### Test 1.3: Exact inference, accepting H0 for linear combinations 
output_CPP_3<-smooth.FEM(observations=data, 
                         covariates = cbind(cov1, cov2),
                         FEMbasis=FEMbasis, lambda=lambda,
                         lambda.selection.criterion='grid', lambda.selection.lossfunction = "GCV", DOF.evaluation = "exact",
                         R_Inference_Data_Object = inf_Obj_3)

# retrieving p_values, H0 is expected to be rejected, hence p_values are expected to be low (e.g. < 0.05) 
output_CPP_3$inference$p_values$wald
output_CPP_3$inference$p_values$speckman

# retrieving confidence intervals, they are expected to contain the real values of beta, i.e. c(2,-1)
output_CPP_3$inference$CI$wald
output_CPP_3$inference$CI$speckman
