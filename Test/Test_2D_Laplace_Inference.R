##########################################
############## TEST SCRIPT ###############
##########################################
### Inference Tests ###
library(fdaPDE)

######## Test 2: c-shaped domain ########
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

# Set inference object
R_Inference_Object=inferenceDataObjectBuilder (test  = 'simultaneous', beta0 = c(1.8,-3.55), interval='one-at-the-time', coeff = matrix(data=c(1.5,0,0,1),nrow = 2, ncol = 2,byrow = TRUE), exact='True', dim=2, level = 0.05)
start=Sys.time()
#### Test 2.1: Without GCV
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda#,R_Inference_Data_Object = R_Inference_Object
                       )
end=Sys.time()
end-start
#### Test 2.2: grid with exact GCV
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', R_Inference_Data_Object = R_Inference_Object)
plot(log10(lambda), output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

#### Test 2.3: grid with stochastic GCV
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV', R_Inference_Data_Object = R_Inference_Object)
plot(log10(lambda), output_CPP$optimization$GCV_vector)
image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

### Test 2.4: Newton exact method with exact  GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis,
                       lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV' , R_Inference_Data_Object = R_Inference_Object)

image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

### Test 2.5: Newton_fd method with exact GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton_fd', DOF.evaluation='exact', lambda.selection.lossfunction='GCV' , R_Inference_Data_Object = R_Inference_Object)

image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta

### Test 2.6: Newton_fd method with stochastic GCV, default initial lambda and tolerance
output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, 
                       lambda.selection.criterion='newton_fd', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV', R_Inference_Data_Object = R_Inference_Object)

image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

output_CPP$solution$beta
