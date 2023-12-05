test_that("SR-PDE unit square second order FE",{
  
options(warn=-1)
foldername = test_path("../data/SR-PDE/test_8/")

x = seq(0,1, length.out = 11)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations, order = 2)

FEMbasis=create.FEM.basis(mesh)

# Test function
a1=1
a2=4
z<-function(p){  
  a1*sin(2*pi*p[,1])*cos(2*pi*p[,2])+a2*sin(3*pi*p[,1])}

# Exact solution (pointwise at nodes)
sol_exact=z(mesh$nodes)

DatiEsatti=z(locations)
ndati = length(DatiEsatti)

# Add error to simulate data
set.seed(7893475)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(ndati, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda= 10^seq(-6,-3,by=0.25)

# Set PDE parameters
PDE_parameters = list(K = matrix(c(1,0,0,4), nrow = 2), b = c(0,0), c = 0)

#### Test 8.1: Without GCV
invisible(capture.output(sol<-smooth.FEM(observations=data, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda, 
                       PDE_parameters=PDE_parameters)))

load(file=paste0(foldername,"/test_8_1.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-sol_ref$fit.FEM$coeff))) < tol, TRUE);

#### Test 8.2: grid with exact GCV
invisible(capture.output(sol<-smooth.FEM(observations=data, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       PDE_parameters=PDE_parameters,
                       lambda.selection.criterion='grid',
                       DOF.evaluation='exact', lambda.selection.lossfunction='GCV')))

load(file=paste0(foldername,"/test_8_2.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-sol_ref$fit.FEM$coeff))) < tol, TRUE);

### Test 8.3: Newton exact method with exact GCV, default initial lambda and tolerance
invisible(capture.output(sol<-smooth.FEM(observations=data, 
                       FEMbasis=FEMbasis, 
                       PDE_parameters=PDE_parameters,
                       lambda.selection.criterion='newton', DOF.evaluation='exact',
                       lambda.selection.lossfunction='GCV')))

load(file=paste0(foldername,"/test_8_3.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-sol_ref$fit.FEM$coeff))) < tol, TRUE);

### Test 8.4: Newton_fd method with exact GCV, default initial lambda and tolerance
invisible(capture.output(sol<-smooth.FEM(observations=data, 
                       FEMbasis=FEMbasis, 
                       PDE_parameters=PDE_parameters,
                       lambda.selection.criterion='newton_fd', DOF.evaluation='exact', 
                       lambda.selection.lossfunction='GCV')))

load(file=paste0(foldername,"/test_8_4.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-sol_ref$fit.FEM$coeff))) < tol, TRUE);

### Test 8.5: Test on inference: generating covariates in order to perform inference
set.seed(509875)
cov1 = rnorm(ndati, mean = 1, sd = 2)
cov2 = sin(locations[,1])

DatiEsatti = DatiEsatti + 2.5*cov1+1.5*cov2

set.seed(7893475)
ran=range(DatiEsatti)
data = DatiEsatti + rnorm(ndati, mean=0, sd=0.05*abs(ran[2]-ran[1]))

### Wald and Speckman CI and p_values for beta; Wald, Sign-flip and ESF p_values for f
inf_obj<-inferenceDataObjectBuilder(test = c(rep('oat',3),rep('sim',3)), interval = c('sim','oat','bonf',rep('none',3)),
                                    component=c(rep('parametric',3), rep('nonparametric',3)),
                                    dim = 2, n_cov = 2, type = c("w","w","s","w","sf","esf"),f0 = fs.test, beta0 = c(2.5,1.5))
invisible(capture.output(sol<-smooth.FEM(observations=data,
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, 
                       PDE_parameters=PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', 
                       lambda.selection.lossfunction='GCV', lambda = lambda,
                       inference.data.object = inf_obj)))

load(file=paste0(foldername,"/test_8_5.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-sol_ref$fit.FEM$coeff))) < tol, TRUE);
expect_equal( max(abs((sol$inference$beta$p_values$wald[[1]]-
                         sol_ref$inference$beta$p_values$wald[[1]]))) < tol, TRUE);
expect_equal( max(abs((sol$inference$beta$p_values$speckman[[1]]-
                         sol_ref$inference$beta$p_values$speckman[[1]]))) < tol, TRUE);
expect_equal( max(abs((sol$inference$f$p_values$wald[[1]]-
                         sol_ref$inference$f$p_values$wald[[1]]))) < tol, TRUE);
})