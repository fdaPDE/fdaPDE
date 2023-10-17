test_that("SR-PDE 2.5D Horseshoe domain",{
  
options(warn=-1)
foldername = test_path("../data/SR-PDE/test_9/")

data(quasicircle2Dareal)
mesh = quasicircle2Dareal$mesh
incidence_matrix = quasicircle2Dareal$incidence_matrix
DatiEsatti = quasicircle2Dareal$data

FEMbasis = create.FEM.basis(mesh)

# Add error to simulate data
set.seed(5839745)
data = DatiEsatti + rnorm(length(DatiEsatti), sd = 0.05*(max(DatiEsatti)-min(DatiEsatti)))

# Set smoothing parameter
lambda = 10^-3

# Set BC (equal to zero)
BC = NULL
BC$BC_indices = which(mesh$nodesmarkers == 1)
BC$BC_values = rep(0,length(BC$BC_indices))

# Set sv-PDE parameters
R = 2.8
K1 = 0.1
K2 = 0.2
beta = 0.5

K_func<-function(points)
{
  output = array(0, c(2, 2, nrow(points)))
  for (i in 1:nrow(points))
    output[,,i] = 10*rbind(c(points[i,2]^2 + K1*points[i,1]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2),
                             (K1-1)*points[i,1]*points[i,2]),
                           c((K1-1)*points[i,1]*points[i,2],
                             points[i,1]^2 + K1*points[i,2]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2)))
  output
}

b_func<-function(points)
{
  output = array(0, c(2, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = 10*beta*c(points[i,1],points[i,2])
  output
}

c_func<-function(points)
{
  rep(c(0), nrow(points))
}

u_func<-function(points)
{
  rep(c(0), nrow(points))
}
PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)

#### Test 9.1: Forcing term = 0  grid without GCV 
invisible(capture.output(sol<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters)))

load(file=paste0(foldername,"/test_9_1.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);

#### Test 9.2: Forcing term = 0  grid with exact  GCV 
invisible(capture.output(sol<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', 
                       lambda.selection.lossfunction='GCV')))

load(file=paste0(foldername,"/test_9_2.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);

#### Test 9.3: Forcing term != 0 without GCV 
# forcing function != 0
u_func<-function(points)
{
  output = array(0, c(1, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = -ifelse((points[i,1]^2+points[i,2]^2)<1,100,0)
  output
}

xgrid=seq(from=-3,to=3,by=0.1)
ygrid=seq(from=-3,to=3,by=0.1)

PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)


invisible(capture.output(sol<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters)))
load(file=paste0(foldername,"/test_9_3.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);


#### Test 9.4: Forcing term != 0 grid with exact GCV 
# forcing function != 0
u_func<-function(points)
{
  output = array(0, c(1, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = -ifelse((points[i,1]^2+points[i,2]^2)<1,100,0)
  output
}

xgrid=seq(from=-3,to=3,by=0.1)
ygrid=seq(from=-3,to=3,by=0.1)

PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)


invisible(capture.output(sol<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', 
                       lambda.selection.lossfunction='GCV')))
load(file=paste0(foldername,"/test_9_4.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);


#### Test 9.5: Forcing term != 0 Newton finite differences with exact GCV 
# forcing function != 0
u_func<-function(points)
{
  output = array(0, c(1, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = -ifelse((points[i,1]^2+points[i,2]^2)<1,100,0)
  output
}

xgrid=seq(from=-3,to=3,by=0.1)
ygrid=seq(from=-3,to=3,by=0.1)

PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)


invisible(capture.output(sol<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='newton_fd', DOF.evaluation='exact',
                       lambda.selection.lossfunction='GCV')))
load(file=paste0(foldername,"/test_9_5.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);

#### Test 9.6: BC != 0      without GCV 
# Add a constat to the data to change true BC
data_backup=data #save a copy of original data
data = data + 5

# Set new value for the BC
BC$BC_values = rep(5,length(BC$BC_indices))

invisible(capture.output(sol<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters)))

load(file=paste0(foldername,"/test_9_6.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);

data=data_backup #restore original data for next tests

#### Test 9.7: BC != 0      grid with exact GCV 
# Add a constat to the data to change true BC
data_backup=data #save a copy of original data
data = data + 5

# Set new value for the BC
BC$BC_values = rep(5,length(BC$BC_indices))

invisible(capture.output(output_CPP<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', 
                       lambda.selection.lossfunction='GCV')))

load(file=paste0(foldername,"/test_9_7.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);

#### Test 9.8 on inference creating covariates and special setting (do not )

# Create covariates
# set.seed(509875)
# cov1 = rnorm(length(DatiEsatti), mean = 2, sd = 1)
# 
# # Add error to simulate data
# set.seed(5839745)
# data = DatiEsatti + 3*cov1 
# data = data + rnorm(length(DatiEsatti), sd = 0.1)
# observations=data
# 
# # Inference obj: separate tests for beta and f. Wald, Speckman, ESF, Enhanced ESF p_values for beta; Wald, Sign-Flip and ESF p_values for f
# inf_beta <- inferenceDataObjectBuilder (test='oat', interval='oat',  dim=2, n_cov=1, type=c('w', 's', 'esf', 'enh-esf'), beta0 = 3, component='parametric', n_flip=1000, f_var=T)
# inf_f <- inferenceDataObjectBuilder(test = 'sim', type=c('w', 'sf', 'esf'), component = 'nonparametric', dim=2, n_cov=1)
# 
# invisible(capture.output(sol<-smooth.FEM(observations=data, 
#                        covariates = cov1,
#                        incidence_matrix = incidence_matrix,
#                        FEMbasis=FEMbasis, 
#                        lambda=lambda,
#                        BC = BC, 
#                        PDE_parameters = PDE_parameters,
#                        lambda.selection.criterion='newton', DOF.evaluation='exact', 
#                        lambda.selection.lossfunction='GCV',
#                        inference.data.object = inf_beta)))
# 
# load(file=paste0(foldername,"/test_9_8.RData"))
# expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);
# expect_equal( max(abs((sol$solution$beta-output_CPP$solution$beta))) < tol, TRUE);
# expect_equal( max(abs((sol$inference$beta$p_values$wald[[1]]-
#                          output_CPP$inference$beta$p_values$wald[[1]]))) < tol, TRUE);
# 
# invisible(capture.output(sol<-smooth.FEM(observations=data, 
#                        covariates = cov1,
#                        incidence_matrix = incidence_matrix,
#                        FEMbasis=FEMbasis, 
#                        lambda=lambda,
#                        BC = BC, 
#                        PDE_parameters = PDE_parameters,
#                        lambda.selection.criterion='newton', DOF.evaluation='exact', 
#                        lambda.selection.lossfunction='GCV',
#                        inference.data.object = inf_f)))
# 
# load(file=paste0(foldername,"/test_9_9.RData"))
# expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);
# expect_equal( max(abs((sol$solution$beta-output_CPP$solution$beta))) < tol, TRUE);
# expect_equal( max(abs((sol$inference$f$p_values$wald[[1]]-
#                          output_CPP$inference$f$p_values$wald[[1]]))) < tol, TRUE);
})