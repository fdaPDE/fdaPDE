
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
invisible(capture.output(sol_ref<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters)))

save(sol_ref, file=paste0(foldername,"/test_9_1.RData"))

#### Test 9.2: Forcing term = 0  grid with exact  GCV 
invisible(capture.output(sol_ref<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', 
                       lambda.selection.lossfunction='GCV')))

save(sol_ref, file=paste0(foldername,"/test_9_2.RData"))

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


invisible(capture.output(sol_ref<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters)))
save(sol_ref, file=paste0(foldername,"/test_9_3.RData"))

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


invisible(capture.output(sol_ref<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', 
                       lambda.selection.lossfunction='GCV')))
save(sol_ref, file=paste0(foldername,"/test_9_4.RData"))

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


invisible(capture.output(sol_ref<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix,
                       FEMbasis=FEMbasis, 
                       BC = BC, 
                       PDE_parameters = PDE_parameters,
                       lambda.selection.criterion='newton_fd', DOF.evaluation='exact',
                       lambda.selection.lossfunction='GCV')))
save(sol_ref, file=paste0(foldername,"/test_9_5.RData"))

#### Test 9.6: BC != 0      without GCV 
# Add a constat to the data to change true BC
data_backup=data #save a copy of original data
data = data + 5

# Set new value for the BC
BC$BC_values = rep(5,length(BC$BC_indices))

invisible(capture.output(sol_ref<-smooth.FEM(observations=data, 
                       incidence_matrix = incidence_matrix, 
                       FEMbasis=FEMbasis, 
                       lambda=lambda,
                       BC = BC, 
                       PDE_parameters = PDE_parameters)))

save(sol_ref, file=paste0(foldername,"/test_9_6.RData"))

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

save(sol_ref, file=paste0(foldername,"/test_9_7.RData"))
