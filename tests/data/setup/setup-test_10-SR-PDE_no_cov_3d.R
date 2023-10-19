
### Second order FEM in 3D ################################
# Compare the accuracy and robustness of first 
# and second order methods in 3D settings.
# Order 2 in 3D is a new functionality.

options(warn=-1)

foldername = test_path("../data/SR-PDE/test_10/")

# Function to generate random points in a sphere
rsphere <- function(n, r = 1.0, surface_only = FALSE) {
  phi       <- runif(n, 0.0, 2.0 * pi)
  cos_theta <- runif(n, -1.0, 1.0)
  sin_theta <- sqrt((1.0-cos_theta)*(1.0+cos_theta))
  radius <- r
  if (surface_only == FALSE) {
    radius <- r * runif(n, 0.0, 1.0)^(1.0/3.0)
  }
  
  x <- radius * sin_theta * cos(phi)
  y <- radius * sin_theta * sin(phi)
  z <- radius * cos_theta
  
  cbind(x, y, z)
}

# Percentage of noise standard deviation on data range
# Try also different values, e.g 0.01 or 1
noisepercent<-0.25

data("sphere3Ddata")

# Build spherical mesh of order 2 
mesh_sphere<-create.mesh.3D(sphere3Ddata$nodes,
                            sphere3Ddata$tetrahedrons,
                            order=2)

FEMbasis <- create.FEM.basis(mesh_sphere)

# Set smoothing parameters
lambda= 10^seq(-6,1,by=.25)

set.seed(5847947)

a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

nnodes<-nrow(mesh_sphere$nodes)

# Evaluate exact solution on mesh nodes
func_evaluation = a1* sin(2*pi*mesh_sphere$nodes[,1]) +  a2* sin(2*pi*mesh_sphere$nodes[,2]) +  a3*sin(2*pi*mesh_sphere$nodes[,3]) + 1

ran=range(func_evaluation)

exact_sol=func_evaluation

# Add noise to exact solution to generate data for the simulation
data= exact_sol + rnorm(nnodes,
                        mean=0,
                        sd=noisepercent * (ran[2]-ran[1]))

# Compute the solution for each lambda
invisible(capture.output(sol_ref <- smooth.FEM(observations = data, 
                         FEMbasis = FEMbasis, 
                         lambda = lambda)))
save(sol_ref, file=paste0(foldername,"/test_10_1.RData"))

### PDE penalization in 3D ################################
# Include PDE parameters and use a penalization term with
# anisotropic diffusion and possibly advection and reaction
# terms. All this is a new functionality.

# Function to generate random points in a sphere
rsphere <- function(n, r = 1.0, surface_only = FALSE) {
  phi       <- runif(n, 0.0, 2.0 * pi)
  cos_theta <- runif(n, -1.0, 1.0)
  sin_theta <- sqrt((1.0-cos_theta)*(1.0+cos_theta))
  radius <- r
  if (surface_only == FALSE) {
    radius <- r * runif(n, 0.0, 1.0)^(1.0/3.0)
  }
  
  x <- radius * sin_theta * cos(phi)
  y <- radius * sin_theta * sin(phi)
  z <- radius * cos_theta
  
  cbind(x, y, z)
}

# Build mesh: Sphere

mesh_sphere<-create.mesh.3D(sphere3Ddata$nodes,
                            sphere3Ddata$tetrahedrons,
                            order=1)

FEMbasis <- create.FEM.basis(mesh_sphere)

set.seed(5847947)

# Exact test function
nnodes = nrow(mesh_sphere$nodes)

# Set smoothing parameter
lambda=10^seq(-6, 1, by=.25)

# Set PDE parameters (in this case they are constant)
PDE_parameters_anys = list(K = diag(c(1,.5,1)), b = c(0,0,0), c = -4*pi^2)

# Evaluate exact solution on mesh nodes
exact_sol =  sin(2*pi*mesh_sphere$nodes[,1]) +  2 * sin(2*pi*mesh_sphere$nodes[,2]) +  sin(2*pi*mesh_sphere$nodes[,3])

# Add noise to generate data - 10% level of noise
data=exact_sol + rnorm(nrow(mesh_sphere$nodes), mean=0, sd=0.10*diff(range(exact_sol)))

# Compute the solution for each lambda
invisible(capture.output(sol_ref <- smooth.FEM(observations = data, PDE_parameters = PDE_parameters_anys,
                         FEMbasis = FEMbasis, lambda = lambda)))
save(sol_ref, file=paste0(foldername,"/test_10_2.RData"))
