### test Evaluator time ! ###

# eval.FEM.time() in FEMEvaluator.R
data(horseshoe2D)
boundary_nodes = horseshoe2D$boundary_nodes
boundary_segments = horseshoe2D$boundary_segments
locations = horseshoe2D$locations

mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
## Create the FEM basis
FEMbasis = create.FEM.basis(mesh)
## Compute the coeff vector evaluating the desired function at the mesh nodes
## In this case we consider the fs.test() function introduced by Wood et al. 2008
time = 1:5
coeff = as.amtrix(rep(fs.test(mesh$nodes[,1], mesh$nodes[,2]),5)*time)
## Create the FEM.time object
FEM_time_function = FEM.time(coeff=coeff, time_mesh=1:5, FEMbasis=FEMbasis, FLAG_PARABOLIC=TRUE)
loc=matrix(c(-0.92,0),ncol=2)

num=3 
loc=matrix(nrow=num,ncol=2)
loc[] = runif(num*2,0,1)

res2D_loc_fdaPDE = eval.FEM.time(FEM_time_function, locations = loc, 
                             time.instants = time)
res2D_loc_new = eval.FEM.time.new(FEM_time_function, locations = loc,
                              time.instants = time)

err2D_loc = abs((res2D_loc_new - res2D_loc_fdaPDE)/res2D_loc_fdaPDE)
err2D_loc
# matrice loc x time
#evaluationsT = matrix(evaluations, nrow=dim(loc)[1],ncol=length(time))
# matrice incidenza 
## Evaluate the mean of the finite element function over the fifth triangle of the mesh
incidence_matrix = matrix(0, ncol = nrow(mesh$triangles))
incidence_matrix[1,5] = 1

res2D_matrix_fdaPDE = fdaPDE::eval.FEM.time(FEM_time_function, incidence_matrix = incidence_matrix,
                                            time.instants = time)
res2D_matrix_new = eval.FEM.time.new(FEM_time_function, incidence_matrix = incidence_matrix,
                                time.instants = time)
err2D_matrix = abs((res2D_matrix_fdaPDE - res2D_matrix_new)/res2D_matrix_fdaPDE)
err2D_matrix

### mesh 2.5D ### 
data(hub2.5D)
hub2.5D.nodes = hub2.5D$hub2.5D.nodes
hub2.5D.triangles = hub2.5D$hub2.5D.triangles

## Create mesh from nodes and connectivity matrix:
mesh = create.mesh.2.5D(nodes = hub2.5D.nodes, triangles = hub2.5D.triangles,order=1)

FEMbasis = fdaPDE::create.FEM.basis(mesh)
time = 1:5
coeff = rep(fs.test.3D(mesh$nodes[,1], mesh$nodes[,2],mesh$nodes[,3]),5)*time
## Create the FEM.time object
FEM_time_function = FEM.time(coeff=coeff, time_mesh=1:5, FEMbasis=FEMbasis, FLAG_PARABOLIC=TRUE)
#Boh coefficienti uguali al variare di T... non saprei che fare ma funziona!
points_ = matrix(nrow=1,ncol=3)
points_[,1:2] = runif(2, min=min(hub2.5D.nodes[,1]), max=max(hub2.5D.nodes[,1]))
points_[,3]   = runif(1, min=0,max=1)
locations = projection.points.2.5D(mesh,points_)

res2.5D_loc_fdaPDE = fdaPDE::eval.FEM.time(FEM_time_function, locations = locations,
                                           time.instants = time, search =  'naive')

res2.5D_loc_new = eval.FEM.time.new(FEM_time_function, locations = locations,
                                    time.instants=time)

err2.5D_loc = (res2.5D_loc_fdaPDE - res2.5D_loc_new)/res2.5D_loc_fdaPDE
err2.5D_loc

## Evaluate the mean of the finite element function over the fifth triangle of the mesh
incidence_matrix = matrix(0, ncol = nrow(mesh$triangles))
incidence_matrix[1,5] = 1

res2.5D_matrix_fdaPDE = fdaPDE::eval.FEM.time(FEM_time_function, incidence_matrix = incidence_matrix,
                                              time.instants = time)
res2.5D_matrix_new = eval.FEM.time.new(FEM_time_function, incidence_matrix = incidence_matrix,
                                          time.instants = time)
err2.5D_matrix = abs((res2.5D_matrix_fdaPDE - res2.5D_matrix_new)/res2.5D_matrix_fdaPDE)
err2.5D_matrix

### mesh3D ###
data(sphere3Ddata)
nodes=sphere3Ddata$nodes
tetrahedrons=sphere3Ddata$tetrahedrons
mesh=create.mesh.3D(nodes,tetrahedrons,order=2)
FEMbasis = fdaPDE::create.FEM.basis(mesh)

time=1:5
coeff = rep(fs.test.3D(mesh$nodes[,1], mesh$nodes[,2],mesh$nodes[,3]),5)*time
## Create the FEM.time object
FEM_time_function = FEM.time(coeff=coeff, time_mesh=1:5, FEMbasis=FEMbasis, FLAG_PARABOLIC=TRUE)

locations = matrix(nrow=1,ncol=3)
locations[,1:3] = runif(3, min=-0.9, max=0.9)

res3D_loc_fdaPDE = fdaPDE::eval.FEM.time(FEM_time_function, locations = locations, search = 'tree',
                                         time.instants = time)
res3D_loc_new = eval.FEM.time.new(FEM_time_function, locations = locations, search = 'tree',
                                  time.instants = time)
err3D_loc = (res3D_loc_fdaPDE - res3D_loc_new)/res3D_loc_fdaPDE
err3D_loc

## Evaluate the mean of the finite element function over the fifth triangle of the mesh
incidence_matrix = matrix(0, ncol = nrow(mesh$tetrahedrons))
incidence_matrix[1,5] = 1

res3D_matrix_fdaPDE = fdaPDE::eval.FEM.time(FEM_time_function, incidence_matrix = incidence_matrix,
                                       time.instants = time)
res3D_matrix_new = eval.FEM.time.new(FEM_time_function, incidence_matrix = incidence_matrix,
                                     time.instants = time)
err3D_matrix = abs((res3D_matrix_fdaPDE - res3D_matrix_new)/res3D_matrix_fdaPDE)
err3D_matrix


