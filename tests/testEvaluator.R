### MESH 1D ###
nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh = create.mesh.1D(nodes,edges,order=1)
plot.mesh.1D(mesh)
points_=matrix(nrow=5,ncol=2)
points_[,1] = runif(5,min=0.25,max=0.75)
points_[,2] = runif(5,min=0.25,max=0.5)

locations=projection.points.1D(mesh,points_)
points(locations[,1],locations[,2], col = 'green')

FEMbasis = fdaPDE::create.FEM.basis(mesh,saveTree = TRUE)
coeff = fs.test(mesh$nodes[,1], mesh$nodes[,2])
FEM = fdaPDE::FEM(coeff, FEMbasis)

res1D_loc_new = eval.FEM_new(FEM,locations,search='naive')
res1D_loc_new

# # # 
incidence_matrix = matrix(0,ncol = nrow(mesh$edges))
incidence_matrix[1,2] = 1

res1D_matrix_new = eval.FEM_new(FEM, incidence_matrix = incidence_matrix,search='naive')
res1D_matrix_new

### MESH 2D ###
data(horseshoe2D)
boundary_nodes = horseshoe2D$boundary_nodes
boundary_segments = horseshoe2D$boundary_segments
locations = horseshoe2D$locations
mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments,order = 2)

FEMbasis = fdaPDE::create.FEM.basis(mesh,saveTree = TRUE)
coeff = fs.test(mesh$nodes[,1], mesh$nodes[,2])
FEMfunction = FEM(coeff, FEMbasis)

res2D_loc_fdaPDE = fdaPDE::eval.FEM(FEMfunction, locations = matrix(c(1, 0.5), ncol = 2))
res2D_loc_new = eval.FEM_new(FEMfunction, locations = matrix(c(1, 0.5), ncol = 2))
err2D_loc = (res2D_loc_fdaPDE - res2D_loc_new)/res2D_loc_fdaPDE
err2D_loc

## Evaluate the mean of the finite element function over the fifth triangle of the mesh
incidence_matrix = matrix(0, ncol = nrow(mesh$triangles))
incidence_matrix[1,5] = 1

res2D_matrix_fdaPDE = fdaPDE::eval.FEM(FEMfunction, incidence_matrix = incidence_matrix)
res2D_matrix_new = eval.FEM_new(FEMfunction, incidence_matrix = incidence_matrix)
err2D_matrix = abs((res2D_matrix_fdaPDE - res2D_matrix_new)/res2D_matrix_fdaPDE)
err2D_matrix

### MESH 2.5D ###
data(hub2.5D)
hub2.5D.nodes = hub2.5D$hub2.5D.nodes
hub2.5D.triangles = hub2.5D$hub2.5D.triangles

## Create mesh from nodes and connectivity matrix:
mesh = create.mesh.2.5D(nodes = hub2.5D.nodes, triangles = hub2.5D.triangles,order=2)

FEMbasis = fdaPDE::create.FEM.basis(mesh)
coeff = fs.test.3D(mesh$nodes[,1], mesh$nodes[,2],mesh$nodes[,3])
FEMfunction = FEM(coeff, FEMbasis)

points_ = matrix(nrow=1,ncol=3)
points_[,1:2] = runif(2, min=min(hub2.5D.nodes[,1]), max=max(hub2.5D.nodes[,1]))
points_[,3]   = runif(1, min=0,max=1)
locations = projection.points.2.5D(mesh,points_)

res2.5D_loc_fdaPDE = fdaPDE::eval.FEM(FEMfunction, locations = locations, search = 'naive')
res2.5D_loc_new = eval.FEM_new(FEMfunction, locations = locations, search = 'naive')
err2.5D_loc = (res2.5D_loc_fdaPDE - res2.5D_loc_new)/res2.5D_loc_fdaPDE
err2.5D_loc

## Evaluate the mean of the finite element function over the fifth triangle of the mesh
incidence_matrix = matrix(0, ncol = nrow(mesh$triangles))
incidence_matrix[1,5] = 1

res2.5D_matrix_fdaPDE = fdaPDE::eval.FEM(FEMfunction, incidence_matrix = incidence_matrix)
res2.5D_matrix_new = eval.FEM_new(FEMfunction, incidence_matrix = incidence_matrix)
err2.5D_matrix = abs((res2.5D_matrix_fdaPDE - res2.5D_matrix_new)/res2.5D_matrix_fdaPDE)
err2.5D_matrix

### MESH 3D ###
data(sphere3Ddata)
nodes=sphere3Ddata$nodes
tetrahedrons=sphere3Ddata$tetrahedrons
mesh=create.mesh.3D(nodes,tetrahedrons,order=2)

FEMbasis = fdaPDE::create.FEM.basis(mesh,TRUE)
coeff = fs.test.3D(mesh$nodes[,1], mesh$nodes[,2],mesh$nodes[,3])
FEMfunction = FEM(coeff, FEMbasis)
locations = matrix(nrow=1,ncol=3)
locations[,1:3] = runif(3, min=-0.9, max=0.9)

res3D_loc_fdaPDE = fdaPDE::eval.FEM(FEMfunction, locations = locations, search = 'tree')
res3D_loc_new = eval.FEM_new(FEMfunction, locations = locations, search = 'tree')
err3D_loc = (res3D_loc_fdaPDE - res3D_loc_new)/res3D_loc_fdaPDE
err3D_loc

## Evaluate the mean of the finite element function over the fifth triangle of the mesh
incidence_matrix = matrix(0, ncol = nrow(mesh$tetrahedrons))
incidence_matrix[1,5] = 1

res3D_matrix_fdaPDE = fdaPDE::eval.FEM(FEMfunction, incidence_matrix = incidence_matrix)
res3D_matrix_new = eval.FEM_new(FEMfunction, incidence_matrix = incidence_matrix)
err3D_matrix = abs((res3D_matrix_fdaPDE - res3D_matrix_new)/res3D_matrix_fdaPDE)
err3D_matrix
