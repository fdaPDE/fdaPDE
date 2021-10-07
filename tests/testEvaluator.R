# MESH 1D
nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh = create.mesh.1D(nodes,edges,order=2)
plot.mesh.1D(mesh)
points_=matrix(nrow=5,ncol=2)
points_[,1] = runif(5,min=0.25,max=0.75)
points_[,2] = runif(5,min=0.25,max=0.5)

locations=projection.points.1D(mesh,points_)
points(locations[,1],locations[,2], col = 'green')

FEMbasis = fdaPDE::create.FEM.basis(mesh)
coeff = fs.test(mesh$nodes[,1], mesh$nodes[,2])
FEM = fdaPDE::FEM(coeff, FEMbasis)

res_point_new = eval.FEM_new(FEM,locations,search='naive')

#una sola regione
incidence_matrix = matrix(0,ncol = nrow(mesh$edges))
incidence_matrix[1,2] = 1
res_region_new = eval.FEM_new(FEM, incidence_matrix = incidence_matrix,search='naive')

# MESH 2D
data(horseshoe2D)
boundary_nodes = horseshoe2D$boundary_nodes
boundary_segments = horseshoe2D$boundary_segments
locations = horseshoe2D$locations
mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)

FEMbasis = fdaPDE::create.FEM.basis(mesh)
coeff = fs.test(mesh$nodes[,1], mesh$nodes[,2])
FEMfunction = FEM(coeff, FEMbasis)

res2D_loc_fdaPDE = fdaPDE::eval.FEM(FEMfunction, locations = matrix(c(1, 0.5), ncol = 2))

res2D_loc_new = eval.FEM_new(FEMfunction, locations = matrix(c(1, 0.5), ncol = 2))

err2D_loc = (res2D_loc_fdaPDE - res2D_loc_new)/res2D_loc_fdaPDE
## Evaluate the mean of the finite element function over the fifth triangle of the mesh
incidence_matrix = matrix(0, ncol = nrow(mesh$triangles))
incidence_matrix[1,5] = 1

res2D_matrix_fdaPDE = fdaPDE::eval.FEM(FEMfunction, incidence_matrix = incidence_matrix)
res2D_matrix_new = eval.FEM_new(FEMfunction, incidence_matrix = incidence_matrix)

err2D_matrix = abs((res2D_matrix_fdaPDE - res2D_matrix_new)/res2D_matrix_fdaPDE)

