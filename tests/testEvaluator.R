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

res_point = eval.FEM_new(FEMfunction,locations,search='naive')

#una sola regione
incidence_matrix = matrix(0,ncol = nrow(mesh$edges))
res_region = eval.FEM_new(FEMfunction, incidence_matrix = incidence_matrix)

