### DE init graph ###
nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh = create.mesh.1.5D(nodes,edges,order=1)
mesh = refine.mesh.1.5D(mesh=mesh,0.01)
plot(mesh)

FEMbasis <- create.FEM.basis(mesh)

set.seed(10)

## Generate data
n <- 40
data_x <- rnorm(n, mean = 0.5,sd = 0.25*0.5)
data_y <- rnorm(n, mean = 0.5,sd = 0.25*0.25)
data <- cbind(data_x, data_y)
data <- projection.points.1.5D(mesh = mesh,locations = data)

lambda = 0.1
sol = DE.heat.FEM(data, FEMbasis, lambda, heatStep=0.1, heatIter=500, init="Heat")
R_plot_graph(FEM(sol$f_init, FEMbasis))

# provo #
nelem = dim(mesh$edges)[1]
incidence_matrix=diag(1,nrow=nelem,ncol=nelem)
integral = sum( eval.FEM(FEM=FEM(sol$f_init,FEMbasis), incidence_matrix = incidence_matrix))
integral
