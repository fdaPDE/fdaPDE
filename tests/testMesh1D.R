
load("~/Scrivania/fdaPDE/data/GraphOrder1.RData")
nodes = GraphOrder1$nodes
edges = GraphOrder1$edges

mesh = create.mesh.1D(nodes, edges)
plot.mesh.1D(mesh)
# refine mesh by length 
ref_mesh = refine.mesh.1D(mesh=mesh, delta=0.1)
plot.mesh.1D(ref_mesh)
########################
#Split mesh1d ORDER 1
n_elem = dim(edges)[1]

refinedMesh <- refine.by.splitting.mesh.1D(mesh = mesh1)
plot.mesh.1D(refinedMesh)

refinedMesh2 <- refine.by.splitting.mesh.1D(mesh = refinedMesh)

plot.mesh.1D(refinedMesh2)

#ORDER 2 
load("~/Scrivania/fdaPDE/data/GraphOrder2.RData")
nodes = GraphOrder2$nodes
edges = GraphOrder2$edges

mesh2 = create.mesh.1D(nodes, edges, order = 2)
plot.mesh.1D(mesh2)
# refine order 2#
ref_mesh.order2 = refine.mesh.1D(mesh2, delta=0.1)
plot.mesh.1D(ref_mesh.order2)

# mesh order2 da grafo ordine1 
load("~/Scrivania/fdaPDE/data/GraphOrder1.RData")
nodes = GraphOrder1$nodes
edges = GraphOrder1$edges
mesh2 = create.mesh.1D(nodes, edges, order = 2)

plot.mesh.1D(mesh2)

#refine.by.splitting order2
mesh2_refined <- refine.by.splitting.mesh.1D(mesh2)
plot.mesh.1D(mesh2_refined)

