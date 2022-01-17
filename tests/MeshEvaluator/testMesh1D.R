
load("~/Scrivania/fdaPDE/data/GraphOrder1.RData")
nodes = GraphOrder1$nodes
edges = GraphOrder1$edges

mesh = create.mesh.1.5D(nodes, edges)
plot.mesh.1.5D(mesh)
# refine mesh by length 
ref_mesh = refine.mesh.1.5D(mesh=mesh, delta=0.1)
plot.mesh.1.5D(ref_mesh)
########################
#Split mesh1d ORDER 1
n_elem = dim(edges)[1]

refinedMesh <- refine.by.splitting.mesh.1.5D(mesh = mesh)
plot.mesh.1.5D(refinedMesh)

refinedMesh2 <- refine.by.splitting.mesh.1.5D(mesh = refinedMesh)

plot.mesh.1.5D(refinedMesh2)

#ORDER 2 
load("~/Scrivania/fdaPDE/data/GraphOrder2.RData")
nodes = GraphOrder2$nodes
edges = GraphOrder2$edges

mesh2 = create.mesh.1.5D(nodes, edges, order = 2)
plot.mesh.1.5D(mesh2)
# refine order 2#
ref_mesh.order2 = refine.mesh.1.5D(mesh2, delta=0.1)
plot.mesh.1.5D(ref_mesh.order2)

# mesh order2 da grafo ordine1 
load("~/Scrivania/fdaPDE/data/GraphOrder1.RData")
nodes = GraphOrder1$nodes
edges = GraphOrder1$edges
mesh2 = create.mesh.1.5D(nodes, edges, order = 2)

plot.mesh.1.5D(mesh2)

#refine.by.splitting order2
mesh2_refined <- refine.by.splitting.mesh.1.5D(mesh2)
plot.mesh.1.5D(mesh2_refined)


### refine on bifurcation ###

nodes=matrix(10*c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh = create.mesh.1.5D(nodes,edges,order=1)
dim(mesh$nodes)
plot(mesh)

mesh_ref = refine.mesh.1.5D(mesh,1)
dim(mesh_ref$nodes)
plot(mesh_ref)

## al contrario ##
nodes = matrix(c(0.75,0.,0.75,0.5,0.5,0.25,0.25,0.25), nrow = 4, byrow=TRUE)
edges=matrix(c(1,3,2,3,3,4),nrow = 3,byrow = TRUE)
mesh = create.mesh.1.5D(nodes,edges,order=1)
plot(mesh)
mesh_ref = refine.mesh.1.5D(mesh,0.01)
plot(mesh_ref)

dijkstra = Dijkstra(mesh_ref,3)



### 

x = seq(from=0., to=10. , length.out=10)
y = x 
nodes = expand.grid(x,y)
mesh.2D = create.mesh.2D(nodes)
plot(mesh.2D)

mesh = create.mesh.1.5D(mesh.2D$nodes, mesh.2D$edges)
plot(mesh, show.nodes = TRUE, pch=16)

mesh.ref = refine.mesh.1.5D(mesh,delta=1)
plot(mesh.ref, show.nodes = TRUE, pch=16)
nrow(mesh.ref$nodes)

#NON VA QUALCOSSA NEL BRICKWALL (?)
library(spatstat)
data("spiders")

vertices = cbind(spiders$domain$vertices$x, spiders$domain$vertices$y)
edges = cbind(spiders$domain$from, spiders$domain$to)
mesh= create.mesh.1.5D(vertices, edges)
plot(mesh,show.nodes = TRUE, pch=16)
vertices

mesh.ref = refine.mesh.1.5D(mesh, 50)

plot(mesh.ref, show.nodes = TRUE, pch=16)

A = vector(mode="integer")
B = vector(mode="integer")
nedges = nrow(mesh.ref$edges)
for(i in 1:nedges ){
  if( mesh.ref$edges[i,1] == 131 ){
    A = append(A, 131)
    B = append(B, mesh.ref$edges[i,2])
    }
  else if( mesh.ref$edges[i,2] == 131 ){
    A = append(A, mesh.ref$edges[i,1])
    B = append(B, 131)
  }
}

cbind(A,B)
points( mesh.ref$nodes[A[1],1], mesh.ref$nodes[A[1],2], col="green",pch=16 )
points( mesh.ref$nodes[B[1],1], mesh.ref$nodes[B[1],2], col="red",pch=16 )

points( mesh.ref$nodes[A[2],1], mesh.ref$nodes[A[2],2], col="green",pch=16 )
points( mesh.ref$nodes[B[2],1], mesh.ref$nodes[B[2],2], col="red",pch=16 )
nrow(mesh.ref$nodes)
