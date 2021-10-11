# plot mesh 1D
plot_mesh_1D<-function(x, ...)
{
  
  if( x$order == 1 ){
    plot(x$nodes, ...)
    segments(x$nodes[x$edges[,1],1], x$nodes[x$edges[,1],2],
           x$nodes[x$edges[,2],1], x$nodes[x$edges[,2],2], ...)
  }
  else{
    plot(x$nodes, ...)
    segments(x$nodes[x$edges[,1],1], x$nodes[x$edges[,1],2],
             x$nodes[x$edges[,2],1], x$nodes[x$edges[,2],2], ...)
  points(x$nodes[x$edges[,3],1], x$nodes[x$edges[,3],2],col="red", ...)
  }
}

load("~/Scrivania/fdaPDE/data/ToyGraph1D.RData")
nodes = ToyGraph$nodes
edges = ToyGraph$egdes

mesh1 = create.mesh.1D(nodes, edges)
plot.mesh.1D(mesh1)

########################
#Split mesh1d ORDER 1
load("~/Scrivania/fdaPDE/data/ToyGraph1D.RData")
nodes = ToyGraph$nodes
edges = ToyGraph$egdes
mesh1 = create.mesh.1D(nodes, edges)
plot.mesh.1D(mesh1)
n_elem = dim(edges)[1]

#ordine1 ok...
for(i in 1:n_elem){
  print( mesh1$neighbors[[i,1]] )
  print( mesh1$neighbors[[i,2]] )
}

refinedMesh <- refine.by.splitting.mesh.1D(mesh = mesh1)
plot.mesh.1D(refinedMesh)
n_elem = dim(refinedMesh$edges)[1]
for(i in 1:n_elem){
  print( refinedMesh$neighbors[[i,1]] )
  print( refinedMesh$neighbors[[i,2]] )
}


#Rifiniamo un'altra volta
refinedMesh2 <- refine.by.splitting.mesh.1D(mesh = refinedMesh)
nodes = refinedMesh2$nodes
edges = refinedMesh2$edges

plot.mesh.1D(refinedMesh2)
#OK
# ORDER 2 (neighbors NON ok)
load("~/Scrivania/fdaPDE/data/ToyGraph1D_Order2.RData")
nodes = ToyGraphOrder2$nodes
edges = ToyGraphOrder2$edges

mesh2 = create.mesh.1D(nodes, edges, order = 2)
plot.mesh.1D(mesh2)

# mesh order2 da setting 1D  (neighbors ok)
load("~/Scrivania/fdaPDE/data/ToyGraph1D.RData")
nodes = ToyGraph$nodes
edges = ToyGraph$egdes
mesh2 = create.mesh.1D(nodes, edges, order = 2)

plot.mesh.1D(mesh2)

#refine.by.splitting order2 (anche qua pare neighbors ok!)
mesh2_refined <- refine.by.splitting.mesh.1D(mesh2)
plot.mesh.1D(mesh2_refined)

