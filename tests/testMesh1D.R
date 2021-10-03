# plot mesh 1D
plot.mesh.1D<-function(x, ...)
{
  
  if( x$order == 1 ){
    x11()
    plot(x$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
    segments(x$nodes[x$edges[,1],1], x$nodes[x$edges[,1],2],
           x$nodes[x$edges[,2],1], x$nodes[x$edges[,2],2], ...)
  }
  else{
    x11()
    plot(x$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
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


refinedMesh <- refine.by.splitting.mesh.1D(mesh = mesh1)
plot.mesh.1D(refinedMesh)

#Rifiniamo un'altra volta
refinedMesh2 <- refine.by.splitting.mesh.1D(mesh = refinedMesh)
nodes = refinedMesh2$nodes
edges = refinedMesh2$edges

plot.mesh.1D(refinedMesh2)
#OK
# ORDER 2
load("~/Scrivania/fdaPDE/data/ToyGraph1D_Order2.RData")
nodes = ToyGraphOrder2$nodes
edges = ToyGraphOrder2$edges

mesh2 = create.mesh.1D(nodes, edges, order = 2)
plot.mesh.1D(mesh2)

# mesh order2 da setting 1D  
load("~/Scrivania/fdaPDE/data/ToyGraph1D.RData")
nodes = ToyGraph$nodes
edges = ToyGraph$egdes
mesh2 = create.mesh.1D(nodes, edges, order = 2)

plot.mesh.1D(mesh2)

#refine.by.splitting order2 
mesh2_refined <- refine.by.splitting.mesh.1D(mesh2)
plot.mesh.1D(mesh2_refined)
