## test projection ##

## mesh.1D ##
nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh_ = create.mesh.1D(nodes,edges,order=2)
plot.mesh.1D(mesh_)

points_=matrix(nrow=5,ncol=2)
points_[,1] = runif(5,min=0.25,max=0.75)
points_[,2] = runif(5,min=0.25,max=0.5)
points( points_[,1],points_[,2],col='blue')

proj = projection.points.1D(mesh_,points_)

points(proj[,1],proj[,2], col = 'green')

## mesh2.5D ##

nodes = matrix(nrow=4,ncol=3)
nodes[,1] = c(0.,0.5,0.5,0.75)
nodes[,2] = c(0.,0.,0.5,0.75)
nodes[,3] = c(0.,0.,0.,0.35)

triangles = matrix(nrow=2,ncol=3)
triangles[1,] = c(1,2,3)
triangles[2,] = c(2,3,4)
mesh_2.5D = create.mesh.2.5D(nodes=nodes,triangles=triangles)

points_ = matrix(nrow=5,ncol=3)
points_[,1] =runif(5,0.,0.75)
points_[,2] =runif(5,0.,0.75)
points_[,3] =runif(5,0.,0.75)

res <- projection.points.2.5D(mesh_2.5D,points_)

#plot (rallenta)
options(rgl.printRglwidget=TRUE)
fdaPDE::plot.mesh.2.5D(mesh_2.5D)
rgl.points(res[,1],res[,2],res[,3],color='green')
rgl.points(points_[,1],points_[,2],points_[,3],color='blue')