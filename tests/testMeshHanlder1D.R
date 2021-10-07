#Mesh Handler 1D test
nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh_ = create.mesh.1D(nodes,edges,order=2)
plot.mesh.1D(mesh_)
###########################################
nodes=matrix(nrow=3,ncol=2)
nodes[,1] = c(0.,0.5,1.)
nodes[,2] = c(0.,0.,0.)
edges=matrix(c(1,2,2,3),nrow = 2,byrow = TRUE)
mesh_ = create.mesh.1D(nodes,edges,order=1)
####################################
ndim = 2
mydim= 1
order =1

######################################
points_ =matrix(nrow=3,ncol=2)
points_[,1] = runif(3, min=0.,max=0.75)
points_[,2] = runif(3, min=0.,max=0.75)
######################################
points_=matrix(nrow=1,ncol=2)
points_[,1] = runif(1,min=0.25,max=0.5)
points_[,2] = runif(1,min=0.1,max=0.375)
######################################
points_=matrix(nrow=3,ncol=2)
points_[1,1] = 0.4 
points_[1,2] = 0.25 + 10*.Machine$double.eps * rnorm(1,0,1)
points_[2,1] = 0.5 + runif(1,min=0.,max=0.25) 
points_[2,2] = points_[2,1] -0.25 +   10*.Machine$double.eps * rnorm(1,0,1)
points_[3,1] = 0.5 + runif(1,min=0.,max=0.25) 
points_[3,2] = -points_[3,1] + 0.75 + 10*.Machine$double.eps * rnorm(1,0,1)



points( points_[,1],points_[,2],col='red')

# Indici in C++ partono da 0!
mesh_$edges = mesh_$edges - 1 

#storage.mode SEMPRE SEMPRE SEMPRE SEMPRE 
storage.mode(ndim) <- "integer"
storage.mode(mydim) <- "integer"
storage.mode(order) <- "integer"
storage.mode(mesh_$edges) <- "integer"
storage.mode(mesh_$nodes) <- "double"

result <- .Call("get_meshHandler", mesh_,order,ndim, mydim, points_ ,PACKAGE = "fdaPDE")
mesh_$edges = mesh_$edges + 1

num_proj = dim(points_)[1]
for (i in 1:num_proj){
  points( result[[2*i+1]][,1],result[[2*i+1]][,2], col='green' )
}

result

########### PROJECTION.h #############
nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh_ = create.mesh.1D(nodes,edges,order=1)
plot.mesh.1D(mesh_)

points_=matrix(nrow=5,ncol=2)
points_[,1] = runif(5,min=0.25,max=0.75)
points_[,2] = runif(5,min=0.25,max=0.5)
points( points_[,1],points_[,2],col='blue')

# Indici in C++ partono da 0!
mesh_$edges = mesh_$edges - 1 
mydim=1
ndim=2
#SEMPRE
storage.mode(mydim)<-"integer"
storage.mode(ndim)<-"integer"
storage.mode(mesh_$order) <- "integer"
storage.mode(mesh_$edges) <- "integer"
storage.mode(mesh_$nodes) <- "double"

res <- .Call("points_projection_1D", mesh_, points_,mydim,ndim)

num_proj = dim(points_)[1]
for (i in 1:num_proj){
  points( res[,1],res[,2], col='green' )
}


### mesh2.5D ###
nodes = matrix(nrow=4,ncol=3)
nodes[,1] = c(0.,0.5,0.5,0.75)
nodes[,2] = c(0.,0.,0.5,0.75)
nodes[,3] = c(0.,0.,0.,0.5)
nodes
triangles = matrix(nrow=2,ncol=3)
triangles[1,] = c(1,2,3)
triangles[2,] = c(2,3,4)
mesh_2.5D = create.mesh.2.5D(nodes=nodes,triangles=triangles)

points_ = matrix(nrow=5,ncol=3)
points_[,1] =runif(5,0.,0.75)
points_[,2] =runif(5,0.,0.75)
points_[,3] =runif(5,0.,0.75)
points_

options(rgl.printRglwidget=TRUE)
plot.mesh.2.5D(mesh_2.5D)
points3d(points_[,1],points_[,2],points_[,3],color='blue')
# Indici in C++ partono da 0!
mesh_2.5D$triangles = mesh_2.5D$triangles - 1 

mydim=2
ndim=3
#SEMPRE
storage.mode(mydim)<-"integer"
storage.mode(ndim)<-"integer"
storage.mode(mesh_2.5D$order) <- "integer"
storage.mode(mesh_2.5D$triangles) <- "integer"
storage.mode(mesh_2.5D$nodes) <- "double"


res <-  .Call("points_projection_1D", mesh_2.5D, points_,mydim,ndim)
res

num_proj = dim(points_)[1]
for (i in 1:num_proj){
  points3d( res[,1],res[,2],res[,3], color='green' )
}

