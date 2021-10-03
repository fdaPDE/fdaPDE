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

##########################
