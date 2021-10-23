###mesh 1D ###
nodes=matrix(c(0.0,0.0,1.0,1.0),nrow=2,ncol=2,byrow=TRUE)
edges=matrix(c(1,2),nrow=1,ncol=2)
mesh=create.mesh.1.5D(nodes,edges,order=2)
plot.mesh.1.5D(mesh)
n = 7
for(i in 1:n){ 
  mesh= refine.by.splitting.mesh.1.5D(mesh)
}

#num_elem=2^n

point_ = matrix(c(0.99,0.99),nrow=1,ncol=2)

mesh$edges=mesh$edges-1
storage.mode(mesh$edges)<-"integer"
storage.mode(mesh$order)<-"integer"

myDim=1
nDim=2
storage.mode(myDim)<-"integer"
storage.mode(nDim)<-"integer"

# res[1] -> Naive Search
# res[2] -> Tree Search
res = .Call("TimingSearch", mesh,mesh$order,myDim,nDim,point_)
res*10^-6 # milliseconds

### mesh2D ###
data(horseshoe2D)
boundary_nodes = horseshoe2D$boundary_nodes
boundary_segments = horseshoe2D$boundary_segments
locations = horseshoe2D$locations
mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments,order = 2)
point_ = matrix(c(1, 0.5), ncol = 2)

mesh$edges=mesh$edges-1
mesh$triangles=mesh$triangles-1
storage.mode(mesh$edges)<-"integer"
storage.mode(mesh$triangles)<-"integer"
storage.mode(mesh$order)<-"integer"

myDim=2
nDim=2

storage.mode(myDim)<-"integer"
storage.mode(nDim)<-"integer"

# res[1] -> Naive Search
# res[2] -> Tree Search
res = .Call("TimingSearch", mesh,mesh$order,myDim,nDim,point_)
res*10^-6 # milliseconds

### mesh2.5D ###
data(hub2.5D)
hub2.5D.nodes = hub2.5D$hub2.5D.nodes
hub2.5D.triangles = hub2.5D$hub2.5D.triangles

## Create mesh from nodes and connectivity matrix:
mesh = create.mesh.2.5D(nodes = hub2.5D.nodes, triangles = hub2.5D.triangles,order=2)
points_ = matrix(nrow=1,ncol=3)
points_[,1:2] = runif(2, min=min(hub2.5D.nodes[,1]), max=max(hub2.5D.nodes[,1]))
points_[,3]   = runif(1, min=0,max=1)
locations = projection.points.2.5D(mesh,points_)

mesh$edges=mesh$edges-1
mesh$triangles=mesh$triangles-1
storage.mode(mesh$edges)<-"integer"
storage.mode(mesh$triangles)<-"integer"
storage.mode(mesh$order)<-"integer"

myDim=2
nDim=3

storage.mode(myDim)<-"integer"
storage.mode(nDim)<-"integer"

# res[1] -> Naive Search
# res[2] -> Tree Search
res = .Call("TimingSearch", mesh,mesh$order,myDim,nDim,locations)
res*10^-6 # milliseconds

### 3D MESH ###
data(sphere3Ddata)
nodes=sphere3Ddata$nodes
tetrahedrons=sphere3Ddata$tetrahedrons
mesh=create.mesh.3D(nodes,tetrahedrons,order=2)
point_ = matrix(nrow=1,ncol=3)
point_[,1:3] = runif(3, min=-0.9, max=0.9)

mesh$tetrahedrons = mesh$tetrahedrons - 1
mesh$faces = mesh$faces - 1
storage.mode(mesh$tetrahedrons) <- "integer"
storage.mode(mesh$faces) <- "integer"
storage.mode(mesh$order)<-"integer"

myDim=3
nDim=3

storage.mode(myDim)<-"integer"
storage.mode(nDim)<-"integer"

# res[1] -> Naive Search
# res[2] -> Tree Search
res = .Call("TimingSearch", mesh,mesh$order,myDim,nDim,point_)
res*10^-6 #milliseconds
