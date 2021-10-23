### isInside ### isInside(mesh,points,search,redundacy) could be used to ensure that points belong to the domain!

nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh = create.mesh.1D(nodes,edges,order=1)
plot.mesh.1D(mesh)
points_=matrix(nrow=5,ncol=2)
points_[,1] = runif(5,min=0.25,max=0.75)
points_[,2] = runif(5,min=0.25,max=0.5)

locations=projection.points.1D(mesh,points_)
points(locations[,1],locations[,2], col = 'green')


res_points = isInside(mesh=mesh,points=points_, search ="naive", redundancy = TRUE)
res_proj   = isInside(mesh=mesh,points=locations,search ="naive",redundancy = TRUE)

### 2D horse-shoe

data(horseshoe2D)

mesh = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)

xeval=runif(1000,-0.5,0.5)
yeval=runif(1000,-0.5,0.5)
points_ = cbind(xeval,yeval)
res = isInside(mesh=mesh,points=points_,search="tree",redundancy = TRUE)

sum( res == FALSE)
dim(points_[res,] )

## 2.5 D

data(hub2.5D)
hub2.5D.nodes = hub2.5D$hub2.5D.nodes
hub2.5D.triangles = hub2.5D$hub2.5D.triangles

## Create mesh from nodes and connectivity matrix:
mesh = create.mesh.2.5D(nodes = hub2.5D.nodes, triangles = hub2.5D.triangles,order=2)

xeval = runif(1000, min=min(hub2.5D.nodes[,1]), max=max(hub2.5D.nodes[,1]))
yeval = runif(1000, min=min(hub2.5D.nodes[,1]), max=max(hub2.5D.nodes[,1]))
zeval   = runif(1000, min=0,max=1)
points_ = cbind(xeval,yeval,zeval)
locations_ = projection.points.2.5D(mesh,points_)

locs = matrix(nrow=1000,ncol=3)
locs[1:500,] =  points_[1:500,]
locs[501:1000,] = locations_[501:1000,]


res =isInside(mesh=mesh,points=locs,search="tree",redundancy = TRUE)
sum(res==FALSE)
