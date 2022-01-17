### eval time 1.5D mesh test ###

nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh = create.mesh.1.5D(nodes,edges)
mesh = refine.mesh.1.5D(mesh=mesh,0.1)
nnodes = dim(mesh$nodes)[1]
plot(mesh)

FEMbasis= create.FEM.basis(mesh=mesh)

time = as.numeric(1:5)
coeff = rep(fs.test(mesh$nodes[,1], mesh$nodes[,2]),5)*time
## Create the FEM.time object
FEM_time_function = FEM.time(coeff=coeff, time_mesh=time, FEMbasis=FEMbasis, FLAG_PARABOLIC=TRUE)

N = 2
points_=matrix(nrow=N,ncol=2)
points_[,1] = runif(N,min=0.25,max=0.75)
points_[,2] = runif(N,min=0.25,max=0.5)
locations=projection.points.1.5D(mesh,points_)
points(locations[,1],locations[,2],col="green",pch=16)

evaluations = eval.FEM.time(FEM_time_function, locations = locations, 
                             time.instants = time)

## Evaluate the mean of the finite element function over the fifth triangle of the mesh
incidence_matrix = matrix(0, ncol = nrow(mesh$edges))
incidence_matrix[1,5] = 1

evaluations = eval.FEM.time(FEM_time_function, incidence_matrix = incidence_matrix, time.instants = time)

### plot FEM.time over graph domains ###
plot.FEM.time(FEM_time_function, time_locations = c(1.2,1.6,1.8))

plot.FEM.time(FEM_time_function, locations = locations)

x = FEM_time_function
N = nrow(x$FEMbasis$mesh$nodes)
time.locations =   
solution <- .Call("eval_FEM_time_nodes",N,x$mesh_time,time_locations = time.locations,x$coeff[,1,1],x$FLAG_PARABOLIC, PACKAGE = "fdaPDE")

plot(x, time_locations = c(1.2,1.4,1.6,1.8,2.0,2.5))


plot(x, locations = locations)

x = FEM_time_function
Nt = 50
t = rep(seq(min(x$mesh_time),max(x$mesh_time),length.out=Nt),times=nrow(locations))
eval_points = cbind(t, rep(locations[,1], each = Nt), rep(locations[,2], each = Nt))


if (!exists("col",inherits=FALSE) || length(col) < nrow(locations)){
  col = rainbow(nrow(locations))
  
}
if (!exists("type",inherits=FALSE)) {
  type = "l"
}
if (!exists("ylab",inherits=FALSE)) {
  ylab = ""
}
if (!exists("xlab",inherits=FALSE)) {
  xlab = "Time"
}
if (!exists("ylim",inherits=FALSE)) {
  ylim = range(eval_sol)
}

eval_sol = eval.FEM.time(FEM.time = x, space.time.locations = eval_points, lambdaT=1, lambdaS=1)
add=FALSE
pdf("prova.pdf")
for (ind in 1:nrow(locations)) {
  if (add == FALSE && ind == 1) {
    plot(t[1:Nt],eval_sol[(ind-1)*Nt+(1:Nt)],type=type,col=col,xlab=xlab,ylab=ylab,ylim=ylim)
  }else{
    points(t[1:Nt],eval_sol[(ind-1)*Nt+(1:Nt)],col=col,type=type)
  }
}

dev.off()


plot.loc.times = function(FEM.time){
  Nt = 50
  t = rep(seq(min(x$mesh_time),max(x$mesh_time),length.out=Nt),times=nrow(locations))
  eval_points = cbind(t, rep(locations[,1], each = Nt), rep(locations[,2], each = Nt))
  
  
  if (!exists("col",inherits=FALSE) || length(col) < nrow(locations)){
    col = rainbow(nrow(locations))
    
  }
  if (!exists("type",inherits=FALSE)) {
    type = "l"
  }
  if (!exists("ylab",inherits=FALSE)) {
    ylab = ""
  }
  if (!exists("xlab",inherits=FALSE)) {
    xlab = "Time"
  }
  if (!exists("ylim",inherits=FALSE)) {
    ylim = range(eval_sol)
  }
  
  eval_sol = eval.FEM.time(FEM.time = x, space.time.locations = eval_points, lambdaT=1, lambdaS=1)
  add=FALSE
  pdf("prova.pdf")
  for (ind in 1:nrow(locations)) {
    if (add == FALSE && ind == 1) {
      plot(t[1:Nt],eval_sol[(ind-1)*Nt+(1:Nt)],type=type,col=col[ind],xlab=xlab,ylab=ylab,ylim=ylim)
    }else{
      points(t[1:Nt],eval_sol[(ind-1)*Nt+(1:Nt)],col=col[ind],type=type)
    }
  }
  
  dev.off()
}
plot.loc.times(x)
