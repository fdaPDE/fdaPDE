
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

plot(mesh.ref)
x=mesh.ref
pdf("/home/aldo/Scrivania/prova.pdf")
for(i in 1:nrow(mesh.ref$edges) ){
  print(paste("EDGE ", i, sep=""))
  #Sys.sleep(2.5) 
  
  plot(x$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n",main=paste("EDGE ",i,sep=""))
  points(vertices, pch=16)
  points( x$nodes[ x$edges[i,1],1], x$nodes[x$edges[i,1],2], pch=16,col="red")
  points( x$nodes[ x$edges[i,2],1], x$nodes[x$edges[i,2],2], pch=16,col="blue")
  segments(x$nodes[x$edges[i,1],1], x$nodes[x$edges[i,1],2],
           x$nodes[x$edges[i,2],1], x$nodes[x$edges[i,2],2])

}
dev.off()

pdf("/home/aldo/Scrivania/provaR.pdf")
x = ret
for(i in 1:nrow(ret$EDGES) ){
  print(paste("EDGE ", i, sep=""))
  #Sys.sleep(2.5) 
  
  plot(x$NODES, xlab="", ylab="", xaxt="n", yaxt="n", bty="n",main=paste("EDGE ",i,sep=""))
  points(vertices, pch=16)
  points( x$NODES[ x$EDGES[i,1],1], x$NODES[x$EDGES[i,1],2], pch=16,col="red")
  points( x$NODES[ x$EDGES[i,2],1], x$NODES[x$EDGES[i,2],2], pch=16,col="blue")
  segments(x$NODES[x$EDGES[i,1],1], x$NODES[x$EDGES[i,1],2],
           x$NODES[x$EDGES[i,2],1], x$NODES[x$EDGES[i,2],2])
  
}
dev.off()


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
###############

refine1D.R <-function(nodes, edges_old, delta){
  
  # stores the length of the i-th old edge;
  lengths = vector(mode="numeric", nrow(edges_old))
  # stores the number of sub edges in which the i-th old edge is splitted into
  num_subs= vector(mode="integer", nrow(edges_old))
  
  num_new_nodes = 0;
  num_tot_edges = 0;
  for(i in 1 : nrow(edges_old) ){
    x0 = nodes[edges_old[i,1],1];
    y0 = nodes[edges_old[i,1],2];
    x1 = nodes[edges_old[i,2],1];
    y1 = nodes[edges_old[i,2],2];
    
    lengths[i] = sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) )
    
    # number of sub_edges
    if( lengths[i] > delta){
        num_subs[i] =  ceiling( lengths[i]/delta )
    }else{
      num_subs[i] = 1
    }
    
    # number of new nodes
    num_new_nodes = num_new_nodes + num_subs[i] - 1;
    # number of edges
    num_tot_edges = num_tot_edges + num_subs[i];
  }
  
  #result = PROTECT(Rf_allocVector(VECSXP,2));
  
  #SET_VECTOR_ELT(result, 1, Rf_allocMatrix(INTSXP, num_tot_edges, 2));
  #RIntegerMatrix edges( VECTOR_ELT(result,1));
  edges = matrix(0,nrow=num_tot_edges,ncol=2)
  
  #SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, num_new_nodes,2));
  nodes_new = matrix(0, nrow=num_new_nodes, ncol=2)
  
  total_nodes = nrow(nodes)
  edges_count = 0
  nodes_new_count = 0
  
  eps = 10 * .Machine$double.eps
  
  for( i in 1:nrow(edges_old) ){
    if( num_subs[i] == 1)
    {
      edges_count = edges_count + 1
      edges[edges_count,1] = edges_old[i,1]; # qua è l'errore di C++ ?
      edges[edges_count,2] = edges_old[i,2]; # NON Sommo 1 ?
    }
    else #// there is at least one new internal node!
      {
        x0 = nodes[edges_old[i,1],1]
        y0 = nodes[edges_old[i,1],2]
        x1 = nodes[edges_old[i,2],1]
        y1 = nodes[edges_old[i,2],2]
        
        cos_ = 0
        sin_ = 0
        
        if( abs(x1-x0) < eps && (y1-y0) > 0  ){
          cos_ =  0.
          sin_ =  1.
        }else if( abs( x1-x0) < eps && (y1-y0) < 0 ){
          cos_ =  0.
          sin_ = -1.
        }else if( abs( y1-y0)<eps && (x1-x0) > 0 ) {
          cos_ =  1.
          sin_ =  0.
        }else if( abs( y1-y0)<eps && (x1-x0) < 0 ){
          cos_ = -1.
          sin_ =  0.
        }else{
           slope = (y1-y0)/(x1-x0)
            if((x1-x0)>0.){
              cos_ = sqrt(1./(1.+slope*slope))
            }else{   
              cos_= -sqrt(1./(1.+slope*slope)) 
            }
           if((y1-y0)>0.){
             sin_ = sqrt(slope*slope/(1.+slope*slope))
           }else{
             sin_ = -sqrt(slope*slope/(1.+slope*slope))
           }
  
        }
        
         delta_ = lengths[i]/num_subs[i];
        # vector storing the global numbers of the new nodes which belong to old current edges
        # If N = num_subs[i] then there are N-1 internal nodes (the newest ones)
        nodes_global_num = vector(mode="integer", length=(num_subs[i] + 1))
        
        #Fixing the first and the last nodes
        nodes_global_num[1] = edges_old[i,1];
        nodes_global_num[ num_subs[i] + 1 ] = edges_old[i,2];
        
        # Computing new nodes coordinates
        for( n in 1:(num_subs[i]-1))
        {
          nodes_global_num[n + 1] = n + total_nodes;
          nodes_new[nodes_new_count + n, 1] = x0 + delta_ * cos_ * n;
          nodes_new[nodes_new_count + n, 2] = y0 + delta_ * sin_ * n;
          }
        
        nodes_new_count = nodes_new_count + num_subs[i] - 1
        total_nodes = total_nodes + num_subs[i] - 1;
        
        for( e in 1:num_subs[i])
        {
            #Indexes in R starts from 1, in C++ from 0, needed transformations!
            edges_count=edges_count+1
            edges[edges_count,1] =  nodes_global_num[e] ;
            edges[edges_count,2] = nodes_global_num[e+1] ;
  
        }
      }
    
  }
  
  NODES = rbind(nodes,nodes_new)
  ret = list( new_nodes = nodes_new, EDGES = edges, NODES=NODES)
  return(ret)
}
