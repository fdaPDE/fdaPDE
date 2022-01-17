euclidean_distance <- function(point.A, point.B){
  
  return( sqrt( (point.A[1] - point.B[1])^2 + (point.A[2] - point.B[2])^2 ))
  
}


Dijkstra <- function( Graph, source)
{
  
  nnodes = nrow(Graph$nodes)
  NODES = 1 : nnodes
  neighbors = matrix(list(), nrow=nnodes, ncol=1)
  
  nedges = nrow(Graph$edges)
  for( edge in 1:nedges){
    left = Graph$edges[edge,1]
    right = Graph$edges[edge,2]
    neighbors[[ left ]] = append(neighbors[[ left ]], right)
    neighbors[[ right ]] = append(neighbors[[ right ]], left)
  }
  
  degree.node = vector(mode="integer", length=nnodes)
  for( node in 1:nnodes){
    degree.node[node] = length(neighbors[[node]])
  }
  distance = vector(mode="numeric", nnodes)
  ordering = order(distance)
  
  previous = vector(mode="integer", nnodes)
  
  distance[] = .Machine$double.xmax #INFINITY
  distance[source] = 0.
  previous[] = -1 #UNDEFINED
  
  while( length(NODES) ){
    current = -1
    
    for( e in 1:nnodes){
      if( is.element( ordering[e], NODES) ){
          current = ordering[e]
          break
      }
    }
   
    # remoing current from NODES
    NODES = NODES[- (which(NODES == current)) ]
   
     # computing distance for each neighbors of current nodes
    for( neigh in neighbors[[current]] ){
      
      dist = distance[current] + euclidean_distance( Graph$nodes[current,], Graph$nodes[neigh, ] )
      
      if( dist < distance[neigh]){
      distance[neigh] = dist
      previous[neigh] = current      
      }
    
    }
    
    ordering = order(distance)
  }
  
  
  ShortestPaths = matrix(list(), nrow=nnodes, ncol=1)
  for( target in 1:nnodes){
    
    e = target

    if( previous[e] != -1 || e == source){
      while( e != source ){
        ShortestPaths[[target]] = append( previous[e], ShortestPaths[[target]])
        e = previous[e]
      }
    }

  }
  return( list(source = source,distance = distance, previous = previous, ShortestPaths = ShortestPaths, degree.nodes = degree.node) )
}

