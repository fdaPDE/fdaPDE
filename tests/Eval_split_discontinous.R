# return an integer "vector" containing 1 if the i-th data is a vertex, 0 otherwise
# summing data is a numeric matrix
is.vertex <- function(Graph, data){

  data = as.matrix(data)
  vertex = Graph$nodes
  nvertex = nrow(Graph$nodes)
  ndata = nrow(data)
  ret = vector(mode="integer", length=ndata)
  for( i in 1:ndata ){
    ret[i] = 0
    for( j in 1: nvertex){
      if( sqrt( (data[i,1]-vertex[j,1])^2 + (data[i,2]-vertex[j,2])^2) < 10 * .Machine$double.eps) {
        ret[i] = j
        break
      }
    }
  }
  
  return(ret)
}
    
# assuming there are no cycles with total length less then 2*bandwidth, 
# i.e I am assuming that for each data point there is only one shortest path 
# connecting the source to this data point with total path length less then bandwidth
equal_split_discontinous <- function(Graph, sigma, dijkstra, x,y){
  
  source = dijkstra$source
  points_ = cbind(x,y)
  
  # index of the edge to which the point belong 
  # (nb. it is useful if i-th data points is not a vertex)
  idx = isInside(Graph, points_)
  
  # 0,            if points_[i] is not a vertex
  # vertex.idx,   if points_[i] is a vertex
  is_vertex = is.vertex(Graph, points_)
  
  Previous = vector(mode="integer", length(x))
  
  # filled with zeros by default
  bandwidth = vector(mode="numeric", length=length(x))
  Dist = vector(mode="numeric", length=length(x))
  delta = vector(mode="numeric", length=length(x))
  weight = vector(mode="numeric", length=length(x))
  
  # points that are vertexes
  idx.vertex = which( is_vertex != 0)
  Dist[ idx.vertex ] = dijkstra$distance[ is_vertex[idx.vertex] ]
  
  bandwidth[ (Dist[idx.vertex] < 5*sigma) ] = 1.0
  delta[idx.vertex] = Dist[idx.vertex]
  # Kernel Density Estimation on a Linear Network (Mc Swiggan, Baddeley, Nair)
  for( i in idx.vertex){
    weight[i] = 2./dijkstra$degree.nodes[ source ]
    current_path = dijkstra$ShortestPaths[[ is_vertex[i] ]]
    if( i != source && length(current_path)>1 ){
      for( prev in 2:length(current_path)){
        weight[i] = weight[i] * 1./( dijkstra$degree.nodes[ current_path[prev] ] - 1 )
      }
    }
  }
  
  # points that not are vertexes
  idx.not.vertex = which(is_vertex==0)

  if( length(idx.not.vertex)>0){
    
    Dist1 = dijkstra$distance[ Graph$edges[idx[idx.not.vertex], 1]]
    Dist2 = dijkstra$distance[ Graph$edges[idx[idx.not.vertex], 2]]
    
    for( i in 1:length(idx.not.vertex)){
      if( Dist1[i] < Dist2[i] ){
        Dist[idx.not.vertex[i]] = Dist1[i]
        Previous[idx.not.vertex[i]] = Graph$edges[idx[idx.not.vertex[i]], 1]
      }else{
        Dist[idx.not.vertex[i]] = Dist2[i]
        Previous[idx.not.vertex[i]] = Graph$edges[ idx[idx.not.vertex[i]] , 2]
      }
    current_length = Dist[idx.not.vertex[i]] + sqrt( (Graph$nodes[Previous[idx.not.vertex[i]], 1]-points_[idx.not.vertex[i], 1])^2 +
                                                     (Graph$nodes[Previous[idx.not.vertex[i]], 2]-points_[idx.not.vertex[i], 2])^2)
    if( current_length < 5*sigma ){
        bandwidth[idx.not.vertex[i]] = 1.  
    }else{
        bandwidth[idx.not.vertex[i]] = 0.0
    }
    
    delta[idx.not.vertex[i]] = current_length
    }
  
  for( i in idx.not.vertex){
      weight[i] = 2./dijkstra$degree.nodes[ source ]
    
      current_path = dijkstra$ShortestPaths[[ Previous[i] ]]
    if(length(current_path)>1){    
      for( prev in 2:length(current_path) ){
        weight[i] = weight[i] * 1./( dijkstra$degree.nodes[current_path[prev]]-1)     
      }
    
      weight[i] = weight[i] * 1./(dijkstra$degree.nodes[ Previous[i] ]-1)
    }
  }
  
  }
    coef = weight * 1/((2*pi)^0.5*sigma) * exp(-delta^2/(2*sigma^2)) * bandwidth  
  
    ret = list( coef = coef, bandwidth=bandwidth, distances=delta, weights = weight )
    return( ret )
    
}