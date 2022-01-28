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

is.inside.edge<-function(points_ ,xA, yA, xB, yB){
  res = vector(mode="logical", length=nrow(points_))
  
  # esclusi estremi
  idx = which( (points_[,1] > min(xA,xB)) &&
               (points_[,1] < max(xA,xB)) &&
               (points_[,2] > min(yA,yB)) &&
               (points_[,2] < max(yA,yB)) )
  
  for( i in idx){
    value = (points_[i,2] - yA)/(yB - yA) - (points_[i,1]- xA)/(xB - xA)
    
    res[i] = value < 10 * .Machine$double.eps
    
  }
  return(res)
}

# equal split discontinous allowing cycle.
equal_split_discontinous.cycle<-function( Graph, sigma, dijkstra, x, y){
  # bandwidth
  h = 5 * sigma
  
  source = dijkstra$source
  coef = vector(mode="numeric", length=length(x))
  
  points_ = cbind(x,y)
  
  # helper data structure 
  neighbors.vertexes = matrix(list(), nrow=nrow(Graph$nodes), ncol=1)
  
  nedges = nrow(Graph$edges)
  for( e in 1:nedges){
    left = Graph$edges[e,1]
    right = Graph$edges[e,2]
    
    neighbors.vertexes[[left]] = append(neighbors.vertexes[[left]], right)
    neighbors.vertexes[[right]] = append(neighbors.vertexes[[right]], left)
    
  }
  
  
  is_vertex = is.vertex(Graph,points_)
  idx = isInside(Graph, points_)
  
  # STEP 1: create an empty queue Q
  
  # Q queue is splitted into four vector to straightforward the implementation
  Q.start = c()
  Q.end = c()
  Q.distance = c()
  Q.weight = c()
  
  # STEP 2 : source is a vertex
  #          I assume that every kernel I consider is centered in a vertex, 
  #          I skip implementation of STEP 2
  
  # STEP 3 : source is a vertex
  
  m = dijkstra$degree.nodes[source]
  for( v in neighbors.vertexes[[source]]){
    Q.start = append(Q.start, source)
    Q.end = append(Q.end, v)
    Q.distance = append(Q.distance, dijkstra$distance[v])
    Q.weight = append(Q.weight, 2./m)
  }
  
  coef[source] = 2./m * 1./((2*pi)^0.5*sigma)
  idx.vertexes = vector(mode="integer", length=0)
  for( i in neighbors.vertexes[[source]] ){
    curr_idx = which(is_vertex == i)
    if( !is.empty(curr_idx))
      idx.vertexes = append(idx.vertexes, curr_idx)
  }
  
  for( i in idx.vertexes)
    coef[i] = 2./m * 1./((2*pi)^0.5*sigma) * exp(- dijkstra$distance[i]^2/(2*sigma^2) )
  
  # STEP 4: filling coef !
  
  while( !is.empty(Q.start)){
    
    start = Q.start[1]
    end = Q.end[1]
    dist = Q.distance[1]
    weight = Q.weight[1]
    
    # removing from queue
    Q.start = Q.start[-1]
    Q.end = Q.end[-1]
    Q.distance = Q.distance[-1]
    Q.weight = Q.weight[-1]
    
    m = dijkstra$degree.nodes[end]
    curr_weight = weight * 1./(m-1)
    for( neigh in neighbors.vertexes[[end]]){
      if(neigh != start){
        idx.no.vertex = which( is.inside.edge(points_, Graph$nodes[neigh,1],
                                                       Graph$nodes[neigh,2],
                                                       Graph$nodes[end,1], 
                                                       Graph$nodes[end,2]) == TRUE )
        
        curr_x = Graph$nodes[neigh,1]
        curr_y = Graph$nodes[neigh,2] 
        if(!is.empty(idx.no.vertex)){
          for( i in idx.no.vertex ){
            curr_dist = dist + sqrt( (curr_x - points_[i,1])^2 + (curr_y - points_[i,2])^2 )
            coef[i] = coef[i] + curr_weight * 1./((2*pi)^0.5*sigma) * exp(-curr_dist^2/(2*sigma^2))  
        
          }
        }
        
        idx.vertex = which( is_vertex == neigh )
        if( !is.empty(idx.vertex) ){
           
          curr_dist = dist + sqrt( (curr_x - Graph$nodes[end,1])^2 +
                                   (curr_y - Graph$nodes[end,2])^2 )
          coef[idx.vertex] = curr_weight * 1/((2*pi)^0.5*sigma) * exp(-curr_dist^2/(2*sigma^2) )+ coef[idx.vertex] # +  coef[idx.vertex] 
          
          if(curr_dist < h ){
            Q.start = append(Q.start, end)
            Q.end = append(Q.end, neigh)
            Q.distance = append(Q.distance, curr_dist) 
            Q.weight =  append(Q.weight, curr_weight)
            }
        
        }
      }
    }
    
  }
  return(coef)
  
}

# connecting the source to this data point with total path length less then bandwidth
equal_split_discontinous.mollifier<- function(Graph, sigma, dijkstra, x,y){
  
  source = dijkstra$source
  points_ = cbind(x,y)
  
  coef = vector(mode="numeric", length=length(x))
  
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
  
  # points that are vertexes
  idx.vertex = which( is_vertex != 0)
  Dist[ idx.vertex ] = dijkstra$distance[ is_vertex[idx.vertex] ]
  
  bandwidth[ (Dist[idx.vertex] < 5*sigma) ] = 1.0
  delta[idx.vertex] = Dist[idx.vertex]
  
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

    
  }
  coef[which(bandwidth == 1.0)] = 1 * (2 * exp( 1.0/( delta[which(bandwidth == 1.0)]^2/(5*sigma)^2  - 1) + 1) - 1  ) 
  
  ret = list( coef = coef, bandwidth=bandwidth, distances=delta)
  return( ret )
  
}




