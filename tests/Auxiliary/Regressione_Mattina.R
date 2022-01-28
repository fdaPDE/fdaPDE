#########################################################################
#################   REGRESSIONE SPAZIALE   ##############################
#########################################################################
library(spatstat)
library(spatstat.linnet)


create.mesh.1D.vertices = function (vertices, edges, delta){
  
  #### D = vettore con le lunghezze di ciascun ramo.
  #### deltas = lunghezza dei sotto-segmenti in ciascun ramo.
  #### subsegments = numero di sotto-segmenti in ciascun ramo.
  edges = as.matrix(edges)
  n_edges = dim(edges)[1]
  D = as.matrix(0,n_edges,1)
  deltas = as.matrix(0,n_edges,1)
  subsegments = as.matrix(0,n_edges,1) 
  
  subel_edge = NULL
  
  for (i in 1:n_edges){
    
    current_edge = t(edges[i,]) # current_edge è un vettore colonna 2x1
    
    node_1 = t(vertices[current_edge[1],])
    node_2 = t(vertices[current_edge[2],])
    
    current_distance = sqrt ( (node_1[1] - node_2[1])^2 + (node_1[2] - node_2[2])^2 )
    
    D[i] = current_distance
    
    approx = round( D[i]/delta ) 
    
    # se la lunghezza del ramo diviso delta è zero, significa che quel ramo è più corto della 
    # lunghezza media che voglio tenere sul network.
    if(approx ==0) { 
      approx = 1
      subsegments[i]=1
      deltas[i] = D[i] 
    }
    else { 
      deltas[i] = D[i] / approx
      subsegments[i] = approx 
    }
    
    subel_edge = rbind(subel_edge, as.matrix(rep(i,approx)))
  }
  
  
  #### connections = lista che contiene, per ogni nodo, i numeri globali dei nodi con cui esso comunica.
  ####               la posizione del nodo in questa matrice definisce la sua numerazione globale.
  #### local_to_global = lista con n_edges elementi. l'elemento i-esimo si riferisce al ramo i-esimo.
  ####                   ogni elemento contiene il vettore dei numeri globali dei nodi contenuti 
  ####                   in tale ramo. ogni vettore considera anche i nodi vertice. 
  ####                   quindi ogni vettore comincia e finisce con un nodo vertice.
  ####                   di conseguenza, i nodi vertice possono essere ripetuti in nodi diversi.
  #### RAMO = vettore che in posizione i contiene il ramo di appartenenza del nodo globale i-esimo.
  #### nodes = matrice che contiene tutti i nodi della mesh.
  ####         la posizione del nodo in questa matrice definisce la sua numerazione globale.
  #### pL = point pattern sul linear network. si basa sulla libreria "spatstat".
  ####      fornisce le caratteristiche principali del network: numero di nodi, vertici, rami.
  connections = list() 
  local_to_global = list()
  
  N_interni = matrix(ncol = 2)
  RAMO = vector()
  count = 1
  
  number_v = dim(vertices)[1]
  
  
  for (i in 1:n_edges){
    
    current_length = deltas[i]
    current_edge = t(edges[i,]) # numeri globali che identificano i vertici estremi del ramo corrente
    node_1 = t(vertices[current_edge[1],]) # coordinate vertice 1
    node_2 = t(vertices[current_edge[2],]) # coordinate vertice 2
    
    local_to_global[[i]] = current_edge[1]
    
    if ( subsegments[i] == 1  ) {
      
      local_to_global[[i]] = c( local_to_global[[i]], current_edge[2] )
      
    }
    
    else {
      
      
      for (j in 1:(subsegments[i]-1) ){  # ciclando in questo intervallo ottengo solo i nodi interni (escludo i vertici)
        
        global_current = number_v + count
        
        local_to_global[[i]] = c( local_to_global[[i]], global_current )
        
        connections[[global_current]] = c(global_current - 1, global_current + 1)
        
        if ( j==1 ){
          connections[[global_current]] = c(current_edge[1], global_current + 1)
          connections[[current_edge[1]]] = c(connections[[current_edge[1]]], global_current)
        }
        
        if ( j==subsegments[i]-1 ){
          connections[[global_current]] = c(global_current - 1, current_edge[2])
          connections[[current_edge[2]]] = c(connections[[current_edge[2]]], global_current)
        }
        
        t = j* current_length / D[i] # coordinata locale 
        current_node = t* node_2 + (1-t)* node_1 
        N_interni = rbind(N_interni, current_node)
        RAMO[count] = i
        count = count+1
      }
      local_to_global[[i]] = c( local_to_global[[i]], current_edge[2] )
    } # graffa dell'if 
  }
  
  N_interni = N_interni[-1,]
  N_interni = as.matrix(N_interni)   # fondamentale, altrimenti non funziona
  
  nodes = rbind(vertices, N_interni)
  n_nodes = dim(nodes)[1]
  
  # pL=lpp(nodes, L)
  
  n_segments = sum(subsegments)
  
  #### costruisco l'analogo della struttura "triangles" nel 2.5D.
  #### la chiamo "elements" e contiene i numeri globali dei nodi che costituiscono ciascun sotto-segmento.
  
  elements = matrix(0, nrow = n_segments, ncol = 2)
  count = 1
  
  for ( i in 1:n_edges ){
    
    v = local_to_global[[i]]
    v = as.matrix(v)
    dim_v = dim(v)[1]
    
    for ( j in 1:(dim_v-1) ){
      elements[count,] = c(v[j],v[j+1])
      count = count +1  }
  }
  
  # ora devo costruire l'oggetto mesh da restituire come output.
  
  output = list(nnodes=n_nodes, nsegments=n_segments, nodes=as.matrix(nodes),
                segments = as.matrix(elements), 
                local_to_global=local_to_global, deltas=deltas,
                subel_edge=subel_edge, D=D)
  
  return(output)
}



###############################################################################

create.FEM.basis.1D = function(mesh) {
  
  nbasis = dim(mesh$nodes)[[1]]
  nodes = mesh$nodes
  segments = mesh$segments
  
  transf_coord = NULL 
  transf_coord$diffx = nodes[segments[,2],1] - nodes[segments[,1],1]
  transf_coord$diffy = nodes[segments[,2],2] - nodes[segments[,1],2]
  
  len = sqrt ( transf_coord$diffx^2 + transf_coord$diffy^2  )
  
  FEMbasis = list(mesh = mesh, nbasis = nbasis, len=len, transf_coord = transf_coord)
  
  return (FEMbasis)
}


###############################################################################

R_mass_1D = function(FEMbasis) {
  
  nodes = FEMbasis$mesh$nodes
  segments = FEMbasis$mesh$segments
  len = FEMbasis$len
  
  nsegments  = dim(segments)[1]
  nnodes  = dim(nodes)[1]
  
  K0 = matrix(0,nrow=nnodes,ncol=nnodes)
  
  K0M = matrix( c( 2,  1,
                   1,  2), ncol=2, nrow=2, byrow=T)/6
  
  for ( i in 1:nsegments ){
    
    ind = segments[i,]
    K0[ind,ind] = K0[ind,ind] + K0M * len[i]
    
  }
  return(K0)
}

##############################################################################

find.boundary = function (m){
  
  return(which(apply(m,1,sum)==1))
  
}

##############################################################################

R_stiff_1D = function(FEMbasis, m) {
  
  nodes = FEMbasis$mesh$nodes
  segments = FEMbasis$mesh$segments
  len = FEMbasis$len
  
  nsegments  = dim(segments)[1]
  nnodes  = dim(nodes)[1]
  
  K1 = matrix(0,nrow=nnodes,ncol=nnodes)
  
  K1M = matrix( c( 1,  -1,
                   -1,  1), ncol=2, nrow=2, byrow=T)
  
  for ( i in 1:nsegments ){
    
    ind = segments[i,]
    K1[ind,ind] = K1[ind,ind] + K1M / len[i]
    
  }
  
  
  
  return (K1)
}

###########################################################################

isBetween = function(a,b,c,epsilon) {
  
  flag = 1
  
  # faccio il prodotto vettoriale tra i vettori ac e ab (vedi quaderno)
  crossproduct = (c[1] - a[1]) * (b[2] - a[2]) - (c[2] - a[2]) * (b[1] - a[1])
  
  if ( abs(crossproduct) > epsilon ){
    flag = 0
  }
  
  dotproduct = (c[1] - a[1]) * (b[1] - a[1]) + (c[2] - a[2])*(b[2] - a[2])
  
  if ( dotproduct < 0 ){
    flag = 0
  }
  
  squaredlengthba = (b[1] - a[1])*(b[1] - a[1]) + (b[2] - a[2])*(b[2] - a[2])
  
  if ( dotproduct > squaredlengthba ){
    flag = 0
  }
  
  return(flag)
}

################################################################################

goLocal = function(vertex_1, vertex_2, my_obs){
  
  len = sqrt ( (vertex_1[2]-vertex_2[2])^2 + (vertex_1[1]-vertex_2[1])^2   )
  
  if ( vertex_2[1] == vertex_1[1] ){
    t = (my_obs[2] - vertex_1[2])/(vertex_2[2] - vertex_1[2])
  }
  else{
    t = (my_obs[1] - vertex_1[1])/(vertex_2[1] - vertex_1[1])
  }
  return(t*len)
}

################################################################################


eval.FEM.1D = function(mesh, edges, vertices, obs, epsilon){
  
  edges = as.matrix(edges)
  n_edges = dim(edges)[1]
  
  deltas = mesh$deltas
  local_to_global = mesh$local_to_global
  nnodes = mesh$nnodes
  
  obs = as.matrix(obs)
  n_obs = dim(obs)[1]
  
  local_coordinate = matrix(0, nrow=nobs, ncol=1)
  
  RAMO_OBS = matrix(0,nrow = n_obs, ncol = 1)
  
  PHI = matrix(0, nrow = n_obs, ncol=nnodes)
  
  for (i in 1:n_edges){
    
    current_length = deltas[i]
    current_edge = t(edges[i,])
    node_1 = t(vertices[current_edge[1],])
    node_2 = t(vertices[current_edge[2],])
    
    for (j in 1:n_obs){
      current_flag = isBetween(node_1, node_2, obs[j,], epsilon)
      if (current_flag==1){ 
        RAMO_OBS[j]=i 
        local_coordinate[j] = goLocal(node_1, node_2, obs[j,])  # coordinata locale t
        inf_local = floor(local_coordinate[j]/current_length)  # numero locale del nodo della mesh inferiore
        sup_local = ceiling(local_coordinate[j]/current_length) # numero locale del nodo della mesh superiore
        inf_global = local_to_global[[i]][inf_local+1]
        sup_global = local_to_global[[i]][sup_local+1]
        
        # a questo punto so a che sottosegmento del ramo i appartiene l'osservazione corrente.
        # in particolare, conosco i nodi globali tra cui è compresa. 
        # se si tratta di un nodo interno, andrò a modificare solo le phi relative a questi due nodi. 
        
        
        PHI[j,inf_global] = (- local_coordinate[j] + sup_local*current_length)/current_length
        PHI[j,sup_global] = ( - inf_local*current_length + local_coordinate[j])/current_length
        
        
      }
    }
  }
  output = list( PHI=PHI, RAMO_OBS=RAMO_OBS, local_coordinate=local_coordinate )
  return (output)
}

###########################################################################################


covariate.FEM.1D = function( W, Z, FEMbasis, PHI, lambda ){
  
  W = as.matrix(W)
  nobs = dim(W)[1]
  
  WW_inv = solve(t(W)%*%W)
  P = W%*% (WW_inv%*%t(W))
  Q = diag(nobs) - P
  
  nnodes = mesh$nnodes
  
  L_mat = t(PHI)%*%Q%*%PHI
  
  R0 = R_mass_1D(FEMbasis)
  R1 = R_stiff_1D(FEMbasis)
  
  row1 = cbind(-L_mat, lambda*t(R1))
  row2 = cbind(-lambda*R1, -lambda*R0)
  
  A = rbind(row1, row2)
  
  b1 = - t(PHI) %*% ( Q%*%Z )
  b2 = matrix(0, nrow = nnodes, ncol=1)
  
  b = rbind(b1,b2)
  
  sol = solve(A, b)
  
  output = sol[1:nnodes,]
  
  
  f_estimate =  PHI%*%output
  beta_estimate = WW_inv%*% t(W) %*% (Z-f_estimate)
  
  output = list ( beta_estimate = beta_estimate, f_coeff = output )
  
  return(output)
  
}


########################################################################################
## aggiunta ##

no.covariate.FEM.1D = function( Z, FEMbasis, PHI, lambda ){
  
  #W = as.matrix(W)
  #nobs = dim(W)[1]
  
  if(is.vector(Z))
    nobs = length(Z)
  else if( is.matrix(Z))
    nobs = dim(Z)[1]
  
  #WW_inv = solve(t(W)%*%W)
  #P = W%*% (WW_inv%*%t(W))
  #Q = diag(nobs) - P
  
  Q = diag(nobs)
  nnodes = mesh$nnodes

  #L_mat = t(PHI)%*%Q%*%PHI
  L_mat = t(PHI)%*%PHI
  
  R0 = R_mass_1D(FEMbasis)
  R1 = R_stiff_1D(FEMbasis)
  
  row1 = cbind(L_mat, -lambda*t(R1))
  row2 = cbind(-lambda*R1, -lambda*R0)
  
  A = rbind(row1, row2)
  
  b1 =  t(PHI) %*% ( Q%*%Z )
  b2 = matrix(0, nrow = nnodes, ncol=1)
  
  b = rbind(b1,b2)
  
  sol = solve(A, b)
  
  output = sol[1:nnodes,]
  
  #f_estimate =  PHI%*%output
  #beta_estimate = WW_inv%*% t(W) %*% (Z-f_estimate)
  
  output = list ( beta_estimate = NULL, f_coeff = output )
  
  return(output)
  
}
