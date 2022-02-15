source("~/Scrivania/fdaPDE/tests/Auxiliary/Regressione_Mattina.R")

source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")
data("spiders")
vertices = cbind(spiders$domain$vertices$x, spiders$domain$vertices$y)
edges = cbind(spiders$domain$from, spiders$domain$to)

M = spiders$domai$m
find.boundary(M)
L = as.linnet(spiders$domain)

#412 nodes -> delta = 40
delta =40
mesh = create.mesh.1D.vertices(vertices, edges, delta)
mesh$nnodes
FEMbasis = create.FEM.basis.1D(mesh)

### fdaPDE - MESH ### 
mesh.fdaPDE = create.mesh.1.5D(mesh$nodes,mesh$segments)
FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)

### CAMPO f ### 
plot(mesh.fdaPDE)
dijkstra.55 = Dijkstra(mesh.fdaPDE,55)
dijkstra.66 = Dijkstra(mesh.fdaPDE,66)
dijkstra.98 = Dijkstra(mesh.fdaPDE, 98)

# sym gaussian 
aux.66 = function(x, y, seg, tp) { 
  
  sigma.66 = 500
  
  res.66 = equal_split_discontinous(mesh.fdaPDE, sigma.66, dijkstra.66, x, y)
  
  res =  res.66$coef 
  return(res)
}
aux.55 = function(x, y, seg, tp) { 
  
  sigma.55 = 500
  
  res.55= equal_split_discontinous(mesh.fdaPDE, sigma.55, dijkstra.55, x, y)
  
  res =  res.55$coef 
  return(res)
}
aux.98 = function(x, y, seg, tp) { 
  
  sigma.98 = 500
  
  res.98= equal_split_discontinous(mesh.fdaPDE, sigma.98, dijkstra.98, x, y)
  
  res =  res.98$coef 
  return(res)
}

########## a-sym gaussian ################
aux.a.sym = function(x,y,tp,seg){
  sigma = 50
  h = 5*sigma
  Graph = mesh.fdaPDE
  dijkstra = dijkstra.46
  source = dijkstra$source
  points_ = cbind(x,y)
  
  coef = vector(mode="numeric", length=nrow(points_))
  is_vertex = is.vertex(Graph, points_)
  idx.vertex = which(is_vertex!=0)
  idx.same.y = which(points_[idx.vertex,2] == Graph$nodes[source,2])
  
  if(!is.empty(idx.same.y)){
    for( i in idx.same.y ){  
      if( dijkstra$distance[idx.vertex[i]] < h ){
      coef[i] = 1./((2*pi)^0.5*sigma) * exp(-dijkstra$distance[idx.vertex[i]]^2/(2*sigma^2))
      }
    }
  }
  
  idx.not.vertex = which(is_vertex==0)
  if(!is.empty(idx.not.vertex)){
    idx.edge = isInside(Graph, points_)
    idx.edge = idx.edge[idx.not.vertex]
    Dist1 = dijkstra$distance[ Graph$edges[idx.edge, 1]]
    Dist2 = dijkstra$distance[ Graph$edges[idx.edge, 2]]
    Dist = vector(mode="numeric", length=length(idx.not.vertex))
    Previous = vector(mode="integer", length(idx.not.vertex))
  
  for( i in 1:length(idx.not.vertex)){
    if( abs( points_[idx.not.vertex[i],2] - Graph$nodes[source,2]) < 10 * .Machine$double.eps){
      
      if( Dist1[i] < Dist2[i] ){
        Dist[i] = Dist1[i]
        Previous[i] = Graph$edges[idx.edge[i], 1]
      }else{
        Dist[i] = Dist2[i]
        Previous[i] = Graph$edges[idx.edge[i], 2]
      }
  
      delta = Dist[i] + sqrt( (Graph$nodes[Previous[i] ,1] - points_[idx.not.vertex[i],1])^2 +
                           (Graph$nodes[Previous[i], 2] - points_[idx.not.vertex[i],2])^2) 
  
      if( delta < h)
        coef[idx.not.vertex[i]] = 1/((2*pi)^0.5 * sigma) * exp(-delta^2/(2*sigma^2))  
  
    }
  }
  }
  return (coef)
}
a.sym.dens = linfun(aux.a.sym, L)
plot(a.sym.dens)
##########################################

# mixture of two gaussian distribution
AUX = function(x,y,seg,tp){
  sym = aux.55(x,y,seg,tp)
  sym.2 = aux.66(x,y,seg,tp)
  sym.3 = aux.98(x,y,seg,tp)
  res = sym + sym.2 + sym.3
  return(res)                    
}

coef.ex = AUX(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2])

R_plot_graph.ggplot2(FEM(coef.ex,FEMbasis.fdaPDE))
density = linfun(AUX,L)

plot(mesh.fdaPDE)
points(mesh.fdaPDE$nodes[46,1], mesh.fdaPDE$nodes[46,2],pch=16, col="red")
points(mesh.fdaPDE$nodes[283,1], mesh.fdaPDE$nodes[283,2],pch=16,col="red")

nobs = c(300,500,700,1000)
N=length(nobs)
M=20

mise = matrix(0,nrow=M,ncol=N)
mise.fdaPDE = matrix(0,nrow=M,ncol=N)
# (||f_{fdaPDE} - f{R} ||_{L2} )^2
err.L2 = matrix(0,nrow=M,ncol=N)

beta_ex = 0.8
rmse.beta = matrix(0,nrow=M,ncol=N)
rmse.fdaPDE.beta = matrix(0,nrow=M,ncol=N)
norms.beta.inf = matrix(0,nrow=M,ncol=N)
norms.beta.2 =matrix(0,nrow=M,ncol=N)

sols = array(0,dim=c(M,N,mesh$nnodes) ) #contiene vettore dei coeff
sols.fdaPDE = array(0,dim=c(M,N,mesh$nnodes) )

times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)
times.fdaPDE.2 = matrix(0,nrow=M,ncol=N)

options(digits.sec=6)
lambda_vector = seq(from=1e-1, to=1000, length.out=10)
lambda = 1000
epsilon=0.0001
for( i in 1:N){
  for( j in 1:M){
    PPP = runiflpp(nobs[i], L)
    points_x = PPP$data$x
    points_y = PPP$data$y
    points = cbind(points_x, points_y)
    
    W = matrix( rnorm(nobs[i], 0, sd=0.5), nrow=nobs[i], ncol=1)
    exact_sol = density(points_x, points_y)
    temp = range(  exact_sol)
    my_range = temp[2]-temp[1]
    EPS = rnorm(nobs[i], 0, sd=0.05*my_range)
    Z = beta_ex*W + exact_sol + EPS
    
    # output.fdaPDE = smooth.FEM( observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W,
    #                             locations=points,
    #                             lambda=lambda_vector, 
    #                             lambda.selection.criterion='grid', 
    #                             DOF.evaluation='exact', 
    #                             lambda.selection.lossfunction='GCV' )
    # 
    # lambda = output.fdaPDE$optimization$lambda_solution
    
    start.fdaPDE = Sys.time()
    output.fdaPDE = smooth.FEM(observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W,
                               locations=points,
                               lambda=lambda)
    end.fdaPDE = Sys.time()
    
    start.time <- Sys.time()
    result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    PHI = result$PHI
    output = covariate.FEM.1D(W, Z, FEMbasis, PHI, lambda)
    end.time  <-Sys.time()
    
    rmse.beta[j,i] = norm(beta_ex[1]- output$beta_estimate[1],type="2")
    rmse.fdaPDE.beta[j,i] = norm(beta_ex[1]- output.fdaPDE$solution$beta[1],type="2")
    
    norms.beta.inf[j,i] = abs( output$beta_estimate[1] - output.fdaPDE$solution$beta[1])
    norms.beta.2[j,i]   = norm(output$beta_estimate[1] - output.fdaPDE$solution$beta[1],type="2")
    
    times[j,i] = difftime(end.time, start.time, units="secs")
    
    times.fdaPDE[j,i] = output.fdaPDE$time
    times.fdaPDE.2[j,i] = difftime(end.fdaPDE, start.fdaPDE, units="secs")
    coef = output$f_coeff
    coef.fdaPDE     = output.fdaPDE$fit.FEM$coeff
    
    sols[j,i,] =  coef
    sols.fdaPDE[j,i,] = coef.fdaPDE
    
    mise[j,i]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
    mise.fdaPDE[j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
    err.L2[j,i] = integrate_f( FEM((coef.fdaPDE-coef)^2,FEMbasis.fdaPDE) )
  }
}

today_ = Sys.Date()
ntest_ = "-test-2"  

file.name = paste("Regression-",today_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

save(nobs,M,N,FEMbasis.fdaPDE,L,
     rmse.beta,
     rmse.fdaPDE.beta,
     mise, mise.fdaPDE, err.L2,
     sols,sols.fdaPDE,
     norms.beta.inf,norms.beta.2,
     times,times.fdaPDE,times.fdaPDE.2, file =save.file)

#####################################################################################
#####################################################################################

source("~/Scrivania/fdaPDE/tests/Auxiliary/Regressione_Mattina.R")

source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")
data("spiders")
vertices = cbind(spiders$domain$vertices$x, spiders$domain$vertices$y)
edges = cbind(spiders$domain$from, spiders$domain$to)

M = spiders$domai$m
find.boundary(M)
L = as.linnet(spiders$domain)

#412 nodes -> delta = 40
delta =40
mesh = create.mesh.1D.vertices(vertices, edges, delta)
mesh$nnodes
FEMbasis = create.FEM.basis.1D(mesh)

### fdaPDE - MESH ### 
mesh.fdaPDE = create.mesh.1.5D(mesh$nodes,mesh$segments)
FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)

dijkstra.55 = Dijkstra(mesh.fdaPDE,55)
dijkstra.66= Dijkstra(mesh.fdaPDE,66)
dijkstra.98 = Dijkstra(mesh.fdaPDE, 98)

# sym gaussian 
aux.55 = function(x, y, seg, tp) { 
  
  sigma.55 = 500
  
  res.55 = equal_split_discontinous(mesh.fdaPDE, sigma.55, dijkstra.55, x, y)
  
  res =  res.55$coef 
  return(res)
}

aux.66 = function(x, y, seg, tp) { 
  
  sigma.66 = 500
  
  res.66= equal_split_discontinous(mesh.fdaPDE, sigma.66, dijkstra.66, x, y)
  
  res =  res.66$coef 
  return(res)
}

aux.98 = function(x, y, seg, tp) { 
  
  sigma.98 = 500
  
  res.98= equal_split_discontinous(mesh.fdaPDE, sigma.98, dijkstra.98, x, y)
  
  res =  res.98$coef 
  return(res)
}

# mixture of three gaussian distributions
AUX = function(x,y,seg,tp){
  sym = aux.98(x,y,seg,tp)
  sym.2 = aux.55(x,y,seg,tp)
  sym.3 = aux.66(x,y,seg,tp)
  res = (sym + sym.2 + sym.3)
  return(res)                    
}

coef.ex = AUX(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2])

nobs = c(300,500,700,1000)
N=length(nobs)
M=20

mise = matrix(0,nrow=M,ncol=N)
mise.fdaPDE = matrix(0,nrow=M,ncol=N)
# (||f_{fdaPDE} - f{R} ||_{L2} )^2
err.L2 = matrix(0,nrow=M,ncol=N)

sols = array(0,dim=c(M,N,mesh$nnodes) ) #contiene vettore dei coeff
sols.fdaPDE = array(0,dim=c(M,N,mesh$nnodes) )

times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)
times.fdaPDE.2 = matrix(0,nrow=M,ncol=N)

options(digits.sec=6)
lambda_vector = seq(from=1e-1, to=1000, length.out=10)
lambda =1000
epsilon=0.0001
for( i in 1:N){
  for( j in 1:M){
    PPP = runiflpp(nobs[i], L)
    points_x = PPP$data$x
    points_y = PPP$data$y
    points = cbind(points_x, points_y)
    
    exact_sol = AUX(points_x, points_y)
    temp = range(exact_sol)
    my_range = temp[2]-temp[1]
    EPS = rnorm(nobs[i], 0, sd=0.05*my_range )
    Z =  exact_sol + EPS
    # 
    # output.fdaPDE = smooth.FEM( observations=Z, FEMbasis=FEMbasis.fdaPDE, 
    #                             covariates=NULL,
    #                             locations=points,
    #                             lambda=lambda_vector, 
    #                             lambda.selection.criterion='grid', 
    #                             DOF.evaluation='exact', 
    #                             lambda.selection.lossfunction='GCV' )
    # 
    #lambda = output.fdaPDE$optimization$lambda_solution
    
    start.fdaPDE = Sys.time()
    output.fdaPDE = smooth.FEM(observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=NULL,
                               locations=points,
                               lambda=lambda)
    end.fdaPDE = Sys.time()
    
    start.time <- Sys.time()
    result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    PHI = result$PHI
    output = no.covariate.FEM.1D(Z, FEMbasis, PHI, lambda)
    end.time  <-Sys.time()
    
    times[j,i] = difftime(end.time, start.time, units="secs")
    
    times.fdaPDE[j,i] = output.fdaPDE$time
    times.fdaPDE.2[j,i] = difftime(end.fdaPDE, start.fdaPDE, units="secs")
    coef = output$f_coeff
    coef.fdaPDE     = output.fdaPDE$fit.FEM$coeff
    
    sols[j,i,] =  coef
    sols.fdaPDE[j,i,] = coef.fdaPDE
    
    mise[j,i]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
    mise.fdaPDE[j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
    err.L2[j,i] = integrate_f( FEM((coef.fdaPDE-coef)^2,FEMbasis.fdaPDE) )
  }
}

today_ = Sys.Date()
ntest_ = "-test-2-no-cov"  

file.name = paste("Regression-",today_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

save(nobs,M,N,
     mise, mise.fdaPDE, err.L2,
     sols,sols.fdaPDE,
     times,times.fdaPDE,times.fdaPDE.2, file =save.file)


R0 = R_mass_1D(FEMbasis)
R1 = R_stiff_1D(FEMbasis)

#################################################################
################## increasing nodes number ######################
#################################################################

source("~/Scrivania/fdaPDE/tests/Auxiliary/Regressione_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")

delta = c(40, 30, 25, 17.5)

nobs =500

N=length(delta)
M=20

mise = matrix(0,nrow=M,ncol=N)
mise.fdaPDE = matrix(0,nrow=M,ncol=N)
# (||f_{fdaPDE} - f{R} ||_{L2} )^2
err.L2 = matrix(0,nrow=M,ncol=N)

beta_ex = 0.8
rmse.beta = matrix(0,nrow=M,ncol=N)
rmse.fdaPDE.beta = matrix(0,nrow=M,ncol=N)
norms.beta.inf = matrix(0,nrow=M,ncol=N)
norms.beta.2 =matrix(0,nrow=M,ncol=N)

times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)
times.fdaPDE.2 = matrix(0,nrow=M,ncol=N)

options(digits.sec=6)
lambda_vector = seq(from=1e-1, to=1000, length.out=10)
epsilon=0.0001

lambda = 1000
for( i in 1:N){
  
  mesh = create.mesh.1D.vertices(vertices, edges, delta[i])
  FEMbasis = create.FEM.basis.1D(mesh)
  nnodes[i] = mesh$nnodes
  
  ### fdaPDE - MESH ### 
  mesh.fdaPDE = create.mesh.1.5D(mesh$nodes,mesh$segments)
  FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)
  dijkstra.283 = Dijkstra(mesh.fdaPDE,283)
  dijkstra.118 = Dijkstra(mesh.fdaPDE,118)
  dijkstra.250 = Dijkstra(mesh.fdaPDE, 250)
  
  # sym gaussian 
  aux.283 = function(x, y, seg, tp) { 
    
    sigma.283 = 500
    
    res.283 = equal_split_discontinous(mesh.fdaPDE, sigma.283, dijkstra.283, x, y)
    
    res =  res.283$coef 
    return(res)
  }
  aux.250 = function(x, y, seg, tp) { 
    
    sigma.250 = 500
    
    res.250= equal_split_discontinous(mesh.fdaPDE, sigma.250, dijkstra.250, x, y)
    
    res =  res.250$coef 
    return(res)
  }
  aux.118 = function(x, y, seg, tp) { 
    
    sigma.118 = 500
    
    res.118= equal_split_discontinous(mesh.fdaPDE, sigma.118, dijkstra.118, x, y)
    
    res =  res.118$coef 
    return(res)
  }
  
  AUX = function(x,y,seg,tp){
    sym = aux.118(x,y,seg,tp)
    sym.2 = aux.250(x,y,seg,tp)
    sym.3 = aux.283(x,y,seg,tp)
    res = (sym + sym.2 + sym.3)
    return(res)                    
  }
  
  density <- linfun(AUX, L)
  coef.ex = AUX(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])
  
  for( j in 1:M){
    PPP = runiflpp(nobs, L)
    points = cbind(PPP$data$x, PPP$data$y)
    
    W = matrix( rnorm(nobs, 0, sd=0.5), nrow=nobs, ncol=1)
    exact_sol = density(points[,1], points[,2])
    temp = range(exact_sol)
    my_range = temp[2]-temp[1]
    EPS = rnorm(nobs, 0, sd=0.05*my_range )
    Z = beta_ex*W + exact_sol + EPS
    
   # output.fdaPDE = smooth.FEM( observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W,
   #                            locations=points,
   #                             lambda=lambda_vector, 
   #                             lambda.selection.criterion='grid', 
   #                             DOF.evaluation='exact', 
   #                             lambda.selection.lossfunction='GCV' )
   #lambda = output.fdaPDE$optimization$lambda_solution
    
    output.fdaPDE = smooth.FEM(observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W,
                               locations=points,
                               lambda=lambda)
    
    #result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    #PHI = result$PHI
    
    start.time <- Sys.time()
    result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    PHI = result$PHI
    output = covariate.FEM.1D(W, Z, FEMbasis, PHI, lambda)
    end.time  <-Sys.time()
    
    times[j,i] = difftime(end.time, start.time, units="secs")
    
    times.fdaPDE[j,i] = output.fdaPDE$time
    times.fdaPDE.2[j,i] = difftime(end.fdaPDE, start.fdaPDE, units="secs")
    coef = output$f_coeff
    coef.fdaPDE     = output.fdaPDE$fit.FEM$coeff
    
    rmse.beta[j,i] = norm(beta_ex- output$beta_estimate[1],type="2")
    rmse.fdaPDE.beta[j,i] = norm(beta_ex- output.fdaPDE$solution$beta[1],type="2")
    
    norms.beta.inf[j,i] = abs( output$beta_estimate[1] - output.fdaPDE$solution$beta[1])
    norms.beta.2[j,i]   = norm(output$beta_estimate[1] - output.fdaPDE$solution$beta[1],type="2")
    
    mise[j,i]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
    mise.fdaPDE[j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
    err.L2[j,i] = integrate_f( FEM((coef.fdaPDE-coef)^2,FEMbasis.fdaPDE) )
  }
}

today_ = Sys.Date()
ntest_ = "-test-2-delta-var" # "-test.1.3" 
#today_ = "2021-11-15"
file.name = paste("Regression-",today_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

save(nobs,M,N, delta, nnodes, 
     rmse.beta, rmse.fdaPDE.beta,
     norms.beta.2, norms.beta.inf,
     mise, mise.fdaPDE, err.L2,
     times,times.fdaPDE, file =save.file)




