### DE test-2  Brickwall
library(spatstat)
library(ggplot2)
source("~/Scrivania/fdaPDE/tests/Auxiliary/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")

data("spiders")
vertices = cbind(spiders$domain$vertices$x, spiders$domain$vertices$y) 
edges = cbind(spiders$domain$from, spiders$domain$to)
L = as.linnet(spiders$domain)
M = spiders$domain$m

mesh = create.mesh.1D.vertices(vertices, edges, delta=40)
nnodes = mesh$nnodes
FEMbasis = create.FEM.basis.1D(mesh) 

mesh.fdaPDE = create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)

plot(mesh.fdaPDE)

#a sym
points(mesh$nodes[66,1],mesh$nodes[66,2],pch=16,col="red")
# sym
points(mesh$nodes[55,1], mesh$nodes[55,2],pch=16,col="blue")
points(mesh$nodes[98,1], mesh$nodes[98,2],pch=16,col="green")

dijkstra.66 = Dijkstra(mesh.fdaPDE,66)
dijkstra.55 = Dijkstra(mesh.fdaPDE,55)
dijkstra.98 = Dijkstra(mesh.fdaPDE,98)

# x = mesh.fdaPDE$nodes[,1]
# y = mesh.fdaPDE$nodes[,2]
# 
# idx.leaf = which( mesh.fdaPDE$nodesmarkers == TRUE)
# res.55 = equal_split_discontinous(mesh.fdaPDE, 90, dijkstra.55, x, y)
# res.55$bandwidth[idx.leaf]
# integrate_f(FEM(res.55$coef, FEMbasis.fdaPDE))
# idx.no.zero.55 = which( res.55$coef != 0.0)
# 
# res.55.cycle = equal_split_discontinous.cycle(mesh.fdaPDE, 90, dijkstra.55, x,y)
# 
# R_plot_graph.ggplot2(FEM(res.55$coef, FEMbasis.fdaPDE))
# res.98 = equal_split_discontinous(mesh.fdaPDE, 50, dijkstra.98, x,y)
# idx.no.zero.98 = which( res.98$coef != 0.0)
# R_plot_graph.ggplot2(FEM(res.98$coef, FEMbasis.fdaPDE))
# plot(mesh.fdaPDE)
# points(mesh.fdaPDE$nodes[idx.no.zero.55,], pch=16,col="green")
# points(mesh.fdaPDE$nodes[idx.no.zero.98,], pch=16,col="green")
# points(mesh.fdaPDE$nodes[idx.no.zero.66,], pch=16,col="green")

################### Density #################
aux.55 = function(x, y, seg, tp) { 
  
  sigma.55 = 90
  
  res.55= equal_split_discontinous(mesh.fdaPDE, sigma.55, dijkstra.55, x, y)
  
  res =  res.55$coef 
  return(res)
}
aux.98 = function(x, y, seg, tp) { 
  
  sigma.98 = 90
  
  res.98= equal_split_discontinous(mesh.fdaPDE, sigma.98, dijkstra.98, x, y)
  res =  res.98$coef 
  
  return(res)
}
aux.66 = function(x,y,tp,seg){
  sigma = 100
  h = 5*sigma
  Graph = mesh.fdaPDE
  dijkstra = dijkstra.66
  source = dijkstra$source
  points_ = cbind(x,y)
  
  coef = vector(mode="numeric", length=nrow(points_))
  is_vertex = is.vertex(Graph, points_)
  idx.vertex = which(is_vertex!=0)
  idx.same.y = which(points_[idx.vertex,2] == Graph$nodes[source,2])
  
  if(!is.empty(idx.same.y)){
    for( i in idx.same.y ){
      if( dijkstra$distance[idx.vertex[i]] < h ){
        coef[i] =  1./((2*pi)^0.5*sigma) * exp(-dijkstra$distance[idx.vertex[i]]^2/(2*sigma^2))
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
AUX = function(x,y,seg,tp){
  sym = aux.55(x,y,seg,tp)
  sym.2 = aux.98(x,y,seg,tp)
  sym.3 = aux.66(x,y,seg,tp)
  res = 1./3. * sym + 1./3. * sym.2 + 1./3. * sym.3
  return(res)                    
}
density <- linfun(AUX, L)
coef.ex = AUX(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2])
R_plot_graph.a.sym.ggplot2(FEM(coef.ex,FEMbasis.fdaPDE))
#############################################

R0 = R_mass_1D(FEMbasis)
R1 = R_stiff_1D(FEMbasis, M)
epsilon = 0.0001 

###
#nobs = c(1000, 1500, 2000, 2500)
nobs = c(300,500,700,1000)
N = length(nobs)

#lambda.350 = c(1e-4, 2.5e-4, 5e-4, 1e-3, 2.5e-3, 5e-3, 1e-2, 5e-2, 1e-1, 1)
#lambda.500= c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 5, 7.5, 10)
lambda = c(1e-4, 2.5e-4, 5e-4, 
           1e-3, 2.5e-3, 5e-3,
           1e-2,5e-2,
           1e-1, 1)
niter = 5000
nfolds=10
eta=0.001

M=20

mise = matrix(0,nrow=M,ncol=N)
mise.fdaPDE = matrix(0,nrow=M,ncol=N)
err.L2 = matrix(0,nrow=M,ncol=N)

times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)

integrals = matrix(0,nrow=M,ncol=N)
integrals.fdaPDE = matrix(0,nrow=M,ncol=N)

sols = array(0,dim=c(M,N,mesh$nnodes) ) #contiene vettore dei coeff
sols.fdaPDE = array(0,dim=c(M,N,mesh$nnodes) )

lambda.opt = matrix(0,nrow=M,ncol=N)

today_ = Sys.Date()
init_ ="-heat" #"-random"#"-null_vector" 
ntest_ = "-test-5-SimplifiedCV"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")


tot.start = Sys.time()
for(i in 1:N){
  for(j in 1:M){
    
    PP = rlpp(nobs[i], density)
    points = cbind(PP$data$x, PP$data$y)
    
    start.fdaPDE <- Sys.time()
    sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=lambda, 
                        nsimulations=niter, 
                        stepProposals = eta,
                        # fvec = exp(cvec),
                        tol1 = 1e-6,
                        nfolds=nfolds, 
                        step_method="Fixed_Step", 
                        direction_method="Gradient",
                        preprocess_method = "SimplifiedCV")
    end.fdaPDE <- Sys.time()
    
    lambda.opt[j,i] = sol.fdaPDE$lambda
    
    start <- Sys.time()
    result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    # NB. R need log(density) init.
    cvec = log( sol.fdaPDE$f_init[, which( abs(lambda - sol.fdaPDE$lambda) <= .Machine$double.eps )] )
    PHI = result$PHI
    sol = my_density.aldo(points, mesh, PHI, FEMbasis, R0, R1, sol.fdaPDE$lambda, eta,
                          cvec, niter, FALSE, tolerance=1e-6 )
    end <- Sys.time()
    
    start.fdaPDE <- Sys.time()
    sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=sol.fdaPDE$lambda, 
                        nsimulations=niter, 
                        stepProposals = eta,
                        fvec = sol.fdaPDE$f_init[, which( abs(lambda - sol.fdaPDE$lambda) <= .Machine$double.eps )],
                        tol1 = 1e-6,
                        nfolds=nfolds, 
                        step_method="Fixed_Step", 
                        direction_method="Gradient")
    end.fdaPDE <- Sys.time()
    tot.end = Sys.time()
    
    coef = exp(sol$density[, sol$iters ])
    coef.fdaPDE     = exp(sol.fdaPDE$g)
    
    sols[j,i,] = coef 
    sols.fdaPDE[j,i,] = coef.fdaPDE
    
    mise[j,i]   = integrate_f( FEM( (coef - coef.ex)^2, FEMbasis.fdaPDE)) 
    mise.fdaPDE[j,i]   = integrate_f( FEM( (coef.fdaPDE - coef.ex)^2, FEMbasis.fdaPDE) )
    err.L2[j,i] = integrate_f(FEM ( (coef.fdaPDE - coef)^2, FEMbasis.fdaPDE) )
    
    times[j,i] = difftime(end, start, units = "mins")
    times.fdaPDE[j,i] = difftime(end.fdaPDE, start.fdaPDE, units = "mins")
    
    integrals[j,i] = integrate_f( FEM(coef,FEMbasis.fdaPDE))
    integrals.fdaPDE[j,i] =   integrate_f(FEM(coef.fdaPDE,FEMbasis.fdaPDE))
    
    lambda.opt[j,i] = sol.fdaPDE$lambda
    
    save(nobs, N, M, coef.ex,niter,
         mise, mise.fdaPDE,err.L2,
         times, times.fdaPDE,
         integrals, integrals.fdaPDE,
         sols, sols.fdaPDE,
         lambda.opt, FEMbasis.fdaPDE,
         file = save.file)  
    }
}
tot.end = Sys.time()
tot.time = difftime(tot.end, tot.start, units="mins")

tot.time
mise
mise.fdaPDE

today_ = Sys.Date()
init_ ="-heat" #"-random"#"-null_vector" 
ntest_ = "-test-5-SimplifiedCV"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")


save(nobs, N, M, coef.ex,niter,
     mise, mise.fdaPDE,err.L2,
     times, times.fdaPDE,
     integrals, integrals.fdaPDE,
     sols, sols.fdaPDE,
     lambda.opt, tot.time, FEMbasis.fdaPDE,
     file = save.file)


mesh.ref = create.mesh.1D.vertices(vertices, edges, delta=15)
mesh.fdaPDE.ref = create.mesh.1.5D(mesh.ref$nodes, mesh.ref$segments)
FEMbasis.fdaPDE.ref = create.FEM.basis(mesh.fdaPDE.ref)

coef.98.ref = eval.FEM( FEM(res.98$coef, FEMbasis.fdaPDE), locations= mesh.fdaPDE.ref$nodes  )

R_plot_graph.ggplot2(FEM(coef.98.ref,FEMbasis.fdaPDE.ref))

######################################################################
################### increasing nodes number ##########################
######################################################################

library(spatstat)
library(ggplot2)
source("~/Scrivania/fdaPDE/tests/Auxiliary/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")

data("spiders")
vertices = cbind(spiders$domain$vertices$x, spiders$domain$vertices$y) 
edges = cbind(spiders$domain$from, spiders$domain$to)
L = as.linnet(spiders$domain)
M = spiders$domain$m

#delta=c(40, 30, 25, 17.5)
#delta=c(75, 65, 45, 25)
delta = c(70, 50, 30, 25)

N = length(delta)
# mesh=create.mesh.1.5D(vertices, edges)
# nnodes = vector(mode="integer", length=N+1) 
# nnodes[1] = nrow(mesh$nodes)
# for(i in 1:N){
#   mesh.ref = refine.mesh.1.5D(mesh, delta[i])
#   nnodes[i+1] = nrow(mesh.ref$nodes)
# }
# nnodes

niter = 5000
nfolds=10
epsilon = 0.0001 

nobs = 500 # 1000
N = length(delta)

#eta= 1e-4 # 1e-3
eta = 5e-4
#eta = 1e-3
nfolds=10
M=20

lambda = c(1e-4, 2.5e-4, 5e-4, 
           1e-3, 2.5e-3, 5e-3,
           1e-2,5e-2,
           1e-1, 1)

####
mise = matrix(0,nrow=M,ncol=N)
mise.fdaPDE = matrix(0,nrow=M,ncol=N)
err.L2 = matrix(0,nrow=M,ncol=N)

norm.l2 = matrix(0, nrow=M, ncol=N)
norm.l2.fdaPDE = matrix(0, nrow=M, ncol=N)

times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)

integrals = matrix(0,nrow=M,ncol=N)
integrals.fdaPDE = matrix(0,nrow=M,ncol=N)

nnodes = vector(mode="integer", length=N)
lambda.opt = matrix(0,nrow=M,ncol=N)

sols = list()
sols.fdaPDE = list()
coefs.ex = list()

tot.start = Sys.time()
for(i in 1:N){
  
  mesh = create.mesh.1D.vertices(vertices, edges, delta[i])
  nodes = mesh$nodes
  nnodes[i] = mesh$nnodes
  FEMbasis = create.FEM.basis.1D(mesh) 
  
  ### fdaPDE ###
  mesh.fdaPDE = create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
  FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)
  
  ################### Density #################
  dijkstra.66 = Dijkstra(mesh.fdaPDE,66)
  dijkstra.55 = Dijkstra(mesh.fdaPDE,55)
  dijkstra.98 = Dijkstra(mesh.fdaPDE,98)
  
  aux.55 = function(x, y, seg, tp) { 
    
    sigma.55 = 90
    
    res.55= equal_split_discontinous(mesh.fdaPDE, sigma.55, dijkstra.55, x, y)
    
    res =  res.55$coef 
    return(res)
  }
  aux.98 = function(x, y, seg, tp) { 
    
    sigma.98 = 90
    
    res.98= equal_split_discontinous(mesh.fdaPDE, sigma.98, dijkstra.98, x, y)
    res =  res.98$coef 
    
    return(res)
  }
  aux.66 = function(x,y,tp,seg){
    sigma = 100
    h = 5*sigma
    Graph = mesh.fdaPDE
    dijkstra = dijkstra.66
    source = dijkstra$source
    points_ = cbind(x,y)
    
    coef = vector(mode="numeric", length=nrow(points_))
    is_vertex = is.vertex(Graph, points_)
    idx.vertex = which(is_vertex!=0)
    idx.same.y = which(points_[idx.vertex,2] == Graph$nodes[source,2])
    
    if(!is.empty(idx.same.y)){
      for( i in idx.same.y ){
        if( dijkstra$distance[idx.vertex[i]] < h ){
          coef[i] =  1./((2*pi)^0.5*sigma) * exp(-dijkstra$distance[idx.vertex[i]]^2/(2*sigma^2))
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
  AUX = function(x,y,seg,tp){
    sym = aux.55(x,y,seg,tp)
    sym.2 = aux.98(x,y,seg,tp)
    sym.3 = aux.66(x,y,seg,tp)
    res = 1./3. * sym + 1./3. * sym.2 + 1./3. * sym.3
    return(res)                    
  }
  density <- linfun(AUX, L)
  coef.ex = AUX(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2])
  #############################################
  
  R0 = R_mass_1D(FEMbasis)
  R1 = R_stiff_1D(FEMbasis, simplenet$m)
  
  sols[[i]] = matrix(0, nrow=nnodes[i], ncol=M)
  sols.fdaPDE[[i]] = matrix(0, nrow=nnodes[i], ncol=M)
  coefs.ex[[i]] = coef.ex
  cvec = matrix(0, nrow=nnodes[i], ncol=1)
  
  for(j in 1:M){
    PPP = rlpp(nobs, density)
    points_x = PPP$data$x
    points_y = PPP$data$y
    points = cbind(points_x, points_y)
    
    #nb. nel tempo di fdaPDE è anche compresa la crossvalidazione.
    start.fdaPDE <- Sys.time()
    sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=lambda, 
                        nsimulations=niter, 
                        stepProposals =eta,
                        fvec = exp(cvec),
                        tol1 = 1e-6,
                        nfolds=nfolds, 
                        step_method="Fixed_Step", 
                        direction_method="Gradient",
                        preprocess_method = "SimplifiedCV")
    end.fdaPDE <- Sys.time()
    difftime(end.fdaPDE, start.fdaPDE, units="mins")
    
    start <- Sys.time()
    result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    # NB. R need log(density) init.
    cvec = log( sol.fdaPDE$f_init[, which( abs(lambda-sol.fdaPDE$lambda)<= .Machine$double.eps ) ])
    PHI = result$PHI
    sol = my_density.aldo(points, mesh, PHI, FEMbasis, R0, R1, 
                          sol.fdaPDE$lambda, eta, cvec, niter, 
                          tolerance=1e-6,
                          Nesterov= FALSE)
    end <- Sys.time()
    
    start.fdaPDE <- Sys.time()
    sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=sol.fdaPDE$lambda, 
                        nsimulations=niter, 
                        stepProposals =eta,
                        fvec = sol.fdaPDE$f_init[, which( abs(lambda-sol.fdaPDE$lambda)<= .Machine$double.eps ) ],
                        tol1 = 1e-6,
                        nfolds=nfolds, 
                        step_method="Fixed_Step", 
                        direction_method="Gradient")
    end.fdaPDE <- Sys.time()
    
    coef = exp(sol$density[,sol$iter])
    coef.fdaPDE     = exp(sol.fdaPDE$g)
    
    sols[[i]][,j] = coef
    sols.fdaPDE[[i]][,j] = coef.fdaPDE
        
    mise[j,i]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
    mise.fdaPDE[j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
    err.L2[j,i] = integrate_f(FEM ( (coef.fdaPDE - coef)^2, FEMbasis.fdaPDE) )
    
    norm.l2[j,i] = norm(coef - coef.ex, type="2")
    norm.l2.fdaPDE[j,i] = norm(coef.fdaPDE - coef.ex, type="2")
    
    times[j,i] = difftime(end, start, units="mins")
    times.fdaPDE[j,i] = difftime(end.fdaPDE, start.fdaPDE, units="mins")
    
    integrals[j,i] = integrate_f( FEM(coef,FEMbasis.fdaPDE))
    integrals.fdaPDE[j,i] =   integrate_f(FEM(coef,FEMbasis.fdaPDE))
    
    lambda.opt[j,i] = sol.fdaPDE$lambda
    
    save(nobs, N, M,niter, nnodes,
         mise, mise.fdaPDE,err.L2, 
         norm.l2, norm.l2.fdaPDE,
         sols, sols.fdaPDE,
         times, times.fdaPDE,
         integrals, integrals.fdaPDE,
         lambda.opt, FEMbasis.fdaPDE,
         file = save.file)
  }
  
}

tot.end = Sys.time()
tot.time = difftime(tot.end, tot.start, units="mins")

today_ = Sys.Date()
init_ = "" # "-fdaPDE""-null_vector"
ntest_ = "test-2-delta-var-500-obs"#"-test-2-delta-var-500-obs"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

# save(nobs, N, M,niter, nnodes,
#      mise, mise.fdaPDE,err.L2,
#      times, times.fdaPDE,
#      integrals, integrals.fdaPDE,
#      lambda.opt, tot.time, FEMbasis.fdaPDE,
#      file = save.file)
# 

######################################################################
##################        fdaPDE ONLY           ######################
################### increasing nodes number ##########################
######################################################################

library(spatstat)
library(ggplot2)
source("~/Scrivania/fdaPDE/tests/Auxiliary/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")

data("spiders")
vertices = cbind(spiders$domain$vertices$x, spiders$domain$vertices$y) 
edges = cbind(spiders$domain$from, spiders$domain$to)
L = as.linnet(spiders$domain)
M = spiders$domain$m

delta=c(75, 65, 45, 25)
#delta=c(40, 30, 25, 17.5)
mesh=create.mesh.1.5D(vertices, edges)
nnodes = c(156, 337, 403, 520, 832)

niter = 500
nfolds=10
epsilon = 0.0001 

nobs = 1000
N = length(delta) + 1 

eta=0.001
nfolds=10
M=10

lambda = c(1e-4, 2.5e-4, 5e-4, 
           1e-3, 2.5e-3, 5e-3,
           1e-2,5e-2,
           1e-1, 1)

lambda = 0.05
####
mise = matrix(0,nrow=M,ncol=N)
mise.fdaPDE = matrix(0,nrow=M,ncol=N)
err.L2 = matrix(0,nrow=M,ncol=N)

times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)

integrals = matrix(0,nrow=M,ncol=N)
integrals.fdaPDE = matrix(0,nrow=M,ncol=N)

lambda.opt = matrix(0,nrow=M,ncol=N)

tot.start = Sys.time()
for(i in 1:N){
  
  if(i == 1){
    mesh.fdaPDE = mesh
  }else{
    mesh.fdaPDE = refine.mesh.1.5D(mesh, delta[1])
  }
    FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)
  
  ################### Density #################
  dijkstra.66 = Dijkstra(mesh.fdaPDE,66)
  dijkstra.55 = Dijkstra(mesh.fdaPDE,55)
  dijkstra.98 = Dijkstra(mesh.fdaPDE,98)
  
  aux.55 = function(x, y, seg, tp) { 
    
    sigma.55 = 90
    
    res.55= equal_split_discontinous(mesh.fdaPDE, sigma.55, dijkstra.55, x, y)
    
    res =  res.55$coef 
    return(res)
  }
  aux.98 = function(x, y, seg, tp) { 
    
    sigma.98 = 90
    
    res.98= equal_split_discontinous(mesh.fdaPDE, sigma.98, dijkstra.98, x, y)
    res =  res.98$coef 
    
    return(res)
  }
  aux.66 = function(x,y,tp,seg){
    sigma = 100
    h = 5*sigma
    Graph = mesh.fdaPDE
    dijkstra = dijkstra.66
    source = dijkstra$source
    points_ = cbind(x,y)
    
    coef = vector(mode="numeric", length=nrow(points_))
    is_vertex = is.vertex(Graph, points_)
    idx.vertex = which(is_vertex!=0)
    idx.same.y = which(points_[idx.vertex,2] == Graph$nodes[source,2])
    
    if(!is.empty(idx.same.y)){
      for( i in idx.same.y ){
        if( dijkstra$distance[idx.vertex[i]] < h ){
          coef[i] =  1./((2*pi)^0.5*sigma) * exp(-dijkstra$distance[idx.vertex[i]]^2/(2*sigma^2))
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
  AUX = function(x,y,seg,tp){
    sym = aux.55(x,y,seg,tp)
    sym.2 = aux.98(x,y,seg,tp)
    sym.3 = aux.66(x,y,seg,tp)
    res = 1./3. * sym + 1./3. * sym.2 + 1./3. * sym.3
    return(res)                    
  }
  density <- linfun(AUX, L)
  coef.ex = AUX(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2])
  #############################################
  
  #R0 = R_mass_1D(FEMbasis)
  #R1 = R_stiff_1D(FEMbasis, simplenet$m)
  
  for(j in 1:M){
    PPP = rlpp(nobs, density)
    points_x = PPP$data$x
    points_y = PPP$data$y
    points = cbind(points_x, points_y)
    
    #nb. nel tempo di fdaPDE è anche compresa la crossvalidazione.
    # start.fdaPDE <- Sys.time()
    # sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=lambda, 
    #                     nsimulations=niter, 
    #                     stepProposals =eta,
    #                     #  fvec = exp(cvec),
    #                     tol1 = 1e-6,
    #                     nfolds=nfolds, 
    #                     step_method="Fixed_Step", 
    #                     direction_method="Gradient",
    #                     preprocess_method = "SimplifiedCV")
    # end.fdaPDE <- Sys.time()
    
    # start <- Sys.time()
    # result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    # # NB. R need log(density) init.
    # cvec = log( sol.fdaPDE$f_init[, which( abs(lambda-sol.fdaPDE$lambda)<= .Machine$double.eps ) ])
    # PHI = result$PHI
    # sol = my_density.aldo(points, mesh, PHI, FEMbasis, R0, R1, 
    #                       sol.fdaPDE$lambda, eta, cvec, niter, 
    #                       tolerance=1e-6,
    #                       Nesterov= FALSE)
    # end <- Sys.time()
    # 
    # start.fdaPDE <- Sys.time()
    # sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=sol.fdaPDE$lambda,
    #                     nsimulations=niter,
    #                     stepProposals =eta,
    #                     fvec = sol.fdaPDE$f_init[, which( abs(lambda-sol.fdaPDE$lambda)<= .Machine$double.eps ) ],
    #                     tol1 = 1e-6,
    #                     nfolds=nfolds,
    #                     step_method="Fixed_Step",
    #                     direction_method="Gradient")
    # end.fdaPDE <- Sys.time()
    # 
    #coef = exp(sol$density[,sol$iter])
    start.fdaPDE <- Sys.time()
    sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=lambda,
                        nsimulations=niter,
                        stepProposals =eta,
                        #fvec = sol.fdaPDE$f_init[, which( abs(lambda-sol.fdaPDE$lambda)<= .Machine$double.eps ) ],
                        tol1 = 1e-6,
                        nfolds=nfolds,
                        step_method="Fixed_Step",
                        direction_method="Gradient")
    end.fdaPDE <- Sys.time()
    
    coef.fdaPDE     = exp(sol.fdaPDE$g)
    
    #mise[j,i]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
    mise.fdaPDE[j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
    #err.L2[j,i] = integrate_f(FEM ( (coef.fdaPDE - coef)^2, FEMbasis.fdaPDE) )
    
    #times[j,i] = difftime(end, start, units="mins")
    times.fdaPDE[j,i] = difftime(end.fdaPDE, start.fdaPDE, units="mins")
  }
  
}

boxplot(mise.fdaPDE)

tot.end = Sys.time()
tot.time = difftime(tot.end, tot.start, units="mins")

today_ = Sys.Date()
init_ = "-heat" # "-fdaPDE""-null_vector"
ntest_ = "-test-2-fdaPDE-ONLY-156-800-nodes-1000obs-parte-2"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

save(nobs, N, M,niter, nnodes,L,density,
     mise.fdaPDE,
     times.fdaPDE,
     tot.time, FEMbasis.fdaPDE,
     file = save.file)
