library(spatstat)
library(ggplot2)
source("~/Scrivania/fdaPDE/tests/Auxiliary/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")

data("brickwall.1.5D")

vertices = nodes
delta=c(0.05, 0.03, 0.015, 0.01)

mesh = create.mesh.1D.vertices(nodes, edges, 0.01)
mesh$nnodes

mesh.fdaPDE = create.mesh.1.5D(mesh$nodes, mesh$segments)
plot(mesh.fdaPDE, show.nodes = TRUE, pch=16)

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

nobs = 300 # 1000
N = length(delta)

#eta= 1e-4 # 1e-3
eta = 5e-4
#eta = 1e-3
nfolds=10
M=20

lambda = seq(from=1e-6, to=1e-4, length.out=10)

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
  dijkstra.7 = Dijkstra(mesh.fdaPDE,7)
  dijkstra.11 = Dijkstra(mesh.fdaPDE,11)
  
  aux.7 = function(x,y,seg, tp){
    sigma = 0.09
    h = 5*sigma
    Graph = mesh.fdaPDE
    dijkstra = dijkstra.7
    source = dijkstra$source
    points_ = cbind(x,y)
    coef = vector(mode="numeric", length=length(x) ) 
    
    idx.ok = which( abs(points_[,2] - Graph$nodes[source,2]) < 10 * .Machine$double.eps )
    for( i in idx.ok){
      
      delta =  abs(points_[i,1]-Graph$nodes[source,1])
      
      if(delta < h ){
        coef[i] = 1/((2*pi)^0.5 * sigma) * exp(-delta^2/(2*sigma^2))
      }
      
    }
    
    return (coef)
  }
  
  aux.11 = function(x,y,seg, tp){
    sigma = 0.09
    h = 5*sigma
    Graph = mesh.fdaPDE
    dijkstra = dijkstra.11
    source = dijkstra$source
    points_ = cbind(x,y)
    coef = vector(mode="numeric", length=length(x) ) 
    
    idx.ok = which( abs(points_[,2] - Graph$nodes[source,2]) < 10 * .Machine$double.eps )
    for( i in idx.ok){
      
      delta =  abs(points_[i,1]-Graph$nodes[source,1])
      
      if(delta < h ){
        coef[i] = 1/((2*pi)^0.5 * sigma) * exp(-delta^2/(2*sigma^2))
      }
      
    }
    
    return (coef)
  }
  
  AUX = function(x,y,seg,tp){
    
    #res = 1/3*( aux.7(x,y,seg,tp) + aux.2(x,y,seg,tp) + aux.3(x,y,seg,tp))
    res = 0.5* ( aux.7(x,y,seg,tp) + aux.11(x,y,seg,tp))
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
  #cvec = matrix(0, nrow=nnodes[i], ncol=1)
  
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
                        #fvec = exp(cvec),
                        tol1 = 1e-5,
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
                          tolerance=1e-5,
                          Nesterov= FALSE)
    end <- Sys.time()
    
    start.fdaPDE <- Sys.time()
    sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=sol.fdaPDE$lambda, 
                        nsimulations=niter, 
                        stepProposals =eta,
                        fvec = sol.fdaPDE$f_init[, which( abs(lambda-sol.fdaPDE$lambda)<= .Machine$double.eps ) ],
                        tol1 = 1e-5,
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
ntest_ = "test-2-brick-aldo-delta-var-300-obs"#"-test-2-delta-var-700-obs"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

######### nobs varying ####################

library(spatstat)
library(ggplot2)
source("~/Scrivania/fdaPDE/tests/Auxiliary/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")

data("brickwall.1.5D")

vertices = nodes
delta=0.015

N = length(delta)

mesh = create.mesh.1D.vertices(vertices, edges, delta)
nodes = mesh$nodes
nnodes = mesh$nnodes
FEMbasis = create.FEM.basis.1D(mesh) 

### fdaPDE ###
mesh.fdaPDE = create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)
plot(mesh.fdaPDE, show.nodes = TRUE, pch=16, asp=1)
################### Density #################
{
dijkstra.7 = Dijkstra(mesh.fdaPDE,7)
dijkstra.11 = Dijkstra(mesh.fdaPDE,11)
#dijkstra.2 = Dijkstra(mesh.fdaPDE,2)
#dijkstra.3 = Dijkstra(mesh.fdaPDE,3)

aux.7 = function(x,y,seg, tp){
  sigma = 0.09
  h = 5*sigma
  Graph = mesh.fdaPDE
  dijkstra = dijkstra.7
  source = dijkstra$source
  points_ = cbind(x,y)
  coef = vector(mode="numeric", length=length(x) ) 
  
  idx.ok = which( abs(points_[,2] - Graph$nodes[source,2]) < 10 * .Machine$double.eps )
  for( i in idx.ok){
    
    delta =  abs(points_[i,1]-Graph$nodes[source,1])
    
    if(delta < h ){
      coef[i] = 1/((2*pi)^0.5 * sigma) * exp(-delta^2/(2*sigma^2))
    }
    
  }
  
  return (coef)
}

aux.2 = function(x,y,seg, tp){
  sigma = 0.09
  h = 5*sigma
  Graph = mesh.fdaPDE
  dijkstra = dijkstra.2
  source = dijkstra$source
  points_ = cbind(x,y)
  coef = vector(mode="numeric", length=length(x) ) 
  
  idx.ok = which( abs(points_[,2] - Graph$nodes[source,2]) < 10 * .Machine$double.eps )
  for( i in idx.ok){
    
    delta =  abs(points_[i,1]-Graph$nodes[source,1])
    
    if(delta < h ){
      coef[i] = 1/((2*pi)^0.5 * sigma) * exp(-delta^2/(2*sigma^2))
    }
    
  }
  
  return (coef)
}
aux.3 = function(x,y,seg, tp){
  sigma = 0.09
  h = 5*sigma
  Graph = mesh.fdaPDE
  dijkstra = dijkstra.3
  source = dijkstra$source
  points_ = cbind(x,y)
  coef = vector(mode="numeric", length=length(x) ) 
  
  idx.ok = which( abs(points_[,2] - Graph$nodes[source,2]) < 10 * .Machine$double.eps )
  for( i in idx.ok){
    
    delta =  abs(points_[i,1]-Graph$nodes[source,1])
    
    if(delta < h ){
      coef[i] = 1/((2*pi)^0.5 * sigma) * exp(-delta^2/(2*sigma^2))
    }
    
  }
  
  return (coef)
}

aux.11 = function(x,y,seg, tp){
  sigma = 0.09
  h = 5*sigma
  Graph = mesh.fdaPDE
  dijkstra = dijkstra.11
  source = dijkstra$source
  points_ = cbind(x,y)
  coef = vector(mode="numeric", length=length(x) ) 
  
  idx.ok = which( abs(points_[,2] - Graph$nodes[source,2]) < 10 * .Machine$double.eps )
  for( i in idx.ok){
    
    delta =  abs(points_[i,1]-Graph$nodes[source,1])
    
    if(delta < h ){
      coef[i] = 1/((2*pi)^0.5 * sigma) * exp(-delta^2/(2*sigma^2))
    }
    
  }
  
  return (coef)
}

AUX = function(x,y,seg,tp){
  
  #res = 1/3*( aux.7(x,y,seg,tp) + aux.2(x,y,seg,tp) + aux.3(x,y,seg,tp))
  res = 0.5* ( aux.7(x,y,seg,tp) + aux.11(x,y,seg,tp))
  return(res)                    
}
density <- linfun(AUX, L)
coef.ex = AUX(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2])
R_plot_graph.a.sym.ggplot2(FEM(coef.ex, FEMbasis.fdaPDE))
}
#############################################

R0 = R_mass_1D(FEMbasis)
R1 = R_stiff_1D(FEMbasis, simplenet$m)
###

niter = 5000
nfolds=10
epsilon = 0.0001 

nobs = c(300, 500, 700, 1000) # 1000
#nobs=c(250,350)
N = length(nobs)

#eta= 1e-4 # 1e-3
eta = 5e-4
#eta = 1e-3
nfolds=10
M=20

lambda = seq(from=1e-6, to=1e-4, length.out=10)

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

lambda.opt = matrix(0,nrow=M,ncol=N)

sols = array(0, dim=c(M,N, nnodes) )
sols.fdaPDE = array(0, dim=c(M,N, nnodes) )


tot.start = Sys.time()
for(i in 1:N){
  for(j in 1:M){
    PPP = rlpp(nobs[i], density)
    points_x = PPP$data$x
    points_y = PPP$data$y
    points = cbind(points_x, points_y)
    
    #nb. nel tempo di fdaPDE è anche compresa la crossvalidazione.
    start.fdaPDE <- Sys.time()
    sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=lambda, 
                        nsimulations=niter, 
                        stepProposals =eta,
                        #fvec = exp(cvec),
                        tol1 = 5e-6,
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
                          tolerance=5e-6,
                          Nesterov= FALSE)
    end <- Sys.time()
    
    start.fdaPDE <- Sys.time()
    sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=sol.fdaPDE$lambda, 
                        nsimulations=niter, 
                        stepProposals =eta,
                        fvec = sol.fdaPDE$f_init[, which( abs(lambda-sol.fdaPDE$lambda)<= .Machine$double.eps ) ],
                        tol1 = 5e-6,
                        nfolds=nfolds, 
                        step_method="Fixed_Step", 
                        direction_method="Gradient")
    end.fdaPDE <- Sys.time()
    
    coef = exp(sol$density[,sol$iter])
    coef.fdaPDE     = exp(sol.fdaPDE$g)
    
    sols[j,i,] = coef
    sols.fdaPDE[j,i,] = coef.fdaPDE
    
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
    
    save(nobs, N, M,niter, nnodes, nobs,
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
ntest_ = "test-2-brick-aldo-nobs-var"#"-test-2-delta-var-500-obs"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")
