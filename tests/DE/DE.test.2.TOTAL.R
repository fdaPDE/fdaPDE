#######################################################
#                                                     #
#                                                     #  
#######################################################

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

###
nobs = c(300,500,700,1000)
N = length(nobs)
delta=c(40, 30, 25, 17.5)
lambda = c(1e-4, 2.5e-4, 5e-4, 
           1e-3, 2.5e-3, 5e-3,
           1e-2,5e-2,
           1e-1, 1)
niter = 500
nfolds=10
eta=0.001
epsilon = 0.0001 

M=20
K=length(delta)

nnodes = vector(mode="numeric", length=N)
FEMbasis.fdaPDEs = list()

mise= list() 
mise.fdaPDE = list()
mise.l2 = list()
mise.fdaPDE.l2 = list() 

times = list()
times.fdaPDE = list()

sols.fdaPDE = list()
sols = list()
coef.exs = list()

err.L2 = list()
err.Inf = list()
err.2 = list()

coef.exs = list( matrix(0,nrow=412, ncol=1),
                 matrix(0,nrow=637,ncol=1),
                 matrix(0,nrow=831,ncol=1),
                 matrix(0,nrow=1143,ncol=1))

lambda = 0.05
####
tot.start = Sys.time()
for( k in 1:K){
  print(paste("######## k = ", k ,"/",K," ########"))
  mise[[k]] =  matrix(0,nrow=M,ncol=N)
  mise.fdaPDE[[k]] = matrix(0,nrow=M,ncol=N)
  err.L2[[k]] = matrix(0,nrow=M,ncol=N)
  
  times[[k]] = matrix(0,nrow=M, ncol=N)
  times.fdaPDE[[k]] = matrix(0,nrow=M,ncol=N)

  lambda.opt [[k]]= matrix(0,nrow=M,ncol=N)
  mise.fdaPDE.l2[[k]] = matrix(0,nrow=M,ncol=N)
  mise.l2[[k]] = matrix(0, nrow=M, ncol=N)
  err.L2[[k]] = matrix(0, nrow=M, ncol=N)
  err.Inf[[k]] = matrix(0, nrow=M, ncol=N)
  err.2[[k]] = matrix(0, nrow=M, ncol=N)
  sols.fdaPDE[[k]] = list(list(), list(), list(), list())
  sols[[k]] = list(list(), list(), list(), list())
  
  for(i in 1:N){
    print(paste("######## i = ", i ,"/",N," ########"))
    #if(i==1){
    mesh = create.mesh.1D.vertices(vertices, edges, delta[i])
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
    
    sols.fdaPDE[[k]][[i]] = matrix(0,nrow= nrow(mesh.fdaPDE$nodes), ncol=M)
    sols[[k]][[i]] = matrix(0,nrow= nrow(mesh.fdaPDE$nodes), ncol=M)
    
    if( k == 1){
      coef.exs[[i]] = coef.ex
      nnodes[i] = mesh$nnodes
      FEMbasis.fdaPDEs[[i]] = FEMbasis.fdaPDE
    }
    
    R0 = R_mass_1D(FEMbasis)
    R1 = R_stiff_1D(FEMbasis, simplenet$m)
    
    for(j in 1:M){
      print(paste("#  k = ", k ,"/",K," ### i = ",i," ### j = ",j," #",sep="") )
      PPP = rlpp(nobs[k], density)
      points_x = PPP$data$x
      points_y = PPP$data$y
      points = cbind(points_x, points_y)
      
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
      # 
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
      start.fdaPDE <- Sys.time()
      sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=lambda, 
                          nsimulations=niter, 
                          stepProposals =eta,
                          tol1 = 1e-6,
                          nfolds=nfolds, 
                          step_method="Fixed_Step", 
                          direction_method="Gradient")
      end.fdaPDE <- Sys.time()
      
      start <- Sys.time()
      result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
      # NB. R need log(density) init.
      PHI = result$PHI
      sol = my_density.aldo(points, mesh, PHI, FEMbasis, R0, R1, 
                            lambda, eta, sol.fdaPDE$f_init, niter, 
                            tolerance=1e-6,
                            Nesterov= FALSE)
      end <- Sys.time()
      
      
      coef = exp(sol$density[,sol$iter])
      coef.fdaPDE     = exp(sol.fdaPDE$g)
      
      sols[[k]][[i]][,j] = coef
      sols.fdaPDE[[k]][[i]][,j] = coef.fdaPDE
      
      mise[[k]][j,i] = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE) )
      mise.fdaPDE[[k]][j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
      
      err.L2[[k]][j,i] = integrate_f(FEM ( (coef.fdaPDE - coef)^2, FEMbasis.fdaPDE) )
      err.2[[k]][j,i] = norm(coef-coef.fdaPDE, type="2")
      err.Inf[[k]][j,i] = max( abs(coef-coef.fdaPDE))
      
      times.fdaPDE[[k]][j,i] = difftime(end.fdaPDE, start.fdaPDE, units="mins")
      times[[k]][j,i] = difftime(end, start, units="mins")
      
      save(nobs, N, M,niter, nnodes, coef.exs, FEMbasis.fdaPDEs,
           mise, mise.fdaPDE,
           err.L2,err.2,err.Inf,L,
           sols, sols.fdaPDE,
           times, times.fdaPDE,
           file = save.file)
      
      }
  }

  
}

tot.end = Sys.time()
tot.time = difftime(tot.end, tot.start, units="mins")

save(nobs, N, M,niter, nnodes, coef.exs, FEMbasis.fdaPDEs,
     mise, mise.fdaPDE,
     err.L2,err.2,err.Inf,L,
     sols, sols.fdaPDE,
     times, times.fdaPDE,
     file = save.file)


today_ = Sys.Date()
init_ = "" # "-fdaPDE""-null_vector"
ntest_ = "-test-2-COMPLETE"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

