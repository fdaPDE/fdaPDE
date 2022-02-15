
library(spatstat)
library(ggplot2)
source("~/Scrivania/fdaPDE/tests/Auxiliary/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")

vertices = cbind(simplenet$vertices$x,simplenet$vertices$y)
edges = cbind(simplenet$from, simplenet$to)
M = simplenet$m
boundary=find.boundary(M)
L = as.linnet(simplenet) 

###
delta=c(0.03, 0.015, 0.0075, 0.00375)

niter = 500
nfolds=10
epsilon = 0.0001 

nobs = c(200, 300, 400, 500) 
N = length(delta)

eta=0.0001
lambda=c(1e-6,5e-5,1e-5,5e-4,2.5e-4,1e-4,7.5e-3,5e-3,1e-3,5e-2)
nfolds=10
M=20
K=length(nobs)

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
tot.time = list()
err.L2 = list()
err.Inf = list()
lambda.opt = list() 

coef.exs = list( matrix(0,nrow=97, ncol=1),
                 matrix(0,nrow=194,ncol=1),
                 matrix(0,nrow=388,ncol=1),
                 matrix(0,nrow=775,ncol=1))
####
for( k in 1:K){
  print(paste("######## k = ", k ,"/",K," ########"))
mise[[k]] =  matrix(0,nrow=M,ncol=N)
mise.fdaPDE[[k]] = matrix(0,nrow=M,ncol=N)
err.L2[[k]] = matrix(0,nrow=M,ncol=N)

times[[k]] = matrix(0,nrow=M, ncol=N)
times.fdaPDE[[k]] = matrix(0,nrow=M,ncol=N)

tot.time[[k]] = 0

lambda.opt [[k]]= matrix(0,nrow=M,ncol=N)
mise.fdaPDE.l2[[k]] = matrix(0,nrow=M,ncol=N)
mise.l2[[k]] = matrix(0, nrow=M, ncol=N)
err.L2[[k]] = matrix(0, nrow=M, ncol=N)
err.Inf[[k]] = matrix(0, nrow=M, ncol=N)
# sols.fdaPDE = list( matrix(0,nrow=97,ncol=M),
#                     matrix(0,nrow=194,ncol=M),
#                     matrix(0,nrow=388,ncol=M),
#                     matrix(0,nrow=775,ncol=M))
# coef.exs = list( matrix(0,nrow=97, ncol=1),
#                  matrix(0,nrow=194,ncol=M),
#                  matrix(0,nrow=388,ncol=M),
#                  matrix(0,nrow=775,ncol=M))

sols.fdaPDE[[k]] = list(list(), list(), list(), list())
sols[[k]] = list(list(), list(), list(), list())

tot.start = Sys.time()
for(i in 1:N){
  print(paste("######## i = ", i ,"/",N," ########"))
    #if(i==1){
    mesh = create.mesh.1D.vertices(vertices, edges, delta[i])
    FEMbasis = create.FEM.basis.1D(mesh) 
  
    ### fdaPDE ###
    mesh.fdaPDE = create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
  #}else{
  #  mesh.fdaPDE = refine.by.splitting.mesh.1.5D(mesh = mesh.fdaPDE)
  #}
  
  FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)
  ### density ### 
  dijkstra.8 = Dijkstra(mesh.fdaPDE,8)
  dijkstra.4 = Dijkstra(mesh.fdaPDE,4)
  aux = function(x, y, seg, tp) { 
    
    sigma.8 = 0.1
    sigma.4 = 0.05
    
    res.8 = equal_split_discontinous(mesh.fdaPDE, sigma.8, dijkstra.8, x, y)
    res.4 = equal_split_discontinous(mesh.fdaPDE, sigma.4, dijkstra.4, x, y) 
    
    res = 0.5 * res.8$coef + 0.5 * res.4$coef
    return(res)
  }
  
  density <- linfun(aux, L)
  coef.ex = aux(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2])
  
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
    PPP = rlpp(nobs[k], density)
    points_x = PPP$data$x
    points_y = PPP$data$y
    points = cbind(points_x, points_y)
    #if(i==4){
    #  eta=1e-4
    #}
    #nb. nel tempo di fdaPDE è anche compresa la crossvalidazione.
    start.fdaPDE <- Sys.time()
    sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=lambda, 
                        nsimulations=niter, 
                        stepProposals =eta,
                        #  fvec = exp(cvec),
                        tol1 = 1e-6,
                        nfolds=nfolds, 
                        step_method="Fixed_Step", 
                        direction_method="Gradient",
                        preprocess_method = "SimplifiedCV")
    end.fdaPDE <- Sys.time()
    
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
    
    sols[[k]][[i]][,j] = coef
    sols.fdaPDE[[k]][[i]][,j] = coef.fdaPDE
    
    mise[[k]][j,i] = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE) )
    mise.fdaPDE[[k]][j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
    
    mise.l2[[k]][j,i] = norm(coef-coef.ex, type="2")
    mise.fdaPDE.l2[[k]][j,i] = norm(coef.fdaPDE-coef.ex, type="2")
    
    err.L2[[k]][j,i] = norm(coef-coef.fdaPDE, type="2")
    err.Inf[[k]][j,i] = max( abs(coef-coef.fdaPDE))
    
    times.fdaPDE[[k]][j,i] = difftime(end.fdaPDE, start.fdaPDE, units="mins")
    times[[k]][j,i] = difftime(end, start, units="mins")
    
    lambda.opt[[k]][j,i] = sol.fdaPDE$lambda
  }
  
  save(nobs, N, M,niter, nnodes, coef.exs, FEMbasis.fdaPDEs,
       mise, mise.fdaPDE,err.L2,err.Inf,
       mise.l2, mise.fdaPDE.l2,
       sols, sols.fdaPDE,
       times, times.fdaPDE,
       lambda.opt, tot.time,
       file = save.file)
}

tot.end = Sys.time()
tot.time[[k]][i] = difftime(tot.end, tot.start, units="mins")

save(nobs, N, M,niter, nnodes, coef.exs, FEMbasis.fdaPDEs,
     mise, mise.fdaPDE,err.L2,
     times, times.fdaPDE,
     lambda.opt, tot.time,
     file = save.file)

}
save(nobs, N, M,niter, nnodes, coef.exs, FEMbasis.fdaPDEs,
     mise, mise.fdaPDE,err.L2,
     times, times.fdaPDE,
     lambda.opt, tot.time,
     file = save.file)


today_ = Sys.Date()
init_ = "" # "-fdaPDE""-null_vector"
ntest_ = "-test-1-COMPLETE"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

save(nobs,N,M,K, niter, nnodes, coef.exs, FEMbasis.fdaPDEs,
     mise, mise.fdaPDE,err.L2,
     times, times.fdaPDE,
     lambda.opt, tot.time,
     file = save.file)

#############################################################################
mise.fine.mesh = matrix(0,nrow=M,ncol=N)
mise.fine.mesh[,N] = mise.fdaPDE[,N]

FEMbasis.fdaPDEs = list( list(), list(), list(), list())
for(i in 1:N){
  if(i==1){
  mesh = create.mesh.1D.vertices(vertices, edges, delta[i])
  nodes = mesh$nodes
  nnodes[i] = mesh$nnodes
  FEMbasis = create.FEM.basis.1D(mesh) 
  
  ### fdaPDE ###
  mesh.fdaPDE = create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
  }
  else{
    mesh.fdaPDE= refine.by.splitting.mesh.1.5D(mesh = mesh.fdaPDE)
  }
  FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)
  ###
  FEMbasis.fdaPDEs[[i]] =  FEMbasis.fdaPDE
  ### density ### 
  dijkstra.8 = Dijkstra(mesh.fdaPDE,8)
  dijkstra.4 = Dijkstra(mesh.fdaPDE,4)
  aux = function(x, y, seg, tp) { 
    
    sigma.8 = 0.1
    sigma.4 = 0.05
    
    res.8 = equal_split_discontinous(mesh.fdaPDE, sigma.8, dijkstra.8, x, y)
    res.4 = equal_split_discontinous(mesh.fdaPDE, sigma.4, dijkstra.4, x, y) 
    
    res = 0.5 * res.8$coef + 0.5 * res.4$coef
    return(res)
  }
  
  density <- linfun(aux, L)
  coef.ex = aux(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2])
  
  coef.exs[[i]] = coef.ex
  }

for(i in 1: (N-1) ){
  for(j in 1:M){
    coef = eval.FEM(FEM(sols.fdaPDE[[i]][,j], FEMbasis.fdaPDEs[[i]]), mesh.fdaPDE$nodes)
    coef.ex = eval.FEM(FEM(coef.exs[[i]], FEMbasis.fdaPDEs[[i]]), mesh.fdaPDE$nodes)
    mise.fine.mesh[j,i] = integrate_f(FEM((coef-coef.ex)^2 , FEMbasis.fdaPDE))
  }
}

boxplot(mise.fine.mesh)
mise.fine.mesh
mise.fdaPDE
