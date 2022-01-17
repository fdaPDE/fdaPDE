### confronto Mattina - fdaPDE 2.0 ###
library(spatstat)
library(ggplot2)
source("~/Scrivania/fdaPDE/tests/Auxiliary/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")

vertices_x = as.matrix(simplenet$vertices$x)
vertices_y = as.matrix(simplenet$vertices$y)
vertices = as.matrix(cbind(vertices_x, vertices_y))
start_coord = cbind(simplenet$lines$ends$x0, simplenet$lines$ends$y0)
end_coord = cbind(simplenet$lines$ends$x1, simplenet$lines$ends$y1)
start = simplenet$from 
end = simplenet$to
edges = cbind(start, end)
M = simplenet$m
boundary=find.boundary(M)
L = as.linnet(simplenet) 

# delta = 0.03 -> nnodes=97
delta=0.015

mesh = create.mesh.1D.vertices(vertices, edges, delta)
nodes = mesh$nodes
nnodes = mesh$nnodes
FEMbasis = create.FEM.basis.1D(mesh) 

### fdaPDE ###
mesh.fdaPDE = create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)
###

### density ### 
dijkstra.8 = Dijkstra(mesh.fdaPDE,8)
dijkstra.37 = Dijkstra(mesh.fdaPDE,37)
aux = function(x, y, seg, tp) { 
  
  sigma.8 = 0.1
  sigma.37 = 0.085
  
  res.8 = equal_split_discontinous(mesh.fdaPDE, sigma.8, dijkstra.8, x, y)
  res.37 = equal_split_discontinous(mesh.fdaPDE, sigma.37, dijkstra.37, x, y) 
  
  res = 0.5 * res.8$coef + 0.5 * res.37$coef
  return(res)
}

density <- linfun(aux, L)
coef.ex = aux(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2])
###
niter = 20000
nfolds=10
R0 = R_mass_1D(FEMbasis)
R1 = R_stiff_1D(FEMbasis, simplenet$m)
epsilon = 0.0001 

###
nobs = c(200,300,400,500)
#nobs = c(50)
N = length(nobs)
#cvec=matrix(0,nrow=mesh$nnodes,ncol=1)

eta=0.001
lambda=c(1e-6,5e-5,1e-5,5e-4,2.5e-4,1e-4,7.5e-3,5e-3,1e-3,5e-2)
nfolds=10
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
                        stepProposals = eta,
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
    cvec = log( sol.fdaPDE$f_init[, which( abs(lambda - sol.fdaPDE$lambda) <= .Machine$double.eps )] )
    PHI = result$PHI
    sol = my_density.aldo(points, mesh, PHI, FEMbasis, R0, R1, 
                          sol.fdaPDE$lambda, eta, cvec, niter, 
                          tolerance=1e-6, Nesterov= FALSE)
    end <- Sys.time()
    
    coef = exp(sol$density[,sol$iter])
    coef.fdaPDE     = exp(sol.fdaPDE$g)
   
    sols[j,i,] = coef 
    sols.fdaPDE[j,i,] = coef.fdaPDE
    
    mise[j,i]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
    mise.fdaPDE[j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
    err.L2[j,i] = integrate_f(FEM ( (coef.fdaPDE - coef)^2, FEMbasis.fdaPDE) )
    
    times[j,i] = difftime(end, start, units="mins")
    times.fdaPDE[j,i] = difftime(end.fdaPDE, start.fdaPDE, units="mins")
    
    integrals[j,i] = integrate_f( FEM(coef,FEMbasis.fdaPDE))
    integrals.fdaPDE[j,i] =   integrate_f(FEM(coef,FEMbasis.fdaPDE))
    
    lambda.opt[j,i] = sol.fdaPDE$lambda
  }
  
  save(nobs, N, M, coef.ex,niter,
       mise, mise.fdaPDE,err.L2,
       times, times.fdaPDE,
       integrals, integrals.fdaPDE,
       sols, sols.fdaPDE,
       lambda.opt, FEMbasis.fdaPDE,
       file = save.file)
}

tot.end = Sys.time()
tot.time = difftime(tot.end, tot.start, units="mins")

today_ = Sys.Date()
init_ = "-heat" # "-fdaPDE""-null_vector"
ntest_ = "-test-1"

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

###############################################################################
####################### increasing nodes number ###############################
###############################################################################

library(spatstat)
library(ggplot2)
source("~/Scrivania/fdaPDE/tests/Auxiliary/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")

vertices_x = as.matrix(simplenet$vertices$x)
vertices_y = as.matrix(simplenet$vertices$y)
vertices = as.matrix(cbind(vertices_x, vertices_y))
start_coord = cbind(simplenet$lines$ends$x0, simplenet$lines$ends$y0)
end_coord = cbind(simplenet$lines$ends$x1, simplenet$lines$ends$y1)
start = simplenet$from 
end = simplenet$to
edges = cbind(start, end)
M = simplenet$m
boundary=find.boundary(M)
L = as.linnet(simplenet) 

###
mesh = create.mesh.1D.vertices(vertices, edges, delta[4])
nnodes[4] = mesh$nnodes
delta=c(0.03, 0.015, 0.0075, 0.00375)

niter = 500
nfolds=10
epsilon = 0.0001 

nobs = 500
N = length(delta)

eta=0.001
lambda=c(1e-6,5e-5,1e-5,5e-4,2.5e-4,1e-4,7.5e-3,5e-3,1e-3,5e-2)
nfolds=10
M=20

####
mise = matrix(0,nrow=M,ncol=N)
mise.fdaPDE = matrix(0,nrow=M,ncol=N)
err.L2 = matrix(0,nrow=M,ncol=N)

times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)

integrals = matrix(0,nrow=M,ncol=N)
integrals.fdaPDE = matrix(0,nrow=M,ncol=N)

nnodes = vector(mode="integer", length=N)
lambda.opt = matrix(0,nrow=M,ncol=N)

tot.start = Sys.time()
for(i in 1:N){
  
  mesh = create.mesh.1D.vertices(vertices, edges, delta[i])
  nodes = mesh$nodes
  nnodes[i] = mesh$nnodes
  FEMbasis = create.FEM.basis.1D(mesh) 
  
  ### fdaPDE ###
  mesh.fdaPDE = create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
  FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)
  ###
  
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
  
  R0 = R_mass_1D(FEMbasis)
  R1 = R_stiff_1D(FEMbasis, simplenet$m)
  
  for(j in 1:M){
    PPP = rlpp(nobs, density)
    points_x = PPP$data$x
    points_y = PPP$data$y
    points = cbind(points_x, points_y)
    if(i==4){
      eta=1e-4
    }
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
    
    mise[j,i]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
    mise.fdaPDE[j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
    err.L2[j,i] = integrate_f(FEM ( (coef.fdaPDE - coef)^2, FEMbasis.fdaPDE) )
    
    times[j,i] = difftime(end, start, units="mins")
    times.fdaPDE[j,i] = difftime(end.fdaPDE, start.fdaPDE, units="mins")
    
    integrals[j,i] = integrate_f( FEM(coef,FEMbasis.fdaPDE))
    integrals.fdaPDE[j,i] =   integrate_f(FEM(coef,FEMbasis.fdaPDE))
    
    lambda.opt[j,i] = sol.fdaPDE$lambda
  }

}

tot.end = Sys.time()
tot.time = difftime(tot.end, tot.start, units="mins")

today_ = Sys.Date()
init_ = "-heat" # "-fdaPDE""-null_vector"
ntest_ = "-test-1-delta-var"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

save(nobs, N, M,niter, nnodes,
     mise, mise.fdaPDE,err.L2,
     times, times.fdaPDE,
     integrals, integrals.fdaPDE,
     lambda.opt, tot.time, FEMbasis.fdaPDE,
     file = save.file)
