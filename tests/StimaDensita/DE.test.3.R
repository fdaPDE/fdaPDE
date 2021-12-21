## DE test 4.1 ## 

source("~/Scrivania/fdaPDE/tests/StimaDensita/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Dijkstra.R")
source("~/Scrivania/fdaPDE/tests/Mesh_Evaluator/isInside.R")
source("tests/R_plot_graph.ggplot2.R")
library(ggplot2)
data("dendrite")

vertices = cbind(dendrite$domain$vertices$x, dendrite$domain$vertices$y) 
edges = cbind(dendrite$domain$from, dendrite$domain$to)
L = as.linnet(dendrite$domain)
M = dendrite$domain$m

mesh =create.mesh.1D.vertices(vertices, edges, 10)
mesh$nnodes

nnodes = mesh$nnodes
FEMbasis = create.FEM.basis.1D(mesh) 

### fdaPDE ###
mesh.fdaPDE = create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)
dim(mesh.fdaPDE$nodes)
###
plot(mesh.fdaPDE)

res.600 = Dijkstra(mesh.fdaPDE, 600)
sigma = 40
AUX = function(x){  1/(sigma*(2*pi)^1.5) * exp(-x^2/(2*sigma^2)) }
coef.600 = AUX(res.600$distance)

#per rlpp
aux.600 = function(x,y,seg,tp){ 
  sigma=40 
  DijkstraDec=res.600
  Graph=mesh.fdaPDE
  
  points_ = cbind(x,y)
  # nb il cbind per renderla ANCHE vettorizzabile
  idx = isInside(Graph, points_)
  Dist1 = DijkstraDec$distance[ Graph$edges[idx, 1]]
  Dist2 = DijkstraDec$distance[ Graph$edges[idx, 2]]
  Dist = vector(mode="numeric", length=length(x))
  Previous = vector(mode="integer", length(x))
  
  for( i in 1:length(x)){
    if( Dist1[i] < Dist2[i] ){
      Dist[i] = Dist1[i]
      Previous[i] = Graph$edges[idx[i], 1]
    }else{
      Dist[i] = Dist2[i]
      Previous[i] = Graph$edges[idx[i], 2]
    }
  }
  
  delta = Dist + sqrt( (Graph$nodes[Previous ,1] - points_[,1])^2 + (Graph$nodes[Previous, 2] - points_[,2])^2) 
  
  return ( 1/(sigma*(2*pi)^1.5) * exp(-delta^2/(2*sigma^2)))
  
}

density <- linfun(aux.600, L)
plot(density)

integral.600 = integrate_f(FEM(coef.600, FEMbasis.fdaPDE))
coef.ex = coef.600 / integral.600 
R_plot_graph.ggplot2(FEM(coef.ex, FEMbasis.fdaPDE))

R0 = R_mass_1D(FEMbasis)
R1 = R_stiff_1D(FEMbasis, M)
epsilon = 0.0001 

###
#nobs = c(1000, 1500, 2000, 2500)
nobs = c(600)
N = length(nobs)

lambda.1 = seq(from=1e-4, to=1, length.out=20)
lambda = seq(from=0.005, to=0.05, length.out=10)


lambda = c(0.01, 0.025, 0.05, 0.075, 0.09, 0.1, 0.25, 0.5, 0.75, 1)
niter = 2500
nfolds=10
eta=0.001
#uniform
#cvec= log( matrix(1e-3,nrow=nnodes,ncol=1))
#random
#cvec= log( 1e-3*matrix( abs(rnorm(nnodes)), nrow=nnodes, ncol=1) )
#integrate_f(FEM(exp(cvec), FEMbasis.fdaPDE))

M=1
M=3
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
    PP = rlpp(nobs[i], density)
    points = cbind(PP$data$x, PP$data$y)
    
    #nb. nel tempo di fdaPDE è anche compresa la crossvalidazione.
    start.fdaPDE <- Sys.time()
    sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=lambda, 
                        nsimulations=niter, 
                        stepProposals = eta,
                       # fvec = exp(cvec),
                        tol1 = 1e-8,
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
    sol = my_density.aldo(points, mesh, PHI, FEMbasis, R0, R1, sol.fdaPDE$lambda, eta,
                          cvec, niter, FALSE, tolerance=1e-8 )
    end <- Sys.time()
    
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
  }
}
tot.end = Sys.time()
tot.time = difftime(tot.end, tot.start, units="mins")


today_ = Sys.Date()
init_ ="-heat" #"-random"#"-null_vector" 
ntest_ = "-test-3-SimplifiedCV"

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

