### DE test-5 (?) ### Brickwall

source("~/Scrivania/fdaPDE/tests/StimaDensita/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Dijkstra.R")
source("~/Scrivania/fdaPDE/tests/Mesh_Evaluator/isInside.R")
source("tests/R_plot_graph.ggplot2.R")
library(ggplot2)

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

plot(mesh.fdaPDE, show.nodes = TRUE, pch=16)
points(mesh.fdaPDE$nodes[25,1], mesh.fdaPDE$nodes[25,2],pch=16, col = "red")

res.25 = Dijkstra(mesh.fdaPDE, source=25)
sigma = 300

aux.25 = function(x,y,seg,tp){ 
  sigma=300 
  DijkstraDec=res.25
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
  
  delta = (Dist + sqrt( (Graph$nodes[Previous ,1] - points_[,1])^2 + (Graph$nodes[Previous, 2] - points_[,2])^2) ) 
  
  return ( as.numeric( abs(y - Graph$nodes[25,2])  <= .Machine$double.eps)* 1/(sigma*(2*pi)^1.5) * exp(-delta^2/(2*sigma^2)) )
  
}

coef.25 = aux.25(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2]) 

density <- linfun(aux.25, L)
plot(density)

integral.25 = integrate_f(FEM(coef.25, FEMbasis.fdaPDE))
coef.ex = coef.25 / integral.25
integrate_f(FEM(coef.ex, FEMbasis.fdaPDE))
R_plot_graph.ggplot2(FEM(coef.ex, FEMbasis.fdaPDE))
R_plot_graph(FEM(coef.ex,FEMbasis.fdaPDE))

R0 = R_mass_1D(FEMbasis)
R1 = R_stiff_1D(FEMbasis, M)
epsilon = 0.0001 

###
#nobs = c(1000, 1500, 2000, 2500)
nobs = c(350,500)
N = length(nobs)
#cvec=matrix(0,nrow=mesh$nnodes,ncol=1)
# 
#lambda = seq(from=1e-4,to=1, length.out=20) 

lambda.350 = c(1e-4, 2.5e-4, 5e-4, 1e-3, 2.5e-3, 5e-3, 1e-2, 5e-2, 1e-1, 1)
#lambda = seq( from=2e-4,to=5e-2, length.out=10)

lambda.500= c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 5, 7.5, 10)

#lambda.1 = seq(from=1e-4, to=1, length.out=20)
#lambda = seq(from=1e-4, to=1, length.out=10)

niter = 5000
nfolds=10
eta=0.001

M=1

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
  if(nobs[i] == 350 )
    lambda = lambda.350
  else
    lambda = lambda.500
  
  for(j in 1:M){
    PP = rlpp(nobs[i], density)
    points = cbind(PP$data$x, PP$data$y)
    # 
    # if(nobs[i]<600){
    #   lambda=lambda.1
    # }else{
    #   lambda=lambda.2
    # }
    # 
    #nb. nel tempo di fdaPDE è anche compresa la crossvalidazione.
    start.fdaPDE <- Sys.time()
    sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=lambda, 
                        nsimulations=niter, 
                        stepProposals = eta,
                        # fvec = exp(cvec),
                        tol1 = 1e-6,
                        nfolds=nfolds, 
                        step_method="Fixed_Step", 
                        direction_method="Gradient",
                        preprocess_method = "RightCV")
    end.fdaPDE <- Sys.time()
    
    start <- Sys.time()
    result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    # NB. R need log(density) init.
    cvec = log( sol.fdaPDE$f_init[, which( abs(lambda - sol.fdaPDE$lambda) <= .Machine$double.eps )] )
    PHI = result$PHI
    sol = my_density.aldo(points, mesh, PHI, FEMbasis, R0, R1, sol.fdaPDE$lambda, eta,
                          cvec, niter, FALSE, tolerance=1e-6 )
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

tot.time
mise
mise.fdaPDE

today_ = Sys.Date()
init_ ="-heat" #"-random"#"-null_vector" 
ntest_ = "-test-5-RightCV"

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


mesh.ref = create.mesh.1D.vertices(vertices, edges, delta=30)
mesh.fdaPDE.ref = create.mesh.1.5D(mesh.ref$nodes, mesh.ref$segments)
FEMbasis.fdaPDE.ref = create.FEM.basis(mesh.fdaPDE.ref)

coef.fdaPDE.ref = eval.FEM( FEM(coef.fdaPDE, FEMbasis.fdaPDE), locations= mesh.fdaPDE.ref$nodes  )

R_plot_graph.ggplot2(FEM(coef.fdaPDE.ref,FEMbasis.fdaPDE.ref))
