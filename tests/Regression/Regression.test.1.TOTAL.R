###########################################################
##          TEST 1 - SIMPLENET - NO COVARIATES           ##
##                   delta varying                       ##
##            nobs = c(200, 300, 400, 500)               ##
##      delta = c(0.03, 0.015, 0.0075, 0.00375)          ##
###########################################################

source("~/Scrivania/fdaPDE/tests/Auxiliary/Regressione_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")

vertices = cbind(simplenet$vertices$x, simplenet$vertices$y)
edges = cbind(simplenet$from, simplenet$to)
M = simplenet$m
find.boundary(M)
L = as.linnet(simplenet)

delta = c(0.03, 0.015, 0.0075, 0.00375)
nnodes = c(97,194,388,775)

nobs =c(200, 300, 400, 500)
K = length(nobs)
N=length(delta)
M=20

coef.exs = list( matrix(0,nrow=97, ncol=1),
                 matrix(0,nrow=194,ncol=1),
                 matrix(0,nrow=388,ncol=1),
                 matrix(0,nrow=775,ncol=1))

mise = list()
mise.fdaPDE = list()
err.L2 = list()
times.fdaPDE = list()
times = list()
FEMbasis.fdaPDEs = list()
sols = list()
sols.fdaPDE = list()

options(digits.sec=6)
lambda_vector = seq(from=1e-6, to=1, length.out=10)
epsilon=0.0001


for( k in 1:K){
  print(paste("######## k = ", k ,"/",K," ########"))
  mise[[k]] = matrix(0,nrow=M,ncol=N)
  mise.fdaPDE[[k]] = matrix(0,nrow=M,ncol=N)
  # (||f_{fdaPDE} - f{R} ||_{L2} )^2
  err.L2[[k]] = matrix(0,nrow=M,ncol=N)
  times[[k]] = matrix(0,nrow=M, ncol=N)
  times.fdaPDE[[k]] = matrix(0,nrow=M,ncol=N)
  
  sols.fdaPDE[[k]] = list(list(), list(), list(), list())
  sols[[k]] = list(list(), list(), list(), list())

  for( i in 1:N){
  print(paste("######## i = ", i ,"/",N," ########"))

    mesh = create.mesh.1D.vertices(vertices, edges, delta[i])
    FEMbasis = create.FEM.basis.1D(mesh)
    
    ### fdaPDE - MESH ### 
    mesh.fdaPDE = create.mesh.1.5D(mesh$nodes,mesh$segments)
    FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)
    idx.leaf = which( mesh.fdaPDE$nodesmarkers==TRUE)
  ### CAMPO f ### 
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
      
    my_dens <- linfun(aux, L)
    coef.ex = aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])
    
    sols.fdaPDE[[k]][[i]] = matrix(0,nrow= nrow(mesh.fdaPDE$nodes), ncol=M)
    sols[[k]][[i]] = matrix(0,nrow= nrow(mesh.fdaPDE$nodes), ncol=M)
    
    if( k == 1){
      coef.exs[[i]] = coef.ex
      FEMbasis.fdaPDEs[[i]] = FEMbasis.fdaPDE
    }
    
    
    for( j in 1:M){
      PPP = runiflpp(nobs[[k]], L)
      points = cbind(PPP$data$x, PPP$data$y)
      
      exact_sol = my_dens(points[,1], points[,2])
      temp = range(exact_sol)
      my_range = temp[2]-temp[1]
      EPS = rnorm(nobs, 0, sd=0.05*my_range )
      Z = exact_sol + EPS
    
      start.fdaPDE = Sys.time()
      output.fdaPDE = smooth.FEM( observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=NULL,
                                locations=points,
                                lambda=lambda_vector, 
                                lambda.selection.criterion='grid', 
                                DOF.evaluation='exact', 
                                lambda.selection.lossfunction='GCV' )
      end.fdaPDE = Sys.time()
    
      lambda = output.fdaPDE$optimization$lambda_solution
    
      output.fdaPDE = smooth.FEM(observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=NULL,
                               locations=points,
                               lambda=lambda)
    
      #result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    #PHI = result$PHI
    
      start.time <- Sys.time()
      result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
      PHI = result$PHI
      output = no.covariate.FEM.1D( Z, FEMbasis, PHI, lambda)
      end.time  <-Sys.time()
    
      times[[k]][j,i] = difftime(end.time, start.time, units="secs")
    
      times.fdaPDE[[k]][j,i] = output.fdaPDE$time
      
      coef = output$f_coeff
      coef.fdaPDE     = output.fdaPDE$fit.FEM$coeff
    
      sols[[k]][[i]][,j] =  coef
      sols.fdaPDE[[k]][[i]][,j] = coef.fdaPDE
    
      mise[[k]][j,i]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
      mise.fdaPDE[[k]][j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
      err.L2[[k]][j,i] = integrate_f( FEM((coef.fdaPDE-coef)^2,FEMbasis.fdaPDE) )
  }
  
  }
  save(nobs,M,N,K, delta, nnodes,
       mise, mise.fdaPDE, err.L2,
       sols, sols.fdaPDE, coef.exs,
       times,times.fdaPDE, file =save.file)
  
}

today_ = Sys.Date()
ntest_ = "-test-1-TOTAL-no-cov" # "-test.1.3" 
#today_ = "2021-11-15"
file.name = paste("Regression-",today_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

save(nobs,M,N,K, delta, nnodes,
     mise, mise.fdaPDE, err.L2,
     sols, sols.fdaPDE, coef.exs,
     times,times.fdaPDE, file =save.file)
