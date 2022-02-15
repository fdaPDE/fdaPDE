###########################################################
##          TEST 2 - BRICKWALL - with COVARIATES         ##
##                   delta varying                       ##
##            nobs = c(300, 500, 700, 1000)              ##
##            delta = c(40, 30, 25, 17.5)                ##
###########################################################

source("~/Scrivania/fdaPDE/tests/Auxiliary/Regressione_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")

data("spiders")
vertices = cbind(spiders$domain$vertices$x, spiders$domain$vertices$y)
edges = cbind(spiders$domain$from, spiders$domain$to)
M = simplenet$m
find.boundary(M)
L = as.linnet(spiders$domain$m)

delta = c(40, 30, 25, 17.5)
nnodes = c(412, 637, 831, 1143) 

nobs =c(300, 500, 700, 1000)
K = length(nobs)
N=length(delta)
M=20

coef.exs = list( matrix(0,nrow=nnodes[1], ncol=1),
                 matrix(0,nrow=nnodes[2],ncol=1),
                 matrix(0,nrow=nnodes[3],ncol=1),
                 matrix(0,nrow=nnodes[4],ncol=1))

mise = list()
mise.fdaPDE = list()
err.L2 = list()
times.fdaPDE = list()
times = list()
rmse = list()
rmse.fdaPDE = list()
FEMbasis.fdaPDEs = list()
sols = list()
sols.fdaPDE = list()

beta_ex = 0.8
options(digits.sec=6)
lambda_vector = seq(from=1e-1, to=1000, length.out=10)
epsilon=0.0001
nnodes = vector(mode="integer",length=N)

for( k in 1:K){
  print(paste("######## k = ", k ,"/",K," ########"))
  
  mise[[k]] = matrix(0,nrow=M,ncol=N)
  mise.fdaPDE[[k]] = matrix(0,nrow=M,ncol=N)
  
  # (||f_{fdaPDE} - f{R} ||_{L2} )^2
  err.L2[[k]] = matrix(0,nrow=M,ncol=N)
  
  times[[k]] = matrix(0,nrow=M, ncol=N)
  times.fdaPDE[[k]] = matrix(0,nrow=M,ncol=N)
  
  rmse[[k]] = matrix(0,nrow=M, ncol=N)
  rmse.fdaPDE[[k]] = matrix(0,nrow=M,ncol=N)
  
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
    
    
    sols.fdaPDE[[k]][[i]] = matrix(0,nrow= nrow(mesh.fdaPDE$nodes), ncol=M)
    sols[[k]][[i]] = matrix(0,nrow= nrow(mesh.fdaPDE$nodes), ncol=M)
    
    if( k == 1){
      coef.exs[[i]] = coef.ex
      FEMbasis.fdaPDEs[[i]] = FEMbasis.fdaPDE
    }
    
    
    for( j in 1:M){
      PPP = runiflpp(nobs[[k]], L)
      points = cbind(PPP$data$x, PPP$data$y)
      
      W = matrix( rnorm(nobs[[k]], 0, sd=0.5), nrow=nobs[[k]], ncol=1)
      exact_sol = density(points[,1], points[,2])
      temp = range(exact_sol)
      my_range = temp[2]-temp[1]
      EPS = rnorm(nobs, 0, sd=0.05*my_range )
      Z = beta_ex*W + exact_sol + EPS
      
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
      
      rmse[[k]][j,i] = norm(beta_ex- output$beta_estimate[1],type="2")
      
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
ntest_ = "-test-2-TOTAL-" # "-test.1.3" 
#today_ = "2021-11-15"
file.name = paste("Regression-",today_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

save(nobs,M,N,K, delta, nnodes,
     mise, mise.fdaPDE, err.L2,
     sols, sols.fdaPDE, coef.exs,
     times,times.fdaPDE, file =save.file)
