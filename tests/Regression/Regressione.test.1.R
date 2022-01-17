### Confronto Mattina - fdaPDE ###
source("~/Scrivania/fdaPDE/tests/Auxiliary/Regressione_Mattina.R")
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
find.boundary(M)
L = as.linnet(simplenet)

delta = 0.03
mesh = create.mesh.1D.vertices(vertices, edges, delta)
FEMbasis = create.FEM.basis.1D(mesh)


### fdaPDE - MESH ### 
mesh.fdaPDE = create.mesh.1.5D(mesh$nodes,mesh$segments)
FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)
idx.leaf = which( mesh.fdaPDE$nodesmarkers==TRUE)

plot(mesh.fdaPDE)
points(mesh.fdaPDE$nodes[8,1], mesh.fdaPDE$nodes[8,2], pch=16)
points(mesh.fdaPDE$nodes[37,1], mesh.fdaPDE$nodes[37,2], pch=16)

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

integral_ex = integrate_f( FEM(aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2]),FEMbasis.fdaPDE) )

my_dens <- linfun(aux, L)

coef.ex = aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])
R_plot_graph.ggplot2(FEM(coef.ex,FEMbasis.fdaPDE))

nobs = c(200,300,400,500)
#nobs = c(200,300,400,500,600,700,800) #test 2
#nobs  = c(100,150,200,250,300,350) #test.??
N=length(nobs)
M=20

beta1 = -0.5 #test-1
beta2 = 0.2  #test-1

#beta1 = 0.5 #test 1.4
#beta2 = -0.2# test 1.4

beta_ex = rbind(beta1,beta2)
# || beta.fdaPDE - beta || 
rmse.beta.1 = matrix(0,nrow=M,ncol=N)
rmse.beta.2 = matrix(0,nrow=M,ncol=N)

mise = matrix(0,nrow=M,ncol=N)
mise.fdaPDE = matrix(0,nrow=M,ncol=N)
# (||f_{fdaPDE} - f{R} ||_{L2} )^2
err.L2 = matrix(0,nrow=M,ncol=N)

rmse.fdaPDE.beta.1 = matrix(0,nrow=M,ncol=N)
rmse.fdaPDE.beta.2 = matrix(0,nrow=M,ncol=N)
norms.beta1.inf = matrix(0,nrow=M,ncol=N)
norms.beta1.2 =matrix(0,nrow=M,ncol=N)
norms.beta2.inf = matrix(0,nrow=M,ncol=N)
norms.beta2.2 =matrix(0,nrow=M,ncol=N)

sols = array(0,dim=c(M,N,mesh$nnodes) ) #contiene vettore dei coeff
sols.fdaPDE = array(0,dim=c(M,N,mesh$nnodes) )

times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)
times.fdaPDE.2 = matrix(0,nrow=M,ncol=N)

options(digits.sec=6)
lambda_vector = c(10^-10,10^-9,10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,10,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10)
epsilon=0.0001
for( i in 1:N){
  for( j in 1:M){
    PPP = runiflpp(nobs[i], L)
    points = cbind(PPP$data$x, PPP$data$y)
  
   ### observations ###
    W1 = rnorm(nobs[i], 0, sd=1.)
    W2 = rnorm(nobs[i], 0, sd=0.5)
 #  W1 = sin(points[,1]) # test 1.4
 #  W2 = cos(points[,2]) # test 1.4
    W = NULL
    W = cbind(W1,W2)
    
    exact_sol = my_dens(points[,1], points[,2])
    temp = range(exact_sol)
    my_range = temp[2]-temp[1]
    EPS = rnorm(nobs[i], 0, sd=0.05*my_range )
    Z = beta1*W1 + beta2*W2 + exact_sol + EPS
    
    start.fdaPDE = Sys.time()
    output.fdaPDE = smooth.FEM( observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W,
                              locations=points,
                              lambda=lambda_vector, 
                              lambda.selection.criterion='grid', 
                              DOF.evaluation='exact', 
                              lambda.selection.lossfunction='GCV' )
    end.fdaPDE = Sys.time()
    
    lambda = output.fdaPDE$optimization$lambda_solution
    
    output.fdaPDE = smooth.FEM(observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W,
                               locations=points,
                               lambda=lambda)
    
    #result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    #PHI = result$PHI
  
    start.time <- Sys.time()
    result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    PHI = result$PHI
    output = covariate.FEM.1D(  W, Z, FEMbasis, PHI, lambda)
    end.time  <-Sys.time()
    
    rmse.beta.1[j,i] = norm(beta_ex[1]- output$beta_estimate[1],type="2")
    rmse.beta.2[j,i] = norm(beta_ex[2]- output$beta_estimate[2],type="2")
    rmse.fdaPDE.beta.1[j,i] = norm(beta_ex[1]- output.fdaPDE$solution$beta[1],type="2")
    rmse.fdaPDE.beta.2[j,i] = norm(beta_ex[2]- output.fdaPDE$solution$beta[2],type="2")

    norms.beta1.inf[j,i] = abs( output$beta_estimate[1] - output.fdaPDE$solution$beta[1])
    norms.beta2.inf[j,i] = abs( output$beta_estimate[2] - output.fdaPDE$solution$beta[2])
    norms.beta1.2[j,i]   = norm(output$beta_estimate[1] - output.fdaPDE$solution$beta[1],type="2")
    norms.beta2.2[j,i]   = norm(output$beta_estimate[2] - output.fdaPDE$solution$beta[2],type="2")
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
ntest_ = "-test-1" # "-test.1.3" 
#today_ = "2021-11-15"
file.name = paste("Regression-",today_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

save(nobs,M,N,
     rmse.beta.1,rmse.beta.2,
     rmse.fdaPDE.beta.1, rmse.fdaPDE.beta.2,
     mise, mise.fdaPDE, err.L2,
     sols,sols.fdaPDE,
     norms.beta1.inf,norms.beta1.2,
     norms.beta2.inf,norms.beta2.2,
     times,times.fdaPDE, file =save.file)

###########################################################
### nobs = 500, delta = c(0.03, 0.015, 0.0075, 0.00375) ###
###########################################################

vertices_x = as.matrix(simplenet$vertices$x)
vertices_y = as.matrix(simplenet$vertices$y)
vertices = as.matrix(cbind(vertices_x, vertices_y))
start_coord = cbind(simplenet$lines$ends$x0, simplenet$lines$ends$y0)
end_coord = cbind(simplenet$lines$ends$x1, simplenet$lines$ends$y1)
start = simplenet$from 
end = simplenet$to
edges = cbind(start, end)
M = simplenet$m
find.boundary(M)
L = as.linnet(simplenet)

delta = c(0.03, 0.015, 0.0075, 0.00375, 0.001875)

coef.ex = aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])

nobs =500

N=length(delta)
M=20

beta1 = -0.5 #test-1
beta2 = 0.2  #test-1

#beta1 = 0.5 #test 1.4
#beta2 = -0.2# test 1.4

beta_ex = rbind(beta1,beta2)
# || beta.fdaPDE - beta || 
rmse.beta.1 = matrix(0,nrow=M,ncol=N)
rmse.beta.2 = matrix(0,nrow=M,ncol=N)

mise = matrix(0,nrow=M,ncol=N)
mise.fdaPDE = matrix(0,nrow=M,ncol=N)
# (||f_{fdaPDE} - f{R} ||_{L2} )^2
err.L2 = matrix(0,nrow=M,ncol=N)

rmse.fdaPDE.beta.1 = matrix(0,nrow=M,ncol=N)
rmse.fdaPDE.beta.2 = matrix(0,nrow=M,ncol=N)
norms.beta1.inf = matrix(0,nrow=M,ncol=N)
norms.beta1.2 =matrix(0,nrow=M,ncol=N)
norms.beta2.inf = matrix(0,nrow=M,ncol=N)
norms.beta2.2 =matrix(0,nrow=M,ncol=N)

times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)
times.fdaPDE.2 = matrix(0,nrow=M,ncol=N)

options(digits.sec=6)
lambda_vector = c(10^-10,10^-9,10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,10,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10)
epsilon=0.0001
nnodes = vector(mode="integer",length=N)

for( i in 1:N){
  
  mesh = create.mesh.1D.vertices(vertices, edges, delta[i])
  FEMbasis = create.FEM.basis.1D(mesh)
  nnodes[i] = mesh$nnodes
  
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
  
  for( j in 1:M){
    PPP = runiflpp(nobs, L)
    points = cbind(PPP$data$x, PPP$data$y)
    
    ### observations ###
    W1 = rnorm(nobs, 0, sd=1.)
    W2 = rnorm(nobs, 0, sd=0.5)
    #  W1 = sin(points[,1]) # test 1.4
    #  W2 = cos(points[,2]) # test 1.4
    W = NULL
    W = cbind(W1,W2)
    
    temp = range(beta1*W1 + beta2*W2 + my_dens(points[,1], points[,2]))
    my_range = temp[2]-temp[1]
    EPS = rnorm(nobs, 0, sd=0.05*my_range )
    Z = beta1*W1 + beta2*W2 + my_dens(points[,1], points[,2]) + EPS
    
    start.fdaPDE = Sys.time()
    output.fdaPDE = smooth.FEM( observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W,
                                locations=points,
                                lambda=lambda_vector, 
                                lambda.selection.criterion='grid', 
                                DOF.evaluation='exact', 
                                lambda.selection.lossfunction='GCV' )
    end.fdaPDE = Sys.time()
    
    lambda = output.fdaPDE$optimization$lambda_solution
    
    output.fdaPDE = smooth.FEM(observations=Z, FEMbasis=FEMbasis.fdaPDE, covariates=W,
                               locations=points,
                               lambda=lambda)
    
    #result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    #PHI = result$PHI
    
    start.time <- Sys.time()
    result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    PHI = result$PHI
    output = covariate.FEM.1D(  W, Z, FEMbasis, PHI, lambda)
    end.time  <-Sys.time()
    
    rmse.beta.1[j,i] = norm(beta_ex[1]- output$beta_estimate[1],type="2")
    rmse.beta.2[j,i] = norm(beta_ex[2]- output$beta_estimate[2],type="2")
    rmse.fdaPDE.beta.1[j,i] = norm(beta_ex[1]- output.fdaPDE$solution$beta[1],type="2")
    rmse.fdaPDE.beta.2[j,i] = norm(beta_ex[2]- output.fdaPDE$solution$beta[2],type="2")
    
    norms.beta1.inf[j,i] = abs( output$beta_estimate[1] - output.fdaPDE$solution$beta[1])
    norms.beta2.inf[j,i] = abs( output$beta_estimate[2] - output.fdaPDE$solution$beta[2])
    norms.beta1.2[j,i]   = norm(output$beta_estimate[1] - output.fdaPDE$solution$beta[1],type="2")
    norms.beta2.2[j,i]   = norm(output$beta_estimate[2] - output.fdaPDE$solution$beta[2],type="2")
    times[j,i] = difftime(end.time, start.time, units="secs")
    
    times.fdaPDE[j,i] = output.fdaPDE$time
    times.fdaPDE.2[j,i] = difftime(end.fdaPDE, start.fdaPDE, units="secs")
    coef = output$f_coeff
    coef.fdaPDE     = output.fdaPDE$fit.FEM$coeff
    
    #sols[j,i,] =  coef
    #sols.fdaPDE[j,i,] = coef.fdaPDE
    
    mise[j,i]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
    mise.fdaPDE[j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
    err.L2[j,i] = integrate_f( FEM((coef.fdaPDE-coef)^2,FEMbasis.fdaPDE) )
  }
}

today_ = Sys.Date()
ntest_ = "-test-1-delta-var" # "-test.1.3" 
#today_ = "2021-11-15"
file.name = paste("Regression-",today_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

save(nobs,M,N, delta, nnodes,
     rmse.beta.1,rmse.beta.2,
     rmse.fdaPDE.beta.1, rmse.fdaPDE.beta.2,
     mise, mise.fdaPDE, err.L2,
     norms.beta1.inf,norms.beta1.2,
     norms.beta2.inf,norms.beta2.2,
     times,times.fdaPDE, file =save.file)


today_= "2022-01-11"
ntest_="-test-1-true-field-mesh"
file.name = paste("Regression-",today_,ntest_,sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")
img.file

MyTheme <- theme(
  axis.text = element_text(size=16),
  axis.title = element_text(size=16),
  title = element_text(size=20),
  legend.text = element_text(size=14),
  legend.key.size = unit(1,"cm") 
)
pdf(img.file, width=8, height=8)
plot(mesh.fdaPDE)
R_plot_graph.ggplot2(FEM(coef.ex,FEMbasis.fdaPDE))
dev.off()

##########################################################################
############### NO COVARIATES ############################################
##########################################################################

source("~/Scrivania/fdaPDE/tests/Regressione/Regressione_Mattina.R")

source("~/Scrivania/fdaPDE/tests/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Mesh_Evaluator/isInside.R")
source("tests/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Dijkstra.R")

vertices_x = as.matrix(simplenet$vertices$x)
vertices_y = as.matrix(simplenet$vertices$y)
vertices = as.matrix(cbind(vertices_x, vertices_y))
start_coord = cbind(simplenet$lines$ends$x0, simplenet$lines$ends$y0)
end_coord = cbind(simplenet$lines$ends$x1, simplenet$lines$ends$y1)
start = simplenet$from 
end = simplenet$to
edges = cbind(start, end)
M = simplenet$m
find.boundary(M)
L = as.linnet(simplenet)

delta = 0.03
mesh = create.mesh.1D.vertices(vertices, edges, delta)
FEMbasis = create.FEM.basis.1D(mesh)

### fdaPDE - MESH ### 
mesh.fdaPDE = create.mesh.1.5D(mesh$nodes,mesh$segments)
FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)
idx.leaf = which( mesh.fdaPDE$nodesmarkers==TRUE)

plot(mesh.fdaPDE)
points(mesh.fdaPDE$nodes[8,1], mesh.fdaPDE$nodes[8,2], pch=16)
points(mesh.fdaPDE$nodes[37,1], mesh.fdaPDE$nodes[37,2], pch=16)

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

integral_ex = integrate_f( FEM(aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2]),FEMbasis.fdaPDE) )

my_dens <- linfun(aux, L)

coef.ex = aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])
R_plot_graph.ggplot2(FEM(coef.ex,FEMbasis.fdaPDE))

nobs = c(200,300,400,500)
#nobs = c(200,300,400,500,600,700,800) #test 2
#nobs  = c(100,150,200,250,300,350) #test.??
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
lambda_vector = c(10^-10,10^-9,10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,10,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10)
lambda_vector = seq(from=1e-6, to=1, length.out=10)
epsilon=0.0001
for( i in 1:N){
  for( j in 1:M){
    PPP = runiflpp(nobs[i], L)
    points = cbind(PPP$data$x, PPP$data$y)
    
    exact_sol = my_dens(points[,1], points[,2])
    temp = range(exact_sol)
    my_range = temp[2]-temp[1]
    EPS = rnorm(nobs[i], 0, sd=0.05*my_range )
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
ntest_ = "-test-1-no-cov" # "-test.1.3" 
#today_ = "2021-11-15"
file.name = paste("Regression-",today_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

save(nobs,M,N,
     mise, mise.fdaPDE, err.L2,
     sols,sols.fdaPDE,
     times,times.fdaPDE, file =save.file)

###### delta varying #####
### nobs = 500, delta = c(0.03, 0.015, 0.0075, 0.00375) ###
###########################################################

vertices_x = as.matrix(simplenet$vertices$x)
vertices_y = as.matrix(simplenet$vertices$y)
vertices = as.matrix(cbind(vertices_x, vertices_y))
start_coord = cbind(simplenet$lines$ends$x0, simplenet$lines$ends$y0)
end_coord = cbind(simplenet$lines$ends$x1, simplenet$lines$ends$y1)
start = simplenet$from 
end = simplenet$to
edges = cbind(start, end)
M = simplenet$m
find.boundary(M)
L = as.linnet(simplenet)

delta = c(0.03, 0.015, 0.0075, 0.00375)

coef.ex = aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])

nobs =500

N=length(delta)
M=20

mise = matrix(0,nrow=M,ncol=N)
mise.fdaPDE = matrix(0,nrow=M,ncol=N)
# (||f_{fdaPDE} - f{R} ||_{L2} )^2
err.L2 = matrix(0,nrow=M,ncol=N)

rmse.fdaPDE.beta.1 = matrix(0,nrow=M,ncol=N)
rmse.fdaPDE.beta.2 = matrix(0,nrow=M,ncol=N)
norms.beta1.inf = matrix(0,nrow=M,ncol=N)
norms.beta1.2 =matrix(0,nrow=M,ncol=N)
norms.beta2.inf = matrix(0,nrow=M,ncol=N)
norms.beta2.2 =matrix(0,nrow=M,ncol=N)

times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)
times.fdaPDE.2 = matrix(0,nrow=M,ncol=N)

options(digits.sec=6)
lambda_vector = c(10^-10,10^-9,10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,10,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10)
lambda_vector = seq(from=1e-6, to=1, length.out=10)
epsilon=0.0001
nnodes = vector(mode="integer",length=N)

for( i in 1:N){
  
  mesh = create.mesh.1D.vertices(vertices, edges, delta[i])
  FEMbasis = create.FEM.basis.1D(mesh)
  nnodes[i] = mesh$nnodes
  
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
  
  for( j in 1:M){
    PPP = runiflpp(nobs, L)
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
    
    times[j,i] = difftime(end.time, start.time, units="secs")
    
    times.fdaPDE[j,i] = output.fdaPDE$time
    times.fdaPDE.2[j,i] = difftime(end.fdaPDE, start.fdaPDE, units="secs")
    coef = output$f_coeff
    coef.fdaPDE     = output.fdaPDE$fit.FEM$coeff
    
    #sols[j,i,] =  coef
    #sols.fdaPDE[j,i,] = coef.fdaPDE
    
    mise[j,i]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
    mise.fdaPDE[j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
    err.L2[j,i] = integrate_f( FEM((coef.fdaPDE-coef)^2,FEMbasis.fdaPDE) )
  }
}

today_ = Sys.Date()
ntest_ = "-test-1-no-cov-delta-var" # "-test.1.3" 
#today_ = "2021-11-15"
file.name = paste("Regression-",today_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

save(nobs,M,N, delta, nnodes,
     mise, mise.fdaPDE, err.L2,
     times,times.fdaPDE, file =save.file)
