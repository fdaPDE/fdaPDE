library(testthat)
library(fdaPDE)

test_that("tSR-PDE - Linear Network", {

options(warn=-1)
foldername <- test_path("../data/tSR-PDE/test_12")

set.seed(0)
eps = 1 / 2
x = c(0., 1)
y = c(0.,eps)
vertices = expand.grid(x,y)
vertices = cbind(vertices[,1], vertices[,2])
edges = matrix(c(1,2,1,3,3,4), nrow=3,ncol=2, byrow=T)

mesh = create.mesh.1.5D(vertices, edges)
mesh = refine.mesh.1.5D(mesh,delta=0.0125)

nnodes=dim(mesh$nodes)[1]
FEMbasis=create.FEM.basis(mesh)
# Locations at nodes
nodesLocations=mesh$nodes

# Exact data - Locations at nodes
nnodes = nrow(mesh$nodes)
TimeNodes = seq(0,pi,length.out = 6)
locations = cbind(rep(TimeNodes,each=nnodes),rep(nodesLocations[,1],length(TimeNodes)),rep(nodesLocations[,2],length(TimeNodes)) )

aux.4 = function(x,y){
  h = 1
  source = 4 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,1] - mesh$nodes[source,1])
    if(delta < h ){
      #coef[i] = 1- 8*delta
      #coef[i] = -4/h^2 * (delta^2 - h^2) - 2
      coef[i] = 3*exp(1/(delta^2 / h^2 -1)+1 ) - 1
      
    }
  }
  
  return(coef)
}
aux.3 = function(x,y){
  
  h = eps
  source = 3 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,1]-mesh$nodes[source,1]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,2] - mesh$nodes[source,2])
    if(delta < h ){
      coef[i] = -1 - 1/h*delta
    }
    
    
  }
  return(coef)
}
aux.1 = function(x,y){
  
  h = 1
  source = 1 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,1] - mesh$nodes[source,1])
    if(delta <= h ){
      coef[i] = -2 - 1/h*delta
    }
  }
  return(coef)
  
}  
func = function(x,y,t){
  
  res = (aux.1(x,y) + aux.3(x,y) + aux.4(x,y)) * cos(t)
  return(res)
}

func_evaluation = func(locations[,2], locations[,3], locations[,1])

lambdaS = 10^seq(-4, -2, by=0.5)
lambdaT = 10^seq(-4, -2, by=0.5)

lambdaS_par=10^seq(-5, -4, 0.25)
lambdaT_par=10^seq(-1, 0, 0.2)

cov1=4*sin(2*pi*locations[,2])*cos(2*pi*locations[,3])
cov2=rnorm(nnodes*length(TimeNodes), mean=3, sd=0.1)*rep(exp(-TimeNodes/length(TimeNodes)),each=nnodes)
W=cbind(cov1,cov2)

# Fix betas
beta_exact=c(1,-0.5)

ran=range(func_evaluation)

data = func_evaluation +rnorm(nnodes*length(TimeNodes),mean=0,sd=0.05*(ran[2]-ran[1]))

ran = range(func_evaluation+ W%*%beta_exact)
datacov=func_evaluation+ W%*%beta_exact +rnorm(nnodes*length(TimeNodes),mean=0,sd=0.05*(ran[2]-ran[1]))

data = matrix(data,nrow(mesh$nodes),length(TimeNodes))
datacov = matrix(datacov,nrow(mesh$nodes),length(TimeNodes))

#########################################SEPARABLE####################################################
invisible(capture.output(sol <- smooth.FEM.time(observations=data,
                         FEMbasis = FEMbasis, time_mesh = TimeNodes, time_locations = TimeNodes,
                         lambdaS = lambdaS, lambdaT = lambdaT,
                         FLAG_PARABOLIC = FALSE)))

load(file=paste0(foldername,"/test_12_1.RData"))
expect_equal( max(abs((sol$fit.FEM.time$coeff-output_CPP$fit.FEM.time$coeff))) < 1e-8, TRUE);

invisible(capture.output(sol <- smooth.FEM.time(observations=datacov, covariates = W,
                            FEMbasis = FEMbasis, time_mesh = TimeNodes,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            FLAG_PARABOLIC = FALSE)))

load(file=paste0(foldername,"/test_12_2.RData"))
expect_equal( max(abs((sol$fit.FEM.time$coeff-output_CPP$fit.FEM.time$coeff))) < 1e-8, TRUE);
expect_equal( max(abs((sol$solution$beta-output_CPP$solution$beta))) < 1e-8, TRUE);
##########################################PARABOLIC####################################################
### MONOLITIC METHOD
invisible(capture.output(sol <- smooth.FEM.time(observations=data,
                         FEMbasis = FEMbasis, time_mesh = TimeNodes, 
                         time_locations = TimeNodes,
                         lambdaS = lambdaS, lambdaT = lambdaT,
                         FLAG_PARABOLIC = TRUE)))

load(file=paste0(foldername,"/test_12_3.RData"))
expect_equal( max(abs((sol$fit.FEM.time$coeff-output_CPP$fit.FEM.time$coeff))) < 1e-8, TRUE);

invisible(capture.output(sol <- smooth.FEM.time(observations=datacov[,2:length(TimeNodes)], 
                                                covariates = W[(1+nrow(mesh$nodes)):(length(TimeNodes)*nrow(mesh$nodes)),],
                            FEMbasis = FEMbasis, time_mesh = TimeNodes,
                            lambdaS = lambdaS, lambdaT = lambdaT,
                            IC=func_evaluation[1:nrow(mesh$nodes)],
                            FLAG_PARABOLIC = TRUE)))

load(file=paste0(foldername,"/test_12_4.RData"))
expect_equal( max(abs((sol$fit.FEM.time$coeff-output_CPP$fit.FEM.time$coeff))) < 1e-8, TRUE);
expect_equal( max(abs((sol$solution$beta-output_CPP$solution$beta))) < 1e-8, TRUE);

### ITERATIVE METHOD
invisible(capture.output(sol <- smooth.FEM.time(observations=data,
                             FEMbasis = FEMbasis, time_mesh = TimeNodes, 
                             time_locations = TimeNodes,
                             lambdaS = lambdaS_par, lambdaT = lambdaT_par,
                             FLAG_PARABOLIC = TRUE, FLAG_ITERATIVE = TRUE)))

load(file=paste0(foldername,"/test_12_5.RData"))
expect_equal( max(abs((sol$solution$beta-output_CPP$solution$beta))) < 1e-8, TRUE);

invisible(capture.output(sol <- smooth.FEM.time(observations=datacov[,2:length(TimeNodes)], 
                                                covariates = W[(1+nrow(mesh$nodes)):(length(TimeNodes)*nrow(mesh$nodes)),],
                                FEMbasis = FEMbasis, time_mesh = TimeNodes,
                                lambdaS = lambdaS, lambdaT = lambdaT,
                                IC=func_evaluation[1:nrow(mesh$nodes)],
                                FLAG_PARABOLIC = TRUE , FLAG_ITERATIVE = TRUE)))

load(file=paste0(foldername,"/test_12_6.RData"))
expect_equal( max(abs((sol$solution$beta-output_CPP$solution$beta))) < 1e-8, TRUE);

})