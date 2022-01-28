# Geometry + Spatial Field 
{
  
  library(spatstat)
  library(purrr)
  source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
  source("tests/Auxiliary/R_plot_graph.ggplot2.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")
  
  FAMILY = "poisson"
  l<-make.link("log")
  link<-l$linkfun
  inv.link<-l$linkinv
  
  vertices = cbind(simplenet$vertices$x, simplenet$vertices$y)
  edges = cbind(simplenet$from , simplenet$to)
  L = as.linnet(simplenet)
  
  delta = 0.03 # 102 nodi, 0.015 199 nodi
  mesh =create.mesh.1.5D(vertices, edges)
  mesh =refine.mesh.1.5D(mesh,delta=0.03)
  
  nnodes = nrow(mesh$nodes)
  FEMbasis = create.FEM.basis(mesh)
  
  
  plot(mesh)
  points(vertices)
  points(mesh$nodes[4,1], mesh$nodes[4,2], pch=16, col="green4")
  points(mesh$nodes[5,1], mesh$nodes[5,2], pch=16, col="green4")
  #points(vertices[6,1], vertices[6,2], pch=16,col="green4")
  points(vertices[8,1], vertices[8,2], pch=16, col="green4")
  points(vertices[10,1], vertices[10,2], pch=16, col="green4")
  
  
  
  ### CAMPO f ### 
  dijkstra.4 = Dijkstra(mesh,4)
  dijkstra.5 = Dijkstra(mesh,5)
  dijkstra.6 = Dijkstra(mesh,6)
  dijkstra.8 = Dijkstra(mesh,8)
  dijkstra.10= Dijkstra(mesh,10)
  
  aux = function(x, y, seg, tp) { 
    
    # sigma.4 = 0.0075
    # sigma.5 = 0.0075
    # sigma.6 = 0.0075
    # sigma.8 = 0.0075
    # sigma.10 = 0.0075
    # 
    sigma.4 = 0.0075 * 2
    sigma.5 = 0.0075 * 2
    sigma.6 = 0.0075 * 2
    sigma.8 = 0.0075 * 2
    sigma.10 = 0.0075 * 2
    
    
    res.4 = equal_split_discontinous(mesh, sigma.4, dijkstra.4, x, y)
    res.5 = equal_split_discontinous(mesh, sigma.5, dijkstra.5, x, y) 
    res.6 = equal_split_discontinous(mesh, sigma.6, dijkstra.6, x, y)
    res.8 = equal_split_discontinous(mesh, sigma.8, dijkstra.8, x, y)
    res.10 = equal_split_discontinous(mesh, sigma.10, dijkstra.10, x, y)
    
    idx.4 = which(res.4$bandwidth==0.0)
    idx.5 = which(res.5$bandwidth==0.0)
    idx.6 = which(res.6$bandwidth==0.0)
    idx.8 = which(res.8$bandwidth==0.0)
    idx.10 = which(res.10$bandwidth==0.0)
    
    res=(res.4$coef + res.5$coef  + res.8$coef + res.10 $coef+ res.6$coef) 
    
    res = res/25
    
    intersect.4.5 = intersect(idx.4, idx.5)
    intersect.8.10 = intersect(idx.10, idx.8)
    intersect.tot = intersect( intersect.4.5, intersect.8.10)
    intersect.tot = intersect(intersect.tot, idx.6)
    
    res[intersect.tot] = -1.0
    return(res)
  }
  
  
  aux.2 = function(x,y){-3/2*sin(2*pi*x)*cos(2*pi*y)+ 2/5*sin(3*pi*x) + 2}
  coef.ex = aux.2(mesh$nodes[,1],mesh$nodes[,2])
  R_plot_graph.ggplot2(FEM(coef.ex,FEMbasis))
} 

# 2D Geometry 
x = seq(from=min(vertices[,1]), to=max(vertices[,1]), length.out=20)
y = seq(from=min(vertices[,2]), to=max(vertices[,2]), length.out=20)
x = seq(from=0, to=1, length.out=20)
y = seq(from=0, to=1, length.out=20)
nodes.2D = expand.grid(x,y)
mesh.2D = create.mesh.2D(nodes.2D)
nnodes.2D = nrow(mesh.2D$nodes)
FEMbasis.2D = create.FEM.basis(mesh.2D)
plot(mesh.2D)
points(mesh$nodes, pch=16, col="red3")
nrow(mesh.2D$nodes)

####
mesh.2D = create.mesh.2D( nodes = rbind(c(0,0), c(0,1), c(1,0), c(1,1)))
mesh.2D = refine.mesh.2D(mesh.2D, maximum_area = 0.001, minimum_angle=20)
FEMbasis.2D = create.FEM.basis(mesh.2D)
plot(mesh.2D)
points(mesh$nodes, pch=16, col="red3")

nrow(mesh.2D$nodes)

###
mesh.2D = create.mesh.2D(vertices)
plot(mesh.2D)
points(vertices, pch=16)
mesh.2D = refine.mesh.2D(mesh.2D, maximum_area = 1e-3)
nnodes.2D = nrow(mesh.2D$nodes)
FEMbasis.2D = create.FEM.basis(mesh.2D)
plot(mesh.2D)
points(mesh$nodes, pch=16, col="red3")
nrow(mesh.2D$nodes)

#nobs = c(400, 600, 800, 1000, 1200, 1400, 1600)
nobs=c(800, 1200, 1600, 2400)
nobs=c(600, 800)
N = length(nobs)
M = 5
times = matrix(0, nrow=M, ncol=N)
mise = matrix(0,nrow=M, ncol=N)
norm.l2 = matrix(0, nrow=M, ncol=N)
times.2D = matrix(0, nrow=M, ncol=N)
mise.2D = matrix(0,nrow=M, ncol=N)
norm.l2.2D = matrix(0,nrow=M, ncol=N)
sols = array(0, dim=c(M,N, nnodes))
sols.2D = array(0, dim=c(M,N, nnodes))
rmse.1 = matrix(0, nrow=M, ncol=N)
rmse.1.2D = matrix(0, nrow=M, ncol=N)
rmse.2 = matrix(0, nrow=M, ncol=N)
rmse.2.2D = matrix(0, nrow=M, ncol=N)

beta_ex = c(0.5,-0.2)

lambda = 10^seq(-5,-3,length.out = 20) # test-1 con e senza covariate
start.tot = Sys.time()
for(i in 1:N){
  for(j in 1:M){
    PP = runiflpp(nobs[i] , L)
    locations = cbind(PP$data$x, PP$data$y)
    
    sol_exact = aux(locations[,1], locations[,2])
    # covariates
    #set.seed(42)
    
    W=matrix(0,nrow=nobs[i], ncol=2)
    #W[,1]= rnorm(nobs[i], mean=2, sd=0.5) # test-3
    #W[,2]= rnorm(nobs[i], mean=0,sd=1)    # test-3
    W[,1]=rbeta(nobs[i],shape1=1.5,shape2=2)
    W[,2]= rbeta(nobs[i], shape1=3, shape2=2)
    #W = sin(2*pi * locations[,1]) * cos(2*pi* locations[,2])
    
    param =  sol_exact + W%*%beta_ex
    mu<-inv.link(param)
    range(mu)
    # sampling response:
    response <- rpois(nobs[i], lambda = mu)
    
    start = Sys.time()
    output_CPP <- fdaPDE::smooth.FEM(location = locations, 
                                     observations = as.numeric(response), 
                                     FEMbasis =FEMbasis, 
                                     covariates = W, #NULL
                                     max.steps=15, 
                                     fam=FAMILY,
                                     mu0=NULL, 
                                     scale.param=NULL,
                                     lambda = lambda, 
                                     lambda.selection.criterion = 'grid', 
                                     DOF.evaluation = 'exact', 
                                     lambda.selection.lossfunction = 'GCV')
    end = Sys.time()
    times[j,i] = difftime(end, start, units="secs")
    
    coef = output_CPP$fit.FEM$coeff[ , output_CPP$optimization$lambda_position]
    mise[j,i] = integrate_f( FEM( (coef.ex - coef)^2, FEMbasis) )
    
    rmse.1[j,i] = norm( beta_ex[1] - output_CPP$solution$beta[1, output_CPP$optimization$lambda_position ], type="2")
    rmse.2[j,i] = norm( beta_ex[2] - output_CPP$solution$beta[2, output_CPP$optimization$lambda_position ], type="2")
    
    sols[j,i,] = coef
    print(paste("lambda pos = ", output_CPP$optimization$lambda_position))
    
    start.2D = Sys.time()
    output_CPP.2D <- fdaPDE::smooth.FEM(location = locations, 
                                        observations = as.numeric(response), 
                                        FEMbasis =FEMbasis.2D, 
                                        covariates = W, # 
                                        max.steps=15, 
                                        fam=FAMILY,
                                        mu0=NULL, 
                                        scale.param=NULL,
                                        lambda = lambda, 
                                        lambda.selection.criterion = 'grid', 
                                        DOF.evaluation = 'exact', 
                                        lambda.selection.lossfunction = 'GCV')
    end.2D = Sys.time()
    
    times.2D[j,i] = difftime(end.2D, start.2D, units="secs")
    coef.2D = eval.FEM( FEM(output_CPP.2D$fit.FEM$coeff[ , output_CPP.2D$optimization$lambda_position],
                            FEMbasis.2D),mesh$nodes) 
    
    rmse.1.2D[j,i] = norm( beta_ex[1] - output_CPP$solution$beta[1, output_CPP.2D$optimization$lambda_position ], type="2")
    rmse.2.2D[j,i] = norm( beta_ex[2] - output_CPP$solution$beta[2, output_CPP.2D$optimization$lambda_position ], type="2")
    
    
    mise.2D[j,i] = integrate_f( FEM( (coef.ex - coef.2D)^2, FEMbasis) )
    sols.2D[j,i,] = coef.2D
    print(paste("lambda.2D pos = ", output_CPP.2D$optimization$lambda_position))
  }
}
end.tot = Sys.time()
tot.time = difftime(end.tot, start.tot, units="mins")
tot.time

today_ = Sys.Date()
ntest_ = "-test-4-with-cov"
path_ = "/home/aldo/Scrivania/fdaPDE-DATA/GSR"

file.name = paste("GSR-",today_,ntest_,sep="")

save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/GSR/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/GSR/",file.name,".pdf",sep="")

save(nobs, N, M, sols, sols.2D,
     mise, mise.2D, rmse.1, rmse.2, rmse.1.2D, rmse.2.2D,
     times, times.2D, file = save.file)

######################################################

library(ggplot2)
library(latex2exp)

MyTheme <- theme(
  axis.text = element_text(size=26),
  axis.title = element_text(size=26),
  title = element_text(size=26),
  plot.title = element_text(hjust = 0.5),
  legend.text = element_text(size=24),
  legend.key.size = unit(1,"cm"),
  legend.key.height = unit(1,"cm"),
  legend.title = element_blank(),
  legend.background = element_rect(fill="white", color="black",
                                   size=c(1,0.5))
)

# legend.position= c(0.825,0.15) # basso dx
# legend.position= c(0.825,0.85) # alto  dx

mise_ = cbind( as.vector(mise), as.vector(mise.2D))
mise_ = as.vector(mise_)
rmse.1_ = cbind(rmse.1, rmse.1.2D)
rmse.1_ = as.vector(rmse.1_)
rmse.2_ = cbind(rmse.2, rmse.2.2D)
rmse.2_ = as.vector(rmse.2_)
times_ = cbind( as.vector(times), as.vector(times.2D)) 
times_ = as.vector(times_)
nobs_ = rep(nobs,each=M)
nobs_ = rep(as.character(nobs_), times=2)
type_ = rep(c("1.5D","2D"),each=(M*N))

data_frame.1 = data.frame(mise_, rmse.1_, rmse.2_, times_,nobs_,type_)

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=mise_,group=interaction(nobs_,type_),fill=type_))+
  scale_x_discrete(limits=as.character(nobs))+
  labs(x="observations",y="",fill="",title="MISE" )+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.85,0.85))+
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=rmse.1_,group=interaction(nobs_,type_),fill=type_))+
  scale_x_discrete(limits=as.character(nobs))+
  labs(x="observations",y="",fill="",title=TeX("$RMSE^{2}$") )+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.85,0.85))+
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=rmse.2_,group=interaction(nobs_,type_),fill=type_))+
  scale_x_discrete(limits=as.character(nobs))+
  labs(x="observations",y="",fill="",title=TeX("$RMSE^{2}$") )+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.85,0.85))+
  MyTheme


##############################################################

{
  library(spatstat)
  library(purrr)
  source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
  source("tests/Auxiliary/R_plot_graph.ggplot2.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")
  
  FAMILY = "poisson"
  l<-make.link("log")
  link<-l$linkfun
  inv.link<-l$linkinv
  
  vertices = cbind(simplenet$vertices$x, simplenet$vertices$y)
  edges = cbind(simplenet$from , simplenet$to)
  L = as.linnet(simplenet)
  
  delta = 0.03 # 102 nodi, 0.015 199 nodi
  mesh =create.mesh.1.5D(vertices, edges)
  mesh =refine.mesh.1.5D(mesh,delta=0.03)
  
  nnodes = nrow(mesh$nodes)
  FEMbasis = create.FEM.basis(mesh)
  
  
  plot(mesh)
  points(vertices)
  points(mesh$nodes[4,1], mesh$nodes[4,2], pch=16, col="green4")
  points(mesh$nodes[5,1], mesh$nodes[5,2], pch=16, col="green4")
  #points(vertices[6,1], vertices[6,2], pch=16,col="green4")
  points(vertices[8,1], vertices[8,2], pch=16, col="green4")
  points(vertices[10,1], vertices[10,2], pch=16, col="green4")
  
  
  
  ### CAMPO f ### 
  dijkstra.4 = Dijkstra(mesh,4)
  dijkstra.5 = Dijkstra(mesh,5)
  #dijkstra.6 = Dijkstra(mesh,6)
  dijkstra.8 = Dijkstra(mesh,8)
  dijkstra.10= Dijkstra(mesh,10)
  
  aux = function(x, y, seg, tp) { 
    
    sigma = 0.03
    
    res.4 = equal_split_discontinous.mollifier(mesh, sigma, dijkstra.4, x, y)
    res.5 = equal_split_discontinous.mollifier(mesh, 0.025, dijkstra.5, x, y) 
    # res.6 = equal_split_discontinous.mollifier(mesh, sigma, dijkstra.6, x, y)
    res.8 = equal_split_discontinous.mollifier(mesh, 0.025, dijkstra.8, x, y)
    res.10 = equal_split_discontinous.mollifier(mesh, sigma, dijkstra.10, x, y)
    
    res= (res.4$coef + res.5$coef  + res.8$coef + res.10 $coef) # + res.6$coef 
    #res = res.4$coef + res.10$coef
    
    idx.4 = which(res.4$bandwidth==0.0)
    idx.5 = which(res.5$bandwidth==0.0)
    #idx.6 = which(res.6$bandwidth==0.0)
    idx.8 = which(res.8$bandwidth==0.0)
    idx.10 = which(res.10$bandwidth==0.0)
    
    intersect.4.5 = intersect(idx.4, idx.5)
    intersect.8.10 = intersect(idx.10, idx.8)
    intersect.tot = intersect( intersect.4.5, intersect.8.10)
    #intersect.tot = intersect(intersect.tot, idx.6)
    #intersect.tot = intersect(idx.4, idx.10)
    
    res[intersect.tot] = -1.0
    return(res)
  }
  
  coef.ex = aux(mesh$nodes[,1],mesh$nodes[,2])
  R_plot_graph.ggplot2.2(FEM(coef.ex,FEMbasis), 
                         color.min  =min(coef.ex),
                         color.max  =max(coef.ex))
} 

# 2D Geometry 
# x = seq(from=min(vertices[,1]), to=max(vertices[,1]), length.out=20)
# y = seq(from=min(vertices[,2]), to=max(vertices[,2]), length.out=20)
# x = seq(from=0, to=1, length.out=20)
# y = seq(from=0, to=1, length.out=20)
# nodes.2D = expand.grid(x,y)
# mesh.2D = create.mesh.2D(nodes.2D)
# nnodes.2D = nrow(mesh.2D$nodes)
# FEMbasis.2D = create.FEM.basis(mesh.2D)
# plot(mesh.2D)
# points(mesh$nodes, pch=16, col="red3")
# nrow(mesh.2D$nodes)

####
mesh.2D = create.mesh.2D( nodes = rbind(c(0,0), c(0,1), c(1,0), c(1,1)))
mesh.2D = refine.mesh.2D(mesh.2D, maximum_area = 0.001, minimum_angle=20)
FEMbasis.2D = create.FEM.basis(mesh.2D)
plot(mesh.2D)
points(mesh$nodes, pch=16, col="red3")
nrow(mesh.2D$nodes)
###
# mesh.2D = create.mesh.2D(vertices)
# plot(mesh.2D)
# points(vertices, pch=16)
# mesh.2D = refine.mesh.2D(mesh.2D, maximum_area = 1e-3)
# nnodes.2D = nrow(mesh.2D$nodes)
# FEMbasis.2D = create.FEM.basis(mesh.2D)
# plot(mesh.2D)
# points(mesh$nodes, pch=16, col="red3")
# nrow(mesh.2D$nodes)

#nobs = c(400, 600, 800, 1000, 1200, 1400, 1600)
nobs=c(400, 600, 800)
N = length(nobs)
M = 20
times = matrix(0, nrow=M, ncol=N)
mise = matrix(0,nrow=M, ncol=N)
norm.l2 = matrix(0, nrow=M, ncol=N)
times.2D = matrix(0, nrow=M, ncol=N)
mise.2D = matrix(0,nrow=M, ncol=N)
norm.l2.2D = matrix(0,nrow=M, ncol=N)
sols = array(0, dim=c(M,N, nnodes))
sols.2D = array(0, dim=c(M,N, nnodes))
rmse.1 = matrix(0, nrow=M, ncol=N)
rmse.1.2D = matrix(0, nrow=M, ncol=N)
rmse.2 = matrix(0, nrow=M, ncol=N)
rmse.2.2D = matrix(0, nrow=M, ncol=N)

beta_ex = c(0.5,-0.2)
lambda.opt = matrix(0,nrow=M, ncol=N)
lambda = 10^seq(-5,-3,length.out = 20) # test-1 con e senza covariate

start.tot = Sys.time()
for(i in 1:N){
  for(j in 1:M){
    PP = runiflpp(nobs[i] , L)
    locations = cbind(PP$data$x, PP$data$y)
    
    sol_exact = aux(locations[,1], locations[,2])
    # covariates
    #set.seed(42)
    
    W=matrix(0,nrow=nobs[i], ncol=2)
    #W[,1]= rnorm(nobs[i], mean=2, sd=0.5) # test-3
    #W[,2]= rnorm(nobs[i], mean=0,sd=1)    # test-3
    W[,1]=rbeta(nobs[i],shape1=1.5,shape2=2)
    W[,2]= rbeta(nobs[i], shape1=3, shape2=2)
    #W = sin(2*pi * locations[,1]) * cos(2*pi* locations[,2])
    
    param =  sol_exact + W%*%beta_ex
    mu<-inv.link(param)
    range(mu)
    # sampling response:
    response <- rpois(nobs[i], lambda = mu)
    
    start = Sys.time()
    output_CPP <- fdaPDE::smooth.FEM(location = locations, 
                                     observations = as.numeric(response), 
                                     FEMbasis =FEMbasis, 
                                     covariates = W, #NULL
                                     max.steps=15, 
                                     fam=FAMILY,
                                     mu0=NULL, 
                                     scale.param=NULL,
                                     lambda = lambda, 
                                     lambda.selection.criterion = 'grid', 
                                     DOF.evaluation = 'exact', 
                                     lambda.selection.lossfunction = 'GCV')
    end = Sys.time()
    times[j,i] = difftime(end, start, units="secs")
    
    coef = output_CPP$fit.FEM$coeff[ , output_CPP$optimization$lambda_position]
    mise[j,i] = integrate_f( FEM( (coef.ex - coef)^2, FEMbasis) )
    
    rmse.1[j,i] = norm( beta_ex[1] - output_CPP$solution$beta[1, output_CPP$optimization$lambda_position ], type="2")
    rmse.2[j,i] = norm( beta_ex[2] - output_CPP$solution$beta[2, output_CPP$optimization$lambda_position ], type="2")
    
    sols[j,i,] = coef
    lambda.opt[j,i] =  output_CPP$optimization$lambda_position
    
    start.2D = Sys.time()
    output_CPP.2D <- fdaPDE::smooth.FEM(location = locations, 
                                        observations = as.numeric(response), 
                                        FEMbasis =FEMbasis.2D, 
                                        covariates = W, # 
                                        max.steps=15, 
                                        fam=FAMILY,
                                        mu0=NULL, 
                                        scale.param=NULL,
                                        lambda = lambda, 
                                        lambda.selection.criterion = 'grid', 
                                        DOF.evaluation = 'exact', 
                                        lambda.selection.lossfunction = 'GCV')
    end.2D = Sys.time()
    
    times.2D[j,i] = difftime(end.2D, start.2D, units="secs")
    coef.2D = eval.FEM( FEM(output_CPP.2D$fit.FEM$coeff[ , output_CPP.2D$optimization$lambda_position],
                            FEMbasis.2D),mesh$nodes) 
    
    rmse.1.2D[j,i] = norm( beta_ex[1] - output_CPP$solution$beta[1, output_CPP.2D$optimization$lambda_position ], type="2")
    rmse.2.2D[j,i] = norm( beta_ex[2] - output_CPP$solution$beta[2, output_CPP.2D$optimization$lambda_position ], type="2")
    
    
    mise.2D[j,i] = integrate_f( FEM( (coef.ex - coef.2D)^2, FEMbasis) )
    sols.2D[j,i,] = coef.2D
    
  }
}
end.tot = Sys.time()
tot.time = difftime(end.tot, start.tot, units="mins")
tot.time

today_ = Sys.Date()
ntest_ = "-test-3-with-cov-22-00"
path_ = "/home/aldo/Scrivania/fdaPDE-DATA/GSR"

file.name = paste("GSR-",today_,ntest_,sep="")

save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/GSR/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/GSR/",file.name,".pdf",sep="")

save(nobs, N, M, sols, sols.2D,
     mise, mise.2D, rmse.1, rmse.2, rmse.1.2D, rmse.2.2D,
     times, times.2D, file = save.file)

##########################

save.file = "/home/aldo/Scrivania/fdaPDE-DATA/GSR/GSR-2022-01-27-test-3-with-cov-22-00.RData"
load(save.file)


{
  PP = runiflpp(nobs[1], L)
  locations = cbind(PP$data$x, PP$data$y)
  W=matrix(0,nrow=nobs[1], ncol=2)
  W[,1]= rnorm(nobs[1], mean=2, sd=0.5)
  W[,2]= rnorm(nobs[1], mean=0,sd=1)  
  #W = sin(2*pi * locations[,1]) * cos(2*pi* locations[,2])
  beta_ex = c(0.5,-0.2)
  sol_exact = aux(locations[,1], locations[,2])
  param =  sol_exact + W%*%beta_ex
  mu<-inv.link(param)
  range(mu)
  # sampling response:
  
  response <- rpois(nobs[1], lambda = mu)
  
  palette = heat.colors(n=(diff(range(response)) + 1) )
  p = palette[response + 1] #discreta
  
  plot(mesh)
  points(locations, pch=16, col=p, cex=1.5)
  
}

{
  mean.coef = matrix(0, nrow=102, ncol=1)
  mean.coef.2D = mean.coef
  
  for(j in 1:M){
    mean.coef =  mean.coef + sols[j,3,]
    mean.coef.2D = mean.coef.2D + sols.2D[j,3,]
  }
  
  mean.coef = mean.coef / M
  mean.coef.2D = mean.coef.2D / M
  
  mesh.ref = refine.mesh.1.5D(mesh, 0.01)
  FEMbasis.ref = create.FEM.basis(mesh.ref)
  
  mean.coef.ref = eval.FEM( FEM(mean.coef, FEMbasis), mesh.ref$nodes)
  mean.coef.2D.ref = eval.FEM( FEM(mean.coef.2D, FEMbasis) , mesh.ref$nodes)
  coef.ex.ref = eval.FEM(FEM(coef.ex, FEMbasis), mesh.ref$nodes)
}

pdf("/home/aldo/Scrivania/fdaPDE-IMG/GSR/GSR-2022-01-27-test-3-with-cov-22-00.pdf")

plot(mesh.2D$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n",type="n")
segments(mesh$nodes[mesh$edges[,1],1], mesh$nodes[mesh$edges[,1],2],
         mesh$nodes[mesh$edges[,2],1], mesh$nodes[mesh$edges[,2],2])
segments(mesh.2D$nodes[mesh.2D$edges[,1],1], mesh.2D$nodes[mesh.2D$edges[,1],2],
         mesh.2D$nodes[mesh.2D$edges[,2],1], mesh.2D$nodes[mesh.2D$edges[,2],2], col=rgb(0,0,0,max=255, alpha=95))
segments(mesh.2D$nodes[mesh.2D$segments[,1],1], mesh.2D$nodes[mesh.2D$segments[,1],2],
         mesh.2D$nodes[mesh.2D$segments[,2],1], mesh.2D$nodes[mesh.2D$segments[,2],2], col="red")
points(locations, pch=16, col=p, cex=1.5)
min.col = min(coef.ex.ref, mean.coef.ref)
max.col = max(coef.ex.ref, mean.coef.ref)

min.col = min(min.col, mean.coef.2D.ref)
max.col = max(max.col, mean.coef.2D)

library(ggplot2)
library(latex2exp)

MyTheme <- theme(
  axis.text = element_text(size=26),
  axis.title = element_text(size=26),
  title = element_text(size=26),
  plot.title = element_text(hjust = 0.5),
  legend.text = element_text(size=24),
  legend.key.size = unit(1,"cm"),
  legend.key.height = unit(1,"cm"),
  legend.title = element_blank(),
  legend.background = element_rect(fill="white", color="black",
                                   size=c(1,0.5))
)

# legend.position= c(0.825,0.15) # basso dx
# legend.position= c(0.825,0.85) # alto  dx

mise_ = cbind( as.vector(mise), as.vector(mise.2D))
mise_ = as.vector(mise_)
rmse.1_ = cbind(rmse.1, rmse.1.2D)
rmse.1_ = as.vector(rmse.1_)
rmse.2_ = cbind(rmse.2, rmse.2.2D)
rmse.2_ = as.vector(rmse.2_)
times_ = cbind( as.vector(times), as.vector(times.2D)) 
times_ = as.vector(times_)
nobs_ = rep(nobs,each=M)
nobs_ = rep(as.character(nobs_), times=2)
type_ = rep(c("1.5D","2D"),each=(M*N))

data_frame.1 = data.frame(mise_, rmse.1_, rmse.2_, times_,nobs_,type_)

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=mise_,group=interaction(nobs_,type_),fill=type_))+
  scale_x_discrete(limits=as.character(nobs))+
  labs(x="observations",y="",fill="",title="MISE" )+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.85,0.85))+
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=rmse.1_,group=interaction(nobs_,type_),fill=type_))+
  scale_x_discrete(limits=as.character(nobs))+
  labs(x="observations",y="",fill="",title=TeX("$RMSE^{2}(\\beta_{1})$") )+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.85,0.85))+
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=rmse.2_,group=interaction(nobs_,type_),fill=type_))+
  scale_x_discrete(limits=as.character(nobs))+
  labs(x="observations",y="",fill="",title=TeX("$RMSE^{2}(\\beta_{2})$") )+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.85,0.85))+
  MyTheme


R_plot_graph.ggplot2(FEM(coef.ex.ref, FEMbasis.ref), 
                     line.size = 1.,
                     title = "True spatial field")

R_plot_graph.ggplot2(FEM(mean.coef.ref, FEMbasis.ref),
                     line.size = 1.,
                     title="1.5D")

R_plot_graph.ggplot2(FEM(mean.coef.2D.ref , FEMbasis.ref),
                     line.size = 1.,
                     title="2D")



R_plot_graph.ggplot2.2(FEM(coef.ex.ref, FEMbasis.ref), 
                       line.size = 1.,
                       title = "True spatial field",
                       color.min =  min.col,
                       color.max = max.col)

R_plot_graph.ggplot2.2(FEM(mean.coef.ref, FEMbasis.ref),
                       line.size = 1.,
                       title="1.5D",
                       color.min =  min.col,
                       color.max = max.col)

R_plot_graph.ggplot2.2(FEM(mean.coef.2D.ref, FEMbasis.ref),
                       line.size = 1.,
                       title="2D",
                       color.min =  min.col,
                       color.max = max.col)

R_plot_graph.R(FEM(coef.ex.ref, FEMbasis.ref), 
               title="True spatial field",
               line.size = 3)
R_plot_graph.R(FEM(mean.coef.ref, FEMbasis.ref), 
               title="1.5D",
               line.size = 3)
R_plot_graph.R(FEM(mean.coef.2D.ref, FEMbasis.ref), 
               title="2D",
               line.size = 3)
dev.off()
