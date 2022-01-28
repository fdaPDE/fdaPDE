# Geometry + Spatial Field 
{
  
  library(spatstat)
  library(purrr)
  source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
  source("tests/Auxiliary/R_plot_graph.ggplot2.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")
  
  FAMILY = "binomial"
  link = function(x) { log(x/(1-x))}
  inv.link = function(x){ exp(x)/(1+exp(x))}

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
    sigma.4 = 0.0075 *2
    sigma.5 = 0.0075 *2
    sigma.6 = 0.0075 *2
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
    
    res = res/30

    intersect.4.5 = intersect(idx.4, idx.5)
    intersect.8.10 = intersect(idx.10, idx.8)
    intersect.tot = intersect( intersect.4.5, intersect.8.10)
    intersect.tot = intersect(intersect.tot, idx.6)

    res[intersect.tot] = -1.0
    return(res)
  }
  
  coef.ex = aux(mesh$nodes[,1],mesh$nodes[,2])
  R_plot_graph.ggplot2(FEM(coef.ex,FEMbasis))
} 

{
  param = coef.ex
  mu<-inv.link(param)
  range(mu)
  # sampling response:
  response <- rbernoulli(nrow(mesh$nodes), p = mu)
  
  
  palette = c("red3", "green3")
  colors = c()
  for(i in 1:nrow(mesh$nodes)){
    if( response[i] == TRUE )
      colors = append(colors, palette[2])
    else
      colors = append(colors, palette[1])
  }
  
  plot(mesh)
  points(mesh$nodes, pch=16, col=colors)

  p = palette[round((response-min(response))/diff(range(response))) + 1] #discreta
  
  plot(mesh)
  points(mesh$nodes, pch=16, col=p)
  }

nobs=c(200, 400, 600, 800)
N = length(nobs)
M = 20
times = matrix(0, nrow=M, ncol=N)
mise = matrix(0,nrow=M, ncol=N)
sols = array(0, dim=c(M,N, nnodes))

lambda = 10^seq(-5,-3,length.out = 20) # test-1 con e senza covariate
lambda=c(1e-10, 1e-8)
i=1
j=1
start.tot = Sys.time()
for(i in 1:N){
  for(j in 1:M){
    PP = runiflpp(nobs[i] , L)
    locations = cbind(PP$data$x, PP$data$y)

    sol_exact = aux(locations[,1], locations[,2])
    # covariates
    #set.seed(42)
    
    #W=matrix(0,nrow=nobs[i], ncol=1)
    #W[,1]= rnorm(nobs[i], mean=2, sd=0.5)
    #W = sin(2*pi * locations[,1]) * cos(2*pi* locations[,2])
    
    param =  sol_exact #+ W * beta_ex
    mu<-inv.link(param)
    range(mu)
    # sampling response:
    response <- rbernoulli(nobs[i], p = mu)
    
    start = Sys.time()
    output_CPP <- fdaPDE::smooth.FEM(location = locations, 
                                     observations = as.numeric(response), 
                                     FEMbasis =FEMbasis, 
                                     covariates = NULL, #W
                                     max.steps=4, 
                                     fam=FAMILY,
                                     mu0=NULL, 
                                     scale.param=NULL,
                                     lambda = lambda , 
                                     lambda.selection.criterion = 'grid', 
                                     DOF.evaluation = 'exact', 
                                     lambda.selection.lossfunction = 'GCV')
    end = Sys.time()
    times[j,i] = difftime(end, start, units="secs")
    
    coef = output_CPP$fit.FEM$coeff[ , output_CPP$optimization$lambda_position]
    mise[j,i] = integrate_f( FEM( (coef.ex - coef)^2, FEMbasis) )
    
    #rmse[j,i] = norm( beta_ex - output_CPP$solution$beta[ output_CPP$optimization$lambda_position ], type="2")
    norm.l2[j,i] = norm(coef.ex- coef, type="2")
    sols[j,i,] = coef
    
    save(nobs, N, M, sols, lambda, coef.ex,i,j,
         mise, 
         times,  file = save.file)
    
    # start.2D = Sys.time()
    # output_CPP.2D <- fdaPDE::smooth.FEM(location = locations, 
    #                                  observations = as.numeric(response), 
    #                                  FEMbasis =FEMbasis.2D, 
    #                                  covariates = W, # 
    #                                  max.steps=15, 
    #                                  fam=FAMILY,
    #                                  mu0=NULL, 
    #                                  scale.param=NULL,
    #                                  lambda = lambda, 
    #                                  lambda.selection.criterion = 'grid', 
    #                                  DOF.evaluation = 'exact', 
    #                                  lambda.selection.lossfunction = 'GCV')
    # end.2D = Sys.time()
    # 
    # times.2D[j,i] = difftime(end.2D, start.2D, units="secs")
    # coef.2D = eval.FEM( FEM(output_CPP.2D$fit.FEM$coeff[ , output_CPP.2D$optimization$lambda_position],
    #                         FEMbasis.2D),mesh$nodes) 
    # 
    # rmse.2D[j,i] = norm( beta_ex - output_CPP$solution$beta[ output_CPP.2D$optimization$lambda_position ], type="2")
    # 
    # mise.2D[j,i] = integrate_f( FEM( (coef.ex - coef.2D)^2, FEMbasis) )
    # sols.2D[j,i,] = coef.2D
    # print(paste("lambda.2D pos = ", output_CPP.2D$optimization$lambda_position))
  }
}
end.tot = Sys.time()
tot.time = difftime(end.tot, start.tot, units="mins")
tot.time

today_ = Sys.Date()
ntest_ = "-test-1-NO-COV-17-00"
path_ = "/home/aldo/Scrivania/fdaPDE-DATA/GSR"

file.name = paste("GSR-",today_,ntest_,sep="")

save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/GSR/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/GSR/",file.name,".pdf",sep="")

save(nobs, N, M, sols, lambda, coef.ex,
     mise, 
     times,  file = save.file)

###########################################################
library(ggplot2)
library(latex2exp)

load("/home/aldo/Scrivania/fdaPDE-DATA/GSR/GSR-2022-01-25-test-1-no-cov-spatial-field.RData")
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

mise_ = as.vector(mise)
times_ = as.vector(times)
nobs_ = rep(as.character(nobs),each=M)

data_frame.1 = data.frame(mise_,times_,nobs_)

{
  PP = runiflpp(nobs[2], L)
  locations = cbind(PP$data$x, PP$data$y)
  param = aux(locations[,1], locations[,2])
  mu<-inv.link(param)
  range(mu)
  # sampling response:
  response <- rbernoulli(nobs[2], p = mu)
  
  palette = c("red3", "green4")
  colors = c()
  for(i in 1:nrow(locations)){
    if( response[i] == TRUE )
      colors = append(colors, palette[2])
    else
      colors = append(colors, palette[1])
  }
  
  plot(mesh)
  points(locations, pch=16, col=colors)
  
  #p = palette[round((response-min(response))/diff(range(response))) + 1] #discreta

}

{
mean.coef = matrix(0, nrow=nnodes, ncol=1)

for(j in 1:M){
  mean.coef =  mean.coef + sols[j,4,]
}

mean.coef = mean.coef / M

mesh.ref = refine.mesh.1.5D(mesh, 0.01)
FEMbasis.ref = create.FEM.basis(mesh.ref)

mean.coef.ref = eval.FEM( FEM(mean.coef, FEMbasis), mesh.ref$nodes)
coef.ex.ref = eval.FEM(FEM(coef.ex, FEMbasis), mesh.ref$nodes)
}

pdf(img.file)

plot(mesh)
points(locations, pch=16, col=colors, cex=2)

min.col = min(coef.ex.ref, mean.coef.ref)
max.col = max(coef.ex.ref, mean.coef.ref)

R_plot_graph.ggplot2(FEM(coef.ex.ref, FEMbasis.ref), 
                     line.size = 1.,
                     title = "True spatial field")

R_plot_graph.ggplot2(FEM(mean.coef.ref, FEMbasis.ref),
                     line.size = 1.,
                     title="fdaPDE")


R_plot_graph.ggplot2.2(FEM(coef.ex.ref, FEMbasis.ref), 
                     line.size = 1.,
                     title = "True spatial field",
                     color.min =  min.col,
                     color.max = max.col)

R_plot_graph.ggplot2.2(FEM(mean.coef.ref, FEMbasis.ref),
                     line.size = 1.,
                     title="fdaPDE",
                     color.min =  min.col,
                     color.max = max.col)

R_plot_graph.R(FEM(coef.ex.ref, FEMbasis.ref), 
               title="True spatial field",
               line.size = 3)
R_plot_graph.R(FEM(mean.coef.ref, FEMbasis.ref), 
               title="fdaPDE",
               line.size = 3)

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=mise_, fill="red"))+
  scale_x_discrete(limits=as.character(nobs))+
  labs(x="observations",y="",title="MISE" )+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")+
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=mise_/20, fill="red"))+
  scale_x_discrete(limits=as.character(nobs))+
  labs(x="observations",y="",title="MISE" )+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")+
  MyTheme


norm.l2 = matrix(0, nrow=M, ncol=N)

for(i in 1:N){
  for(j in 1:M)
    norm.l2[j,i] = norm(coef.ex - sols[j,i,], type="2")
}

################################################################################
# Test 1 True field sum of five "Mollifier"

{
  library(spatstat)
  library(purrr)
  source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
  source("tests/Auxiliary/R_plot_graph.ggplot2.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")
  
  FAMILY = "binomial"
  link = function(x) { log(x/(1-x))}
  inv.link = function(x){ exp(x)/(1+exp(x))}
  
  vertices = cbind(simplenet$vertices$x, simplenet$vertices$y)
  edges = cbind(simplenet$from , simplenet$to)
  L = as.linnet(simplenet)
  
  delta = 0.03 # 102 nodi, 0.015 199 nodi
  mesh =create.mesh.1.5D(vertices, edges)
  mesh =refine.mesh.1.5D(mesh,delta=0.015)
  
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
    #res.5 = equal_split_discontinous.mollifier(mesh, 0.025, dijkstra.5, x, y) 
   # res.6 = equal_split_discontinous.mollifier(mesh, sigma, dijkstra.6, x, y)
    #res.8 = equal_split_discontinous.mollifier(mesh, 0.025, dijkstra.8, x, y)
    res.10 = equal_split_discontinous.mollifier(mesh, sigma, dijkstra.10, x, y)
    
    #res= (res.4$coef + res.5$coef  + res.8$coef + res.10 $coef) # + res.6$coef 
    res = res.4$coef + res.10$coef
    
    idx.4 = which(res.4$bandwidth==0.0)
    #idx.5 = which(res.5$bandwidth==0.0)
    #idx.6 = which(res.6$bandwidth==0.0)
    #idx.8 = which(res.8$bandwidth==0.0)
    idx.10 = which(res.10$bandwidth==0.0)
    
    #intersect.4.5 = intersect(idx.4, idx.5)
    #intersect.8.10 = intersect(idx.10, idx.8)
    #intersect.tot = intersect( intersect.4.5, intersect.8.10)
    #intersect.tot = intersect(intersect.tot, idx.6)
    intersect.tot = intersect(idx.4, idx.10)
    
    res[intersect.tot] = -1.0
    return(res)
  }
  
  #coef.ex = aux(mesh$nodes[,1],mesh$nodes[,2])
  R_plot_graph.ggplot2.2(FEM(coef.ex,FEMbasis), 
                         color.min  =min(coef.ex),
                         color.max  =max(coef.ex))
} 

{
  param = coef.ex
  mu<-inv.link(param)
  range(mu)
  # sampling response:
  response <- rbernoulli(nrow(mesh$nodes), p = mu)
  
  
  palette = c("red3", "green3")
  colors = c()
  for(i in 1:nrow(mesh$nodes)){
    if( response[i] == TRUE )
      colors = append(colors, palette[2])
    else
      colors = append(colors, palette[1])
  }
  
  plot(mesh)
  points(mesh$nodes, pch=16, col=colors)
  
  p = palette[round((response-min(response))/diff(range(response))) + 1] #discreta
  
  plot(mesh)
  points(mesh$nodes, pch=16, col=p)
}

nobs=c(200, 400, 600, 800)
N = length(nobs)
M = 20
times = matrix(0, nrow=M, ncol=N)
mise = matrix(0,nrow=M, ncol=N)
sols = array(0, dim=c(M,N, nnodes))
norm.l2 = matrix(0, nrow=M, ncol=N)
lambda = 10^seq(-7,-6,length.out = 20) # test-1 con e senza covariate
lambda.opt = matrix(0, nrow=M, ncol=N)
start.tot = Sys.time()
for(i in 1:N){
  for(j in 1:M){
    PP = runiflpp(nobs[i] , L)
    locations = cbind(PP$data$x, PP$data$y)
    
    sol_exact = aux(locations[,1], locations[,2])
    
    param =  sol_exact #+ W * beta_ex
    mu<-inv.link(param)
    range(mu)
    # sampling response:
    response <- rbernoulli(nobs[i], p = mu)
    
    start = Sys.time()
    output_CPP <- fdaPDE::smooth.FEM(location = locations, 
                                     observations = as.numeric(response), 
                                     FEMbasis =FEMbasis, 
                                     covariates = NULL, #W
                                     max.steps=500, 
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
    lambda.opt[j,i] = output_CPP$optimization$lambda_position
    
    norm.l2[j,i] = norm(coef.ex- coef, type="2")
    sols[j,i,] = coef
  }
}
end.tot = Sys.time()
tot.time = difftime(end.tot, start.tot, units="mins")
tot.time

today_ = Sys.Date()
ntest_ = "-test-1-NO-COV-23-00"
path_ = "/home/aldo/Scrivania/fdaPDE-DATA/GSR"

file.name = paste("GSR-",today_,ntest_,sep="")

save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/GSR/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/GSR/",file.name,".pdf",sep="")

save(nobs, N, M, sols, lambda, coef.ex,
     mise, 
     times,  file = save.file)
###########################################################

today_ = "2022-01-25"
ntest_ = "-test-1-NO-COV-23-00"
path_ = "/home/aldo/Scrivania/fdaPDE-DATA/GSR"

file.name = paste("GSR-",today_,ntest_,sep="")

save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/GSR/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/GSR/",file.name,".pdf",sep="")
load(save.file)

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

mise_ = as.vector(mise[,2:N])
times_ = as.vector(times[,2:N])
nobs_ = rep(as.character(nobs[2:N]),each=M)

data_frame.1 = data.frame(mise_,times_,nobs_)

{
  PP = runiflpp(nobs[2], L)
  locations = cbind(PP$data$x, PP$data$y)
  param = aux(locations[,1], locations[,2])
  mu<-inv.link(param)
  range(mu)
  # sampling response:
  response <- rbernoulli(nobs[2], p = mu)
  
  palette = c("red3", "green4")
  colors = c()
  for(i in 1:nrow(locations)){
    if( response[i] == TRUE )
      colors = append(colors, palette[2])
    else
      colors = append(colors, palette[1])
  }
  
  plot(mesh)
  points(locations, pch=16, col=colors)
  
  #p = palette[round((response-min(response))/diff(range(response))) + 1] #discreta
  
}

{
  mean.coef = matrix(0, nrow=199, ncol=1)
  
  for(j in 1:M){
    mean.coef =  mean.coef + sols[j,4,]
  }
  
  mean.coef = mean.coef / M
  
  mesh.ref = refine.mesh.1.5D(mesh, 0.01)
  FEMbasis.ref = create.FEM.basis(mesh.ref)
  
  mean.coef.ref = eval.FEM( FEM(mean.coef, FEMbasis), mesh.ref$nodes)
  coef.ex.ref = eval.FEM(FEM(coef.ex, FEMbasis), mesh.ref$nodes)
}

pdf(img.file)
plot(mesh)
points(locations, pch=16, col=colors, cex=1.5)
min.col = min(coef.ex.ref, mean.coef.ref)
max.col = max(coef.ex.ref, mean.coef.ref)

R_plot_graph.ggplot2(FEM(coef.ex.ref, FEMbasis.ref), 
                     line.size = 1.,
                     title = "True spatial field")

R_plot_graph.ggplot2(FEM(mean.coef.ref, FEMbasis.ref),
                     line.size = 1.,
                     title="fdaPDE")


R_plot_graph.ggplot2.2(FEM(coef.ex.ref, FEMbasis.ref), 
                       line.size = 1.,
                       title = "True spatial field",
                       color.min =  min.col,
                       color.max = max.col)

R_plot_graph.ggplot2.2(FEM(mean.coef.ref, FEMbasis.ref),
                       line.size = 1.,
                       title="fdaPDE",
                       color.min =  min.col,
                       color.max = max.col)

R_plot_graph.R(FEM(coef.ex.ref, FEMbasis.ref), 
               title="True spatial field",
               line.size = 3)
R_plot_graph.R(FEM(mean.coef.ref, FEMbasis.ref), 
               title="fdaPDE",
               line.size = 3)

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=mise_, fill="red"))+
  scale_x_discrete(limits=as.character(nobs[2:N]))+
  labs(x="observations",y="",title="MISE" )+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")+
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=mise_/20, fill="red"))+
  scale_x_discrete(limits=as.character(nobs[2:N]))+
  labs(x="observations",y="",title="MISE" )+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")+
  MyTheme

dev.off()
