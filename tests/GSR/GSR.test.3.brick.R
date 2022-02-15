
  
  library(spatstat)
  library(purrr)
  source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
  source("tests/Auxiliary/R_plot_graph.ggplot2.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")
  # scala
  {
  eps = 0.05
  x = c(0., 1, 0.05, 0.10, 0.45, 0.5, 0.55, 0.90, 0.95)
  y = c(0.,eps)
  vertices = expand.grid(x,y)
  vertices = cbind(vertices[,1], vertices[,2])

  edges = matrix(c(1,3,
                   3,4,
                   4,5,
                   5,6,
                   6,7,
                   7,8,
                   8,9,
                   9,2,
                   10,12,
                   12,13,
                   13,14,
                   14,15,
                   15,16,
                   16,17,
                   17,18,
                   18,11,
                   3,12,
                   4,13,
                   5,14,
                   6,15,
                   7,16,
                   8,17,
                   9,18), nrow=23, ncol=2, byrow=T)
  m = matrix(FALSE, nrow=nrow(vertices),ncol=nrow(vertices))
  for(e in 1:nrow(edges) ){
    m[edges[e,1], edges[e,2]] = TRUE
    m[edges[e,2], edges[e,1]] = TRUE
    
  }
  vertices.L = ppp(vertices[,1], vertices[,2])
  L = linnet(vertices.L, m) 
  PP = runiflpp(100, L)
  }
  
  # "C-shaped #
  {
    eps = 1 / 2
    x = c(0., 1)
    y = c(0.,eps)
    vertices = expand.grid(x,y)
    vertices = cbind(vertices[,1], vertices[,2])
    
    edges = matrix(c(1,2,1,3,3,4), nrow=3,ncol=2, byrow=T)
    m = matrix(FALSE, nrow=nrow(vertices),ncol=nrow(vertices))
    for(e in 1:nrow(edges) ){
      m[edges[e,1], edges[e,2]] = TRUE
      m[edges[e,2], edges[e,1]] = TRUE
      
    }
    vertices.L = ppp(vertices[,1], vertices[,2])
    L = linnet(vertices.L, m) 
    PP = runiflpp(100, L)
    
  }
  
  mesh = create.mesh.1.5D(vertices, edges)
  mesh = refine.mesh.1.5D(mesh,delta=0.0125)
  
  plot(mesh,show.nodes = T)
  points(PP$data$x, PP$data$y, pch=16)
  nnodes = nrow(mesh$nodes)
  
  FAMILY = "poisson"
  l<-make.link("log")
  link<-l$linkfun
  inv.link<-l$linkinv
 
  FEMbasis = create.FEM.basis(mesh)
  
  ### CAMPO f ### 
   # Scala
   { 
  aux.11  = function(x,y){
    h = 0.5
    source = 11 
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
      else{
        coef[i] = -1.
      }
    }

    return(coef)
  }
# 
#   aux.11 = function(x,y){
# 
#     h = 0.05
#     source = 11
#     points_ = cbind(x,y)
#     idx.ok = which( abs(points_[,1]-mesh$nodes[source,1]) < 10 * .Machine$double.eps )
#     coef = vector(mode="numeric", length=length(x))
# 
#     for(i in idx.ok){
#       delta = abs(points_[i,2] - mesh$nodes[source,2])
#       if(delta < h ){
#         coef[i] = -2.5/h^2 * (delta^2 - h^2)  - 1
#       }
#     }
# 
#     return(coef)
# 
#   }
#   
  aux.1  = function(x,y){
    h = 0.5
    source = 1
    points_ = cbind(x,y)
    idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
    coef = vector(mode="numeric", length=length(x))
    
    for(i in idx.ok){
      delta = abs(points_[i,1] - mesh$nodes[source,1])
      if(delta < h ){
        #coef[i] = 1- 8*delta
        coef[i] = -4/h^2 * (delta^2 - h^2) - 2
        coef[i] = 3*exp(1/(delta^2 / h^2 -1)+1 ) - 1
      }
      else{
        coef[i] = -1.
      }
    }
  
    return(coef)
  }
  
  AUX = function(x,y){
    res = aux.1(x,y) + aux.11(x,y) #+ aux.11(x,y)
    return (res)
  }
   
   }
   
  # C-shaped
   {
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
   
     AUX = function(x,y){
       
       res = aux.1(x,y) + aux.3(x,y) + aux.4(x,y)
       return(res)
     }
     coef.1 = aux.1(mesh$nodes[,1], mesh$nodes[,2])
     coef.2 = aux.3(mesh$nodes[,1], mesh$nodes[,2])
     coef.4 = aux.4(mesh$nodes[,1], mesh$nodes[,2])
     R_plot_graph.a.sym.ggplot2(FEM(coef.4,FEMbasis))
     
     coef.ex = AUX(mesh$nodes[,1], mesh$nodes[,2])
     R_plot_graph.a.sym.ggplot2(FEM(coef.ex,FEMbasis))
    
     }
  coef.ex = AUX(mesh$nodes[,1],mesh$nodes[,2])
  R_plot_graph.a.sym.ggplot2(FEM(coef.ex,FEMbasis))
  R_plot_graph.ggplot2.2(FEM(coef.ex,FEMbasis),ratio=2.5)
 

# 2D Geometry 
#x = seq(from=min(vertices[,1])-0.02, to=max(vertices[,1])+0.02, length.out=50)
#y = seq(from=min(vertices[,2])-0.02, to=max(vertices[,2])+0.02, length.out=10)
x=c(min(vertices[,1])-0.02, max(vertices[,1])+0.02)
y=c(min(vertices[,2])-0.02, max(vertices[,2])+0.02)  
nodes.2D = expand.grid(x,y)
mesh.2D = create.mesh.2D(nodes.2D)
mesh.2D = refine.mesh.2D(mesh.2D, maximum_area = 0.001*eps/0.5, minimum_angle=20)
nnodes.2D = nrow(mesh.2D$nodes)
FEMbasis.2D = create.FEM.basis(mesh.2D)
plot(mesh.2D,asp=1)
points(mesh$nodes, pch=16, col="red3")
nrow(mesh.2D$nodes)

nobs=c(300, 500, 700, 1000)
N = length(nobs)
M = 20
nnodes = nrow(mesh$nodes)
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
lambda = 10^seq(-3.5,-2,length.out=10)
lambda = 10^seq(-2,0,length.out = 10) 
start.tot = Sys.time()
for(i in 1:N){
  for(j in 1:M){
    PP = runiflpp(nobs[i] , L)
    locations = cbind(PP$data$x, PP$data$y)
    
    sol_exact = AUX(locations[,1], locations[,2])
    # covariates
    #set.seed(42)
    
    W=matrix(0,nrow=nobs[i], ncol=2)
    W[,1]= rnorm(nobs[i], mean=2, sd=0.5) # test-3
    W[,2]= rnorm(nobs[i], mean=0,sd=1)    # test-3
    # W[,1]=rbeta(nobs[i],shape1=1.5,shape2=2)
    # W[,2]= rbeta(nobs[i], shape1=3, shape2=2)
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
    
    norm.l2[j,i] = norm(coef-coef.ex, type="2")
    
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
    
    norm.l2.2D[j,i] = norm(coef.2D-coef.ex, type = "2")
    
    mise.2D[j,i] = integrate_f( FEM( (coef.ex - coef.2D)^2, FEMbasis) )
    sols.2D[j,i,] = coef.2D
    print(paste("lambda.2D pos = ", output_CPP.2D$optimization$lambda_position))
  }
}
end.tot = Sys.time()
tot.time = difftime(end.tot, start.tot, units="mins")
tot.time

today_ = Sys.Date()
ntest_ = "-test-2-C-shape-DELTA-0-5"
path_ = "/home/aldo/Scrivania/fdaPDE-DATA/GSR"

file.name = paste("GSR-",today_,ntest_,sep="")

save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/GSR/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/GSR/",file.name,".pdf",sep="")

save(nobs, N, M, sols, sols.2D,
     mise, mise.2D, rmse.1, rmse.2, rmse.1.2D, rmse.2.2D,
     norm.l2, norm.l2.2D,
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

mise_ = cbind( as.vector(mise[,2:4]), as.vector(mise.2D[,2:4]))
mise_ = as.vector(mise_)
rmse.1_ = cbind(rmse.1[,2:4], rmse.1.2D[,2:4])
rmse.1_ = as.vector(rmse.1_)
rmse.2_ = cbind(rmse.2[,2:4], rmse.2.2D[,2:4])
rmse.2_ = as.vector(rmse.2_)
times_ = cbind( as.vector(times[,2:4]), as.vector(times.2D[,2:4])) 
times_ = as.vector(times_)
nobs_ = rep(nobs[2:4],each=M)
nobs_ = rep(as.character(nobs_), times=2)
type_ = rep(c("1.5D","2D"),each=(M*(N-1)) )

data_frame.1 = data.frame(mise_, rmse.1_, rmse.2_, times_,nobs_,type_)

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=mise_,group=interaction(nobs_,type_),fill=type_))+
  scale_x_discrete(limits=as.character(nobs[2:4]))+
  labs(x="observations",y="",fill="",title="MISE" )+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.875,0.85))+
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=rmse.1_,group=interaction(nobs_,type_),fill=type_))+
  scale_x_discrete(limits=as.character(nobs))+
  labs(x="observations",y="",fill="",title=TeX("$RMSE^{2}$") )+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.875,0.85))+
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=rmse.2_,group=interaction(nobs_,type_),fill=type_))+
  scale_x_discrete(limits=as.character(nobs))+
  labs(x="observations",y="",fill="",title=TeX("$RMSE^{2}$") )+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.875,0.85))+
  MyTheme

{
  N_ = 500
  PP = runiflpp(N_, L)
  locations = cbind(PP$data$x, PP$data$y)
  W=matrix(0,nrow=N_, ncol=2)
  # W[,1]=rbeta(150,shape1=1.5,shape2=2)
  # W[,2]= rbeta(150, shape1=3, shape2=2)
  W[,1]= rnorm(N_, mean=2, sd=0.5) # test-3
  W[,2]= rnorm(N_, mean=0,sd=1)    # test-3
  beta_ex = c(0.5,-0.2)
  sol_exact = AUX(locations[,1], locations[,2])
  param =  sol_exact + W%*%beta_ex
  mu<-inv.link(param)
  range(mu)
  # sampling response:
  
  response <- rpois(N_, lambda = mu)
  
  palette = heat.colors(n=(diff(range(response)) + 1) )
  p = palette[response + 1] #discreta

  plot(mesh)
  points(locations, pch=16, col=p, cex=1.5)
  
  palette = jet.col(n=(diff(range(response)) + 1))
  p = palette[response + 1] #discreta

  TOT = 7
  palette = jet.col(n= floor((diff(range(response)))/TOT) )
  p = palette[floor(response/TOT) + 1] #discreta
  
  plot(mesh)
  points(locations, pch=16, col=p, cex=1.5)
}

{
  mean.coef = matrix(0, nrow=nnodes, ncol=1)
  mean.coef.2D = mean.coef
  i_ = 4
  
  for(j in 1:M){
    mean.coef =  mean.coef + sols[j,i_,]
    mean.coef.2D = mean.coef.2D + sols.2D[j,i_,]
  }
  
  mean.coef = mean.coef / M
  mean.coef.2D = mean.coef.2D / M
  
  mesh.ref = refine.mesh.1.5D(mesh, 0.01)
  FEMbasis.ref = create.FEM.basis(mesh.ref)
  
  mean.coef.ref = eval.FEM( FEM(mean.coef, FEMbasis), mesh.ref$nodes)
  mean.coef.2D.ref = eval.FEM( FEM(mean.coef.2D, FEMbasis) , mesh.ref$nodes)
  coef.ex.ref = eval.FEM(FEM(coef.ex, FEMbasis), mesh.ref$nodes)
}

pdf("/home/aldo/Scrivania/fdaPDE-IMG/GSR/GSR-2022-02-09-test-2-C-shape-DELTA-0-25.pdf")

plot(mesh.2D$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n",type="n")
segments(mesh$nodes[mesh$edges[,1],1], mesh$nodes[mesh$edges[,1],2],
         mesh$nodes[mesh$edges[,2],1], mesh$nodes[mesh$edges[,2],2], cex=5)
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

mise_ = cbind( as.vector(mise[,2:4]), as.vector(mise.2D[,2:4]))
mise_ = as.vector(mise_)
rmse.1_ = cbind(rmse.1[,2:4], rmse.1.2D[,2:4])
rmse.1_ = as.vector(rmse.1_)
rmse.2_ = cbind(rmse.2[,2:4], rmse.2.2D[,2:4])
rmse.2_ = as.vector(rmse.2_)
times_ = cbind( as.vector(times[,2:4]), as.vector(times.2D[,2:4])) 
times_ = as.vector(times_)
nobs_ = rep(nobs[2:4],each=M)
nobs_ = rep(as.character(nobs_), times=2)
type_ = rep(c("1.5D","2D"),each=(M*(N-1)) )

data_frame.1 = data.frame(mise_, rmse.1_, rmse.2_, times_,nobs_,type_)

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=mise_,group=interaction(nobs_,type_),fill=type_))+
  scale_x_discrete(limits=as.character(nobs[2:4]))+
  labs(x="observations",y="",fill="",title="MISE" )+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.875,0.85))+
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=rmse.1_,group=interaction(nobs_,type_),fill=type_))+
  scale_x_discrete(limits=as.character(nobs[2:4]))+
  labs(x="observations",y="",fill="",title=TeX("$RMSE^{2}(\\beta_{1})$") )+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.875,0.85))+
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=rmse.2_,group=interaction(nobs_,type_),fill=type_))+
  scale_x_discrete(limits=as.character(nobs[2:4]))+
  labs(x="observations",y="",fill="",title=TeX("$RMSE^{2}(\\beta_{2})$") )+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.875,0.85))+
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
                       color.max = max.col, 
                       ratio = 1 )

R_plot_graph.ggplot2.2(FEM(mean.coef.ref, FEMbasis.ref),
                       line.size = 1.,
                       title="1.5D",
                       color.min =  min.col,
                       color.max = max.col,
                       ratio = 1)

R_plot_graph.ggplot2.2(FEM(mean.coef.2D.ref, FEMbasis.ref),
                       line.size = 1.,
                       title="2D",
                       color.min =  min.col,
                       color.max = max.col,
                       ratio = 1)

R_plot_graph.R(FEM(coef.ex.ref, FEMbasis.ref), 
               title="True spatial field",
               line.size = 3)
R_plot_graph.R(FEM(mean.coef.ref, FEMbasis.ref), 
               title="1.5D",
               line.size = 3)
R_plot_graph.R(FEM(mean.coef.2D.ref, FEMbasis.ref), 
               title="2D",
               line.size = 3)

x=vector(mode="double")
y=vector(mode="double")
grp.nodes=vector(mode="integer")

num_edges= dim(mesh$edges)[1]
for(e in 1:num_edges){
  x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
  y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
  grp.nodes = append(grp.nodes, rep(e,times=2))
}

data=data.frame(x,y,grp.nodes)

p = jet.col(n=128, alpha=0.8)  
TOT = 7
p = jet.col(n= 12 )


data=data.frame(x,y,grp.nodes)

x.2D=vector(mode="double")
y.2D=vector(mode="double")
grp.nodes.2D=vector(mode="integer")

num_edges= dim(mesh.2D$edges)[1]
for(e in 1:num_edges){
  x.2D = append(x.2D, c(mesh.2D$nodes[ mesh.2D$edges[e,1], 1], mesh.2D$nodes[ mesh.2D$edges[e,2], 1]))
  y.2D = append(y.2D, c(mesh.2D$nodes[ mesh.2D$edges[e,1], 2], mesh.2D$nodes[ mesh.2D$edges[e,2], 2]))
  grp.nodes.2D = append(grp.nodes.2D, rep(e,times=2))
}
data.2D = data.frame(x.2D, y.2D, grp.nodes.2D)

x.2D=vector(mode="double")
y.2D=vector(mode="double")
grp.nodes.2D=vector(mode="integer")
num_segments= dim(mesh.2D$segments)[1]
for(e in 1:num_segments){
  x.2D = append(x.2D, c(mesh.2D$nodes[ mesh.2D$segments[e,1], 1], mesh.2D$nodes[ mesh.2D$segments[e,2], 1]))
  y.2D = append(y.2D, c(mesh.2D$nodes[ mesh.2D$segments[e,1], 2], mesh.2D$nodes[ mesh.2D$segments[e,2], 2]))
  grp.nodes.2D = append(grp.nodes.2D, rep(e,times=2))
}
data.2D.segments = data.frame(x.2D, y.2D, grp.nodes.2D)

num_points = nrow(locations)
coef.points = mu

x.points = locations[,1]
y.points = locations[,2]
data.points = data.frame(x.points, y.points, coef.points)

ggplot() + 
  geom_point(data=data, aes(x=x,y=y,group=grp.nodes), 
             alpha=0.0) + 
  geom_line(data=data, aes(x=x,y=y,group=grp.nodes), 
            size=1, alpha = 0.8) +
  geom_line(data=data.2D, aes(x=x.2D, y=y.2D, group=grp.nodes.2D), 
            size=0.5, alpha=0.1) +
  geom_line(data=data.2D.segments, aes(x=x.2D, y=y.2D, group=grp.nodes.2D), 
            size=1, color="red", alpha=1) +
  geom_point(data=data.points,aes(x=x.points,y=y.points, color=coef.points),
             size=3) +
  labs(x="",y="",color="", title="") + 
  scale_color_gradientn(colours=p)+ 
  coord_fixed(ratio=1) + 
  theme_void() +
  theme(plot.title = element_text(hjust=0.5),
        title = element_text(size=26),
        legend.title = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size=20),
        legend.key.size = unit(1,"cm"),
        legend.key.width = unit(0.5,"cm"))

dev.off()


pdf("/home/aldo/Scrivania/Caricare/GSR-test-2-domain.pdf")
plot(mesh, show.nodes = F, lwd=3, asp=5)
dev.off()
