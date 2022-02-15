################ ST-Regression #####################
#                                                  #
####################################################

library(spatstat)
library(purrr)
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")

# C shaped linear network
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
nnodes = nrow(mesh$nodes)
FEMbasis = create.FEM.basis(mesh)

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
  
  AUX = function(x,y,t){
    
    res = (aux.1(x,y) + aux.3(x,y) + aux.4(x,y)) * cos(t)
    return(res)
  }
}

NumTimePoints=6
TimePoints=seq(0,pi,length.out = NumTimePoints)

SpacePoints=mesh$nodes
SpaceTimePoints=cbind(rep(mesh$nodes[,1],NumTimePoints),rep(mesh$nodes[,2],NumTimePoints),rep(TimePoints,each=nrow(mesh$nodes)))

coef.ex = matrix(0, nrow= nrow(mesh$nodes), ncol=NumTimePoints)
for(i in 1:NumTimePoints){
  coef.ex[,i] = AUX(mesh$nodes[,1],
                    mesh$nodes[,2],
                    TimePoints[i])
  
}
coef.0 = AUX(mesh$nodes[,1], mesh$nodes[,2],t=0)
R_plot_graph.ggplot2(FEM(coef.0, FEMbasis))
min.col = min(coef.ex)
max.col = max(coef.ex)

# True Space-Time field
for( i in 1:NumTimePoints){
  
  print( R_plot_graph.ggplot2.2(FEM(coef.ex[,i], FEMbasis),
                                line.size = 1,
                                title = paste("t = ",TimePoints[i],sep=""),
                                color.min = min.col,
                                color.max = max.col))
}

# smooth.FEM.time<-function(locations = NULL, time_locations = NULL, observations, FEMbasis, time_mesh=NULL,
#                           covariates = NULL, PDE_parameters = NULL,  BC = NULL,
#                           incidence_matrix = NULL, areal.data.avg = TRUE,
#                           FLAG_MASS = FALSE, FLAG_PARABOLIC = FALSE,
#                           FLAG_ITERATIVE = FALSE, threshold = 10^(-4), 
#                           max.steps = 50, IC = NULL,
#                           search = "tree", bary.locations = NULL,
#                           lambda.selection.criterion = "grid", DOF.evaluation = NULL, 
#                           lambda.selection.lossfunction = NULL,
#                           lambdaS = NULL, lambdaT = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, 
#                           DOF.matrix = NULL, GCV.inflation.factor = 1, lambda.optimization.tolerance = 0.05)

# Exact solution (pointwise at nodes)
sol_exact=AUX(SpaceTimePoints[,1],
              SpaceTimePoints[,2],
              SpaceTimePoints[,3])

# Add error to simulate data
M = 20 
sols = list()
sols.PARABOLIC = list()
rmse = matrix(0, nrow=M, ncol=1)
rmse.PARABOLIC = matrix(0, nrow=M, ncol=1)
mise = matrix(0, nrow=M, ncol=1)
mise.PARABOLIC = matrix(0, nrow=M, ncol=1)
mean.coef = matrix(0, nrow= nnodes * NumTimePoints, ncol=1)
mean.coef.PARABOLIC = matrix(0, nrow= nnodes * NumTimePoints, ncol=1)
times = matrix(0, nrow=M, ncol=1)
times.PARABOLIC =   matrix(0, nrow=M, ncol=1)

start.tot = Sys.time()  
for( j in 1:M){
ran = range(sol_exact)
data = sol_exact + rnorm(length(sol_exact), mean=0, sd=0.05*abs(ran[2]-ran[1]))
observations=matrix(data,nrow(SpacePoints),NumTimePoints)

lambdaS = 10^(-5:-3)
lambdaT = 10^(-5:-3)

start = Sys.time()
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                            lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
end = Sys.time()

start.PARABOLIC = Sys.time()
output_CPP.PARABOLIC<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT, FLAG_PARABOLIC = TRUE,
                            lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
end.PARABOLIC = Sys.time()

coef.estimated.PARABOLIC = eval.FEM.time( output_CPP.PARABOLIC$fit.FEM.time, locations=mesh$nodes, time.instants =  TimePoints,
                                          lambdaS = output_CPP$bestlambda[1],
                                          lambdaT = output_CPP$bestlambda[2])
  
coef.estimated = eval.FEM.time( output_CPP$fit.FEM.time, locations=mesh$nodes, time.instants =  TimePoints,
                                lambdaS = output_CPP$bestlambda[1],
                                lambdaT = output_CPP$bestlambda[2])

mean.coef = mean.coef + coef.estimated / M
mean.coef.PARABOLIC = mean.coef.PARABOLIC + coef.estimated.PARABOLIC / M

times[j] = difftime(end, start, units="mins")
times.PARABOLIC[j] = difftime(end.PARABOLIC, start.PARABOLIC, units="mins")

sols[[i]] = output_CPP
sols.PARABOLIC[[i]] = output_CPP.PARABOLIC

diff = coef.estimated - as.vector(coef.ex)
diff.PARABOLIC = coef.estimated.PARABOLIC - as.vector(coef.ex)

tmp = matrix(0, nrow=NumTimePoints, ncol=1)
tmp.PARABOLIC = matrix(0,nrow=NumTimePoints, ncol=1)

tmp.mise = matrix(0, nrow=NumTimePoints, ncol=1)
tmp.mise.PARABOLIC = matrix(0, nrow=NumTimePoints, ncol=1)

for( i in 1:NumTimePoints){
  
  tmp[i] = norm( diff[(1 + (i-1)*nnodes) : ( nnodes + (i-1)*nnodes)], type="2") 
  tmp.PARABOLIC[i] = norm( diff.PARABOLIC[(1 + (i-1)*nnodes) : ( nnodes + (i-1)*nnodes)], type="2")
  
  tmp.mise[i] = integrate_f(FEM(diff[(1 + (i-1)*nnodes) : ( nnodes + (i-1)*nnodes)]^2,FEMbasis))
  tmp.mise.PARABOLIC[i] = integrate_f(FEM(diff.PARABOLIC[(1 + (i-1)*nnodes) : ( nnodes + (i-1)*nnodes)]^2, FEMbasis))
  
  }

rmse[j] = max(tmp) 
rmse.PARABOLIC[j] = max(tmp.PARABOLIC)

mise[j] = max(tmp.mise)
mise.PARABOLIC[j] = max(tmp.mise.PARABOLIC)

save(rmse, sols, sols.PARABOLIC,
     coef.ex, mean.coef, mean.coef.PARABOLIC,
     times, times.PARABOLIC,
     mise, mise.PARABOLIC,
     file=save.file)
}
end.tot = Sys.time()
tot.time  = difftime(end.tot, start.tot, units="mins")
tot.time 

today_ = Sys.Date()
ntest_ ="-test-1"
path_ = "/home/aldo/Scrivania/fdaPDE-DATA/TS-Regression"

file.name = paste("TS-Regression-",today_,ntest_,sep="")

save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/TS-Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/TS-Regression/",file.name,".pdf",sep="")

save(rmse, sols, sols.PARABOLIC,
     coef.ex, mean.coef, mean.coef.PARABOLIC,
     times, times.PARABOLIC,
     mise, mise.PARABOLIC,
     file=save.file)

################################################################################

for( i in 1:NumTimePoints){

  print( norm( diff[(1 + (i-1)*nnodes) : ( nnodes + (i-1)*nnodes)], type="2") )
}

#coef.estimated = matrix(0, nrow=102, ncol=1)
min.col = min(coef.estimated, coef.ex)
max.col = max(coef.estimated, coef.ex)


for( i in 1:NumTimePoints){
  
  print( R_plot_graph.ggplot2.2(FEM(coef.estimated[(1 + (i-1)*nnodes):( nnodes + (i-1)*nnodes ),1], 
                                    FEMbasis),
                                paste("t = ",TimePoints[i],sep=""),
                                line.size = 1,
                                color.min = min.col,
                                color.max = max.col))
}

# True Space-Time field
for( i in 1:NumTimePoints){
  
  print( R_plot_graph.ggplot2.2(FEM(coef.ex[,i], FEMbasis),
                                line.size = 1,
                                title = paste("t = ",TimePoints[i],sep=""),
                                color.min = min.col,
                                color.max = max.col))
}


nrow(mesh$nodes)

#####################################################

save.file = "/home/aldo/Scrivania/fdaPDE-DATA/TS-Regression/TS-Regression-2022-01-27-test-1.RData"
load(save.file)

# mean time-space field
{
coef.estimated = matrix(0, nrow=length(coef.ex) , ncol=1 )

for( i in 1:M){
  output_CPP = sols[[i]]
  coef.estimated = coef.estimated + eval.FEM.time( output_CPP$fit.FEM.time, locations=mesh$nodes, time.instants =  TimePoints,
                                lambdaS = output_CPP$bestlambda[1],
                                lambdaT = output_CPP$bestlambda[2])
}

coef.estimated = coef.estimated / M
}

mean.coef = matrix(0, nrow=199, ncol=1)

for(j in 1:M){
  mean.coef =  mean.coef + sols[j,4,]
}

mean.coef = mean.coef / M

mesh.ref = refine.mesh.1.5D(mesh, 0.01)
FEMbasis.ref = create.FEM.basis(mesh.ref)

mean.coef.ref = eval.FEM( FEM(mean.coef, FEMbasis), mesh.ref$nodes)
coef.ex.ref = eval.FEM(FEM(coef.ex, FEMbasis), mesh.ref$nodes)



min.col = min(coef.estimated, coef.ex)
max.col = max(coef.estimated, coef.ex)
# plot true field

pdf(img.file)
for( i in 1:NumTimePoints){
  
  print( R_plot_graph.ggplot2.2(FEM(coef.ex[,i], FEMbasis),
                                line.size = 1,
                                title = paste("t = ",TimePoints[i],sep=""),
                                color.min = min.col,
                                color.max = max.col))
}

# plot mean time-space field
for( i in 1:NumTimePoints){
  
  print( R_plot_graph.ggplot2.2(FEM(coef.estimated[(1 + (i-1)*nnodes):( nnodes + (i-1)*nnodes ),1], 
                                    FEMbasis),
                                paste("t = ",TimePoints[i],sep=""),
                                line.size = 1,
                                color.min = min.col,
                                color.max = max.col))
}

for( i in 1:NumTimePoints){
  
   R_plot_graph.R(FEM(coef.ex[,i], FEMbasis),
                                line.size = 3,
                                title = paste("t = ",TimePoints[i],sep=""))

R_plot_graph.R(FEM(coef.estimated[(1 + (i-1)*nnodes):( nnodes + (i-1)*nnodes ),1], FEMbasis),
                       line.size = 3,
                       title = paste("t = ",TimePoints[i],sep=""))

}

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


data_frame.1 = data.frame(rmse_ = rmse, nobs_ = rep(as.character(102), times=M))
ggplot(data_frame.1)+
  geom_boxplot(aes(x=nobs_, y=as.vector(rmse_), fill="red"))+
  labs(x="ST-PDE", y="", title="RMSE")+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none", 
        axis.text.x= element_blank())+
  MyTheme


dev.off()


min.col = min(coef.ex)
max.col = max(coef.ex)


ggplot_list = list()

for( i in 1:length(TimePoints)){
  
  if(i == 1)
    title_ = TeX(sprintf("$t = 0"))
  else if(i==2)
    title_ = TeX(sprintf("$t = \\pi / 5$"))
  else if(i != 6)
    title_ = TeX(sprintf("$t =  %d \\pi / 5$",i-1))
  else
    title_ = TeX(sprintf("$t =  \\pi$"))
    
  
  #print(
    ggplot_list[[i]] <- R_plot_graph.ggplot2.2(FEM(coef.ex[,i], FEMbasis),
                                line.size = 1,
                                title = title_,
                                color.min = min.col,
                                color.max = max.col,
                                ratio=1) 
    #) 
}


library(gridExtra)

pdf("/home/aldo/Scrivania/Caricare/TS-Regression-test-1.pdf",
    height=9,width = 9)
grid.arrange(grobs = ggplot_list, 
             layout_matrix = rbind(c(1,2),
                                   c(3,4),
                                   c(5,6)) )



grid.arrange(grobs = ggplot_list)

mise_ = cbind(mise, mise.PARABOLIC)
mise_ = as.vector(mise_)
type_ = rep(c("separable", "parabolic"), each=20) 
nobs_ = rep(1, times=40)

data_frame.1 = data.frame(nobs_, mise_, type_)

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

ggplot(data_frame.1)+
  geom_boxplot(aes(x=nobs_, y=mise_, fill=type_))+
  labs(x="", y="", title="MISE")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.85,0.875),
        axis.text.x= element_blank())+
  MyTheme
dev.off()
