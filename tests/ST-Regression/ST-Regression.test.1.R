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

vertices = cbind(simplenet$vertices$x, simplenet$vertices$y)
edges = cbind(simplenet$from , simplenet$to)
L = as.linnet(simplenet)

delta = 0.03 # 102 nodi, 0.015 199 nodi
mesh =create.mesh.1.5D(vertices, edges)
mesh =refine.mesh.1.5D(mesh,delta=0.03)

nnodes = nrow(mesh$nodes)
FEMbasis = create.FEM.basis(mesh)

source = 4

aux = function(x, y, t, x0=mesh$nodes[source,1], y0=mesh$nodes[source,2]){
  sigma = 0.03
  dist.2 = (x-x0)^2 + (y-y0)^2
  
  value =  exp( 1.0 /( dist.2/(5*sigma)^2 - 1 ) + 1 ) * exp(-t) *( sqrt(dist.2) < 5*sigma) 
}

NumTimePoints=5
TimePoints=seq(0,1,length.out = NumTimePoints)

SpacePoints=mesh$nodes
SpaceTimePoints=cbind(rep(mesh$nodes[,1],NumTimePoints),rep(mesh$nodes[,2],NumTimePoints),rep(TimePoints,each=nrow(mesh$nodes)))

coef.ex = matrix(0, nrow= nrow(mesh$nodes), ncol=NumTimePoints)
for(i in 1:NumTimePoints){
  coef.ex[,i] = aux(mesh$nodes[,1],
                    mesh$nodes[,2],
                    TimePoints[i])
  
}
coef.0 = aux(mesh$nodes[,1], mesh$nodes[,2],t=0)
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

# Exact solution (pointwise at nodes)
sol_exact=aux(SpaceTimePoints[,1],
              SpaceTimePoints[,2],
              SpaceTimePoints[,3])

# Add error to simulate data
M = 20 
sols = list()
rmse = matrix(0, nrow=M, ncol=1)
mean.coef = matrix(0, nrow=102*NumTimePoints, ncol=1)

start = Sys.time()  
for( j in 1:M){
ran = range(sol_exact)
data = sol_exact + rnorm(length(sol_exact), mean=0, sd=0.05*abs(ran[2]-ran[1]))
observations=matrix(data,nrow(SpacePoints),NumTimePoints)

lambdaS = 10^(-5:-3)
lambdaT = 10^(-5:-3)
output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                            lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

coef.estimated = eval.FEM.time( output_CPP$fit.FEM.time, locations=mesh$nodes, time.instants =  TimePoints,
                                lambdaS = output_CPP$bestlambda[1],
                                lambdaT = output_CPP$bestlambda[2])
mean.coef = mean.coef + coef.estimated / M

sols[[i]] = output_CPP
diff = coef.estimated - as.vector(coef.ex)

tmp = matrix(0, nrow=NumTimePoints, ncol=1)
for( i in 1:NumTimePoints){
  
  tmp[i] = norm( diff[(1 + (i-1)*nnodes) : ( nnodes + (i-1)*nnodes)], type="2") 
}

rmse[j] = max(tmp) 

}
end = Sys.time()
tot.time  = difftime(end, start, units="mins")
tot.time 

today_ = Sys.Date()
ntest_ ="-test-1-18-00"
path_ = "/home/aldo/Scrivania/fdaPDE-DATA/TS-Regression"

file.name = paste("TS-Regression-",today_,ntest_,sep="")

save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/TS-Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/TS-Regression/",file.name,".pdf",sep="")

save(rmse, sols, coef.ex, mean.coef,
     tot.time, file=save.file)

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
