# loading data
date_ = "2021-12-10"
ntest_ = "-test-2" 

file.name = paste("Regression-",date_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

load(save.file)

# postprocessing #
library(ggplot2)
library(latex2exp)

MyTheme <- theme(
  axis.text = element_text(size=16),
  axis.title = element_text(size=16),
  title = element_text(size=20),
  legend.text = element_text(size=14),
  legend.key.size = unit(1,"cm") 
)

rmse.beta_ = cbind( as.vector(rmse.beta), as.vector(rmse.fdaPDE.beta)) 
rmse.beta_ = as.vector(rmse.beta_)
mise_ = cbind( as.vector(mise), as.vector(mise.fdaPDE))
mise_ = as.vector(mise_)
times_ = cbind( as.vector(times), as.vector(times.fdaPDE)) 
times_ = as.vector(times_)
nobs_ = rep(nobs,each=M)
nobs_ = rep(as.character(nobs_), times=2) 
type_ = rep(c("R","fdaPDE"),each=(M*N))

err.L2_ = as.vector(err.L2)

# data_frame boxplots #
data_frame.1 = data.frame(rmse.beta_, mise_, times_, nobs_, type_)

pdf(img.file)

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=rmse.beta_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="",y="",fill="",title=TeX( "$RMSE(\\beta)^{2}$") )+
  scale_x_discrete(limits=c("300","500","700","1000"))+
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=mise_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title=TeX( "MISE"))+
  scale_x_discrete(limits=c("300","500","700","1000"))+
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme

ggplot(data_frame.1) + 
  geom_boxplot(aes(x=nobs_,y=times_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title="TIME [s]") +
  scale_x_discrete(limits=c("300","500","700","1000"))+
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme

rmse.beta_ = as.vector(rmse.beta)
rmse.fdaPDE.beta_ = as.vector(rmse.fdaPDE.beta)
mise_ = as.vector(mise)
mise.fdaPDE_ = as.vector(mise.fdaPDE)
times_ = as.vector(times)
times.fdaPDE_ = as.vector(times.fdaPDE)
err.L2_ = as.vector(err.L2)
nobs_ = rep(as.character(nobs),each=M)

data_frame.2 = data.frame(nobs_,
                          rmse.beta_,
                          rmse.fdaPDE.beta_,
                          mise_,mise.fdaPDE_,err.L2_,
                          times_,times.fdaPDE_)

ggplot(data_frame.2)+
  geom_boxplot(aes(x=nobs_,y=err.L2_,group=nobs_))+
  labs(x="observations",y="",fill="",
       title= TeX("\ $||f_{R} - f_{fdaPDE}||_{L_{2}}^2$") ) +
  scale_x_discrete(limits=c("300","500","700","1000"))+
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme

for(i in 1:N){
  data_frame.tmp = subset(data_frame.2, nobs_ %in% c(nobs[i]))
  
  x.min = min(data_frame.tmp["rmse.beta_"], data_frame.tmp["rmse.fdaPDE.beta_"])
  x.max = max(data_frame.tmp["rmse.beta_"], data_frame.tmp["rmse.fdaPDE.beta_"])
  
  print( ggplot(data_frame.tmp, aes(x=rmse.beta_,y=rmse.fdaPDE.beta_)) +
           geom_point() +
           geom_segment(aes(x=x.min, y=x.min, xend=x.max, yend=x.max),color="red",alpha=0.01)+  
           labs(x="R", y="fdaPDE", 
                title= TeX("Scatter  $RMSE(\\beta)^{2}$"),
                subtitle= paste("nobs = ",nobs[i],sep=""))+
           theme(plot.title = element_text(hjust=0.5),
                 plot.subtitle = element_text(hjust=0.5))+
           MyTheme )
  
  x.min = min(data_frame.tmp["mise_"],data_frame.tmp["mise.fdaPDE_"])
  x.max = max(data_frame.tmp["mise_"],data_frame.tmp["mise.fdaPDE_"])
  
  print(
    ggplot(data_frame.tmp, aes(x=mise_,y=mise.fdaPDE_)) +
      geom_point() +
      geom_segment(aes(x=x.min, y=x.min, xend=x.max, yend=x.max),color="red",alpha=0.01) +
      labs(x="R", y="fdaPDE", 
           title= TeX("Scatter $MISE$"),
           subtitle= paste("nobs = ",nobs[i],sep="") )+
      theme(plot.title = element_text(hjust=0.5),
            plot.subtitle = element_text(hjust=0.5))+
      MyTheme)
  
}

FEM.ex = FEM(coef.ex, FEMbasis.fdaPDE)
FEM.fdaPDE = FEM( sols.fdaPDE[ which( mise.fdaPDE[,4] <= min(mise.fdaPDE[,4])), 4,], FEMbasis.fdaPDE)
FEM.R = FEM( sols[ which( mise[,4] <= min(mise[,4])), 4,], FEMbasis.fdaPDE )
R_plot_graph.ggplot2(FEM.ex)
R_plot_graph.ggplot2(FEM.fdaPDE)
R_plot_graph.ggplot2(FEM.R)

dev.off()



########################################################################
##### Regression PostProcessing Test-2 delta-var #######################
########################################################################
library(latex2exp)

date_ = "2022-01-12-"
ntest_ = "test-2-delta-var" # "test-2"

MyTheme <- theme(
  axis.text = element_text(size=16),
  axis.title = element_text(size=16),
  title = element_text(size=20),
  legend.text = element_text(size=14),
  legend.key.size = unit(1,"cm") 
)


file.name = paste("Regression-",date_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

load(save.file)

rmse.beta_ = cbind( as.vector(rmse.beta), as.vector(rmse.fdaPDE.beta)) 
rmse.beta_ = as.vector(rmse.beta_)
mise_ = cbind( as.vector(mise), as.vector(mise.fdaPDE))
mise_ = as.vector(mise_)
times_ = cbind( as.vector(times), as.vector(times.fdaPDE)) 
times_ = as.vector(times_)
delta_ = rep(delta,each=M)
delta_ = rep(as.character(delta_),times=2)
nodes_ = rep(nnodes,each=M)
nodes_ = rep(as.character(nodes_), times=2 )
type_ = rep(c("R","fdaPDE"),each=(M*N))

err.L2_ = as.vector(err.L2)

# data_frame boxplots #
data_frame.1 = data.frame(rmse.beta_, 
                          mise_,times_,delta_,type_, nodes_)

pdf(img.file)

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nodes_,y=rmse.beta_,group=interaction(delta_,type_),fill=type_))+
  labs(x="nodes",y="",fill="",title=TeX( "$RMSE(\\beta_{1})^{2}$") )+
  scale_x_discrete(limits=c(as.character(nnodes)))+
  theme(plot.title = element_text(hjust = 0.5)) +
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nodes_,y=mise_,group=interaction(delta_,type_),fill=type_))+
  labs(x="nodes",y="",fill="",title=TeX( "MISE"))+
  scale_x_discrete(limits= as.character(nnodes))+
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme

ggplot(data_frame.1) + 
  geom_boxplot(aes(x=delta_,y=times_,group=interaction(delta_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title="TIME [s]") +
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme

avg.times = vector(mode="numeric", length=N)
avg.fdaPDE.times = vector(mode="numeric", length=N)
for( i in 1:N){
  avg.times[i] = mean(times[,i])
  avg.fdaPDE.times[i] = mean(times.fdaPDE[,i])
}

avg.times
avg.fdaPDE.times
avg.times_ = cbind(avg.times, avg.fdaPDE.times)
avg.times_ = as.vector(avg.times_)
type_ = rep(c("R","fdaPDE"), each=N)
nnodes_ = rep(as.character(nnodes), times = 2)

data.frame.3 = data.frame(avg.times_, nnodes_, type_)

ggplot(data=data.frame.3, aes(x = nnodes_, y=log2(avg.times_), group=type_, fill=type_ ))+
  geom_point(aes(shape=type_,color=type_))+
  geom_line(aes(color=type_), linetype="dashed")+
  scale_x_discrete(limits=as.character(nnodes)) + 
  labs(x=TeX("$nodes$"),y=TeX("$log_{2}(times)$"), title="Average running time")+
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())+
  MyTheme

dev.off()

