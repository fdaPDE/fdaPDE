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
  scale_x_discrete(limits=c("300","500","700","1000","1500"))+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=mise_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title=TeX( "MISE"))+
  scale_x_discrete(limits=c("300","500","700","1000","1500"))+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data_frame.1) + 
  geom_boxplot(aes(x=nobs_,y=times_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title="TIME [s]") +
  scale_x_discrete(limits=c("300","500","700","1000","1500"))+
  theme(plot.title = element_text(hjust = 0.5))

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
  scale_x_discrete(limits=c("300","500","700","1000","1500"))+
  theme(plot.title = element_text(hjust = 0.5))

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
                 plot.subtitle = element_text(hjust=0.5))  )
  
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
            plot.subtitle = element_text(hjust=0.5)) )
  
}

FEM.ex = FEM(coef.ex, FEMbasis.fdaPDE)
FEM.fdaPDE = FEM( sols.fdaPDE[ which( mise.fdaPDE[,4] <= min(mise.fdaPDE[,4])), 4,], FEMbasis.fdaPDE)
FEM.R = FEM( sols[ which( mise[,4] <= min(mise[,4])), 4,], FEMbasis.fdaPDE )
R_plot_graph.ggplot2(FEM.ex)
R_plot_graph.ggplot2(FEM.fdaPDE)
R_plot_graph.ggplot2(FEM.R)

dev.off()
