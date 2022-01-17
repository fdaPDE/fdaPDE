### DE post processing ### 
library(ggplot2)
library(latex2exp)

MyTheme <- theme(
  axis.text = element_text(size=16),
  axis.title = element_text(size=16),
  title = element_text(size=20),
  legend.text = element_text(size=14),
  legend.key.size = unit(1,"cm") 
)

# load data # DE-2021-12-17-heat-test-3-RIGHT_CV
date_ = "2021-12-19" # "2021-12-17"
init_ = "-heat" #"-fdaPDE" "-null_vector"
ntest_ = "-test-5-RightCV"


file.name = paste("DE-",date_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

load(save.file)
##

mise_ =  cbind(as.vector(mise) , as.vector(mise.fdaPDE))
mise_= as.vector(mise_)
times_ = cbind( as.vector(times),as.vector(times.fdaPDE))
times_= as.vector(times_)
nobs_ = rep(nobs,each=M)
nobs_ = rep(as.character(nobs_),times=2)
type_ = rep(c("R","fdaPDE"),each=(M*N))

# data frame per boxplots #
data_frame.1 = data.frame(nobs_,
                          mise_,times_,type_)


#data frame per scatterplots#
pdf(img.file)

ggplot(data_frame.1)+
  geom_boxplot(aes(x=nobs_,y=mise_, group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y= "",fill="",title="MISE")+
  scale_x_discrete(limits=as.character(nobs) ) +
  MyTheme +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data_frame.1)+
  geom_boxplot(aes(x=nobs_,y=times_, group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title="Running Time [mins]",subtitle="(fdaPDE + CV)")+
  scale_x_discrete(limits=as.character(nobs) ) +
  MyTheme +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5))

mise_ = as.vector(mise)
mise.fdaPDE_ = as.vector(mise.fdaPDE)
times_ = as.vector(times)
times.fdaPDE_ = as.vector(times.fdaPDE)
err.L2_ = as.vector(err.L2)
nobs_ = rep(as.character(nobs),each=M)
# data frame per scatterplots #
data_frame.2 = data.frame(nobs_,
                          mise_,mise.fdaPDE_,err.L2_,
                          times_,times.fdaPDE_)

ggplot(data_frame.2)+
  geom_boxplot(aes(x=nobs_,y=err.L2_,group=nobs_))+
  labs(x="observations",y="",fill="",
       title=TeX(  "$||f_{R} - f_{fdaPDE}||_{L_{2}}^{2}$")) +
  scale_x_discrete(limits=as.character(nobs) ) +
  MyTheme +
  theme(plot.title = element_text(hjust = 0.5))

for( i in 1:N){
  data_frame.tmp = subset(data_frame.2, nobs_ %in% c(nobs[i]))
  
  x.min = min(data_frame.tmp["mise_"],data_frame.tmp["mise_"])
  x.max = max(data_frame.tmp["mise_"],data_frame.tmp["mise_"])
  
  print(  ggplot(data_frame.tmp,aes(x=mise_,y=mise.fdaPDE_))+
            geom_point() +
            geom_segment(aes(x=x.min, y=x.min, xend=x.max, yend=x.max),color="red",alpha=0.01)+
            labs(x="R", y="fdaPDE", 
                 title= TeX("Scatter $ MISE$"), 
                 subtitle= paste("nobs = ",nobs[i],sep="") ) +
            MyTheme +
            theme(plot.title = element_text(hjust=0.5), 
                  plot.subtitle = element_text(hjust=0.5) )    )
  
}

for(i in 1:N){
  title.mise = paste("MISE, nobs =",as.character(nobs[i]),seq="")
  plot(1:M,mise[,i], pch=16, col="red", xlab="",ylab="",
       main=title.mise ,ylim=c(min(mise[,i],mise.fdaPDE[,i]),max(mise[,i],mise.fdaPDE[,i])) )
  points(1:M,mise.fdaPDE[,i], pch=13, col ="blue")
  segments(1:M,mise[,i],1:M,mise.fdaPDE[,i], lty=3)
  legend("topright",legend=c("Mattina","fdaPDE"),col=c("red","blue"),pch=c(16,13) )
}

#ref mesh

coef.fda = sols.fdaPDE[1,4,] 
coef.R = sols[1,4,]

R_plot_graph.a.sym.ggplot2(FEM(coef.ex,FEMbasis.fdaPDE))
R_plot_graph.a.sym.ggplot2(FEM(coef.fda,FEMbasis.fdaPDE))
R_plot_graph.a.sym.ggplot2(FEM(coef.R, FEMbasis.fdaPDE))
###

# point pattern #

PP = rlpp(300, density)
plot(mesh.fdaPDE)
points(PP$data$x, PP$data$y, pch=16, col="darkgreen",cex=0.75)
dev.off()

################################################################################
######################## increasing nodes number ###############################
################################################################################

MyTheme <- theme(
  axis.text = element_text(size=16),
  axis.title = element_text(size=16),
  title = element_text(size=20),
  legend.text = element_text(size=14),
  legend.key.size = unit(1,"cm") 
)

# load data # DE-2021-12-17-heat-test-3-RIGHT_CV
date_ = "" # "2021-12-17"
init_ = "-heat" #"-fdaPDE" "-null_vector"
ntest_ = "-test-1-delta-var"

file.name = paste("DE-",date_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

load(save.file)
##

mise_ =  cbind(as.vector(mise) , as.vector(mise.fdaPDE))
mise_= as.vector(mise_)
times_ = cbind( as.vector(times),as.vector(times.fdaPDE))
times_= as.vector(times_)
nnodes_ = rep(nnodes,each=M)
nnodes_ = rep(as.character(nnodes_),times=2)
type_ = rep(c("R","fdaPDE"),each=(M*N))

# data frame per boxplots #
data_frame.1 = data.frame(nnodes_,
                          mise_,times_,type_)


#data frame per scatterplots#
pdf(img.file)

ggplot(data_frame.1)+
  geom_boxplot(aes(x=nnodes_,y=mise_, group=interaction(nnodes_,type_),fill=type_))+
  labs(x="nodes",y= "",fill="",title="MISE")+
  scale_x_discrete(limits=as.character(nnodes))+
  MyTheme +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data_frame.1)+
  geom_boxplot(aes(x=nnodes_,y=times_, group=interaction(nnodes_,type_),fill=type_))+
  labs(x="nodes",y="",fill="",title="Running Time [mins]",subtitle="(fdaPDE + CV)")+
  scale_x_discrete(limits=as.character(nnodes))+
  MyTheme +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5))

mise_ = as.vector(mise)
mise.fdaPDE_ = as.vector(mise.fdaPDE)
times_ = as.vector(times)
times.fdaPDE_ = as.vector(times.fdaPDE)
err.L2_ = as.vector(err.L2)
nnodes_ = rep(as.character(nnodes), each=M)
# data frame per scatterplots #
data_frame.2 = data.frame(nnodes_,
                          mise_,mise.fdaPDE_,err.L2_,
                          times_,times.fdaPDE_)

ggplot(data_frame.2)+
  geom_boxplot(aes(x=nnodes_,y=err.L2_,group=nnodes_))+
  labs(x="nodes",y="",fill="",
       title=TeX(  "$||f_{R} - f_{fdaPDE}||_{L_{2}}^{2}$")) +
  MyTheme +
  scale_x_discrete(limits=as.character(nnodes))+
  theme(plot.title = element_text(hjust = 0.5))


avg.times = vector(mode="numeric", length=N)
avg.fdaPDE.times = vector(mode="numeric", length=N)
for( i in 1:N){
  avg.times[i] = mean(times[,i])
  avg.fdaPDE.times[i] = mean(times.fdaPDE[,i])
}

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
