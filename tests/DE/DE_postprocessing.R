### DE post processing ### 
library(ggplot2)
library(latex2exp)
library(scales)

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


# load data # DE-2021-12-17-heat-test-3-RIGHT_CV
date_ = "2022-01-16" # "2022-01-14"
init_ = "-heat" # 
ntest_ =  "-test-5-SimplifiedCV"# "-test-1" "-test-5-SimplifiedCV"


file.name = paste("DE-",date_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/","DE-test-2-parte-1-500-obs",".pdf",sep="")
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
  theme(legend.position= c(0.825,0.85))

ggplot(data_frame.1)+
  geom_boxplot(aes(x=nobs_,y=times_, group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title="Running Time [mins]",subtitle="(fdaPDE + CV)")+
  scale_x_discrete(limits=as.character(nobs) ) +
  MyTheme +
  theme(legend.position= c(0.825,0.85), plot.subtitle=element_text(hjust=0.5))

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

dev.off()
#ref mesh # test 1
# coef.fda = matrix(0,nrow=97,ncol=1)
# coef.R = matrix(0,nrow=97,ncol=1)
# 
# for(j in 1:M){
#   coef.fda = coef.fda + sols.fdaPDE[j,4,]
#   coef.R = coef.R + sols[j,4,]
# }
# 
# coef.fda = coef.fda / M
# coef.R = coef.R / M
# 
# 
# R_plot_graph.a.sym.ggplot2(FEM(coef.ex,FEMbasis.fdaPDE))
# R_plot_graph.a.sym.ggplot2(FEM(coef.fda,FEMbasis.fdaPDE))
# R_plot_graph.a.sym.ggplot2(FEM(coef.R, FEMbasis.fdaPDE))
# ######
# vertices = cbind(simplenet$vertices$x, simplenet$vertices$y)
# edges = cbind(simplenet$from, simplenet$to)
# #####
# delta=0.015
# mesh = create.mesh.1D.vertices(vertices, edges, delta)
# 
# coef.fda.ref = eval.FEM(FEM(coef.fdaPDE, FEMbasis.fdaPDE), mesh$nodes)
# coef.R.ref = eval.FEM(FEM(coef.R, FEMbasis.fdaPDE), mesh$nodes)
# coef.ex.ref = eval.FEM(FEM(coef.ex, FEMbasis.fdaPDE), mesh$nodes)
# 
# mesh.fdaPDE.ref = create.mesh.1.5D(mesh$nodes, mesh$segments)
# FEMbasis.fdaPDE.ref = create.FEM.basis(mesh=mesh.fdaPDE.ref)
###
###ref mesh # test 2
coef.fda = matrix(0,nrow=412,ncol=1)
coef.R = matrix(0,nrow=412,ncol=1)

for(j in 1:M){
  coef.fda = coef.fda + sols.fdaPDE[j,4,]
  coef.R = coef.R + sols[j,4,]
}

coef.fda = coef.fda / M
coef.R = coef.R / M


R_plot_graph.a.sym.ggplot2(FEM(coef.ex,FEMbasis.fdaPDE))
R_plot_graph.a.sym.ggplot2(FEM(coef.fda,FEMbasis.fdaPDE))
R_plot_graph.a.sym.ggplot2(FEM(coef.R, FEMbasis.fdaPDE))
######
data("spiders")
vertices = cbind(spiders$domain$vertices$x, spiders$domain$vertices$y) 
edges = cbind(spiders$domain$from, spiders$domain$to)
#####
delta=30
mesh = create.mesh.1D.vertices(vertices, edges, delta)

coef.fda.ref = eval.FEM(FEM(coef.fda, FEMbasis.fdaPDE), mesh$nodes)
coef.R.ref = eval.FEM(FEM(coef.R, FEMbasis.fdaPDE), mesh$nodes)
coef.ex.ref = eval.FEM(FEM(coef.ex, FEMbasis.fdaPDE), mesh$nodes)

mesh.fdaPDE.ref = create.mesh.1.5D(mesh$nodes, mesh$segments)
FEMbasis.fdaPDE.ref = create.FEM.basis(mesh=mesh.fdaPDE.ref)
###


# point pattern #

PP = rlpp(300, density)
plot(mesh.fdaPDE)
points(PP$data$x, PP$data$y, pch=16, col="darkgreen",cex=0.75)


################################################################################
######################## increasing nodes number ###############################
################################################################################
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

# load data # DE-2021-12-17-heat-test-3-RIGHT_CV
date_ = "2022-01-16" # "2021-12-16"
init_ = "-heat" #"-fdaPDE" "-null_vector"
ntest_ = "-test-1-delta-var"

file.name = paste("DE-",date_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

load(save.file)
##

img.file= "/home/aldo/Scrivania/Caricare/DE-test-1-300-obs-parte-2.pdf"
mise_ =  cbind(as.vector(mise) , as.vector(mise.fdaPDE)) #[,2:4]
mise_= as.vector(mise_)
times_ = cbind( as.vector(times),as.vector(times.fdaPDE)) #[,2:4]
times_= as.vector(times_)
nnodes_ = rep(nnodes,each=M) #[,2:4]
nnodes_ = rep(as.character(nnodes_),times=2)
type_ = rep(c("R","fdaPDE"),each=(M*(N)) ) #(N-1)

# data frame per boxplots #
data_frame.1 = data.frame(nnodes_,
                          mise_,times_,type_)


#data frame per scatterplots#
pdf(img.file)

ggplot(data_frame.1)+
  geom_boxplot(aes(x=nnodes_,y=mise_, group=interaction(nnodes_,type_),fill=type_))+
  labs(x="nodes",y= "",fill="",title="MISE")+
  scale_x_discrete(limits=as.character(nnodes))+ #[2:4]
  MyTheme +
  theme(legend.position= c(0.825,0.85))

ggplot(data_frame.1)+
  geom_boxplot(aes(x=nnodes_,y=times_, group=interaction(nnodes_,type_),fill=type_))+
  labs(x="nodes",y="",fill="",title="Running Time [mins]",subtitle="(fdaPDE + CV)")+
  scale_x_discrete(limits=as.character(nnodes))+ #[2:4]
  MyTheme +
  theme(legend.position= c(0.825,0.85))

mise_ = as.vector(mise) #[,2:4]
mise.fdaPDE_ = as.vector(mise.fdaPDE) #[2:4]
times_ = as.vector(times) # [,2:4]
times.fdaPDE_ = as.vector(times.fdaPDE) #[,2:4]
err.L2_ = as.vector(err.L2) #[,2:4]
nnodes_ = rep(as.character(nnodes), each=M) #[2:4]
# data frame per scatterplots #
data_frame.2 = data.frame(nnodes_,
                          mise_,mise.fdaPDE_,err.L2_,
                          times_,times.fdaPDE_)

ggplot(data_frame.2)+
  geom_boxplot(aes(x=nnodes_,y=err.L2_,group=nnodes_))+
  labs(x="nodes",y="",fill="",
       title=TeX(  "$||f_{R} - f_{fdaPDE}||_{L_{2}}^{2}$")) +
  MyTheme +
  scale_x_discrete(limits=as.character(nnodes))+ #[2:4]
  theme(plot.title = element_text(hjust = 0.5))


avg.times = vector(mode="numeric", length=(N)) #(N-1)
avg.fdaPDE.times = vector(mode="numeric", length=(N)) #(N-1)
for( i in 1:(N)){ # (N-1)
  avg.times[i] = mean(times[,i])
  avg.fdaPDE.times[i] = mean(times.fdaPDE[,i])
}

avg.times_ = cbind(avg.times, avg.fdaPDE.times)
avg.times_ = as.vector(avg.times_)
type_ = rep(c("R","fdaPDE"), each=(N)) #(N-1)
nnodes_ = rep(nnodes, times = 2) #[2:4]
data.frame.3 = data.frame(avg.times_, nnodes_, type_)

library(scales)
# EXPS = seq(from=min(log2(nnodes[2:4])),to=max(log2(nnodes[2:4])), length.out=3)
# x.1 = 2^EXPS
# y.2 = 2^( EXPS - 10)
# y.3 = 2^(2* (EXPS -10))
# x_ = rep(x.1, times = 2)
# y_ = c(y.2, y.2)

ggplot(data=data.frame.3, aes(x = nnodes_, y=avg.times_, group=type_, fill=type_ ))+
  #geom_line(aes(x=nnodes_, y=nnodes_) , linetype="dashed", size=0.5, alpha=0.8)+
  geom_line(aes(x=nnodes_, y=nnodes_/2^12) , linetype="dashed", size=0.5, alpha=0.8)+
  geom_line(aes(x=nnodes_, y=nnodes_^2/2^16) , linetype="dashed", size=0.5, alpha=0.8)+
  geom_point(aes(shape=type_,color=type_),size=2)+
  geom_line(aes(color=type_), linetype="dashed",size=1)+
  scale_x_continuous(trans="log2", 
                     breaks=trans_breaks("log2",function(x) 2^x),
                     labels=trans_format("log2", math_format(.x) ) )+
  scale_y_continuous(trans="log2", 
                     breaks=trans_breaks("log2",function(x) 2^x),
                     labels=trans_format("log2", math_format(.x) ) )+
  labs(x=TeX("$log_{2}(nodes)$"),y=TeX("$log_{2}(times)$"), title="Average running time")+
  theme(legend.position= c(0.175,0.85))+
  MyTheme

dev.off()


################################################################################
########################      fdaPDE ONLY        ###############################
######################## increasing nodes number ###############################
################################################################################
library(scales)
library(latex2exp)

MyTheme <- theme(
  axis.text = element_text(size=24),
  axis.title = element_text(size=26),
  title = element_text(size=26),
  legend.text = element_text(size=26),
  legend.key.size = unit(1,"cm") 
)

# load data # DE-2021-12-17-heat-test-3-RIGHT_CV
date_ = "2022-01-16" # "2021-12-16"
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

mise_ = as.vector(mise.fdaPDE)
times_ = as.vector(times.fdaPDE)
nnodes_ = rep(nnodes, each=M)
nnodes_=as.character(nnodes_)
type_ = rep(c("fdaPDE"), times=M*N)

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

#mise_ = as.vector(mise)
mise.fdaPDE_ = as.vector(mise.fdaPDE)
#times_ = as.vector(times)
times.fdaPDE_ = as.vector(times.fdaPDE)
#err.L2_ = as.vector(err.L2)
nnodes_ = rep(as.character(nnodes), each=M)
# data frame per scatterplots #
data_frame.2 = data.frame(nnodes_,
                          mise.fdaPDE_,
                          times.fdaPDE_)

#avg.times = vector(mode="numeric", length=N)
avg.fdaPDE.times = vector(mode="numeric", length=N)
for( i in 1:N){
 # avg.times[i] = mean(times[,i])
  avg.fdaPDE.times[i] = mean(times.fdaPDE[,i])
}

avg.times_ = avg.fdaPDE.times

type_ = rep(c("fdaPDE"), each=N)
nnodes_ = nnodes
data.frame.3 = data.frame(avg.times_, nnodes_, type_)

EXPS = seq(from=min(log2(nnodes)),to=max(log2(nnodes)), length.out=5)
x.1 = 2^EXPS
y.1 = 2^EXPS
y.2 = 2^( EXPS - 10)
y.3 = 2^(EXPS - 1)
y.4 = 2^(EXPS+ 1)
y.3 = 2^(2* (EXPS -10))
x_ = rep(x.1, times = 1)
y_ = rep(y.2, times = 1)


ggplot(data=data.frame.3, aes(x = nnodes_, y=avg.times_, group=type_, fill=type_ ))+
  geom_line(aes(x=x_, y=y_) , linetype="dashed", size=0.5, alpha=0.8)+
  geom_point(aes(shape=type_,color=type_),size=2)+
  geom_line(aes(color=type_), linetype="dashed",size=1)+
  scale_x_continuous(trans="log2", 
                     breaks=trans_breaks("log2",function(x) 2^x),
                     labels=trans_format("log2", scales::math_format(2^.x) ) )+
  scale_y_continuous(trans="log2", 
                     breaks=trans_breaks("log2",function(x) 2^x),
                     labels=trans_format("log2", scales::math_format(2^.x ) ) )+
  labs(x=TeX("$nodes$"),y=TeX("$log_{2}(times)$"), title="Average running time")+
  theme( legend.position= c(0.125,0.85))+
  MyTheme

dev.off()


