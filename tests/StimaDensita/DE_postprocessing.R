### DE post processing ### 
library(ggplot2)
library(latex2exp)

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
  labs(x="observations",y= "MISE",fill="",title="f")+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data_frame.1)+
  geom_boxplot(aes(x=nobs_,y=times_, group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title="Running Time [mins]",subtitle="(fdaPDE + CV)")+
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
source("~/Scrivania/fdaPDE/tests/R_plot_graph.ggplot2.R")
  
#ref mesh
mesh.ref = create.mesh.1D.vertices(vertices, edges, delta=30)
mesh.fdaPDE.ref = create.mesh.1.5D(mesh.ref$nodes, mesh.ref$segments)
FEMbasis.fdaPDE.ref = create.FEM.basis(mesh.fdaPDE.ref)

coef.25.ref = aux.25(mesh.fdaPDE.ref$nodes[,1], mesh.fdaPDE.ref$nodes[,2])
integral.25.ref = integrate_f(FEM(coef.25.ref, FEMbasis.fdaPDE.ref))
coef.ex.ref = coef.25.ref / integral.25.ref

coef.fdaPDE = sols.fdaPDE[ which(mise.fdaPDE[,2] <= min(mise.fdaPDE[,2])), 2, ]
coef.R = sols[which(mise[,2] <= min(mise[,2])), 2,]

coef.fdaPDE.ref = eval.FEM( FEM(coef.fdaPDE, FEMbasis.fdaPDE), locations= mesh.fdaPDE.ref$nodes  )
coef.R.ref = eval.FEM(FEM(coef.R, FEMbasis.fdaPDE), locations= mesh.fdaPDE.ref$nodes) 

R_plot_graph.ggplot2(FEM(coef.ex.ref, FEMbasis.fdaPDE.ref))
R_plot_graph.ggplot2(FEM(coef.fdaPDE.ref,FEMbasis.fdaPDE.ref))
R_plot_graph.ggplot2(FEM(coef.R.ref, FEMbasis.fdaPDE.ref))

dev.off()
