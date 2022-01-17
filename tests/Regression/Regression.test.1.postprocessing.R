# NB questo script solo per post processing di confronto ! test.1.3 e test.1.4

library(ggplot2)
library(latex2exp)


## loading data ##
date_ = "2022-01-11-"
ntest_ = "test-1" # "test-2"

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
###

## post - processing ## 
rmse.beta.1_ = cbind( as.vector(rmse.beta.1), as.vector(rmse.fdaPDE.beta.1)) 
rmse.beta.1_ = as.vector(rmse.beta.1_)
rmse.beta.2_ = cbind( as.vector(rmse.beta.2), as.vector(rmse.fdaPDE.beta.2)) 
rmse.beta.2_ = as.vector(rmse.beta.2_)
mise_ = cbind( as.vector(mise), as.vector(mise.fdaPDE))
mise_ = as.vector(mise_)
times_ = cbind( as.vector(times), as.vector(times.fdaPDE)) 
times_ = as.vector(times_)
nobs_ = rep(nobs,each=M)
nobs_ = rep(as.character(nobs_), times=2)
type_ = rep(c("R","fdaPDE"),each=(M*N))

err.L2_ = as.vector(err.L2)

# data_frame boxplots #
data_frame.1 = data.frame(rmse.beta.1_, rmse.beta.2_,
                          mise_,times_,nobs_,type_)

pdf(img.file)

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=rmse.beta.1_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title=TeX( "$RMSE(\\beta_{1})^{2}$") )+
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=rmse.beta.2_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title=TeX( "$RMSE( \\beta_{2})^{2}$") )+
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=mise_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title=TeX( "MISE"))+
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme

ggplot(data_frame.1) + 
  geom_boxplot(aes(x=nobs_,y=times_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title="TIME [s]") +
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme


rmse.beta.1_ = as.vector(rmse.beta.1)
rmse.beta.2_ = as.vector(rmse.beta.2)
rmse.fdaPDE.beta.1_ = as.vector(rmse.fdaPDE.beta.1)
rmse.fdaPDE.beta.2_ = as.vector(rmse.fdaPDE.beta.2)
mise_ = as.vector(mise)
mise.fdaPDE_ = as.vector(mise.fdaPDE)
times_ = as.vector(times)
times.fdaPDE_ = as.vector(times.fdaPDE)
err.L2_ = as.vector(err.L2)
nobs_ = rep(nobs,each=M)
# data frame per scatterplots #
data_frame.2 = data.frame(nobs_,
                          rmse.beta.1_, rmse.beta.2_,
                          rmse.fdaPDE.beta.1_, rmse.fdaPDE.beta.2_,
                          mise_,mise.fdaPDE_,err.L2_,
                          times_,times.fdaPDE_)

ggplot(data_frame.2)+
  geom_boxplot(aes(x=nobs_,y=err.L2_,group=nobs_))+
  labs(x="observations",y="",fill="",
       title= TeX("Scatter \ $||f_{R} - f_{fdaPDE}||_{L_{2}}^2$") ) +
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme


for(i in 1:N){
  data_frame.tmp = subset(data_frame.2, nobs_ %in% c(nobs[i]))
  
  x.min = min(data_frame.tmp["rmse.beta.1_"], data_frame.tmp["rmse.fdaPDE.beta.1_"])
  x.max = max(data_frame.tmp["rmse.beta.1_"], data_frame.tmp["rmse.fdaPDE.beta.1_"])
  
  print( ggplot(data_frame.tmp, aes(x=rmse.beta.1_,y=rmse.fdaPDE.beta.1_)) +
           geom_point() +
           geom_segment(aes(x=x.min, y=x.min, xend=x.max, yend=x.max),color="red",alpha=0.01)+  
           labs(x="R", y="fdaPDE", 
                title= TeX("Scatter  $RMSE(\\beta_{1})^{2}$"),
                subtitle= paste("nobs = ",nobs[i],sep=""))+
           theme(plot.title = element_text(hjust=0.5),
                 plot.subtitle = element_text(hjust=0.5))+
           MyTheme )
    
  
  x.min = min(data_frame.tmp["rmse.beta.2_"], data_frame.tmp["rmse.fdaPDE.beta.2_"])
  x.max = max(data_frame.tmp["rmse.beta.2_"], data_frame.tmp["rmse.fdaPDE.beta.2_"])
  
  print(  ggplot(data_frame.tmp, aes(x=rmse.beta.2_,y=rmse.fdaPDE.beta.2_)) +
            geom_point() +
            geom_segment(aes(x=x.min, y=x.min, xend=x.max, yend=x.max),color="red",alpha=0.01) +
            labs(x="R", y="fdaPDE", 
                 title= TeX("Scatter  $RMSE(\\beta_{2})^{2}$"),
                 subtitle= paste("nobs = ",nobs[i],sep=""))+
            theme(plot.title = element_text(hjust=0.5),
                  plot.subtitle = element_text(hjust=0.5))+
            MyTheme)
  
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


dev.off()


## no ggplot2 ## 
img.file.no.ggplot2 = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,"-no-ggplot2.pdf",sep="")
pdf(img.file.no.ggplot2)
for( i in 1:N){
  title.norm.2 = paste("Norm l^2, nobs =",as.character( nobs[i]),sep=" ")
  title.norm.inf = paste("Norm l^inf, nobs =",as.character( nobs[i]),sep=" ")
  #x11()
  plot(1:M,rmse.beta.2[,i],pch=16,col="red",main=paste("RMSE"," beta1, nobs =",as.character( nobs[i]),sep=" "),xlab="",ylab="")
  points(1:M,rmse.fdaPDE.beta.2[,i],pch=13,col="blue")
  segments(1:M,rmse.beta.2[,i],1:M,rmse.fdaPDE.beta.2[,i], lty=3)
  legend("topright", legend=c("Mattina","fdaPDE"),pch=c(16,13),col=c("red","blue"))
  
  plot(1:M, norms.beta1.inf[,i], pch=16, col="red", main=paste("Beta1",title.norm.inf,sep=" "),xlab="",ylab="")
  plot(1:M, norms.beta1.2[,i],pch = 16, col="red", main =paste("Beta1",title.norm.2,sep=" "),xlab="",ylab="")
  
  plot(1:M,rmse.beta.2[,i],pch=16,col="red",main=paste("RMSE"," beta2, nobs =",as.character( nobs[i]),sep=" "),xlab="",ylab="")
  points(1:M,rmse.fdaPDE.beta.2[,i],pch=13,col="blue")
  segments(1:M,rmse.beta.2[,i],1:M,rmse.fdaPDE.beta.2[,i], lty=3)
  legend("topright", legend=c("Mattina","fdaPDE"),pch=c(16,13),col=c("red","blue"))
  plot(1:M, norms.beta2.inf[,i], pch=16, col="red", main=paste("Beta2",title.norm.inf,sep=" "),xlab="",ylab="")
  plot(1:M, norms.beta2.2[,i],   pch=16, col="red", main =paste("Beta2",title.norm.2,sep=" "),xlab="",ylab="")
  
}
boxplot(rmse.beta.1[,1],rmse.fdaPDE.beta.1[,1],
        rmse.beta.1[,2],rmse.fdaPDE.beta.1[,2],
        rmse.beta.1[,3],rmse.fdaPDE.beta.1[,3],
        rmse.beta.1[,4],rmse.fdaPDE.beta.1[,4],
        #main="beta1",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11),
        rmse.beta.1[,5],rmse.fdaPDE.beta.1[,5],
        rmse.beta.1[,6],rmse.fdaPDE.beta.1[,6],
        rmse.beta.1[,7],rmse.fdaPDE.beta.1[,7],
        main="beta1",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20),
        col=c("red","blue"))
legend("topright", legend=c("Mattina","fdaPDE"), col=c("red","blue"),pch=c(15,15) )

boxplot(rmse.beta.2[,1],rmse.fdaPDE.beta.2[,1],
        rmse.beta.2[,2],rmse.fdaPDE.beta.2[,2],
        rmse.beta.2[,3],rmse.fdaPDE.beta.2[,3],
        rmse.beta.2[,4],rmse.fdaPDE.beta.2[,4],
        #main="beta2",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11),
        rmse.beta.2[,5],rmse.fdaPDE.beta.2[,5],
        rmse.beta.2[,6],rmse.fdaPDE.beta.2[,6],
        rmse.beta.2[,7],rmse.fdaPDE.beta.2[,7],
        main="beta2",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20),
        col=c("red","blue"))
legend("topright", legend=c("Mattina","fdaPDE"), col=c("red","blue"),pch=c(15,15) )

boxplot(mise[,1],mise.fdaPDE[,1],
        mise[,2],mise.fdaPDE[,2],
        mise[,3],mise.fdaPDE[,3],
        mise[,4],mise.fdaPDE[,4],
        #main="beta1",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11),
        mise[,5],mise.fdaPDE[,5],
        mise[,6],mise.fdaPDE[,6],
        mise[,7],mise.fdaPDE[,7],
        main="MISE",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20),
        col=c("red","blue"))
legend("topright", legend=c("Mattina","fdaPDE"), col=c("red","blue"),pch=c(15,15) )

boxplot(err.L2[,1],
        err.L2[,2],
        err.L2[,3],
        err.L2[,4],
        #main="beta1",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11),
        err.L2[,5],
        err.L2[,6],
        err.L2[,7],
        main="L2 norm",names=nobs,at=c(1,2,3,4,5,6,7))



time.plot = matrix(0,nrow=1,ncol=N)
time.fdaPDE.plot = matrix(0,nrow=1,ncol=N)

for(i in 1:N){
  time.plot[i] = mean(times[,i])
  time.fdaPDE.plot[i] = mean(times.fdaPDE[,i])
}

#x11()
plot(nobs, time.plot,type="b",col="red",pch=15 ,xlab = "nobs",ylab = "times [s]", main="average time." ,ylim=c(min(time.fdaPDE.plot),max(time.plot)))
points(nobs, time.fdaPDE.plot,type="b",pch = 16, col="blue")
legend("topleft",legend=c("eval+algo Mattina", "fdaPDE"), 
       col=c("red","blue"),pch=c(15,16),cex = 0.8,text.font=4 )


#x11()
boxplot(times[,1],times.fdaPDE[,1],
        times[,2],times.fdaPDE[,2],
        times[,3],times.fdaPDE[,3],
        times[,4],times.fdaPDE[,4],
        #main="time",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11),
        times[,5],times.fdaPDE[,5],
        times[,6],times.fdaPDE[,6],
        times[,7],times.fdaPDE[,7],
        main="times",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20),
        col=c("red","blue"))
legend("topright", legend=c("eval+algo Mattina","fdaPDE"), col=c("red","blue"),pch=c(15,15) )

dev.off()

########################################################################
##### Regression PostProcessing Test-1 delta-var #######################
########################################################################
date_ = "2022-01-11-"
ntest_ = "test-1-delta-var" # "test-2"

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

rmse.beta.1_ = cbind( as.vector(rmse.beta.1), as.vector(rmse.fdaPDE.beta.1)) 
rmse.beta.1_ = as.vector(rmse.beta.1_)
rmse.beta.2_ = cbind( as.vector(rmse.beta.2), as.vector(rmse.fdaPDE.beta.2)) 
rmse.beta.2_ = as.vector(rmse.beta.2_)
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
data_frame.1 = data.frame(rmse.beta.1_, rmse.beta.2_,
                          mise_,times_,delta_,type_, nodes_)

data_frame.1 = data.frame(
                          mise_,times_,delta_,type_, nodes_)

pdf(img.file)

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nodes_,y=rmse.beta.1_,group=interaction(delta_,type_),fill=type_))+
  labs(x="nodes",y="",fill="",title=TeX( "$RMSE(\\beta_{1})^{2}$") )+
  scale_x_discrete(limits=c(as.character(nnodes)))+
  theme(plot.title = element_text(hjust = 0.5)) +
  MyTheme

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nodes_,y=rmse.beta.2_,group=interaction(delta_,type_),fill=type_))+
  labs(x="nodes",y="",fill="",title=TeX( "$RMSE( \\beta_{2})^{2}$") )+
  scale_x_discrete(limits= as.character(nnodes))+
  theme(plot.title = element_text(hjust = 0.5))+
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

avg.mise = vector(mode="numeric", length=N)
avg.fdaPDE.mise = vector(mode="numeric", length=N)
for(i in 1:N){
  avg.mise[i] = mean(mise[,i])
  avg.fdaPDE.mise[i] = mean(mise.fdaPDE[,i])
}


dev.off()

plot(nnodes, avg.times, log="y", col="red", pch=16,
     ylim=c(min(min(avg.times),min(avg.fdaPDE.times)),max(max(avg.times),max(avg.fdaPDE.times)))) 
points(nnodes, avg.fdaPDE.times, log="y", col="blue", pch=16)

y.min = min(min(log2(avg.times)),min(log2(avg.fdaPDE.times)))
y.max = max(max(log2(avg.times)),max(log2(avg.fdaPDE.times)))
plot(nnodes, log2(avg.times), col="red",
     ylim=c(y.min,y.max) ) 
points(nnodes, log2(avg.fdaPDE.times), col="blue")


avg.times
avg.fdaPDE.times
avg.times_ = cbind(log2(avg.times), log2(avg.fdaPDE.times))
avg.times_ = as.vector(avg.times_)
type_ = rep(c("R","fdaPDE"), each=N)
nnodes_ = rep(as.character(nnodes), times = 2)

data.frame.3 = data.frame(avg.times_, nnodes_, type_)

ggplot(data=data.frame.3, aes(x = nnodes_, y=avg.times_, group=type_, fill=type_ ))+
  geom_point(aes(x = nnodes_, y=avg.times_, group=type_, fill=type_ ,shape=type_,color=type_))+
  geom_line(aes(color=type_), linetype="dashed")+
  scale_x_discrete(limits=as.character(log2(nnodes)), 
                   labels = paste(format(log2(nnodes),nsmall=4), "..", sep="")) + 
  labs(x=TeX("$log_{2}(nodes)$"),y=TeX("$log_{2}(times)$"), title="Average running time")+
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())

##########################################################
################# NO COVARIATES ##########################
##########################################################
library(ggplot2)
library(latex2exp)

date_ = "2022-01-12-"
ntest_ = "test-1-no-cov" 

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
###

## post - processing ## 
mise_ = cbind( as.vector(mise), as.vector(mise.fdaPDE))
mise_ = as.vector(mise_)
times_ = cbind( as.vector(times), as.vector(times.fdaPDE)) 
times_ = as.vector(times_)
nobs_ = rep(nobs,each=M)
nobs_ = rep(as.character(nobs_), times=2)
type_ = rep(c("R","fdaPDE"),each=(M*N))

err.L2_ = as.vector(err.L2)

# data_frame boxplots #
data_frame.1 = data.frame(mise_,times_,nobs_,type_)

pdf(img.file)

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=mise_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title=TeX( "MISE"))+
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme

ggplot(data_frame.1) + 
  geom_boxplot(aes(x=nobs_,y=times_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title="TIME [s]") +
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme

mise_ = as.vector(mise)
mise.fdaPDE_ = as.vector(mise.fdaPDE)
times_ = as.vector(times)
times.fdaPDE_ = as.vector(times.fdaPDE)
err.L2_ = as.vector(err.L2)
nobs_ = rep(nobs,each=M)
# data frame per scatterplots #
data_frame.2 = data.frame(nobs_,
                          mise_,mise.fdaPDE_,err.L2_,
                          times_,times.fdaPDE_)

ggplot(data_frame.2)+
  geom_boxplot(aes(x=nobs_,y=err.L2_,group=nobs_))+
  labs(x="observations",y="",fill="",
       title= TeX("Scatter \ $||f_{R} - f_{fdaPDE}||_{L_{2}}^2$") ) +
  theme(plot.title = element_text(hjust = 0.5))+
  MyTheme


for(i in 1:N){
  data_frame.tmp = subset(data_frame.2, nobs_ %in% c(nobs[i]))
  
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


dev.off()
