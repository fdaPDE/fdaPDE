# NB questo script solo per post processing di confronto ! test.1.3 e test.1.4

library(ggplot2)
library(latex2exp)


## loading data ##
date_ = ""
ntest_ = "test-1" # "test-2"

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
  labs(x="",y="",fill="",title=TeX( "$RMSE(\\beta_{1})^{2}$") )+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=rmse.beta.2_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="",y="",fill="",title=TeX( "$RMSE( \\beta_{1})^{2}$") )+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data_frame.1) +
  geom_boxplot(aes(x=nobs_,y=mise_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title=TeX( "MISE"))+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data_frame.1) + 
  geom_boxplot(aes(x=nobs_,y=times_,group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="",fill="",title="TIME [s]") +
  theme(plot.title = element_text(hjust = 0.5))


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
  theme(plot.title = element_text(hjust = 0.5))


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
                 plot.subtitle = element_text(hjust=0.5))  )
  
  x.min = min(data_frame.tmp["rmse.beta.2_"], data_frame.tmp["rmse.fdaPDE.beta.2_"])
  x.max = max(data_frame.tmp["rmse.beta.2_"], data_frame.tmp["rmse.fdaPDE.beta.2_"])
  
  print(  ggplot(data_frame.tmp, aes(x=rmse.beta.2_,y=rmse.fdaPDE.beta.2_)) +
            geom_point() +
            geom_segment(aes(x=x.min, y=x.min, xend=x.max, yend=x.max),color="red",alpha=0.01) +
            labs(x="R", y="fdaPDE", 
                 title= TeX("Scatter  $RMSE(\\beta_{2})^{2}$"),
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
