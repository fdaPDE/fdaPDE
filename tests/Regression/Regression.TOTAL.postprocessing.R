###########################################################
##            POSTPROCESSING REGRESSION TOTAL            ##
##          TEST 1 - SIMPLENET - NO COVARIATES           ##
##                   delta varying                       ##
##            nobs = c(200, 300, 400, 500)               ##
##      delta = c(0.03, 0.015, 0.0075, 0.00375)          ##
###########################################################

###########################################################
##############      loading data               ############
###########################################################
date_ = "2022-01-19"
ntest_ = "-test-1-no-cov-part2-COMPLETE"
file.name = paste("Regression-",date_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/","Regression-test-1-no-cov-part-2",".pdf",sep="")
load(save.file)
###########################################################

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
# legend.position= c(0.175,0.85) # alto sx
data_frame.1 = list()
data_frame.2 = list()
data_frame.3 = list()
nnodes = c(97, 194, 388, 775)
NNODES_ = nnodes_ = rep(nnodes, times = 2)
### Filling data frames ###
for(k in 1:K){
  mise_ = cbind( as.vector(mise[[k]]), as.vector(mise.fdaPDE[[k]]))
  mise_ = as.vector(mise_)
  times_ = cbind( as.vector(times[[k]]), as.vector(times.fdaPDE[[k]])) 
  times_ = as.vector(times_)
  delta_ = rep(delta,each=M)
  delta_ = rep(as.character(delta_),times=2)
  nodes_ = rep(nnodes,each=M)
  nodes_ = rep(as.character(nodes_), times=2 )
  type_ = rep(c("R","fdaPDE"),each=(M*N))
  
  err.L2_ = as.vector(err.L2[[k]])
  data_frame.1[[k]] = data.frame(mise_, times_, delta_, type_, nodes_)
  
  avg.times = vector(mode="numeric", length=N)
  avg.fdaPDE.times = vector(mode="numeric", length=N)
  for( i in 1:N){
    avg.times[i] = mean(times[[k]][,i])
    avg.fdaPDE.times[i] = mean(times.fdaPDE[[k]][,i])
  }
  
  avg.times_ = cbind(avg.times, avg.fdaPDE.times)
  avg.times_ = as.vector(avg.times_)
  type_ = rep(c("R","fdaPDE"), each=N)
  nnodes_ = rep(nnodes, times = 2)
  
  data_frame.3[[k]] = data.frame(avg.times_, nnodes_, type_)
}

# boxplots of L2, boxplots norm l-2, boxplots of l-inf
err.2 = list()
err.Inf = list()
for(k in 1:K){
  err.2[[k]] = matrix(0,nrow=M,ncol=N)
  err.Inf[[k]] = matrix(0, nrow=M, ncol=N)
  
  for(i in 1:N){
    for( j in 1:M){
      err.2[[k]][j,i] = norm( sols[[k]][[i]][,j] - sols.fdaPDE[[k]][[i]][,j], type="2")
      err.Inf[[k]][j,i] = max(abs( sols[[k]][[i]][,j] - sols.fdaPDE[[k]][[i]][,j]) )
    }
  }
}
for(k in 1:K){
  mise_ = as.vector(mise[[k]])
  mise.fdaPDE_ = as.vector(mise.fdaPDE[[k]])
  err.L2_ = as.vector(err.L2[[k]])
  err.2_ = as.vector(err.2[[k]])
  err.Inf_ = as.vector(err.Inf[[k]])
  #nobs_ = rep(nobs,each=M)
  nnodes_ = rep(as.character(nnodes), each=M)
  
  data_frame.2[[k]] = data.frame(nnodes_,
                                 mise_,mise.fdaPDE_,
                                 err.L2_, err.2, err.Inf_)
  
}
###########################
nnodes = c(97, 194, 388, 775)
# EXPS = seq(from=min(log2(nnodes)), to=max(log2(nnodes)), length.out=N)
# x.1 = 2^EXPS
# y.1 = 2^( EXPS - 10 )
# x_ = rep(x.1, times=2)
# y_ = rep(y.1, times=2)

######### IMGs ###########
img.file="/home/aldo/Scrivania/Caricare/Regression-test-1-no-cov-part-2"
pdf(img.file)
for(k in 1:K){
  
  # boxplots of MISE 
  print(ggplot(data_frame.1[[k]]) +
          geom_boxplot(aes(x=nodes_,y=mise_,group=interaction(delta_,type_),fill=type_))+
          labs(x="nodes",y="",fill="",title=paste("MISE (n=",nobs[k],")",sep=""))+
          scale_x_discrete(limits= as.character(nnodes))+
          theme(legend.position= c(0.825,0.85))+
          MyTheme
  )
  
  # Average running time in log2 log2 scale
  print(ggplot(data=data_frame.3[[k]], aes(x = nnodes_, y=avg.times_, group=type_, fill=type_ ))+
          #geom_line(aes(x=x_, y=y_) , linetype="dashed", size=0.5, alpha=0.8)+
          geom_line(aes(x=NNODES_, y=NNODES_ / 2^14) , linetype="dashed", size=0.5, alpha=0.8)+
          geom_line(aes(x=NNODES_, y=NNODES_^2 / 2^18) , linetype="dashed", size=0.5, alpha=0.8)+
          geom_point(aes(shape=type_,color=type_),size=4)+
          geom_line(aes(color=type_), linetype="dashed",size=2)+
          scale_x_continuous(trans="log2", 
                             breaks=trans_breaks("log2",function(x) 2^x),
                             labels=trans_format("log2", scales::math_format(.x) ) )+
          scale_y_continuous(trans="log2", 
                             breaks=trans_breaks("log2",function(x) 2^x),
                             labels=trans_format("log2", scales::math_format(.x ) ) )+
          labs(x=TeX("$log_{2}(nodes)$"),
               y=TeX("$log_{2}(times)$"), 
               title=paste("Average running time (n=",nobs[k],")",sep="") )+
          theme(legend.position= c(0.175,0.85))+
          MyTheme
  )
  
  ### Boxplots L2 norm - l_2 norm - l_inf norm ###
  print(ggplot(data_frame.2[[k]])+
          geom_boxplot(aes(x=nnodes_,y=err.L2_,group=nnodes_))+
          labs(x="nodes",y="",fill="",
               title= paste("L2 error (n=",nobs[k],")",sep="") ) +
          scale_x_discrete(limits=as.character(nnodes))+
          theme(plot.title = element_text(hjust = 0.5))+
          MyTheme
  )
  
  print(ggplot(data_frame.2[[k]])+
          geom_boxplot(aes(x=nnodes_,y=err.2_,group=nnodes_))+
          labs(x="nodes",y="",fill="",
               title= paste("l2 error (n=",nobs[k],")",sep="") ) +
          scale_x_discrete(limits=as.character(nnodes))+
          theme(plot.title = element_text(hjust = 0.5))+
          MyTheme
  )
  
  print(ggplot(data_frame.2[[k]])+
          geom_boxplot(aes(x=nnodes_,y=err.Inf_,group=nnodes_))+
          labs(x="nodes",y="",fill="",
               title= paste("l Inf error (n=",nobs[k],")",sep="") ) +
          scale_x_discrete(limits=as.character(nnodes))+
          theme(plot.title = element_text(hjust = 0.5))+
          MyTheme
  )
}
dev.off()

##########################

###########################################################
##            POSTPROCESSING REGRESSION TOTAL            ##
##          TEST 2 - BRICKWALL - with COVARIATES         ##
##            nobs = c(300, 500, 700, 1000)              ##
##            delta = c(40, 30, 25, 17.5)                ##
###########################################################

###########################################################
##############      loading data               ############
###########################################################
date_ = "2022-01-20-"
ntest_ = "test-2-part2-COMPLETE-2"
file.name = paste("Regression-",date_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")

img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/","Regression-test-2-part-2",".pdf",sep="")
img.file = "/home/aldo/Scrivania/Caricare/Regression-test-2-part-2.pdf"
load(save.file)
###########################################################

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

data_frame.1 = list()
data_frame.2 = list()
data_frame.3 = list()

for(k in 1:K){
  mise_ = cbind( as.vector(mise[[k]]), as.vector(mise.fdaPDE[[k]]))
  mise_ = as.vector(mise_)
  rmse_ = cbind(as.vector(rmse[[k]]), as.vector(rmse.fdaPDE[[k]]) )
  rmse_ = as.vector(rmse_)
  times_ = cbind( as.vector(times[[k]]), as.vector(times.fdaPDE[[k]])) 
  times_ = as.vector(times_)
  delta_ = rep(delta,each=M)
  delta_ = rep(as.character(delta_),times=2)
  nodes_ = rep(nnodes,each=M)
  nodes_ = rep(as.character(nodes_), times=2 )
  type_ = rep(c("R","fdaPDE"),each=(M*N))
  
  err.L2_ = as.vector(err.L2[[k]])
  data_frame.1[[k]] = data.frame(mise_, rmse_, times_, delta_, type_, nodes_)
  
  avg.times = vector(mode="numeric", length=N)
  avg.fdaPDE.times = vector(mode="numeric", length=N)
  for( i in 1:N){
    avg.times[i] = mean(times[[k]][,i])
    avg.fdaPDE.times[i] = mean(times.fdaPDE[[k]][,i])
  }
  
  avg.times_ = cbind(avg.times, avg.fdaPDE.times)
  avg.times_ = as.vector(avg.times_)
  type_ = rep(c("R","fdaPDE"), each=N)
  nnodes_ = rep(nnodes, times = 2)
  
  data_frame.3[[k]] = data.frame(avg.times_, nnodes_, type_)
  
}

# boxplots of L2, boxplots norm l-2, boxplots of l-inf
err.2 = list()
err.Inf = list()
for(k in 1:K){
  err.2[[k]] = matrix(0,nrow=M,ncol=N)
  err.Inf[[k]] = matrix(0, nrow=M, ncol=N)
  
  for(i in 1:N){
    for( j in 1:M){
      err.2[[k]][j,i] = norm( sols[[k]][[i]][,j] - sols.fdaPDE[[k]][[i]][,j], type="2")
      err.Inf[[k]][j,i] = max(abs( sols[[k]][[i]][,j] - sols.fdaPDE[[k]][[i]][,j]) )
    }
  }
}
for(k in 1:K){
  mise_ = as.vector(mise[[k]])
  mise.fdaPDE_ = as.vector(mise.fdaPDE[[k]])
  err.L2_ = as.vector(err.L2[[k]])
  err.2_ = as.vector(err.2[[k]])
  err.Inf_ = as.vector(err.Inf[[k]])
  #nobs_ = rep(nobs,each=M)
  nnodes_ = rep(as.character(nnodes), each=M)
  
  data_frame.2[[k]] = data.frame(nnodes_,
                                 mise_,mise.fdaPDE_,
                                 err.L2_, err.2, err.Inf_)
  
}
###########################

# 
# EXPS = seq(from=min(log2(nnodes)),to=max(log2(nnodes)), length.out=4)
# x.1 = 2^EXPS
# y.1 = 2^( EXPS - 10 )
# x_ = rep(x.1, times=2)
# y_ = rep(y.1, times=2)
NNODES_ = rep(nnodes, times=2)

######## IMGs ############
pdf(img.file)
for(k in 1:K){
  
  #boxplots of RMSE
  print(ggplot(data_frame.1[[k]]) +
          geom_boxplot(aes(x=nodes_,y=rmse_,group=interaction(delta_,type_),fill=type_))+
          labs(x="nodes",y="",fill="",title=paste("RMSE (n=",nobs[k],")",sep=""))+
          scale_x_discrete(limits= as.character(nnodes))+
          theme(legend.position= c(0.825,0.85))+
          MyTheme
  )
  
  # boxplots of MISE 
  print(ggplot(data_frame.1[[k]]) +
          geom_boxplot(aes(x=nodes_,y=mise_,group=interaction(delta_,type_),fill=type_))+
          labs(x="nodes",y="",fill="",title=paste("MISE (n=",nobs[k],")",sep=""))+
          scale_x_discrete(limits= as.character(nnodes))+
          theme(legend.position= c(0.825,0.85))+
          MyTheme
  )
  
  # Average running time in log2 log2 scale
  print(ggplot(data=data_frame.3[[k]], aes(x = nnodes_, y=avg.times_, group=type_, fill=type_ ))+
          #geom_line(aes(x=x_, y=y_) , linetype="dashed", size=0.5, alpha=0.8)+
          geom_line(aes(x=NNODES_, y=NNODES_ / 2^10) , linetype="dashed", size=0.5, alpha=0.8)+
          geom_line(aes(x=NNODES_, y=NNODES_^2 / 2^16) , linetype="dashed", size=0.5, alpha=0.8)+
          geom_point(aes(shape=type_,color=type_), size=4)+
          geom_line(aes(color=type_), linetype="dashed", size=2)+
          scale_x_continuous(trans="log2", 
                             breaks=trans_breaks("log2",function(x) 2^x),
                             labels=trans_format("log2", scales::math_format(.x) ) )+
          scale_y_continuous(trans="log2", 
                             breaks=trans_breaks("log2",function(x) 2^x),
                             labels=trans_format("log2", scales::math_format(.x ) ) )+
          labs(x=TeX("$log_{2}(nodes)$"),
               y=TeX("$log_{2}(times)$"), 
               title=paste("Average running time (n=",nobs[k],")",sep="") )+
          theme(legend.position= c(0.175,0.85))+
          MyTheme
  )
  
  ### Boxplots L2 norm - l_2 norm - l_inf norm ###
  print(ggplot(data_frame.2[[k]])+
          geom_boxplot(aes(x=nnodes_,y=err.L2_,group=nnodes_))+
          labs(x="nodes",y="",fill="",
               title= paste("L2 error (n=",nobs[k],")",sep="") ) +
          scale_x_discrete(limits=as.character(nnodes))+
          theme(plot.title = element_text(hjust = 0.5))+
          MyTheme
  )
  
  print(ggplot(data_frame.2[[k]])+
          geom_boxplot(aes(x=nnodes_,y=err.2_,group=nnodes_))+
          labs(x="nodes",y="",fill="",
               title= paste("l2 error (n=",nobs[k],")",sep="") ) +
          scale_x_discrete(limits=as.character(nnodes))+
          theme(plot.title = element_text(hjust = 0.5))+
          MyTheme
  )
  
  print(ggplot(data_frame.2[[k]])+
          geom_boxplot(aes(x=nnodes_,y=err.Inf_,group=nnodes_))+
          labs(x="nodes",y="",fill="",
               title= paste("l Inf error (n=",nobs[k],")",sep="") ) +
          scale_x_discrete(limits=as.character(nnodes))+
          theme(plot.title = element_text(hjust = 0.5))+
          MyTheme
  )
}
dev.off()
##########################