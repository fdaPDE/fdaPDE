################################################################################
######################## increasing nodes number ###############################
################################################################################
library(scales)
library(ggplot2)
library(latex2exp)
# load data # DE-2021-12-17-heat-test-3-RIGHT_CV
date_ = "2022-01-20" # "2021-12-16"
init_ = "" #"-fdaPDE" "-null_vector"
ntest_ = "-test-1-COMPLETE"

file.name = paste("DE-",date_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

load(save.file)
##

data_frame.1 = list()
data_frame.2 = list()
data_frame.3 = list()
K = 4
for( k in 1:K){

mise_ =  cbind(as.vector(mise[[k]]) , as.vector(mise.fdaPDE[[k]]))
mise_= as.vector(mise_)
times_ = cbind( as.vector(times[[k]]),as.vector(times.fdaPDE[[k]]))
times_= as.vector(times_)
nnodes_ = rep(nnodes,each=M)
nnodes_ = rep(as.character(nnodes_),times=2)
type_ = rep(c("R","fdaPDE"),each=(M*N))

# data frame per boxplots #
data_frame.1[[k]] = data.frame(nnodes_,
                          mise_,times_,type_)

mise_ = as.vector(mise[[k]])
mise.fdaPDE_ = as.vector(mise.fdaPDE[[k]])
times_ = as.vector(times[[k]])
times.fdaPDE_ = as.vector(times.fdaPDE[[k]])
nnodes_ = rep(as.character(nnodes), each=M)
err.L2_ = as.vector(err.L2[[k]])

# data frame per scatterplots #
data_frame.2[[k]] = data.frame(nnodes_,
                          mise_,mise.fdaPDE_, err.L2_,
                          times_,times.fdaPDE_)

# data frame per tempi 
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

img.file = paste("/home/aldo/Scrivania/Caricare/","DE-test-1-part-2",".pdf",sep="")

# IMGs #
{
pdf(img.file) #  width = 7, height = 7

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

for( k in 1:K){
  sub.title = paste("n = ",nobs[k],sep="")
 print(ggplot(data_frame.1[[k]])+
        geom_boxplot(aes(x=nnodes_,y=mise_, group=interaction(nnodes_,type_),fill=type_))+
         labs(x="nodes",y= "",
              title=paste("MISE (n=",nobs[k],")",sep=""))+
         scale_x_discrete(limits=as.character(nnodes))+
         MyTheme +
         theme(legend.position= c(0.825,0.85))
)

EXPS = seq(from=min(log2(nnodes)),to=max(log2(nnodes)), length.out=4)
x.1 = 2^EXPS
y.2 = 2^( EXPS - 10)
y.3 = 2^(2* (EXPS -10))
x_ = rep(x.1, times = 2)
y_ = c(y.2, y.2)
NNODES_ = rep(nnodes, times=2)
print(ggplot(data=data_frame.3[[k]], aes(x = nnodes_, y=avg.times_, group=type_, fill=type_ ))+
        #geom_line(aes(x=x_, y=y_) , linetype="dashed", size=0.5, alpha=0.8)+
        geom_line(aes(x=NNODES_, y=NNODES_/2^13) , linetype="dashed", size=0.5, alpha=0.8)+
        geom_line(aes(x=NNODES_, y=NNODES_^2/2^16) , linetype="dashed", size=0.5, alpha=0.8)+
        geom_point(aes(shape=type_,color=type_),size=4)+
        geom_line(aes(color=type_), linetype="dashed",size=2)+
        scale_x_continuous(trans="log2", 
                           breaks=trans_breaks("log2",function(x) 2^x),
                           labels=trans_format("log2", math_format(.x) ) )+
        scale_y_continuous(trans="log2", 
                           breaks=trans_breaks("log2",function(x) 2^x),
                           labels=trans_format("log2", math_format(.x) ) )+ #scales::math_format(2^.x )
        labs(x=TeX("$log_{2}(nodes)$"),
             y=TeX("$log_{2}(times)$"), 
             title=paste("Average running time (n=",nobs[k],")",sep="") )+
        MyTheme + 
        theme( legend.position= c(0.175,0.85) )
        
)
}

dev.off()
}
