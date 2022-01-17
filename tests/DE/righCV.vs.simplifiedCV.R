# RightCV vs SimplifiedCV #
# RightCV Gradient descent with fixed step vs BFGS with fixed step #

# DE test 2 #

source("~/Scrivania/fdaPDE/Auxiliary/integrate_f.R")
library(spatstat)
data("dendrite")

vertices = cbind(dendrite$domain$vertices$x, dendrite$domain$vertices$y) 
edges = cbind(dendrite$domain$from, dendrite$domain$to)
L = as.linnet(dendrite$domain)
M = dendrite$domain$m

mesh = create.mesh.1.5D(nodes=vertices,edges=edges)
FEMbasis = create.FEM.basis(mesh=mesh)

plot(mesh)
points(mesh$nodes[600,1],mesh$nodes[600,2],pch=16)
points(mesh$nodes[70,1],mesh$nodes[70,2],pch=16)
points(mesh$nodes[10,1],mesh$nodes[10,2],pch=16)
# density
mu1 = mesh$nodes[70,]
mu2 = mesh$nodes[600,]
mu3 = mesh$nodes[10,]
sigmax = 20
sigmay = 20
aux = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], x3=mu3[1], y3=mu3[2], sigma_x=sigmax, sigma_y=sigmay) 
{   ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))+ 
    ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x3)^2 + (1/(2*sigma_y^2))*(y-y3)^2 ))}


integral_ex = integrate_f( FEM(aux(mesh$nodes[,1],mesh$nodes[,2]),FEMbasis) )

aux = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], x3=mu3[1], y3=mu3[2], sigma_x=sigmax, sigma_y=sigmay) 
{   1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))+
    1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x3)^2 + (1/(2*sigma_y^2))*(y-y3)^2 ))}
my_dens <- linfun(aux, L) # FUNZIONE SU UN NETWORK

coef.ex = aux(mesh$nodes[,1],mesh$nodes[,2])
plot(my_dens)

niter = 25
nfolds=5
 
nnodes = dim(mesh$nodes)[1]
###
nobs = c(100,150,200,250,500)
#nobs = c(150)
N = length(nobs)
#cvec=matrix(0,nrow=mesh$nnodes,ncol=1)
# 
lambda = c(1e-4,1e-3,1e-2,1e-1,1) 

M=30
mise.RightCV.Grad = matrix(0,nrow=M,ncol=N)
mise.RightCV.BFGS = matrix(0,nrow=M,ncol=N)
mise.SimplifiedCV.Grad = matrix(0,nrow=M,ncol=N) # RightCV vs Smplified con Gradient (?)
mise.SimplifiedCV.BFGS = matrix(0,nrow=M,ncol=N)

err.L2.RightCV.Grad.vs.SimplifiedCV.Grad= matrix(0,nrow=M,ncol=N)
err.L2.RightCV.BFGS.vs.SimplifiedCV.BFGS = matrix(0,nrow=M,ncol=N)

times.RightCV.Grad = matrix(0,nrow=M, ncol=N)
times.RightCV.BFGS = matrix(0,nrow=M,ncol=N)
times.SimplifiedCV.Grad = matrix(0,nrow=M,ncol=N)
times.SimplifiedCV.BFGS = matrix(0,nrow=M,ncol=N)

sols.RightCV.Grad = array(0,dim=c(M,N,nnodes) ) #contiene vettore dei coeff
sols.RightCV.BFGS = array(0,dim=c(M,N,nnodes) )
sols.SimplifiedCV.Grad = array(0,dim=c(M,N,nnodes) )
sols.SimplifiedCV.BFGS = array(0,dim=c(M,N,nnodes) )

lambda.RightCV.Grad = matrix(0,nrow=M,ncol=N)
lambda.RightCV.BFGS = matrix(0,nrow=M,ncol=N)
lambda.SimplifiedCV.Grad = matrix(0,nrow=M,ncol=N)
lambda.SimplifiedCV.BFGS = matrix(0,nrow=M,ncol=N)

tot.start = Sys.time()
for(i in 1:N){
  for(j in 1:M){
    print(paste("##### nobs = ",nobs[i]," #####",sep=""))
    print(paste("### ",j ,"/M ###",sep=""))
    PPP = rlpp(nobs[i], my_dens)
    points_x = PPP$data$x
    points_y = PPP$data$y
    points = NULL 
    points = cbind(points_x, points_y)
    points = as.matrix(points)
    colnames(points) = NULL
    
    # RightCV Gradient
    print("# RightCV - Gradient #")
    start <- Sys.time()
    sol = DE.FEM(data=points,FEMbasis=FEMbasis,lambda=lambda, 
                        nsimulations=niter, 
                        #fvec = exp(cvec),
                        tol1=1e-5,
                        nfolds=nfolds, 
                        step_method="Fixed_Step", 
                        direction_method="Gradient",
                        preprocess_method = "RightCV")
    end <- Sys.time()
    coef = exp(sol$g)
    sols.RightCV.Grad[j,i,] = coef 
    mise.RightCV.Grad[j,i]   = integrate_f( FEM( (coef - coef.ex)^2, FEMbasis)) 
    times.RightCV.Grad[j,i] = difftime(end, start, units="secs")
    lambda.RightCV.Grad[j,i] = sol$lambda
    
    # RightCV BFGS
    print("# RightCV - BFGS#")
    start <- Sys.time()
    sol = DE.FEM(data=points,FEMbasis=FEMbasis,lambda=lambda, 
                 nsimulations=niter, 
                 #fvec = exp(cvec),
                 tol1 = 1e-5,
                 nfolds=nfolds, 
                 step_method="Fixed_Step", 
                 direction_method="BFGS",
                 preprocess_method = "RightCV")
    end <- Sys.time()
    coef = exp(sol$g)
    sols.RightCV.BFGS[j,i,] = coef 
    mise.RightCV.BFGS[j,i]   = integrate_f( FEM( (coef - coef.ex)^2, FEMbasis)) 
    times.RightCV.BFGS[j,i] = difftime(end, start, units="secs")
    lambda.RightCV.BFGS[j,i] = sol$lambda
    
    # SimplifiedCV Gradient
    print("# SimplifiedCV - Gradient #")
    start <- Sys.time()
    sol = DE.FEM(data=points,FEMbasis=FEMbasis,lambda=lambda, 
                 nsimulations=niter, 
                 #fvec = exp(cvec),
                 tol1 = 1e-5,
                 nfolds=nfolds, 
                 step_method="Fixed_Step", 
                 direction_method="Gradient",
                 preprocess_method = "SimplifiedCV")
    end <- Sys.time()
    coef = exp(sol$g)
    sols.SimplifiedCV.Grad[j,i,] = coef 
    mise.SimplifiedCV.Grad[j,i]   = integrate_f( FEM( (coef - coef.ex)^2, FEMbasis)) 
    times.SimplifiedCV.Grad[j,i] = difftime(end, start, units="secs")
    lambda.SimplifiedCV.Grad[j,i] = sol$lambda
    
    # SimplifiedCV BFGS
    print("# SimplifiedCV - BFGS #")
    start <- Sys.time()
    sol = DE.FEM(data=points,FEMbasis=FEMbasis,lambda=lambda, 
                 nsimulations=niter, 
                 #fvec = exp(cvec),
                 tol1 = 1e-5,
                 nfolds=nfolds, 
                 step_method="Fixed_Step", 
                 direction_method="BFGS",
                 preprocess_method = "SimplifiedCV")
    end <- Sys.time()
    coef = exp(sol$g)
    sols.SimplifiedCV.BFGS[j,i,] = coef 
    mise.SimplifiedCV.BFGS[j,i]   = integrate_f( FEM( (coef - coef.ex)^2, FEMbasis)) 
    times.SimplifiedCV.BFGS[j,i] = difftime(end, start, units="secs")
    lambda.SimplifiedCV.BFGS[j,i] = sol$lambda
    
    
    err.L2.RightCV.Grad.vs.SimplifiedCV.Grad[j,i] = 
            integrate_f( FEM( (sols.RightCV.Grad[j,i,] - sols.SimplifiedCV.Grad[j,i,])^2, FEMbasis))
                    
    err.L2.RightCV.BFGS.vs.SimplifiedCV.BFGS[j,i] = 
            integrate_f( FEM( (sols.RightCV.BFGS[j,i,] - sols.SimplifiedCV.BFGS[j,i,])^2, FEMbasis))
  
  }
}
tot.end = Sys.time()
tot.time = difftime(tot.end,tot.start, units="secs")


today_ = Sys.Date()
init_ = "-fdaPDE" # "-null_vector"
ntest_ = "-test.2.3"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")


save(nobs, N, M, coef.ex,niter,
     mise.RightCV.Grad, mise.RightCV.BFGS,
     mise.SimplifiedCV.Grad, mise.SimplifiedCV.BFGS,
     err.L2.RightCV.Grad.vs.SimplifiedCV.Grad,
     err.L2.RightCV.BFGS.vs.SimplifiedCV.BFGS,
     times.RightCV.Grad, times.RightCV.BFGS,
     times.SimplifiedCV.Grad, times.SimplifiedCV.BFGS,
     sols.RightCV.Grad, sols.RightCV.BFGS,
     sols.SimplifiedCV.Grad, sols.SimplifiedCV.BFGS,
     lambda.RightCV.Grad, lambda.RightCV.BFGS,
     lambda.SimplifiedCV.Grad, lambda.SimplifiedCV.BFGS,
     file = save.file)

# post-processing #
library(ggplot2)
library(latex2exp)

date_ = "2021-12-08"
init_ = "-fdaPDE" # "-null_vector"
ntest_ = "-test.2.3"

file.name = paste("DE-",date_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

load(save.file)

mise.RightCV.Grad_ = as.vector(mise.RightCV.Grad)
mise.RightCV.BFGS_ = as.vector(mise.RightCV.BFGS)
mise.SimplifiedCV.Grad_ = as.vector(mise.SimplifiedCV.Grad)
mise.SimplifiedCV.BFGS_ = as.vector(mise.SimplifiedCV.BFGS)
mise_ = cbind(mise.RightCV.Grad_,
              mise.RightCV.BFGS_,
              mise.SimplifiedCV.Grad_,
              mise.SimplifiedCV.BFGS_)
mise_ = as.vector(mise_)

times.RightCV.Grad_ = as.vector(times.RightCV.Grad)
times.RightCV.BFGS_ = as.vector(times.RightCV.BFGS)
times.SimplifiedCV.Grad_ = as.vector(times.SimplifiedCV.Grad)
times.SimplifiedCV.BFGS_ = as.vector(times.SimplifiedCV.BFGS)
times_ = cbind(times.RightCV.Grad_,
               times.RightCV.BFGS_,
               times.SimplifiedCV.Grad_,
               times.SimplifiedCV.BFGS_)
times_ = as.vector(times_)

err.L2.RightCV.Grad.vs.SimplifiedCV.Grad_  =
                as.vector(err.L2.RightCV.Grad.vs.SimplifiedCV.Grad)
err.L2.RightCV.BFGS.vs.SimplifiedCV.BFGS_ =
                as.vector(err.L2.RightCV.BFGS.vs.SimplifiedCV.BFGS)

nobs_ = rep(nobs,  each=M)
nobs_ = rep(as.character(nobs_), times=4)

CV_.1 = rep(c("RightCV"),times=2*M*N)
CV_.2 = rep(c("SimplifiedCV"), times=2*M*N)
CV_ = cbind( CV_.1, CV_.2)
CV_ = as.vector( CV_)

Direction_ = rep(c("GD","BFGS"), each = M*N)
Direction_ = rep(Direction_, times=2)

data.frame.1 = data.frame(mise_, times_, nobs_, CV_, Direction_)

pdf(img.file)
ggplot( data.frame.1)+
  geom_boxplot(aes(x=nobs_, y=mise_, group=interaction(nobs_, CV_, Direction_), fill=interaction(CV_ ,Direction_)) )+
  labs(x="observations",y= "MISE",fill="",title="f")+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1) )

sub.data.frame.1 = subset(data.frame.1, Direction_ %in% c("GD"))

ggplot(sub.data.frame.1)+
  geom_boxplot(aes(x=nobs_, y=mise_, group=interaction(nobs_, CV_), fill=CV_ ))+
  labs(x="observations",y= "MISE",fill="",title="f")+
  theme(plot.title = element_text(hjust = 0.5))

ggplot( sub.data.frame.1)+
  geom_boxplot(aes(x=nobs_, y=times_, group=interaction(nobs_, CV_), fill=CV_))+
  labs(x="observations",y= "",fill="",title="TIME [s]")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1) )

for(i in 1:N){
  data.frame.tmp = subset(data.frame.1, nobs_ %in% c(nobs[i]))
  print(ggplot(data.frame.tmp)+
        geom_boxplot(aes(x=interaction(CV_, Direction_), y=times_ ))+
        labs(x="",y= "",fill="",title="Running Time [s]", subtitle=paste("nobs = ",nobs[i],sep=""))+
        theme(plot.title = element_text(hjust = 0.5), 
              plot.subtitle = element_text(hjust=0.5),
              axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1) ) )
  
    
}

nobs_ = rep(nobs, each=M)

data.frame.2 = data.frame(nobs_,
                    mise.RightCV.Grad_, mise.SimplifiedCV.Grad_,
                    mise.RightCV.BFGS_, mise.SimplifiedCV.BFGS_,
                    times.RightCV.Grad_, times.SimplifiedCV.Grad_,
                    times.RightCV.BFGS_, times.SimplifiedCV.BFGS_,
                    err.L2.RightCV.Grad.vs.SimplifiedCV.Grad_, 
                    err.L2.RightCV.BFGS.vs.SimplifiedCV.BFGS_)
ggplot(data.frame.2)+
  geom_boxplot(aes(x=nobs_,y=err.L2.RightCV.Grad.vs.SimplifiedCV.Grad_, group=nobs_))+
  labs(x="observations", y="", 
       title= TeX(" $||f_{RightCV} \ - \ f_{SimplCV} ||_{L_{2}^{2}$"),
       subtitle="GD")+
  theme(  plot.title = element_text(hjust=0.5), 
          plot.subtitle = element_text(hjust=0.5) ) 

ggplot(data.frame.2)+
  geom_boxplot(aes(x=nobs_,y=err.L2.RightCV.BFGS.vs.SimplifiedCV.BFGS_, group=nobs_))+
  labs(x="observations", y="", 
       title= TeX(" $||f_{RightCV} \ - \ f_{SimplCV} ||_{L_{2}^{2}$"),
       subtitle="BFGS")+
  theme(  plot.title = element_text(hjust=0.5), 
          plot.subtitle = element_text(hjust=0.5) ) 

for(i in 1:N){
  data.frame.tmp = subset(data.frame.2, nobs_ %in% c(nobs[i]))
  
  x.min = min( data.frame.tmp["mise.RightCV.Grad_"], data.frame.tmp["mise.SimplifiedCV.Grad_"])
  x.max = max( data.frame.tmp["mise.RightCV.Grad_"], data.frame.tmp["mise.SimplifiedCV.Grad_"])
  print( ggplot(data.frame.tmp)+
        geom_point( aes(x=mise.RightCV.Grad_, y= mise.SimplifiedCV.Grad_))+
        geom_segment(aes(x=x.min, y=x.min, xend=x.max, yend=x.max),color="red",alpha=0.01)+
        labs(x="RightCV", y="Simplified", 
        title= TeX("Scatter $MISE$ - GD"), 
        subtitle= paste("nobs = ",nobs[i],sep="") ) +
        theme(plot.title = element_text(hjust=0.5), 
             plot.subtitle = element_text(hjust=0.5) )    )
  
  x.min = min( data.frame.tmp["mise.RightCV.BFGS_"], data.frame.tmp["mise.SimplifiedCV.BFGS_"])
  x.max = max( data.frame.tmp["mise.RightCV.BFGS_"], data.frame.tmp["mise.SimplifiedCV.BFGS_"])
  print( ggplot(data.frame.tmp)+
           geom_point( aes(x=mise.RightCV.BFGS_, y= mise.SimplifiedCV.BFGS_))+
           geom_segment(aes(x=x.min, y=x.min, xend=x.max, yend=x.max),color="red",alpha=0.01)+
           labs(x="RightCV", y="Simplified", 
                title= TeX("Scatter $MISE$ - BFGS"), 
                subtitle= paste("nobs = ",nobs[i],sep="") ) +
           theme(plot.title = element_text(hjust=0.5), 
                 plot.subtitle = element_text(hjust=0.5) )    )
  
}


reps_ = rep(1:M, times = 2)
type_ = rep(c("RightCV","SimplifiedCV"), each=M)
for(i in 1:N){
  lambda.GD_ = cbind(log10(lambda.RightCV.Grad[,i]), log10(lambda.SimplifiedCV.Grad[,i]))
  lambda.GD_ = as.vector(lambda.GD_)
  data.frame.3 = data.frame(lambda.GD_, reps_, type_)
  print(
    ggplot(data.frame.3, aes(x=reps_,y=lambda.GD_, shape=type_, color=type_))+
    geom_point(size=2.5)+
    scale_shape_manual(values=c(16,13))+
    scale_color_manual(values=c("red","blue"))+
    labs(x="repetions",y=TeX("$ log_{10}\\lambda $"),
       title="Smoothing parameter",
       subtitle = paste("nobs = ",nobs[i],sep="") )+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5), 
        plot.subtitle = element_text(hjust=0.5) )  
)
}

boxplot(times.RightCV.Grad[,5],times.SimplifiedCV.Grad[,5],
        times.RightCV.BFGS[,5],times.SimplifiedCV.BFGS[,5],
        names=c("RighCV Grad","SimplfiedCV Grad","RightCV BFGS","SimplifiedCV BFGS"),
        main = "Running Time [s]")


### densities ### 
source("~/Scrivania/fdaPDE/tests/R_plot_graph.ggplot2.R")

idxs.min = which( abs(mise.RightCV.BFGS[,5] - min(mise.RightCV.BFGS[,5])) <= .Machine$double.eps)

R_plot_graph.ggplot2(FEM( sols.RightCV.BFGS[idxs.min,5,], FEMbasis ))
R_plot_graph.ggplot2(FEM( coef.ex, FEMbasis))

plot(lambda.RightCV.Grad[1:30,5],pch=16,col="red",ylab="lambdas",xlab="")
points(lambda.SimplifiedCV.Grad[1:30,5],pch = 13,col = "blue")

plot(lambda.RightCV.Grad[1:30,5],pch=16,col="red",ylab="lambdas",xlab="")
points(lambda.SimplifiedCV.Grad[1:30,5],pch = 13,col = "blue")

dev.off()