### confronto Mattina - fdaPDE 2.0 ###

source("~/Scrivania/fdaPDE/tests/StimaDensita/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/integrate_f.R")

vertices_x = as.matrix(simplenet$vertices$x)
vertices_y = as.matrix(simplenet$vertices$y)
vertices = as.matrix(cbind(vertices_x, vertices_y))
start_coord = cbind(simplenet$lines$ends$x0, simplenet$lines$ends$y0)
end_coord = cbind(simplenet$lines$ends$x1, simplenet$lines$ends$y1)
start = simplenet$from 
end = simplenet$to
edges = cbind(start, end)
M = simplenet$m
boundary=find.boundary(M)
L = as.linnet(simplenet) 

# delta = 0.03 -> nnodes=97
delta=0.03

mesh = create.mesh.1D.vertices(vertices, edges, delta)
nodes = mesh$nodes
nnodes = mesh$nnodes
FEMbasis = create.FEM.basis.1D(mesh) 

### fdaPDE ###
mesh.fdaPDE = create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)
###

mu1 = vertices[6,]
mu2 = vertices[8,]
sigmax = 0.1
sigmay = 0.1
aux = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], sigma_x=sigmax, sigma_y=sigmay) 
{ ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))}


integral_ex = integrate_f( FEM(aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2]),FEMbasis.fdaPDE) )

aux = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], sigma_x=sigmax, sigma_y=sigmay) 
{ 1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))}
my_dens <- linfun(aux, L) # FUNZIONE SU UN NETWORK

coef.ex = aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])

###
niter = 1000
nfolds=5
R0 = R_mass_1D(FEMbasis)
R1 = R_stiff_1D(FEMbasis, M)
epsilon = 0.0001 

###
nobs = c(100,200,300,400,500)
#nobs = c(50)
N = length(nobs)
#cvec=matrix(0,nrow=mesh$nnodes,ncol=1)

lambda = c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1) 
eta=0.001

M=30
mise = matrix(0,nrow=M,ncol=N)
mise.fdaPDE = matrix(0,nrow=M,ncol=N)
err.L2 = matrix(0,nrow=M,ncol=N)

times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)

integrals = matrix(0,nrow=M,ncol=N)
integrals.fdaPDE = matrix(0,nrow=M,ncol=N)

sols = array(0,dim=c(M,N,mesh$nnodes) ) #contiene vettore dei coeff
sols.fdaPDE = array(0,dim=c(M,N,mesh$nnodes) )

lambda.opt = matrix(0,nrow=M,ncol=N)

tot.start = Sys.time()
for(i in 1:N){
  for(j in 1:M){
    PPP = rlpp(nobs[i], my_dens)
    points_x = PPP$data$x
    points_y = PPP$data$y
    points = NULL 
    points = cbind(points_x, points_y)
    points = as.matrix(points)
    colnames(points) = NULL
    
    #nb. nel tempo di fdaPDE è anche compresa la crossvalidazione.
    start.fdaPDE <- Sys.time()
    sol.fdaPDE = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=lambda, 
                        nsimulations=niter, 
                        stepProposals = eta,
                      #  fvec = exp(cvec),
                        tol1 = 1e-16,
                        nfolds=nfolds, 
                        step_method="Fixed_Step", 
                        direction_method="Gradient",
                        preprocess_method = "RightCV")
    end.fdaPDE <- Sys.time()
    
    start <- Sys.time()
    result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
    # NB. R need log(density) init.
    cvec = log( sol.fdaPDE$f_init[, which( abs(lambda - sol.fdaPDE$lambda) <= .Machine$double.eps )] )
    PHI = result$PHI
    sol = my_density(points, mesh, PHI, FEMbasis, R0, R1, sol.fdaPDE$lambda, eta, cvec, niter, FALSE)
    end <- Sys.time()
    
    coef = exp(sol$density[,niter])
    coef.fdaPDE     = exp(sol.fdaPDE$g)
   
    sols[j,i,] = coef 
    sols.fdaPDE[j,i,] = coef.fdaPDE
    
    mise[j,i]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
    mise.fdaPDE[j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
    err.L2[j,i] = integrate_f(FEM ( (coef.fdaPDE - coef)^2, FEMbasis.fdaPDE) )
    
    times[j,i] = end - start
    times.fdaPDE[j,i] = end.fdaPDE - start.fdaPDE
    
    integrals[j,i] = integrate_f( FEM(coef,FEMbasis.fdaPDE))
    integrals.fdaPDE[j,i] =   integrate_f(FEM(coef,FEMbasis.fdaPDE))
    
    lambda.opt[j,i] = sol.fdaPDE$lambda
    }
}
tot.end = Sys.time()
tot.time = tot.end-tot.start

today_ = Sys.Date()
init_ = "-fdaPDE" # "-fdaPDE""-null_vector"
ntest_ = "-test.2.2"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

save(nobs, N, M, coef.ex,niter,
     mise, mise.fdaPDE,err.L2,
     times, times.fdaPDE,
     integrals, integrals.fdaPDE,
     sols, sols.fdaPDE,
     lambda.opt, tot.time, FEMbasis.fdaPDE,
     file = save.file)


### post-processing ###
#today_ = "2021-12-04"
today_ = Sys.Date()
init_ = "-fdaPDE" # "-null_vector"
ntest_ = "-test.2.1"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/DE/",file.name,".pdf",sep="")

load(save.file)

### ggplot2 ###
library(ggplot2)
library(latex2exp)

##!!! crea mesh -> FEMbasis !!!##

err.L2 = matrix(0,nrow=M,ncol=N)
for( i in 1:N){
  for(j in 1:M){
    
    err.L2[j,i] = integrate_f( FEM( (sols[j,i,]-sols.fdaPDE[j,i,])^2, FEMbasis.fdaPDE ) )
  }
}


mise_ = cbind(as.vector(mise) , as.vector(mise.fdaPDE))
mise_= as.vector(mise_)
times_ = cbind( as.vector(times),as.vector(times.fdaPDE))
times_= as.vector(times_)
nobs_ = rep(nobs,each=M)
nobs_ = rep(nobs_,times=2)
type_ = rep(c("R","fdaPDE"),each=(M*N))

# data frame per boxplots #
data_frame.1 = data.frame(nobs_,
                          mise_,times_,type_)


#data frame per scatterplots#
pdf(img.file)

ggplot(data_frame.1)+
  geom_boxplot(aes(x=nobs_,y=mise_, group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y=bquote( NORMA[L[2]] ^2 ),fill="",title="f")+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data_frame.1)+
  geom_boxplot(aes(x=nobs_,y=times_, group=interaction(nobs_,type_),fill=type_))+
  labs(x="observations",y="time [s]",fill="",title="TIME")+
  theme(plot.title = element_text(hjust = 0.5))

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
  labs(x="observations",y="",fill="",title=bquote(  (f[R] - f[fdaPDE])[L[2]]^2 ) ) +
  theme(plot.title = element_text(hjust = 0.5))

#data_frame.50 = subset(data_frame.2, nobs_ %in% c(50))
#data_frame.100 = subset(data_frame.2, nobs_ %in% c(100))
#data_frame.250 = subset(data_frame.2, nobs_ %in% c(250))
#data_frame.500 = subset(data_frame.2, nobs_ %in% c(250))

for( i in 1:N){
  data_frame.tmp = subset(data_frame.2, nobs_ %in% c(nobs[i]))
    
    x.min = min(data_frame.tmp["mise_"],data_frame.tmp["mise_"])
    x.max = max(data_frame.tmp["mise_"],data_frame.tmp["mise_"])
    
    print(  ggplot(data_frame.tmp,aes(x=mise_,y=mise.fdaPDE_))+
              geom_point() +
              geom_segment(aes(x=x.min, y=x.min, xend=x.max, yend=x.max),color="red",alpha=0.01)+
              labs(x="R", y="fdaPDE", 
                    title= TeX("Scatter $||f - f_{ex}||^{L_{2}}^2$"), 
                    subtitle= paste("nobs = ",nobs[i],sep="") ) +
              theme(plot.title = element_text(hjust=0.5), 
                    plot.subtitle = element_text(hjust=0.5) )    )
  
}

# ggplot(data_frame.50,aes(x=mise_,y=mise.fdaPDE_))+
#  geom_point()+
#  geom_smooth(method="lm",se=TRUE)+
#   labs(x="R", y="fdaPDE", title= bquote("Scatterplot "*(f[R] - f[fdaPDE])[L[2]]^2* " (nobs = 50)" ) )+
#   theme(plot.title = element_text(hjust=0.5))
# 
# ggplot(data_frame.100,aes(x=mise_,y=mise.fdaPDE_))+
#   geom_point()+
#   geom_smooth(method="lm",se=TRUE)+
#   labs(x="R", y="fdaPDE", title= bquote("Scatterplot "*(f[R] - f[fdaPDE])[L[2]]^2* " (nobs = 100)" ) )+
#   theme(plot.title = element_text(hjust=0.5))
# 
# ggplot(data_frame.250,aes(x=mise_,y=mise.fdaPDE_))+
#   geom_point()+
#   geom_smooth(method="lm",se=TRUE,level=0.95)+
#   labs(x="R", y="fdaPDE", title= bquote("Scatterplot "*(f[R] - f[fdaPDE])[L[2]]^2* " (nobs = 250)" ) )+
#   theme(plot.title = element_text(hjust=0.5))
# 
# ggplot(data_frame.500, aes(x=mise_,y=mise.fdaPDE_)) +
#   geom_point()+
#   geom_smooth(method="lm",se=TRUE,level=0.95)+
#   labs(x="R", y="fdaPDE", title= bquote("Scatterplot "*(f[R] - f[fdaPDE])[L[2]]^2* " (nobs = 500)" ) )+
#   theme(plot.title = element_text(hjust=0.5))

for(i in 1:N){
  title.mise = paste("MISE^2, nobs =",as.character(nobs[i]),seq="")
  plot(1:M,mise[,i], pch=16, col="red", xlab="",ylab="",
       main=title.mise ,ylim=c(min(mise[,i],mise.fdaPDE[,i]),max(mise[,i],mise[,i])) )
  points(1:M,mise.fdaPDE[,i], pch=13, col ="blue")
  segments(1:M,mise[,i],1:M,mise.fdaPDE[,i], lty=3)
  legend("topright",legend=c("Mattina","fdaPDE"),col=c("red","blue"),pch=c(16,13) )
}

dev.off()


M
pdf(img.file)
for(i in 1:N){
  title.mise = paste("MISE^2, nobs =",as.character(nobs[i]),seq="")
  plot(1:M,mise[,i], pch=16, col="red", xlab="",ylab="",
       main=title.mise ,ylim=c(min(mise[,i],mise.fdaPDE[,i]),max(mise[,i],mise[,i])) )
  points(1:M,mise.fdaPDE[,i], pch=13, col ="blue")
  segments(1:M,mise[,i],1:M,mise.fdaPDE[,i], lty=3)
  legend("topright",legend=c("Mattina","fdaPDE"),col=c("red","blue"),pch=c(16,13) )
}

boxplot(mise[,1],mise.fdaPDE[,1],
        mise[,2],mise.fdaPDE[,2],
        mise[,3],mise.fdaPDE[,3],
        mise[,4],mise.fdaPDE[,4],
        main="MISE^2",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11),
        col=c("red","blue"),xlab = "nobs")
legend("topright", legend=c("Mattina","fdaPDE"), col=c("red","blue"),pch=c(15,15) )

time.plot = matrix(0,nrow=1,ncol=N)
time.fdaPDE.plot = matrix(0,nrow=1,ncol=N)

for(i in 1:N){
  time.plot[i] = mean(times[,i])
  time.fdaPDE.plot[i] = mean(times.fdaPDE[,i])
}

plot(nobs, time.plot,type="b",col="red",pch=15 ,xlab = "nobs",ylab = "times [s]", main="average time." ,ylim=c(min(time.fdaPDE.plot),max(time.plot)))
points(nobs, time.fdaPDE.plot,type="b",pch = 16, col="blue")
legend("topleft",legend=c("Mattina", "fdaPDE (+CV)"), 
       col=c("red","blue"),pch=c(15,16),cex = 0.8,text.font=4 )

boxplot(times[,1],times.fdaPDE[,1],
        times[,2],times.fdaPDE[,2],
       times[,3],times.fdaPDE[,3],
       times[,4],times.fdaPDE[,4],
        main="time",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11), ylab="times [s]",xlab = "nobs",
        col=c("red","blue"))
legend("topright", legend=c("eval+algo Mattina","fdaPDE"), col=c("red","blue"),pch=c(15,15) )

dev.off()
