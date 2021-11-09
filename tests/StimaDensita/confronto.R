### Confronto Mattina - fdaPDE ###

# APPLICAZIONE A SIMPLENET
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
x11()
plot(my_dens)


coef.ex = aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])

x11()
R_plot_graph.ggplot2( FEM(coef.ex, FEMbasis.fdaPDE) )

R_plot_graph(FEM(coef.ex,FEMbasis.fdaPDE))

###
niter = 1000
nfolds=10
R0 = R_mass_1D(FEMbasis)
R1 = R_stiff_1D(FEMbasis, M)
cvec=matrix(0,nrow=mesh$nnodes,ncol=1)
epsilon = 0.0001 

###
nobs = c(50,100,250,500)
N = length(nobs)
times = matrix(0,nrow=N,ncol=2)
lambda = c(0.0001, 0.001, 0.01, 0.1, 1) 
eta=0.001

M=30
mise = matrix(0,nrow=M,ncol=N)
mise.fdaPDE = matrix(0,nrow=M,ncol=N)
times = matrix(0,nrow=M, ncol=N)
times.fdaPDE = matrix(0,nrow=M,ncol=N)
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
                    nfolds=nfolds, 
                    step_method="Fixed_Step", 
                    direction_method="Gradient",
                    preprocess_method = "RightCV")
  end.fdaPDE <- Sys.time()
  
  start <- Sys.time()
  result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
  PHI = result$PHI
  sol = my_density(points, mesh, PHI, FEMbasis, R0, R1, sol.fdaPDE$lambda, eta, cvec, niter, FALSE)
  end <- Sys.time()
  
  coef = exp(sol$density[,niter])/sol$int_tot
  coef.fdaPDE     = exp(sol.fdaPDE$g)
  integral.fdaPDE = integrate_f( FEM(coef.fdaPDE,FEMbasis.fdaPDE)) 
  coef.fdaPDE = 1./integral.fdaPDE * coef.fdaPDE
  
  mise[j,i]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
  mise.fdaPDE[j,i]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
  times[j,i] = end - start
  times.fdaPDE[j,i] = end.fdaPDE - start.fdaPDE
  }
}
tot.end = Sys.time()
tot.time = tot.end-tot.start

x11()
boxplot(mise[,1],mise.fdaPDE[,1],
        mise[,2],mise.fdaPDE[,2],
        mise[,3],mise.fdaPDE[,3],
        mise[,4],mise.fdaPDE[,4],
        main="MISE",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11),
        col=c("red","blue"))
legend("topright", legend=c("Mattina","fdaPDE"), col=c("red","blue"),pch=c(15,15) )

time.plot = matrix(0,nrow=1,ncol=N)
time.fdaPDE.plot = matrix(0,nrow=1,ncol=N)

for(i in 1:N){
  time.plot[i] = mean(times[,i])
  time.fdaPDE.plot[i] = mean(times.fdaPDE[,i])
}



x11()
plot(nobs, time.plot,type="b",col="red",pch=15 ,xlab = "nobs",ylab = "times", main="mean time." ,ylim=c(min(time.fdaPDE.plot),max(time.plot)))
points(nobs, time.fdaPDE.plot,type="b",pch = 16, col="blue")
legend("topleft",legend=c("Mattina", "fdaPDE (+CV)"), 
       col=c("red","blue"),pch=c(15,16),cex = 0.8,text.font=4 )


x11()
boxplot(times[,1],times.fdaPDE[,1],
        times[,2],times.fdaPDE[,2],
        times[,3],times.fdaPDE[,3],
        times[,4],times.fdaPDE[,4],
        #main="time",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11),
        col=c("red","blue"))
legend("topright", legend=c("eval+algo Mattina","fdaPDE"), col=c("red","blue"),pch=c(15,15) )


### cleaned ####

out.times.fdaPDE.1 = boxplot(times.fdaPDE[,1],plot=FALSE)$out
out.times.fdaPDE.2 = boxplot(times.fdaPDE[,2],plot=FALSE)$out
out.times.fdaPDE.3 = boxplot(times.fdaPDE[,3],plot=FALSE)$out
out.times.fdaPDE.4 = boxplot(times.fdaPDE[,4],plot=FALSE)$out

times.fdaPDE.1 = times.fdaPDE[- which(times.fdaPDE[,1] >= out.times.fdaPDE.1),1]
times.fdaPDE.2 = times.fdaPDE[- which(times.fdaPDE[,2] >= out.times.fdaPDE.2),2]
times.fdaPDE.3 = times.fdaPDE[- which(times.fdaPDE[,3] >= out.times.fdaPDE.3),3]
times.fdaPDE.4 = times.fdaPDE[- which(times.fdaPDE[,4] >= out.times.fdaPDE.4),4]

x11()
boxplot(times[,1],times.fdaPDE.1,
        times[,2],times.fdaPDE.2,
        times[,3],times.fdaPDE.3,
        times[,4],times.fdaPDE.4,
        main="time",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11),
        col=c("red","blue"))
legend("center", legend=c("eval+algo Mattina","fdaPDE + CV"), col=c("red","blue"),pch=c(15,15) )

#### ####

mise.1 = mise[-which(is.nan(mise[,1])),1]
mise.2 = mise[-which(is.nan(mise[,2])),2]
mise.3 = mise[-which(is.nan(mise[,3])),3]
mise.4 = mise[-which(is.nan(mise[,4])),4]

mise.fdaPDE.1 = mise.fdaPDE[-which(is.nan(mise.fdaPDE[,1])),1]
mise.fdaPDE.2 = mise.fdaPDE[-which(is.nan(mise.fdaPDE[,2])),2]
mise.fdaPDE.3 = mise.fdaPDE[-which(is.nan(mise.fdaPDE[,3])),3]
mise.fdaPDE.4 = mise.fdaPDE[-which(is.nan(mise.fdaPDE[,4])),4]
mise.fdaPDE.3 = mise.fdaPDE[,3]

x11()
boxplot(mise.1,mise.fdaPDE.1,
        mise.2,mise.fdaPDE.2,
        mise.3,mise.fdaPDE.3,
        mise.4,mise.fdaPDE.4,
        main="MISE",names=rep(nobs,each=2),at=c(1,2,4,5,7,8,10,11),
        col=c("red","blue"),ylim=c(0.,3.))
legend("topright", legend=c("Mattina","fdaPDE"), col=c("red","blue"),pch=c(15,15) )

out.1 = boxplot(mise.1,plot=FALSE)$out
out.2 = boxplot(mise.2,plot=FALSE)$out
out.3 = boxplot(mise.3,plot=FALSE)$out
out.4 = boxplot(mise.4,plot=FALSE)$out

out.fdaPDE.1 = boxplot(mise.fdaPDE.1,plot=FALSE)$out
out.fdaPDE.2 = boxplot(mise.fdaPDE.2,plot=FALSE)$out
out.fdaPDE.3 = boxplot(mise.fdaPDE.3,plot=FALSE)$out
out.fdaPDE.4 = boxplot(mise.fdaPDE.4,plot=FALSE)$out

mise.fdaPDE.1 = mise.fdaPDE[-which(is.nan(mise.fdaPDE[,1])),1]
mise.fdaPDE.2 = mise.fdaPDE[-which(is.nan(mise.fdaPDE[,2])),2]
mise.fdaPDE.3 = mise.fdaPDE[-which(is.nan(mise.fdaPDE[,3])),3]
mise.fdaPDE.4 = mise.fdaPDE[-which(is.nan(mise.fdaPDE[,4])),4]
