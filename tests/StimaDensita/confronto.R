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
plot(my_dens)

coef.ex = aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])
x11()
R_plot_graph.ggplot2( FEM(coef.ex, FEMbasis.fdaPDE) )


###
niter = 20000
nfolds=10
R0 = R_mass_1D(FEMbasis)
R1 = R_stiff_1D(FEMbasis, M)
cvec=matrix(0,nrow=mesh$nnodes,ncol=1)
epsilon = 0.0001 

###
nobs = c(50,100,250,500)
N = length(nobs)
err.l2 = matrix(0,nrow=N,ncol=1)
err =matrix(0,nrow=N,ncol=2)
times = matrix(0,nrow=N,ncol=2)
lambda = c(0.0001, 0.001, 0.01, 0.1, 1) 
eta=0.001

for(i in 1:N){
  PPP = rlpp(nobs[1], my_dens)
  points_x = PPP$data$x
  points_y = PPP$data$y
  points = NULL 
  points = cbind(points_x, points_y)
  points = as.matrix(points)
  colnames(points) = NULL
  
  result = eval.FEM.1D( mesh, edges, vertices, points, epsilon )
  PHI = result$PHI
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
  sol = my_density(points, mesh, PHI, FEMbasis, R0, R1, sol.fdaPDE$lambda, eta, cvec, niter, FALSE)
  end <- Sys.time()
  
  coef = exp(sol$density[,niter])/sol$int_tot
  coef.fdaPDE     = exp(sol.fdaPDE$g)
  integral.fdaPDE = integrate_f( FEM(coef.fdaPDE,FEMbasis.fdaPDE)) 
  coef.fdaPDE = 1./integral.fdaPDE * coef.fdaPDE
  
  err.l2[i]  = norm(coef-coef.fdaPDE,type="2")
  
  err[i,1]   = integrate_f( FEM( (coef-coef.ex)^2, FEMbasis.fdaPDE)) 
  err[i,2]   = integrate_f( FEM( (coef.fdaPDE-coef.ex)^2, FEMbasis.fdaPDE) )
  times[i,1] = end - start
  times[i,2] = end.fdaPDE - start.fdaPDE
}

x11()
plot(1:N, times[,1],type="b",pch = 15,col="red", xlab='',ylab='times', main="time (min)",
     ylim=c(min(times),max(times)),asp=1)
points(1:N, times[,2],type="b",pch = 16, col="blue")
legend("topleft",legend=c( "Mattina R","fdaPDE R"), 
       col=c("red","blue"),pch=c(15,16),cex = 0.8,text.font=4 )

x11()
plot(1:N,err.l2,type="b", pch=15 ,col="red",xlab='',ylab='f_m - f_fda', main="l2 norm",
     ylim=c(min(err.l2),max(err.l2)))

x11()
plot(1:N,err[,1],type="b",pch=15,col="red",xlab='',ylab='f - f.ex', main="L2 norm",ylim=c(min(err),max(err)))
points(1:N,err[,2],type="b",pch=16,col="blue")
legend("topleft",legend=c( "Mattina ","fdaPDE "), 
       col=c("red","blue"),pch=c(15,16),cex = 0.8,text.font=4 )

R_plot_graph.ggplot2(FEM(coef.fdaPDE,FEMbasis.fdaPDE) ) 
integrate_f( FEM(coef.fdaPDE, FEMbasis.fdaPDE))
