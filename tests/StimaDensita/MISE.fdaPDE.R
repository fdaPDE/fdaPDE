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
delta=0.03

mesh = create.mesh.1D.vertices(vertices, edges, delta)

### fdaPDE ###
mesh.fdaPDE = create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)
###

mu11 = vertices[6,]
mu22 = vertices[8,]
mu1 = mu11
mu2 = mu22
sigmax = 0.1
sigmay = 0.1
w1 = 0.5
w2 = 0.5

aux = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], sigma_x=sigmax, sigma_y=sigmay) 
{ ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))}

integral_ex = integrate_f( FEM( aux(mesh$nodes[,1],mesh$nodes[,2]), FEMbasis=FEMbasis.fdaPDE) )  
aux_n = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], sigma_x=sigmax, sigma_y=sigmay) 
{ 1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))}


my_dens <- linfun(aux_n, L) # FUNZIONE SU UN NETWORK
plot(my_dens)
coef.ex = aux_(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])
R_plot_graph.ggplot2( FEM(coef.ex,FEMbasis.fdaPDE) )

nobs=c(50,100,250,500)
N = length(nobs)
M = 30
mise.l2 = matrix(0,nrow=M,ncol=N)
mise.int = matrix(0,nrow=M,ncol=N)
integrals = matrix(0,nrow=M,ncol=N)

for( i in 1:N){
  for(j in 1:M){
    PPP = rlpp(nobs[i], my_dens)
    points_x = PPP$data$x
    points_y = PPP$data$y
    points = NULL 
    points = cbind(points_x, points_y)
    points = as.matrix(points)
    colnames(points) = NULL

    lambda = c(0.00001,0.0001, 0.001, 0.01, 0.1, 1) 
    nfolds=10
    sol = DE.FEM(data=points,FEMbasis=FEMbasis.fdaPDE,lambda=lambda, 
                    nfolds=nfolds, 
                    step_method="Fixed_Step", 
                    direction_method="Gradient",
                    preprocess_method = "RightCV")
      coef = exp(sol$g)
      integrals[j,i] = integrate_f( FEM(coef,FEMbasis.fdaPDE) )
      coef = 1./integrals[j,i] * coef
      mise.int[j,i]  = integrate_f( FEM( (coef-coef.ex)^2 , FEMbasis.fdaPDE))
      mise.l2[j,i]   = norm(coef-coef.ex,type="2")
  }
}

save(mise.l2,mise.int,integrals,file="mise_03_11.RData")


boxplot(mise)
boxplot(mise.int)
boxplot(mise.int.adj)
