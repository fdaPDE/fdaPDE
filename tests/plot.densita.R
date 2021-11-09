library(ggplot2)
library(plotrix)
source("tests/StimaDensita/Densita_Mattina.R")
source("tests/integrate_f.R")
source("tests/R_plot_graph.R")

### Bifurc Domain ### 
nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh = create.mesh.1.5D(nodes,edges,order=2)
mesh = refine.mesh.1.5D(mesh=mesh,0.1)
plot(mesh)

FEMbasis=create.FEM.basis(mesh)
f=function(x,y){ sin(x)*cos(y) + sin(x*y) }
coef.ex = f(mesh$nodes[,1],mesh$nodes[,2])
FEM = FEM(coef.ex,FEMbasis)

R_plot_graph(FEM)

x11()
R_plot_graph.ggplot2(FEM(coef.ex,FEMbasis))

mesh.ref = refine.mesh.1.5D(mesh,0.1/10)
coef.ref = f(mesh.ref$nodes[,1],mesh.ref$nodes[,2])
FEM.ref  = FEM( coef.ref, create.FEM.basis(mesh.ref))
x11()
R_plot_graph.ggplot2(FEM.ref)

#####

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
coef.ex = aux_n(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])
x11()
plot(my_dens)

R_plot_graph(FEM(coef.ex, FEMbasis.fdaPDE) )

x11()
R_plot_graph.ggplot2(FEM(coef.ex, FEMbasis.fdaPDE))

# mesh refinement (or eval.FEM( ref.nodes ... ) ) 
mesh.ref = refine.mesh.1.5D(mesh.fdaPDE,delta/6)
coef.ref = aux_n(mesh.ref$nodes[,1],mesh.ref$nodes[,2])

x11()
R_plot_graph.ggplot2(FEM(coef.ref, create.FEM.basis(mesh.ref) ))

### 2.5D ###

data(hub2.5D)
mesh <- create.mesh.2.5D(nodes = hub2.5D$hub2.5D.nodes,triangles = hub2.5D$hub2.5D.triangles,order=2)
FEMbasis <- create.FEM.basis(mesh)

loc = mesh$nodes
nloc = dim(loc)[1]
nloc
# 2.5D random field (function f)
a1 = rnorm(1,mean = 0, sd = 1)
a2 = rnorm(1,mean = 0, sd = 1)
a3 = rnorm(1,mean = 0, sd = 1)

sol_exact = numeric(nloc)
for (i in 0:(nloc-1)){
  sol_exact[i+1] = a1* sin(2*pi*loc[i+1,1]) +  a2* sin(2*pi*loc[i+1,2]) +  a3*sin(2*pi*loc[i+1,3]) + 7
}

plot( FEM(sol_exact,FEMbasis))

