library(ggplot2)
library(plotrix)
source("tests/StimaDensita/Densita_Mattina.R")
source("tests/integrate_f.R")
source("tests/R_plot_graph.R")

# simplenet #
vertices = as.matrix(cbind(simplenet$vertices$x, simplenet$vertices$x) )
start_coord = cbind(simplenet$lines$ends$x0, simplenet$lines$ends$y0)
end_coord = cbind(simplenet$lines$ends$x1, simplenet$lines$ends$y1)
edges = cbind(simplenet$from,simplenet$to)
L = as.linnet(simplenet)
delta=0.03

mesh.fdaPDE = create.mesh.1.5D(nodes=vertices,edges=edges)
mesh.fdaPDE = refine.mesh.1.5D(mesh.fdaPDE, delta)

FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)

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

integral_ex = integrate_f( FEM( aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2]), FEMbasis=FEMbasis.fdaPDE) )  
aux_n = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], sigma_x=sigmax, sigma_y=sigmay) 
{ 1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))}

my_dens = linfun(aux_n,L)
plot(my_dens)
coef.ex = aux_n(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])

AUX = function(x){  1/(2*pi)^1.5 * exp(-x^2/(2*sigmax)) }

res = Dijkstra(mesh.fdaPDE, source=6)

coef.AUX = AUX( res$distance ) 
R_plot_graph(FEM(coef.AUX, FEMbasis.fdaPDE))
plot(mesh.fdaPDE)

# dendrite #
library(spatstat)
data("dendrite")

vertices = cbind(dendrite$domain$vertices$x, dendrite$domain$vertices$y) 
edges = cbind(dendrite$domain$from, dendrite$domain$to)
L = as.linnet(dendrite$domain)
M = dendrite$domain$m

### fdaPDE ###
mesh.fdaPDE = create.mesh.1.5D(nodes=vertices,edges=edges)
FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)
###

plot(mesh.fdaPDE)
points(mesh.fdaPDE$nodes[600,1],mesh.fdaPDE$nodes[600,2],pch=16)
points(mesh.fdaPDE$nodes[70,1],mesh.fdaPDE$nodes[70,2],pch=16)
points(mesh.fdaPDE$nodes[10,1],mesh.fdaPDE$nodes[10,2],pch=16)
# density
mu1 = mesh.fdaPDE$nodes[70,]
mu2 = mesh.fdaPDE$nodes[600,]
mu3 = mesh.fdaPDE$nodes[10,]
sigmax = 20
sigmay = 20
aux = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], x3=mu3[1], y3=mu3[2], sigma_x=sigmax, sigma_y=sigmay) 
{   ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))+ 
    ( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x3)^2 + (1/(2*sigma_y^2))*(y-y3)^2 ))}


integral_ex = integrate_f( FEM(aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2]),FEMbasis.fdaPDE) )

aux_n = function(x, y, seg, tp, x1=mu1[1], y1=mu1[2], x2=mu2[1], y2=mu2[2], x3=mu3[1], y3=mu3[2], sigma_x=sigmax, sigma_y=sigmay) 
{   1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x1)^2 + (1/(2*sigma_y^2))*(y-y1)^2 ))+
    1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x2)^2 + (1/(2*sigma_y^2))*(y-y2)^2 ))+
    1./integral_ex*( 1/(2*pi*sigma_x*sigma_y)  )*exp(- ( (1/(2*sigma_x^2))*(x-x3)^2 + (1/(2*sigma_y^2))*(y-y3)^2 ))}
my_dens <- linfun(aux, L) # FUNZIONE SU UN NETWORK

coef.ex = aux_n(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])

mesh.ref.2 = refine.mesh.1.5D(mesh.fdaPDE,20)
coef.ref.2 = aux_n(mesh.ref.2$nodes[,1],mesh.ref$nodes[,2])
R_plot_graph.ggplot2(FEM(coef.ref.2, create.FEM.basis(mesh.ref.2) ))


#### 


source("~/Scrivania/fdaPDE/tests/StimaDensita/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Mesh_Evaluator/isInside.R")
source("tests/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Dijkstra.R")
library(ggplot2)

idx.leaf = which(mesh.fdaPDE$nodesmarkers==TRUE)
plot(mesh.fdaPDE)
points(mesh.fdaPDE$nodes[idx.leaf,1], mesh.fdaPDE$nodes[idx.leaf,2], pch=16)
points(mesh.fdaPDE$nodes[250,1], mesh.fdaPDE$nodes[250,2], pch=16,col="red")

dijkstra.130 = Dijkstra(mesh.fdaPDE,130)
coef.130 = equal_split_discontinous(mesh.fdaPDE, 12.5, dijkstra.130, mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])

dijkstra.250 = Dijkstra(mesh.fdaPDE, 250)
coef.250 = equal_split_discontinous(mesh.fdaPDE, 12.5, dijkstra.250, mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])

R_plot_graph.ggplot2(FEM(coef$coef, FEMbasis.fdaPDE))

integrate_f(FEM(coef.250$coef, FEMbasis.fdaPDE))

aux <- function(x,y,seg,tp){
    
    sigma = 12.5
    res.130 = equal_split_discontinous(mesh.fdaPDE, sigma, dijkstra.130, x, y)
    res.250 = equal_split_discontinous(mesh.fdaPDE, sigma, dijkstra.250, x, y) 
    
    res = 0.5 * res.130$coef + 0.5 * res.250$coef 
    
    return(res)
}

integrate_f(FEM(aux(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2]), FEMbasis.fdaPDE ))
my_dens = linfun(aux, L)
plot(FEM(my_dens(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2]), FEMbasis.fdaPDE) )#questo ok

plot(my_dens)

PP = rlpp(100, my_dens)
plot(mesh.fdaPDE)
points(PP$data$x, PP$data$y,pch=16)
