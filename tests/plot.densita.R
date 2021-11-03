### plot densità ###
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
coef.ex = eval.FEM( FEM(aux_n(mesh$nodes[,1],mesh$nodes[,2]),FEMbasis.fdaPDE), locations=mesh$nodes)
x11()
plot(my_dens)

R_plot_graph(FEM(coef.ex, FEMbasis.fdaPDE) )
x11()
R_plot_graph.ggplot2(FEM(coef.ex, FEMbasis.fdaPDE))

# utilizzare un delta 10 volte più piccolo? 0.03 -> 0.003 (anche 5 volte più piccolo è sufficiente)
mesh.ref = refine.mesh.1.5D(mesh.fdaPDE,delta/5)
coef.ref = aux_n(mesh.ref$nodes[,1],mesh.ref$nodes[,2])

x11()
R_plot_graph.ggplot2(FEM(coef.ref, create.FEM.basis(mesh.ref) ))

M=max(coef.ref)
m=min(coef.ref)



min.v =  -floor((m-floor(m))*10)*1e-1 + floor(m)
max.v =  floor((ceiling(M)-M)*10)*1e-1+ floor(M)
