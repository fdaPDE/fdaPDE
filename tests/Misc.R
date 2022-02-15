#### Regression Test 1 Simplenet - NO COVARIATE ### -----------------------
# GEOMETRY + DENSITY ------------------------------------------------------
{
source("~/Scrivania/fdaPDE/tests/Auxiliary/Regressione_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")
vertices_x = as.matrix(simplenet$vertices$x)
vertices_y = as.matrix(simplenet$vertices$y)
vertices = as.matrix(cbind(vertices_x, vertices_y))
start_coord = cbind(simplenet$lines$ends$x0, simplenet$lines$ends$y0)
end_coord = cbind(simplenet$lines$ends$x1, simplenet$lines$ends$y1)
start = simplenet$from 
end = simplenet$to
edges = cbind(start, end)
L = as.linnet(simplenet)

delta = 0.03
mesh = create.mesh.1D.vertices(vertices, edges, delta)
FEMbasis = create.FEM.basis.1D(mesh)

### fdaPDE - MESH ### 
mesh.fdaPDE = create.mesh.1.5D(mesh$nodes,mesh$segments)
FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)
idx.leaf = which( mesh.fdaPDE$nodesmarkers==TRUE)

plot(mesh.fdaPDE)
points(mesh.fdaPDE$nodes[8,1], mesh.fdaPDE$nodes[8,2], pch=16)
points(mesh.fdaPDE$nodes[37,1], mesh.fdaPDE$nodes[37,2], pch=16)

### CAMPO f ### 
dijkstra.8 = Dijkstra(mesh.fdaPDE,8)
dijkstra.37 = Dijkstra(mesh.fdaPDE,37)
aux = function(x, y, seg, tp) { 
  
  sigma.8 = 0.1
  sigma.37 = 0.085
  
  res.8 = equal_split_discontinous(mesh.fdaPDE, sigma.8, dijkstra.8, x, y)
  res.37 = equal_split_discontinous(mesh.fdaPDE, sigma.37, dijkstra.37, x, y) 
  
  res = 0.5 * res.8$coef + 0.5 * res.37$coef
  return(res)
}

integral_ex = integrate_f( FEM(aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2]),FEMbasis.fdaPDE) )

my_dens <- linfun(aux, L)

coef.ex = aux(mesh.fdaPDE$nodes[,1],mesh.fdaPDE$nodes[,2])
R_plot_graph.ggplot2(FEM(coef.ex,FEMbasis.fdaPDE))
} 

# LOADING DATA ------------------------------------------------------------
{
today_ = "2022-01-12"
ntest_ = "-test-1-no-cov"
file.name = paste("Regression-",today_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
load(save.file)
img.file = paste("/home/aldo/Scrivania/fdaPDE-IMG/Regression/",file.name,".pdf",sep="")
nnodes = nrow(mesh.fdaPDE$nodes)
# mean coef 500 obs
mean.coef.fdaPDE = matrix(0, nrow=nnodes, ncol=1)
mean.coef.R = mean.coef.fdaPDE

for(j in 1:M){
  mean.coef.fdaPDE =  mean.coef.fdaPDE + sols.fdaPDE[j,4,]
  mean.coef.R = mean.coef.R + sols[j,4,]
}

mean.coef.fdaPDE = mean.coef.fdaPDE / M
mean.coef.R = mean.coef.R / M

mesh.ref = refine.mesh.1.5D(mesh.fdaPDE, 0.01)
FEMbasis.ref = create.FEM.basis(mesh.ref)

mean.coef.fdaPDE.ref = eval.FEM( FEM(mean.coef.fdaPDE, FEMbasis.fdaPDE), mesh.ref$nodes)
mean.coef.R.ref = eval.FEM( FEM(mean.coef.R, FEMbasis.fdaPDE), mesh.ref$nodes)
coef.ex.ref = eval.FEM(FEM(coef.ex, FEMbasis.fdaPDE), mesh.ref$nodes)
}

# IMGs #
img.file = "/home/aldo/Scrivania/fdaPDE-IMG/Regression/Regression-TEST-1-MISC.pdf"
pdf(img.file)
plot(mesh.fdaPDE)
R_plot_graph.ggplot2(FEM(coef.ex.ref, FEMbasis.ref),
                     title="True spatial field",
                     line.size = 1)
R_plot_graph.ggplot2(FEM(mean.coef.fdaPDE.ref, FEMbasis.ref), 
                     title="fdaPDE",
                     line.size = 1)
R_plot_graph.ggplot2(FEM(mean.coef.R.ref, FEMbasis.ref),
                     title="R",
                     line.size=1)
R_plot_graph.R(FEM(coef.ex.ref, FEMbasis.ref), 
               title="True spatial field",
               line.size = 3)
R_plot_graph.R(FEM(mean.coef.fdaPDE.ref, FEMbasis.ref), 
               title="fdaPDE",
               line.size = 3)
R_plot_graph.R(FEM(mean.coef.R.ref, FEMbasis.ref), 
               title="R",
               line.size = 3)

min.col = min(coef.ex.ref, mean.coef.fdaPDE.ref)
min.col = min(min.col, mean.coef.R.ref)

max.col = max(coef.ex.ref, mean.coef.fdaPDE.ref)
max.col = max(max.col, mean.coef.R.ref)

R_plot_graph.ggplot2.2(FEM(coef.ex.ref, FEMbasis.ref),
                       title="True spatial field",
                       line.size = 1,
                       ratio = 1,
                       color.min = min.col,
                       color.max = max.col)
R_plot_graph.ggplot2.2(FEM(mean.coef.fdaPDE.ref, FEMbasis.ref),
                       title="fdaPDE",
                       line.size = 1,
                       ratio = 1,
                       color.min = min.col,
                       color.max = max.col)
R_plot_graph.ggplot2.2(FEM(mean.coef.R.ref, FEMbasis.ref),
                       title="R",
                       line.size = 1,
                       ratio = 1,
                       color.min = min.col,
                       color.max = max.col)


dev.off()



#### Regression Test 2 Brickwall - with COVARIATES #### -------------------
# GEOMETRY + DENSITY ----------------------------------------------------
{
source("~/Scrivania/fdaPDE/tests/Auxiliary/Regressione_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")
data("spiders")
vertices = cbind(spiders$domain$vertices$x, spiders$domain$vertices$y)
edges = cbind(spiders$domain$from, spiders$domain$to)

L = as.linnet(spiders$domain)

#412 nodes -> delta = 40
delta =40
mesh = create.mesh.1D.vertices(vertices, edges, delta)
mesh$nnodes
FEMbasis = create.FEM.basis.1D(mesh)

### fdaPDE - MESH ### 
mesh.fdaPDE = create.mesh.1.5D(mesh$nodes,mesh$segments)
FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)

dijkstra.55 = Dijkstra(mesh.fdaPDE,55)
dijkstra.66= Dijkstra(mesh.fdaPDE,66)
dijkstra.98 = Dijkstra(mesh.fdaPDE, 98)

# sym gaussian 
aux.55 = function(x, y, seg, tp) { 
  
  sigma.55 = 500
  
  res.55 = equal_split_discontinous(mesh.fdaPDE, sigma.55, dijkstra.55, x, y)
  
  res =  res.55$coef 
  return(res)
}

aux.66 = function(x, y, seg, tp) { 
  
  sigma.66 = 500
  
  res.66= equal_split_discontinous(mesh.fdaPDE, sigma.66, dijkstra.66, x, y)
  
  res =  res.66$coef 
  return(res)
}

aux.98 = function(x, y, seg, tp) { 
  
  sigma.98 = 500
  
  res.98= equal_split_discontinous(mesh.fdaPDE, sigma.98, dijkstra.98, x, y)
  
  res =  res.98$coef 
  return(res)
}

# mixture of three gaussian distributions
AUX = function(x,y,seg,tp){
  sym = aux.98(x,y,seg,tp)
  sym.2 = aux.55(x,y,seg,tp)
  sym.3 = aux.66(x,y,seg,tp)
  res = (sym + sym.2 + sym.3)
  return(res)                    
}

coef.ex = AUX(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2])
}

# LOADING DATA ------------------------------------------------------------
{
today_ = "2022-01-20"
ntest_ = "-test-2"
file.name = paste("Regression-",today_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/Regression/",file.name,".RData",sep="")
load(save.file)
nnodes = nrow(mesh.fdaPDE$nodes)
# mean coef 1000 obs
mean.coef.fdaPDE = matrix(0, nrow=nnodes, ncol=1)
mean.coef.R = mean.coef.fdaPDE

for(j in 1:M){
  mean.coef.fdaPDE =  mean.coef.fdaPDE + sols.fdaPDE[j,4,]
  mean.coef.R = mean.coef.R + sols[j,4,]
}

mean.coef.fdaPDE = mean.coef.fdaPDE / M
mean.coef.R = mean.coef.R / M

mesh.ref = refine.mesh.1.5D(mesh.fdaPDE, 20)
FEMbasis.ref = create.FEM.basis(mesh.ref)

mean.coef.fdaPDE.ref = eval.FEM( FEM(mean.coef.fdaPDE, FEMbasis.fdaPDE), mesh.ref$nodes)
mean.coef.R.ref = eval.FEM( FEM(mean.coef.R, FEMbasis.fdaPDE), mesh.ref$nodes)
coef.ex.ref = eval.FEM(FEM(coef.ex, FEMbasis.fdaPDE), mesh.ref$nodes)

PP = runiflpp(300, L)
W = rnorm(300, mean=0, sd = 0.5)
}

img.file = "/home/aldo/Scrivania/fdaPDE-IMG/Regression/Regression-TEST-2-MISC.pdf"
pdf(img.file)

plot(mesh.fdaPDE)
p <- heat.colors(n=48,alpha=1)
points(PP$data$x, PP$data$y,pch=16, cex=2, col=p)

R_plot_graph.ggplot2(FEM(coef.ex.ref, FEMbasis.ref),
                     title="True spatial field",
                     line.size = 1)
R_plot_graph.ggplot2(FEM(mean.coef.fdaPDE.ref, FEMbasis.ref), 
                     title="fdaPDE",
                     line.size = 1)
R_plot_graph.ggplot2(FEM(mean.coef.R.ref, FEMbasis.ref),
                     title="R",
                     line.size=1)

R_plot_graph.R(FEM(coef.ex.ref, FEMbasis.ref), 
               title="True spatial field",
               line.size = 3)
R_plot_graph.R(FEM(mean.coef.fdaPDE.ref, FEMbasis.ref), 
               title="fdaPDE",
               line.size = 3)
R_plot_graph.R(FEM(mean.coef.R.ref, FEMbasis.ref), 
               title="R",
               line.size = 3)

min.col = min(coef.ex.ref, mean.coef.fdaPDE.ref)
min.col = min(min.col, mean.coef.R.ref)

max.col = max(coef.ex.ref, mean.coef.fdaPDE.ref)
max.col = max(max.col, mean.coef.R.ref)

R_plot_graph.ggplot2.2(FEM(coef.ex.ref, FEMbasis.ref),
                       title="True spatial field",
                       line.size = 1,
                       ratio = 1,
                       color.min = min.col,
                       color.max = max.col,
                       a.sym = T)
R_plot_graph.ggplot2.2(FEM(mean.coef.fdaPDE.ref, FEMbasis.ref),
                       title="fdaPDE",
                       line.size = 1,
                       ratio = 1,
                       color.min = min.col,
                       color.max = max.col,
                       a.sym = T)
R_plot_graph.ggplot2.2(FEM(mean.coef.R.ref, FEMbasis.ref),
                       title="R",
                       line.size = 1,
                       ratio = 1,
                       color.min = min.col,
                       color.max = max.col,
                       a.sym = T)


  x=vector(mode="double")
  y=vector(mode="double")
  grp.nodes=vector(mode="integer")
  
  num_edges= dim(mesh.fdaPDE$edges)[1]
  for(e in 1:num_edges){
    x = append(x, c(mesh.fdaPDE$nodes[ mesh.fdaPDE$edges[e,1], 1], mesh.fdaPDE$nodes[ mesh.fdaPDE$edges[e,2], 1]))
    y = append(y, c(mesh.fdaPDE$nodes[ mesh.fdaPDE$edges[e,1], 2], mesh.fdaPDE$nodes[ mesh.fdaPDE$edges[e,2], 2]))
    grp.nodes = append(grp.nodes, rep(e,times=2))
  }
  
  coef.points = W
  x.points = PP$data$x
  y.points = PP$data$y
  data.points = data.frame(x.points, y.points, coef.points)
  
  MyTheme <- theme(
    axis.text = element_text(size=24),
    axis.title = element_text(size=26),
    title = element_text(size=26),
    legend.text = element_text(size=20),
    legend.key.size = unit(1,"cm"),
    legend.position = "right"
  )
  
  data=data.frame(x,y,grp.nodes)
  p = jet.col(n=128, alpha=0.8) 
  ggplot() + 
    geom_point(data=data, aes(x=x,y=y,group=grp.nodes),alpha=0.0) +
    geom_point(data=data.points,aes(x=x.points,y=y.points, color=coef.points), size=3)+
    geom_line(data=data, aes(x=x,y=y,group=grp.nodes),size=0.5)+
    labs(x="",y="",color="", title="Covariate") +  
    coord_fixed(ratio=1) + 
    scale_color_gradientn(colours=p)+
    theme_void() +
    theme(plot.title = element_text(hjust=0.5),
          title = element_text(size=26),
          legend.title = element_blank(),
          axis.title = element_blank(),
          legend.text = element_text(size=20),
          legend.key.size = unit(1,"cm"),
          legend.key.width = unit(0.5,"cm"))


dev.off()
  

#### DE Test 1 SIMPLENET -------------------------------------------------
# GEOMETRY + DENSITY -----------------------------------------------------
{
library(spatstat)
library(ggplot2)
source("~/Scrivania/fdaPDE/tests/Auxiliary/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")
vertices_x = as.matrix(simplenet$vertices$x)
vertices_y = as.matrix(simplenet$vertices$y)
vertices = as.matrix(cbind(vertices_x, vertices_y))
start_coord = cbind(simplenet$lines$ends$x0, simplenet$lines$ends$y0)
end_coord = cbind(simplenet$lines$ends$x1, simplenet$lines$ends$y1)
start = simplenet$from 
end = simplenet$to
edges = cbind(start, end)
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

### density ### 
dijkstra.8 = Dijkstra(mesh.fdaPDE,8)
dijkstra.37 = Dijkstra(mesh.fdaPDE,37)
aux = function(x, y, seg, tp) { 
  
  sigma.8 = 0.1
  sigma.37 = 0.085
  
  res.8 = equal_split_discontinous(mesh.fdaPDE, sigma.8, dijkstra.8, x, y)
  res.37 = equal_split_discontinous(mesh.fdaPDE, sigma.37, dijkstra.37, x, y) 
  
  res = 0.5 * res.8$coef + 0.5 * res.37$coef
  return(res)
}

density <- linfun(aux, L)
coef.ex = aux(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2])
 
}

# LOADING DATA ------------------------------------------------------------
{
today_ = "2022-01-14"
init_ = "-heat" # "-fdaPDE""-null_vector"
ntest_ = "-test-1"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
load(save.file)

M = 20
coef.fdaPDE = matrix(0, nrow=nnodes, ncol=1)
coef.R = coef.fdaPDE

for(j in 1:M){
  coef.fdaPDE =  coef.fdaPDE + sols.fdaPDE[j,4,]
  coef.R = coef.R + sols[j,4,]
}

coef.fdaPDE = coef.fdaPDE / M
coef.R = coef.R / M
mesh.ref = refine.mesh.1.5D(mesh.fdaPDE, 0.01)
FEMbasis.ref = create.FEM.basis(mesh.ref)

coef.fdaPDE.ref = eval.FEM( FEM(coef.fdaPDE, FEMbasis.fdaPDE), mesh.ref$nodes)
coef.R.ref = eval.FEM( FEM(coef.R, FEMbasis.fdaPDE), mesh.ref$nodes)
coef.ex.ref = eval.FEM(FEM(coef.ex, FEMbasis.fdaPDE), mesh.ref$nodes)

PP = rlpp(200, density)

}

# IMGs #
img.file = "/home/aldo/Scrivania/fdaPDE-IMG/DE/DE-TEST-1-MISC.pdf"

pdf(img.file)
plot(mesh.fdaPDE)
points(PP$data$x, PP$data$y, pch=16, cex=1.5, col="darkgreen")
R_plot_graph.ggplot2(FEM(coef.ex.ref, FEMbasis.ref),
                     title="Density",
                     line.size = 1)
R_plot_graph.ggplot2(FEM(coef.fdaPDE.ref, FEMbasis.ref), 
                     title="fdaPDE",
                     line.size = 1)
R_plot_graph.ggplot2(FEM(coef.R.ref, FEMbasis.ref),
                     title="R",
                     line.size=1)
R_plot_graph.R(FEM =FEM(coef.ex.ref, FEMbasis.ref), 
               title="Density", 
               line.size = 3)
R_plot_graph.R(FEM =FEM(coef.fdaPDE.ref, FEMbasis.ref), 
               title="fdaPDE", 
               line.size = 3)
R_plot_graph.R(FEM =FEM(coef.R.ref, FEMbasis.ref), 
               title="R", 
               line.size = 3)

min.col = min(coef.ex.ref, coef.fdaPDE.ref)
min.col = min(min.col, coef.R.ref)

max.col = max(coef.ex.ref, coef.fdaPDE.ref)
max.col = max(max.col, coef.R.ref)

R_plot_graph.ggplot2.2(FEM(coef.ex.ref, FEMbasis.ref),
                     title="Density",
                     line.size = 1,
                     color.min = min.col,
                     color.max = max.col)
R_plot_graph.ggplot2.2(FEM(coef.fdaPDE.ref, FEMbasis.ref), 
                     title="fdaPDE",
                     line.size = 1,
                     color.min = min.col,
                     color.max = max.col)
R_plot_graph.ggplot2.2(FEM(coef.R.ref, FEMbasis.ref),
                     title="R",
                     line.size=1,
                     color.min = min.col,
                     color.max = max.col)

dev.off()

img.file = "/home/aldo/Scrivania/prova.pdf"

pdf(img.file)
plot(mesh.fdaPDE)
points(PP$data$x, PP$data$y, pch=16, cex=1.5, col="darkgreen")
R_plot_graph.ggplot2.2(FEM(coef.ex.ref, FEMbasis.ref),
                     title="Density",
                     line.size = 1)
R_plot_graph.ggplot2.2(FEM(coef.fdaPDE.ref, FEMbasis.ref), 
                     title="fdaPDE",
                     line.size = 1)
R_plot_graph.ggplot2.2(FEM(coef.R.ref, FEMbasis.ref),
                     title="R",
                     line.size=1)

dev.off()

#### DE Test 2 BRICKWALL -------------------------------------------------
# GEOMETRY + DENSITY -----------------------------------------------------
{
source("~/Scrivania/fdaPDE/tests/Auxiliary/Densita_Mattina.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
source("tests/Auxiliary/R_plot_graph.ggplot2.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")

data("spiders")
vertices = cbind(spiders$domain$vertices$x, spiders$domain$vertices$y) 
edges = cbind(spiders$domain$from, spiders$domain$to)
L = as.linnet(spiders$domain)

mesh = create.mesh.1D.vertices(vertices, edges, delta=40)
nnodes = mesh$nnodes
FEMbasis = create.FEM.basis.1D(mesh) 

mesh.fdaPDE = create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)

dijkstra.66 = Dijkstra(mesh.fdaPDE,66)
dijkstra.55 = Dijkstra(mesh.fdaPDE,55)
dijkstra.98 = Dijkstra(mesh.fdaPDE,98)

################### Density 
aux.55 = function(x, y, seg, tp) { 
  
  sigma.55 = 90
  
  res.55= equal_split_discontinous(mesh.fdaPDE, sigma.55, dijkstra.55, x, y)
  
  res =  res.55$coef 
  return(res)
}
aux.98 = function(x, y, seg, tp) { 
  
  sigma.98 = 90
  
  res.98= equal_split_discontinous(mesh.fdaPDE, sigma.98, dijkstra.98, x, y)
  res =  res.98$coef 
  
  return(res)
}
aux.66 = function(x,y,tp,seg){
  sigma = 100
  h = 5*sigma
  Graph = mesh.fdaPDE
  dijkstra = dijkstra.66
  source = dijkstra$source
  points_ = cbind(x,y)
  
  coef = vector(mode="numeric", length=nrow(points_))
  is_vertex = is.vertex(Graph, points_)
  idx.vertex = which(is_vertex!=0)
  idx.same.y = which(points_[idx.vertex,2] == Graph$nodes[source,2])
  
  if(!is.empty(idx.same.y)){
    for( i in idx.same.y ){
      if( dijkstra$distance[idx.vertex[i]] < h ){
        coef[i] =  1./((2*pi)^0.5*sigma) * exp(-dijkstra$distance[idx.vertex[i]]^2/(2*sigma^2))
      }
    }
  }
  
  idx.not.vertex = which(is_vertex==0)
  if(!is.empty(idx.not.vertex)){
    idx.edge = isInside(Graph, points_)
    idx.edge = idx.edge[idx.not.vertex]
    Dist1 = dijkstra$distance[ Graph$edges[idx.edge, 1]]
    Dist2 = dijkstra$distance[ Graph$edges[idx.edge, 2]]
    Dist = vector(mode="numeric", length=length(idx.not.vertex))
    Previous = vector(mode="integer", length(idx.not.vertex))
    
    for( i in 1:length(idx.not.vertex)){
      if( abs( points_[idx.not.vertex[i],2] - Graph$nodes[source,2]) < 10 * .Machine$double.eps){
        
        if( Dist1[i] < Dist2[i] ){
          Dist[i] = Dist1[i]
          Previous[i] = Graph$edges[idx.edge[i], 1]
        }else{
          Dist[i] = Dist2[i]
          Previous[i] = Graph$edges[idx.edge[i], 2]
        }
        
        delta = Dist[i] + sqrt( (Graph$nodes[Previous[i] ,1] - points_[idx.not.vertex[i],1])^2 +
                                  (Graph$nodes[Previous[i], 2] - points_[idx.not.vertex[i],2])^2) 
        
        if( delta < h)
          coef[idx.not.vertex[i]] = 1/((2*pi)^0.5 * sigma) * exp(-delta^2/(2*sigma^2))  
        
      }
    }
  }
  return (coef)
}
AUX = function(x,y,seg,tp){
  sym = aux.55(x,y,seg,tp)
  sym.2 = aux.98(x,y,seg,tp)
  sym.3 = aux.66(x,y,seg,tp)
  res = 1./3. * sym + 1./3. * sym.2 + 1./3. * sym.3
  return(res)                    
}
density <- linfun(AUX, L)
coef.ex = AUX(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2])
}


# LOADING DATA ------------------------------------------------------------
{
  
today_ = "2022-01-16"
init_ ="-heat" #"-random"#"-null_vector" 
ntest_ = "-test-5-SimplifiedCV"

file.name = paste("DE-",today_,init_,ntest_,sep="")
save.file = paste("/home/aldo/Scrivania/fdaPDE-DATA/DE/",file.name,".RData",sep="")
load(save.file)

M = 20
coef.fdaPDE = matrix(0, nrow=nnodes, ncol=1)
coef.R = coef.fdaPDE

for(j in 1:M){
  coef.fdaPDE =  coef.fdaPDE + sols.fdaPDE[j,4,]
  coef.R = coef.R + sols[j,4,]
}

coef.fdaPDE = coef.fdaPDE / M
coef.R = coef.R / M


mesh.ref = refine.mesh.1.5D(mesh.fdaPDE, 20)
FEMbasis.ref = create.FEM.basis(mesh.ref)

coef.fdaPDE.ref = eval.FEM( FEM(coef.fdaPDE, FEMbasis.fdaPDE), mesh.ref$nodes)
coef.R.ref = eval.FEM( FEM(coef.R, FEMbasis.fdaPDE), mesh.ref$nodes)
coef.ex.ref = eval.FEM(FEM(coef.ex, FEMbasis.fdaPDE), mesh.ref$nodes)

PP = rlpp(300, density)
}

# IMGs #
img.file = "/home/aldo/Scrivania/fdaPDE-IMG/DE/DE-TEST-2-MISC.pdf"

pdf(img.file)
plot(mesh.fdaPDE)
points(PP$data$x, PP$data$y, pch=16, cex=1.5, col="darkgreen")
R_plot_graph.ggplot2(FEM(coef.ex.ref, FEMbasis.ref),
                     title="Density",
                     line.size = 1)
R_plot_graph.ggplot2(FEM(coef.fdaPDE.ref, FEMbasis.ref), 
                     title="fdaPDE",
                     line.size = 1)
R_plot_graph.ggplot2(FEM(coef.R.ref, FEMbasis.ref),
                     title="R",
                     line.size=1)
R_plot_graph.R(FEM =FEM(coef.ex.ref, FEMbasis.ref), 
               title="Density", 
               line.size = 3)
R_plot_graph.R(FEM =FEM(coef.fdaPDE.ref, FEMbasis.ref), 
               title="fdaPDE", 
               line.size = 3)
R_plot_graph.R(FEM =FEM(coef.R.ref, FEMbasis.ref), 
               title="R", 
               line.size = 3)
dev.off()

#### DE Test 2 BRICKWALL ALDO -------------------------------------------------
# GEOMETRY + DENSITY -----------------------------------------------------
{library(ggplot2)
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Densita_Mattina.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
  source("tests/Auxiliary/R_plot_graph.ggplot2.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")
  
  data("brickwall.1.5D")
  
  vertices = nodes
  delta=0.03
  
  N = length(delta)
  
  mesh = create.mesh.1D.vertices(vertices, edges, delta)
  nodes = mesh$nodes
  nnodes = mesh$nnodes
  FEMbasis = create.FEM.basis.1D(mesh) 
  
  ### fdaPDE ###
  mesh.fdaPDE = create.mesh.1.5D(nodes=mesh$nodes,edges=mesh$segments)
  FEMbasis.fdaPDE = create.FEM.basis(mesh=mesh.fdaPDE)
  plot(mesh.fdaPDE, show.nodes = TRUE, pch=16)
  ################### Density #################
  dijkstra.7 = Dijkstra(mesh.fdaPDE,7)
  
  aux.7 = function(x,y,seg, tp){
    sigma = 0.09
    h = 5*sigma
    Graph = mesh.fdaPDE
    dijkstra = dijkstra.7
    source = dijkstra$source
    points_ = cbind(x,y)
    coef = vector(mode="numeric", length=length(x) ) 
    
    idx.ok = which( abs(points_[,2] - Graph$nodes[source,2]) < 10 * .Machine$double.eps )
    for( i in idx.ok){
      
      delta =  abs(points_[i,1]-Graph$nodes[source,1])
      
      if(delta < h ){
        coef[i] = 1/((2*pi)^0.5 * sigma) * exp(-delta^2/(2*sigma^2))
      }
      
    }
    
    return (coef)
  }
  
  aux.11 = function(x,y,seg, tp){
    sigma = 0.09
    h = 5*sigma
    Graph = mesh.fdaPDE
    dijkstra = dijkstra.11
    source = dijkstra$source
    points_ = cbind(x,y)
    coef = vector(mode="numeric", length=length(x) ) 
    
    idx.ok = which( abs(points_[,2] - Graph$nodes[source,2]) < 10 * .Machine$double.eps )
    for( i in idx.ok){
      
      delta =  abs(points_[i,1]-Graph$nodes[source,1])
      
      if(delta < h ){
        coef[i] = 1/((2*pi)^0.5 * sigma) * exp(-delta^2/(2*sigma^2))
      }
      
    }
    
    return (coef)
  }
  
  AUX = function(x,y,seg,tp){
  
    res = 0.5* ( aux.7(x,y,seg,tp) + aux.11(x,y,seg,tp))
    return(res)                    
  }
  density <- linfun(AUX, L)
  coef.ex = AUX(mesh.fdaPDE$nodes[,1], mesh.fdaPDE$nodes[,2])
  
  R_plot_graph.ggplot2(FEM(coef.ex,FEMbasis.fdaPDE))
}


# LOADING DATA ------------------------------------------------------------
{
  
  load("/home/aldo/Scrivania/fdaPDE-DATA/DE/DE-2022-02-04test-2-brick-aldo-nobs-var.RData")
  
  M = 20
  coef.fdaPDE = matrix(0, nrow=nnodes, ncol=1)
  coef.R = coef.fdaPDE
  
  for(j in 1:M){
    coef.fdaPDE =  coef.fdaPDE + sols.fdaPDE[j,4,]
    coef.R = coef.R + sols[j,4,]
  }
  
  coef.fdaPDE = coef.fdaPDE / M
  coef.R = coef.R / M
  
  PP = rlpp(50, density)
}

# IMGs #
img.file = "/home/aldo/Scrivania/fdaPDE-IMG/DE/DE-TEST-2-aldo-brick-MISC.pdf"

min.col = min(coef.ex, coef.fdaPDE)
min.col = min(min.col, coef.R)
max.col = max(coef.ex, coef.fdaPDE)
max.col = max(max.col, coef.R)

pdf(img.file)
plot(mesh.fdaPDE,asp=1)
points(PP$data$x, PP$data$y, pch=16, cex=1.5, col="darkgreen")



R_plot_graph.ggplot2.2(FEM(coef.ex, FEMbasis.fdaPDE),
                     title="Density",
                     color.min = min.col,
                     color.max = max.col,
                     line.size = 1,
                     ratio = 1,
                     a.sym = TRUE)
R_plot_graph.ggplot2.2(FEM(coef.fdaPDE, FEMbasis.fdaPDE), 
                     title="fdaPDE",
                     color.min = min.col,
                     color.max = max.col,
                     line.size = 1,
                     ratio = 1,
                     a.sym = TRUE)

R_plot_graph.ggplot2.2(FEM(coef.R, FEMbasis.fdaPDE),
                     title="R",
                     color.min = min.col,
                     color.max = max.col,
                     line.size=1,
                     ratio = 1,
                     a.sym = TRUE)
R_plot_graph.R(FEM =FEM(coef.ex, FEMbasis.fdaPDE), 
               title="Density", 
               line.size = 3)
R_plot_graph.R(FEM =FEM(coef.fdaPDE, FEMbasis.fdaPDE), 
               title="fdaPDE", 
               line.size = 3)
R_plot_graph.R(FEM =FEM(coef.R, FEMbasis.fdaPDE), 
               title="R", 
               line.size = 3)
dev.off()

