# TESTS FOR MESH C 2D AREAL DATA
rm(list=ls())
graphics.off()

library(fdaPDE)
library(purrr)
library(pracma)
library(ggplot2)

# FAMILY CHOICE
FAMILY = "poisson"

l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

# BETA
beta1= 5

# lambda 
lambda = 10^seq(-8,3,length.out = 20)

# scale param 
scale.param = 1


#data(AREAL_MeshC)
data(horseshoe2D)
mesh = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
mesh = refine.mesh.2D(mesh, maximum_area = 0.044, minimum_angle = 30)

nodes = mesh$nodes

#load("../data/meshC_areal.RData")

# x11()
plot(mesh, lwd=3, cex = 1.9)
nnodes = dim(mesh$nodes)[1]
FEMbasis = fdaPDE::create.FEM.basis(mesh)

# 2D random field (function f) ------------------------------
a_fun <- function(p){
  
  if(p[1]>= 0 && p[2] > 0){
    pi/4 + p[1]
  }else{
    if(p[1]>= 0 && p[2] <= 0){
      -pi/4 -p[1]
    }else{
      if(p[1] < 0){
        -0.5*atan(p[2]/p[1])
      }
    }
  }
}

d_fun <- function(p){
  
  if(p[1]>= 0 && p[2] > 0){
    -0.5 + p[2]
  }else{
    if(p[1]>= 0 && p[2] <= 0){
      -0.5 - p[2]
    }else{
      if(p[1] < 0){
        sqrt(p[1]^2 + p[2]^2) - 0.5
      }
    }
  }
}

z <- function(p){
  a_fun(p) + d_fun(p)^2
}

list_of_ones = c(239, 244, 596, 602, 898, 900, 1256, 1258, 1504, 1517,
                 1838, 1839, 2172, 2173, 2497, 2498, 2796, 2848, 3125,
                 3128, 3482, 3491, 3857, 3859, 4191, 4212, 4542, 4543,
                 4892, 4985, 5343, 5345, 5698, 5706, 6030, 6053, 6402,
                 6404, 6763, 6775)
tri = mesh$triangles
incidence_matrix = zeros(dim(tri)[1],length(lambda))
incidence_matrix[list_of_ones]=1
dim(incidence_matrix)
dim(tri)
integration_nodes = data.frame(matrix(nrow = dim(incidence_matrix)[2],ncol = 11))
names(integration_nodes) = c("T1","T2","N1","N2","N3","N4", "xmin", "xmax", "ymin", "ymax", "label")
nodi_plot = NULL

for(i in 1:dim(incidence_matrix)[2]){

  tri_used = which(incidence_matrix[,i] == 1)
  integration_nodes$T1[i] = tri_used[1]
  integration_nodes$T2[i] = tri_used[2]
  nodes_used = unique(c(tri[tri_used[1],],tri[tri_used[2],]))
  integration_nodes$N1[i] = nodes_used[1]
  integration_nodes$N2[i] = nodes_used[2]
  integration_nodes$N3[i] = nodes_used[3]
  integration_nodes$N4[i] = nodes_used[4]
  integration_nodes$label[i] = i
  xvec = c(nodes[integration_nodes$N1[i],1],nodes[integration_nodes$N2[i],1],nodes[integration_nodes$N3[i],1],nodes[integration_nodes$N4[i],1])
  yvec = c(nodes[integration_nodes$N1[i],2],nodes[integration_nodes$N2[i],2],nodes[integration_nodes$N3[i],2],nodes[integration_nodes$N4[i],2])
  integration_nodes$xmin[i] = min(xvec)
  integration_nodes$xmax[i] = max(xvec)
  integration_nodes$ymin[i] = min(yvec)
  integration_nodes$ymax[i] = max(yvec)
  
  nodi_tmp = rbind(mesh$nodes[integration_nodes$N1[i],],mesh$nodes[integration_nodes$N2[i],],mesh$nodes[integration_nodes$N3[i],],mesh$nodes[integration_nodes$N4[i],])
  nodi_tmp = cbind(nodi_tmp,rep(i,4))
  nodi_plot = rbind(nodi_plot,nodi_tmp)
  
}

nodi_plot <- as.data.frame(nodi_plot)
names(nodi_plot) = c("X","Y","Tri")
nodi_plot$Tri = as.factor(nodi_plot$Tri)

# plot mesh areale ---------------------------------------------

nnodes = dim(mesh$nodes)[1]
plot_values = numeric(nnodes)

for(i in 1:dim(incidence_matrix)[2]){

    plot_values[integration_nodes$N1[i]] = 15
    plot_values[integration_nodes$N2[i]] = 15
    plot_values[integration_nodes$N3[i]] = 15
    plot_values[integration_nodes$N4[i]] = 15
    

}

image(FEM(plot_values,FEMbasis))

# ---------------------------------------------------------------
nnodes = dim(mesh$nodes)[1]
sol_exact <- numeric(nnodes)
for(i in 1:nnodes){
  sol_exact[i] <- z(mesh$nodes[i,])
}

z_xy <- function(x,y){z(c(x,y))}
sol_integrand = matrix(nrow = 20 , ncol = 1)
for(i in 1:20){
  res = integral2(fun = z_xy, xmin = integration_nodes$xmin[i], xmax = integration_nodes$xmax[i] , ymin = integration_nodes$ymin[i], ymax = integration_nodes$ymax[i])
  sol_integrand[i] = res$Q
}

range(sol_integrand)

set.seed(2020)
# covariates ---------------------------
desmat=matrix(0,nrow=length(sol_integrand),ncol=1)
desmat[,1]=rbeta(length(sol_integrand),shape1=2,shape2=2)  # sampling covariates from beta distr.


theta = sol_integrand + desmat%*%beta1
range(theta)

mu = inv.link(theta)
range(mu)

# sampling response:
response <- rpois(length(mu), lambda = mu)


output_CPP_new <- fdaPDE::smooth.FEM(observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                                  lambda = lambda, max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                                  incidence_matrix = t(incidence_matrix), areal.data.avg = FALSE,
                                  lambda.selection.criterion="grid", DOF.evaluation='exact')



plot(log10(lambda),output_CPP$GCV)
best_fit <- output_CPP$bestlambda
image(FEM(output_CPP$fit.FEM$coeff[,best_fit], FEMbasis))
image(FEM(sol_integrand,FEMbasis))
t(sol_integrand)




