if(!dir.exists(test_path("../data/GSR-PDE")))
  dir.create(test_path("../data/GSR-PDE"))

if(!dir.exists(test_path("../data/GSR-PDE/test_23"))){
  dir.create(test_path("../data/GSR-PDE/test_23"))  

options(warn=-1)
foldername <- test_path("../data/GSR-PDE/test_23/")

# family
FAMILY = "poisson"

l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

# beta
beta1 = 5

# lambda
lambda = 10^seq(-8,3,length.out = 20)

# scale param
scale.param = 1

# mesh data
load(test_path("../data/domain/horseshoe2D_areal.RData"))
nodes = mesh$nodes

nnodes = dim(mesh$nodes)[1]
FEMbasis = create.FEM.basis(mesh)

# 2D random field (function f)
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

tri = mesh$triangles
dim(incidence_matrix)
dim(tri)
integration_nodes = data.frame(matrix(nrow = dim(incidence_matrix)[2],ncol = 11))
names(integration_nodes) = c("T1","T2","N1","N2","N3","N4", "xmin", "xmax", "ymin", "ymax", "label")
nodi_plot = NULL

for(i in 1:dim(incidence_matrix)[2])
{
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

# plot areal mesh
nnodes = dim(mesh$nodes)[1]
plot_values = numeric(nnodes)

for(i in 1:dim(incidence_matrix)[2])
{
  plot_values[integration_nodes$N1[i]] = 15
  plot_values[integration_nodes$N2[i]] = 15
  plot_values[integration_nodes$N3[i]] = 15
  plot_values[integration_nodes$N4[i]] = 15
}

#image(FEM(plot_values,FEMbasis))

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

# covariates
desmat=matrix(0,nrow=length(sol_integrand),ncol=1)
desmat[,1]=rbeta(length(sol_integrand),shape1=2,shape2=2)  # sampling covariates from beta distr.

theta = sol_integrand + desmat%*%beta1

mu = inv.link(theta)

# sampling response
response <- rpois(length(mu), lambda = mu)

#### Test 23.1: Without GCV
invisible(capture.output(sol_ref <- smooth.FEM(observations = as.numeric(response), FEMbasis = FEMbasis, covariates = desmat,
                         max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                         incidence_matrix = t(incidence_matrix), areal.data.avg = FALSE,
                         lambda = lambda)))
save(sol_ref, file=paste0(foldername,"/test_23_1.RData"))

#### Test 23.2: grid with exact GCV
invisible(capture.output(sol_ref <- smooth.FEM(observations = as.numeric(response), FEMbasis = FEMbasis, covariates = desmat,
                         max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL,
                         incidence_matrix = t(incidence_matrix), areal.data.avg = FALSE,
                         lambda = lambda, lambda.selection.criterion = 'grid', 
                         DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')))
save(sol_ref, file=paste0(foldername,"/test_23_2.RData"))
}
