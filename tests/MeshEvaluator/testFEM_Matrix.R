# Test FEM matrix (Expr Templates)
#NB) funzioni in smoothing_CPP.R

nodes=matrix(c(0.25,0.25,0.5,0.25,0.75,0.5,0.75,0.), nrow = 4, byrow=TRUE)
edges=matrix(c(1,2,2,3,2,4),nrow = 3,byrow = TRUE)
mesh_ = create.mesh.1.5D(nodes,edges,order=2)
plot.mesh.1.5D(mesh_)

FEMbasis <- create.FEM.basis(mesh = mesh_, saveTree = FALSE)

Stiff <- fdaPDE::CPP_get.FEM.Stiff.Matrix(FEMbasis = FEMbasis)
Mass  <- fdaPDE::CPP_get.FEM.Mass.Matrix(FEMbasis = FEMbasis)

###############################

nodes=matrix(nrow=2,ncol=2)
nodes[,1] = c(0.,1.0)
nodes[,2] = c(0.,0.)

edges=matrix(c(1,2),nrow = 1,byrow = TRUE)
mesh_ = create.mesh.1.5D(nodes,edges,order=2)
plot.mesh.1.5D(mesh_)

FEMbasis <- create.FEM.basis(mesh = mesh_, saveTree = FALSE)

#NB utilizzi formula quadratura ordine 2 quando ORDER=1 (matrice massa!)
#   utilizzi formula  --        ordine 4 quando ORDER=2  
Stiff <- fdaPDE::CPP_get.FEM.Stiff.Matrix(FEMbasis = FEMbasis)
Mass  <- fdaPDE::CPP_get.FEM.Mass.Matrix(FEMbasis = FEMbasis)

Stiff
Mass

# simplenet -------------------------------------------
library(spatstat)
source("~/Scrivania/fdaPDE/tests/Auxiliary/Regressione_Mattina.R")
vertices_x = as.matrix(simplenet$vertices$x)
vertices_y = as.matrix(simplenet$vertices$y)
vertices = as.matrix(cbind(vertices_x, vertices_y))
start_coord = cbind(simplenet$lines$ends$x0, simplenet$lines$ends$y0)
end_coord = cbind(simplenet$lines$ends$x1, simplenet$lines$ends$y1)
start = simplenet$from 
end = simplenet$to
edges = cbind(start, end)
M = simplenet$m
find.boundary(M)
L = as.linnet(simplenet)

delta = 0.03
mesh = create.mesh.1D.vertices(vertices, edges, delta)
FEMbasis = create.FEM.basis.1D(mesh)


### fdaPDE - MESH ### 
mesh.fdaPDE = create.mesh.1.5D(mesh$nodes,mesh$segments)
FEMbasis.fdaPDE = fdaPDE::create.FEM.basis(mesh.fdaPDE)


# end mesh

R0 = R_mass_1D(FEMbasis)
R1 = R_stiff_1D(FEMbasis)

R0.fdaPDE = fdaPDE::CPP_get.FEM.Mass.Matrix(FEMbasis = FEMbasis.fdaPDE)
R1.fdaPDE = fdaPDE::CPP_get.FEM.Stiff.Matrix(FEMbasis = FEMbasis.fdaPDE)

err.R0 = abs(R0- R0.fdaPDE)
err.R1 = abs(R1 - R1.fdaPDE)

which(err.R0 > 10 * .Machine$double.eps)
which(err.R1 > 10 * .Machine$double.eps) # errore dell'ordine di 10^-14 e inferiore.
