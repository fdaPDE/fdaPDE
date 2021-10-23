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

