if(!dir.exists(test_path("../data/fPCA")))
  dir.create(test_path("../data/fPCA"))

if(!dir.exists(test_path("../data/fPCA/test_31"))){
  dir.create(test_path("../data/fPCA/test_31"))  

  options(warn=-1)
  foldername <- test_path("../data/fPCA/test_31/")
  
data(sphere3Ddata)
mesh = create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)

nnodes = dim(mesh$nodes)[1]

FEMbasis = create.FEM.basis(mesh)

eigenfunc1=FEM(sin(2*pi*mesh$nodes[,1]), FEMbasis)
eigenfunc2=FEM(cos(2*pi*mesh$nodes[,1]), FEMbasis)

truedatarange<-max(c(eigenfunc1$coeff,eigenfunc2$coeff))-min(c(eigenfunc1$coeff,eigenfunc2$coeff))

set.seed(5847947)

nSamples=50

sd_score1<-0.1
sd_score2<-0.05
sd_error<-0.05

score1<-rnorm(n=nSamples,sd=sd_score1*truedatarange)
score2<-rnorm(n=nSamples,sd=sd_score2*truedatarange)

datamatrix.pointwise.exact<-matrix(score1)%*%t(matrix(eigenfunc1$coeff))+matrix(score2)%*%t(matrix(eigenfunc2$coeff))
dm.pointwise.centred.exact<-datamatrix.pointwise.exact-matrix(apply(datamatrix.pointwise.exact,2,mean),ncol=ncol(datamatrix.pointwise.exact),nrow=nrow(datamatrix.pointwise.exact),byrow=TRUE)

error<-rnorm(n=nSamples*nnodes,sd=sd_error*truedatarange)

datamatrix.pointwise<-datamatrix.pointwise.exact+error
dm.pointwise.centred<-datamatrix.pointwise-matrix(apply(datamatrix.pointwise,2,mean),ncol=ncol(datamatrix.pointwise),nrow=nrow(datamatrix.pointwise),byrow=TRUE)

#### Test 31.1: Without GCV
lambda = 10^-2
invisible(capture.output(sol_ref <-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,
                        lambda=lambda,
                        nPC=2)))
save(sol_ref, file=paste0(foldername,"/test_31_1.RData"))

#### Test 31.2: K-Fold validation - Exact
lambda = 10^c(-6,-5,-4,-3,-2)

invisible(capture.output(sol_ref <-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=2,
                        validation='KFold',GCVmethod='Exact')))
save(sol_ref, file=paste0(foldername,"/test_31_2.RData"))

#### Test 31.3: GCV validation - Exact
invisible(capture.output(sol_ref <-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=2,
                        validation='GCV',GCVmethod='Exact')))
save(sol_ref, file=paste0(foldername,"/test_31_3.RData"))
}
