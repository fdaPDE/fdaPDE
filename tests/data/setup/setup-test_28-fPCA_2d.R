if(!dir.exists(test_path("../data/fPCA")))
  dir.create(test_path("../data/fPCA"))

if(!dir.exists(test_path("../data/fPCA/test_28"))){
  dir.create(test_path("../data/fPCA/test_28"))  
  
  options(warn=-1)
  foldername <- test_path("../data/fPCA/test_28/")

x = seq(0,1, length.out = 20)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)

nnodes = dim(mesh$nodes)[1]

FEMbasis = create.FEM.basis(mesh)

eigenfunc1=FEM(sin(2*pi*mesh$nodes[,1]), FEMbasis)
eigenfunc2=FEM(cos(2*pi*mesh$nodes[,1]), FEMbasis)
eigenfunc3=FEM(cos(2*pi*mesh$nodes[,2]), FEMbasis)

truedatarange<-max(c(eigenfunc1$coeff,eigenfunc2$coeff,eigenfunc3$coeff))-min(c(eigenfunc1$coeff,eigenfunc2$coeff,eigenfunc3$coeff))
truecoeff<-cbind(eigenfunc1$coeff,eigenfunc2$coeff,eigenfunc3$coeff)

set.seed(5847947)

nSamples=50

sd_score1<-0.1
sd_score2<-0.05
sd_score3<-0.03 
sd_error<-0.05

score1<-rnorm(n=nSamples,sd=sd_score1*truedatarange)
score2<-rnorm(n=nSamples,sd=sd_score2*truedatarange)
score3<-rnorm(n=nSamples,sd=sd_score3*truedatarange)

datamatrix.pointwise.exact<-matrix(score1)%*%t(matrix(eigenfunc1$coeff))+matrix(score2)%*%t(matrix(eigenfunc2$coeff))+matrix(score3)%*%t(matrix(eigenfunc3$coeff))
dm.pointwise.centred.exact<-datamatrix.pointwise.exact-matrix(apply(datamatrix.pointwise.exact,2,mean),ncol=ncol(datamatrix.pointwise.exact),nrow=nrow(datamatrix.pointwise.exact),byrow=TRUE)

nnodes=dim(mesh$nodes)[1]

error<-rnorm(n=nSamples*nnodes,sd=sd_error*truedatarange)

datamatrix.pointwise<-datamatrix.pointwise.exact+error
dm.pointwise.centred<-datamatrix.pointwise-matrix(apply(datamatrix.pointwise,2,mean),ncol=ncol(datamatrix.pointwise),nrow=nrow(datamatrix.pointwise),byrow=TRUE)

#### Test 28.1: Without GCV
lambda = 10^-2
invisible(capture.output(sol_ref <- FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,
                        lambda=lambda,
                        nPC=3)))
save(sol_ref, file=paste0(foldername,"/test_28_1.RData"))

#### Test 28.2: K-Fold validation - Exact
lambda = 10^c(-6,-5,-4,-3,-2)

invisible(capture.output(sol_ref <-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                        validation='KFold',GCVmethod='Exact')))
save(sol_ref, file=paste0(foldername,"/test_28_2.RData"))

#### Test 28.3: GCV validation - Exact
invisible(capture.output(sol_ref <-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                        validation='GCV',GCVmethod='Exact')))
save(sol_ref, file=paste0(foldername,"/test_28_3.RData"))

}
