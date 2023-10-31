if(!dir.exists(test_path("../data/fPCA")))
  dir.create(test_path("../data/fPCA"))

if(!dir.exists(test_path("../data/fPCA/test_29"))){
  dir.create(test_path("../data/fPCA/test_29"))  
  
  options(warn=-1)
  foldername <- test_path("../data/fPCA/test_29/")

x = seq(0,1, length.out = 21)
y = x
vertex = expand.grid(x,y)

mesh = create.mesh.2D(vertex,order = 2)

FEMbasis = create.FEM.basis(mesh)

x = seq(0,1, length.out = 11)
y = x
locations = expand.grid(x,y)

ndata = dim(locations)[1]

eigenfunc1=sin(2*pi*locations[,1])
eigenfunc2=cos(2*pi*locations[,1])
eigenfunc3=cos(2*pi*locations[,2])

truedatarange<-max(c(eigenfunc1,eigenfunc2,eigenfunc3))-min(c(eigenfunc1,eigenfunc2,eigenfunc3))

set.seed(5847947)

nSamples=50

sd_score1<-0.1
sd_score2<-0.05
sd_score3<-0.03 
sd_error<-0.05

score1<-rnorm(n=nSamples,sd=sd_score1*truedatarange)
score2<-rnorm(n=nSamples,sd=sd_score2*truedatarange)
score3<-rnorm(n=nSamples,sd=sd_score3*truedatarange)

datamatrix.pointwise.exact<-matrix(score1)%*%t(matrix(eigenfunc1))+matrix(score2)%*%t(matrix(eigenfunc2))+matrix(score3)%*%t(matrix(eigenfunc3))
dm.pointwise.centred.exact<-datamatrix.pointwise.exact-matrix(apply(datamatrix.pointwise.exact,2,mean),ncol=ncol(datamatrix.pointwise.exact),nrow=nrow(datamatrix.pointwise.exact),byrow=TRUE)

error<-rnorm(n=nSamples*ndata,sd=sd_error*truedatarange)

datamatrix.pointwise<-datamatrix.pointwise.exact+error
dm.pointwise.centred<-datamatrix.pointwise-matrix(apply(datamatrix.pointwise,2,mean),ncol=ncol(datamatrix.pointwise),nrow=nrow(datamatrix.pointwise),byrow=TRUE)

#### Test 29.1: Without GCV
lambda = 10^-2
invisible(capture.output(sol_ref <-FPCA.FEM(locations = locations,
                        datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,
                        lambda=lambda,
                        nPC=3)))
save(sol_ref, file=paste0(foldername,"/test_29_1.RData"))

#### Test 29.2: K-Fold validation - Exact
lambda = 10^c(-6,-5,-4,-3,-2)
invisible(capture.output(sol_ref <-FPCA.FEM(locations = locations,
                        datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                        validation='KFold',GCVmethod='Exact')))
save(sol_ref, file=paste0(foldername,"/test_29_2.RData"))

#### Test 29.3: GCV validation - Exact
invisible(capture.output(sol_ref <-FPCA.FEM(locations = locations,
                        datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                        validation='GCV',GCVmethod='Exact')))
save(sol_ref, file=paste0(foldername,"/test_29_3.RData"))
}
