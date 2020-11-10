##########################################
############## TEST SCRIPT ###############
################ for FPCA ################
##########################################

library(fdaPDE)

#### Test 1: 2D square domain ####
#            locations = nodes 
#            order FE = 1
rm(list=ls())
graphics.off()

x = seq(0,1, length.out = 60)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)
plot(mesh)

nnodes = dim(mesh$nodes)[1]

FEMbasis = create.FEM.basis(mesh)

eigenfunc1=FEM(sin(2*pi*mesh$nodes[,1]), FEMbasis)
eigenfunc2=FEM(cos(2*pi*mesh$nodes[,1]), FEMbasis)
eigenfunc3=FEM(cos(2*pi*mesh$nodes[,2]), FEMbasis)

plot(eigenfunc1)
plot(eigenfunc2)
plot(eigenfunc3)

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

#### Test 1.1: Without GCV
lambda = 10^-2
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,
                        lambda=lambda,
                        nPC=3)

plot(sol.pointwise$loadings.FEM)

#### Test 1.2: K-Fold validation - Stochastic
lambda = 10^c(-6,-5,-4,-3,-2)

sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                        validation='KFold',GCVmethod='Stochastic')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda

#### Test 1.3: K-Fold validation - Exact
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                        validation='KFold',GCVmethod='Exact')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda

#### Test 1.4: GCV validation - Stochastic
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                        validation='GCV',GCVmethod='Stochastic')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda

#### Test 1.5: GCV validation - Exact
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                        validation='GCV',GCVmethod='Exact')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda


#### Test 2: 2D square domain ####
#            locations != nodes
#            order FE = 2
rm(list=ls())
graphics.off()

x = seq(0,1, length.out = 31)
y = x
vertex = expand.grid(x,y)

mesh = create.mesh.2D(vertex,order = 2)

FEMbasis = create.FEM.basis(mesh)

x = seq(0,1, length.out = 41)
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

#### Test 2.1: Without GCV
lambda = 10^-2
sol.pointwise<-FPCA.FEM(locations = locations,
                        datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,
                        lambda=lambda,
                        nPC=3)

plot(sol.pointwise$loadings.FEM)

#### Test 2.2: K-Fold validation - Stochastic
lambda = 10^c(-6,-5,-4,-3,-2)

sol.pointwise<-FPCA.FEM(locations = locations,
                        datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                        validation='KFold',GCVmethod='Stochastic')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda

#### Test 2.3: K-Fold validation - Exact
sol.pointwise<-FPCA.FEM(locations = locations,
                        datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                        validation='KFold',GCVmethod='Exact')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda

#### Test 2.4: GCV validation - Stochastic
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                        validation='GCV',GCVmethod='Stochastic')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda

#### Test 2.5: GCV validation - Exact
sol.pointwise<-FPCA.FEM(locations = locations,
                        datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                        validation='GCV',GCVmethod='Exact')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda


#### Test 3: 2.5D sphere domain ####
#            locations = nodes 
#            order FE = 1
rm(list=ls())
graphics.off()

data(sphere2.5D)
mesh<-create.mesh.2.5D(sphere2.5D$nodes, sphere2.5D$triangles)

nnodes = dim(mesh$nodes)[1]

FEMbasis = create.FEM.basis(mesh)

eigenfunc1=FEM(sin(2*pi*mesh$nodes[,1]), FEMbasis)
eigenfunc2=FEM(cos(2*pi*mesh$nodes[,1]), FEMbasis)

plot(eigenfunc1)
plot(eigenfunc2)

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

#### Test 3.1: Without GCV
lambda = 10^-2
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,
                        lambda=lambda,
                        nPC=2)

plot(sol.pointwise$loadings.FEM)

#### Test 3.2: K-Fold validation - Stochastic
lambda = 10^c(-6,-5,-4,-3,-2)

sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=2,
                        validation='KFold',GCVmethod='Stochastic')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda

#### Test 3.3: K-Fold validation - Exact
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=2,
                        validation='KFold',GCVmethod='Exact')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda

#### Test 3.4: GCV validation - Stochastic
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=2,
                        validation='GCV',GCVmethod='Stochastic')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda

#### Test 3.5: GCV validation - Exact
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                        validation='GCV',GCVmethod='Exact')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda


#### Test 4: 3D sphere domain ####
#            locations = nodes 
#            order FE = 1
rm(list=ls())
graphics.off()

data(sphere3Ddata)
mesh = create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)
plot(mesh)

nnodes = dim(mesh$nodes)[1]

FEMbasis = create.FEM.basis(mesh)

eigenfunc1=FEM(sin(2*pi*mesh$nodes[,1]), FEMbasis)
eigenfunc2=FEM(cos(2*pi*mesh$nodes[,1]), FEMbasis)

plot(eigenfunc1)
plot(eigenfunc2)

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

#### Test 1.1: Without GCV
lambda = 10^-2
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,
                        lambda=lambda,
                        nPC=2)

plot(sol.pointwise$loadings.FEM)

#### Test 1.2: K-Fold validation - Stochastic
lambda = 10^c(-6,-5,-4,-3,-2)

sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=2,
                        validation='KFold',GCVmethod='Stochastic')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda

#### Test 1.3: K-Fold validation - Exact
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=2,
                        validation='KFold',GCVmethod='Exact')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda

#### Test 1.4: GCV validation - Stochastic
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=2,
                        validation='GCV',GCVmethod='Stochastic')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda

#### Test 1.5: GCV validation - Exact
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,
                        FEMbasis=FEMbasis,lambda=lambda,nPC=2,
                        validation='GCV',GCVmethod='Exact')

plot(sol.pointwise$loadings.FEM)
sol.pointwise$lambda