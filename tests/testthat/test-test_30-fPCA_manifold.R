test_that("fPCA - manifold",{
  options(warn=-1)
  foldername <- test_path("../data/fPCA/test_30/")
  
  load(test_path("../data/domain/sphere2.5D.RData"))
  mesh<-create.mesh.2.5D(sphere2.5D$nodes, sphere2.5D$triangles)
  
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
  
  #### Test 30.1: Without GCV
  lambda = 10^-2
  invisible(capture.output(sol <-FPCA.FEM(datamatrix=dm.pointwise.centred,
                                              FEMbasis=FEMbasis,
                                              lambda=lambda,
                                              nPC=2)))
  load(file=paste0(foldername,"/test_30_1.RData"))
  expect_equal( max(abs((sol$loadings.FEM$coeff-sol_ref$loadings.FEM$coeff))) < tol, TRUE);
  expect_equal( max(abs((sol$scores-sol_ref$scores))) < tol, TRUE);
  
  #### Test 30.2: K-Fold validation - Exact
  invisible(capture.output(sol <-FPCA.FEM(datamatrix=dm.pointwise.centred,
                                              FEMbasis=FEMbasis,lambda=lambda,nPC=2,
                                              validation='KFold',GCVmethod='Exact')))
  load(file=paste0(foldername,"/test_30_2.RData"))
  expect_equal( max(abs((sol$loadings.FEM$coeff-sol_ref$loadings.FEM$coeff))) < tol, TRUE);
  expect_equal( max(abs((sol$scores-sol_ref$scores))) < tol, TRUE);
  
  #### Test 30.3: GCV validation - Exact
  invisible(capture.output(sol <-FPCA.FEM(datamatrix=dm.pointwise.centred,
                                              FEMbasis=FEMbasis,lambda=lambda,nPC=3,
                                              validation='GCV',GCVmethod='Exact')))
  load(file=paste0(foldername,"/test_30_3.RData"))
  expect_equal( max(abs((sol$loadings.FEM$coeff-sol_ref$loadings.FEM$coeff))) < tol, TRUE);
  expect_equal( max(abs((sol$scores-sol_ref$scores))) < tol, TRUE);
})