test_that("SR-PDE Unit Square",{
  
  foldername = test_path("../data/SR-PDE/test_6/")  
  x = seq(0,1, length.out = 20)
  y = x
  locations = expand.grid(x,y)
  
  mesh = create.mesh.2D(locations)
  
  nnodes = dim(mesh$nodes)[1]
  
  FEMbasis = create.FEM.basis(mesh)
  
  # Test function
  f = function(x, y, z = 1)
  {
    coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
    sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
  }
  
  # Exact solution (pointwise at nodes)
  sol_exact = f(mesh$nodes[,1], mesh$nodes[,2])
  
  # Add error to simulate data
  set.seed(7893475)
  ran = range(sol_exact)
  data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  # Set smoothing parameter
  lambda = 10^seq(-6,-3,length=20)
  
  options(warn=-1)  
  #### Test 6.1: grid with exact GCV
  invisible(capture.output(sol<-smooth.FEM(observations=data, 
                                           FEMbasis=FEMbasis, 
                                           lambda=lambda, 
                                           lambda.selection.criterion='grid', 
                                           DOF.evaluation='exact', 
                                           lambda.selection.lossfunction='GCV')))
  load(file = paste0(foldername, "/test_6_1.RData"))
  expect_equal( max(abs((sol$fit.FEM$coeff-sol_ref$fit.FEM$coeff))) < tol, TRUE);
  
  ### Test 6.2: Newton exact method with exact GCV, default initial lambda and tolerance
  invisible(capture.output(sol<-smooth.FEM(observations=data, FEMbasis=FEMbasis, 
                                           lambda.selection.criterion='newton',
                                           DOF.evaluation='exact', lambda.selection.lossfunction='GCV')))
  
  load(file = paste0(foldername, "/test_6_2.RData"))
  expect_equal( max(abs((sol$fit.FEM$coeff-sol_ref$fit.FEM$coeff))) < tol, TRUE);
})
