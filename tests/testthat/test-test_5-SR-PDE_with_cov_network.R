test_that("Spatial Regression (with covariates) - Linear Network", {
  
  foldername <- test_path("../data/SR-PDE/test_5/")
  
  eps = 1 / 2
  x = c(0., 1)
  y = c(0.,eps)
  vertices = expand.grid(x,y)
  vertices = cbind(vertices[,1], vertices[,2])
  edges = matrix(c(1,2,1,3,3,4), nrow=3,ncol=2, byrow=T)
  
  mesh = create.mesh.1.5D(vertices, edges)
  mesh = refine.mesh.1.5D(mesh,delta=0.0125)
  
  nnodes=dim(mesh$nodes)[1]
  FEMbasis=create.FEM.basis(mesh)
  
  # Exact solution (pointwise at nodes)
  aux.4 = function(x,y){
    h = 1
    source = 4 
    points_ = cbind(x,y)
    idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
    coef = vector(mode="numeric", length=length(x))
    
    for(i in idx.ok){
      delta = abs(points_[i,1] - mesh$nodes[source,1])
      if(delta < h ){
        coef[i] = 3*exp(1/(delta^2 / h^2 -1)+1 ) - 1
        
      }
    }
    
    return(coef)
  }
  aux.3 = function(x,y){
    
    h = eps
    source = 3 
    points_ = cbind(x,y)
    idx.ok = which( abs(points_[,1]-mesh$nodes[source,1]) < 10 * .Machine$double.eps )
    coef = vector(mode="numeric", length=length(x))
    
    for(i in idx.ok){
      delta = abs(points_[i,2] - mesh$nodes[source,2])
      if(delta < h ){
        coef[i] = -1 - 1/h*delta
      }
      
      
    }
    return(coef)
  }
  aux.1 = function(x,y){
    
    h = 1
    source = 1 
    points_ = cbind(x,y)
    idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
    coef = vector(mode="numeric", length=length(x))
    
    for(i in idx.ok){
      delta = abs(points_[i,1] - mesh$nodes[source,1])
      if(delta <= h ){
        coef[i] = -2 - 1/h*delta
      }
    }
    return(coef)
    
  }
  
  f= function(x,y){
    
    res = aux.1(x,y) + aux.3(x,y) + aux.4(x,y)
    return(res)
  }
  
  sol_exact=f(mesh$nodes[,1],mesh$nodes[,2])
  
  # Generate covariate and data
  locations = mesh$nodes
  
  set.seed(32)
  cov1 = rnorm(nnodes, mean=2,sd=1)
  DatiEsatti = f(locations[,1],locations[,2]) + cov1
  
  # Add error to simulate data
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  # Set smoothing parameter
  lambda = 10^seq(-4,-2,by=0.25)
  
  #### Test 2.1: Without GCV
  # invisible(capture.output(sol<-smooth.FEM(observations=data, 
  #                        locations = NULL, 
  #                        covariates = cov1,
  #                        FEMbasis=FEMbasis, 
  #                        lambda=lambda[1])))
  # load(file=paste0(foldername,"/test_5_1.RData"))
  # expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);
  # expect_equal( max(abs((sol$solution$beta-output_CPP$solution$beta))) < tol, TRUE);
  
  #### Test 2.2: grid with exact GCV
  invisible(capture.output(sol<-smooth.FEM(observations=data, locations = NULL,
                         covariates = cov1,
                         FEMbasis=FEMbasis, lambda=lambda,
                         lambda.selection.criterion='grid', 
                         DOF.evaluation='exact', lambda.selection.lossfunction='GCV')))
  
  load(file=paste0(foldername,"/test_5_2.RData"))
  expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);
  expect_equal( max(abs((sol$solution$beta-output_CPP$solution$beta))) < tol, TRUE);
  
  # #### Test 2.3: grid with stochastic GCV
  # invisible(capture.output(sol<-smooth.FEM(observations=data, locations = NULL, 
  #                        covariates = cov1,
  #                        FEMbasis=FEMbasis, lambda=lambda,
  #                        lambda.selection.criterion='grid', 
  #                        DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')))
  # 
  # load(file=paste0(foldername,"/test_5_3.RData"))
  # expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);
  # expect_equal( max(abs((sol$solution$beta-output_CPP$solution$beta))) < tol, TRUE);
   
  #### Test 5.4: Newton method with exact GCV, default initial lambda and tolerance
  invisible(capture.output(sol<-smooth.FEM(observations=data, locations = NULL, 
                         covariates = cov1,
                         FEMbasis=FEMbasis, 
                         lambda.selection.criterion='newton',
                         DOF.evaluation='exact', lambda.selection.lossfunction='GCV')))
  
  
  load(file=paste0(foldername,"/test_5_4.RData"))
  expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);
  expect_equal( max(abs((sol$solution$beta-output_CPP$solution$beta))) < tol, TRUE);
  
  #### Test 5.5: Newton_fd method with exact GCV, default initial lambda and tolerance
  invisible(capture.output(sol<-smooth.FEM(observations=data, locations = NULL, 
                         covariates = cov1,
                         FEMbasis=FEMbasis, 
                         lambda.selection.criterion='newton_fd', DOF.evaluation='exact',  
                         lambda.selection.lossfunction='GCV')))
  
  
  load(file=paste0(foldername,"/test_5_5.RData"))
  expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);
  expect_equal( max(abs((sol$solution$beta-output_CPP$solution$beta))) < tol, TRUE);
  
  #### Test 5.7: Inference on beta and f, hypothesis testing and CI, Wald and Enhanced-ESF implementation
  inf_obj <- inferenceDataObjectBuilder(test = c("sim", "oat"), interval = "oat", component = c("both", "parametric"), type = c("w", "enh-esf"), dim = 2, n_cov = 1, beta0 = 1, f0 = f, n_flip = 150000, locations_by_nodes = T)
  
  invisible(capture.output(sol<-smooth.FEM(observations=data, locations = NULL, 
                         covariates = cov1,
                         FEMbasis=FEMbasis, 
                         lambda.selection.criterion='newton', DOF.evaluation='exact', 
                         lambda.selection.lossfunction='GCV', inference.data.object = inf_obj)))
  
  
  load(file=paste0(foldername,"/test_5_7.RData"))
  expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);
  expect_equal( max(abs((sol$solution$beta-output_CPP$solution$beta))) < tol, TRUE);
  expect_equal( max(abs((sol$inference$beta$p_values$wald[[1]]-
                        output_CPP$inference$beta$p_values$wald[[1]]))) < tol, TRUE);
  expect_equal( max(abs((sol$inference$f$p_values$wald[[1]]-
                        output_CPP$inference$f$p_values$wald[[1]]))) < tol, TRUE);
})
