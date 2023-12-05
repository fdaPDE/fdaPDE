
if(!dir.exists(test_path("../data/SR-PDE")))
  dir.create(test_path("../data/SR-PDE"))

if(!dir.exists(test_path("../data/SR-PDE/test_7"))){
  dir.create(test_path("../data/SR-PDE/test_7"))

  foldername = test_path("../data/SR-PDE/test_7/")  
  
  data(horseshoe2D)
  
  mesh = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
  locations = refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
  mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
  
  FEMbasis = create.FEM.basis(mesh)
  
  ndata = nrow(locations)
  
  # Create covariates
  set.seed(509875)
  cov1 = rnorm(ndata, mean = 1, sd = 2)
  cov2 = sin(locations[,1])
  
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test(locations[,1], locations[,2]) + 2*cov1 - cov2
  
  # Add error to simulate data
  set.seed(543663)
  ran = range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  # Set smoothing parameter
  lambda = 10^seq(-3,3,by=0.25)
  
  options(warn=-1)
  #### Test 7.1: grid with exact GCV
  invisible(capture.output(sol_ref<-smooth.FEM(locations = locations, observations=data, 
                                           covariates = cbind(cov1, cov2),
                                           FEMbasis=FEMbasis, lambda=lambda,
                                           lambda.selection.criterion='grid',
                                           DOF.evaluation='exact', 
                                           lambda.selection.lossfunction='GCV')))
  
  save(sol_ref,file=paste0(foldername,"/test_7_1.RData"))
  
  ### Test 7.2: Newton exact method with exact  GCV, default initial lambda and tolerance
  invisible(capture.output(sol_ref<-smooth.FEM(locations = locations, observations=data, 
                                           covariates = cbind(cov1, cov2),
                                           FEMbasis=FEMbasis,
                                           lambda.selection.criterion='newton', 
                                           DOF.evaluation='exact', 
                                           lambda.selection.lossfunction='GCV')))
  
  #save(output_CPP,file=paste0(foldername,"/test_7_2.RData"))
  save(sol_ref,file=paste0(foldername,"/test_7_2.RData"))
  
  ### Test 7.3: Inference on beta, hypothesis testing, Wald, Speckman, ESF, and enhanced ESF p_values
  inf_obj<-inferenceDataObjectBuilder(test = c("sim", rep("oat",3)), dim = 2, n_cov = 2, type = c("w", "s", "esf", "enh-esf"), beta0 = c(2,-1))
  invisible(capture.output(sol_ref<-smooth.FEM(locations = locations, observations=data, 
                                           covariates = cbind(cov1, cov2),
                                           FEMbasis=FEMbasis,
                                           lambda=lambda,
                                           lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', 
                                           inference.data.object=inf_obj)))

  save(sol_ref,file=paste0(foldername,"/test_7_3.RData"))
   
  ### Test 7.4: Inference on beta, hypothesis testing and confidence intervals of linear combinations, Wald and Speckman p_values  
  inf_obj<-inferenceDataObjectBuilder(test = "oat", interval = "oat", dim = 2, n_cov = 2, type = c("w", "s"), coeff = matrix(data = c(1,1,1,-1), nrow = 2, byrow = T))
  invisible(capture.output(sol_ref<-smooth.FEM(locations = locations, observations=data, 
                                           covariates = cbind(cov1, cov2),
                                           FEMbasis=FEMbasis, 
                                           lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', 
                                           inference.data.object=inf_obj)))
  
  save(sol_ref,file=paste0(foldername,"/test_7_4.RData"))
  
  ### Test 7.5: Inference on f, hypothesis testing, equality to f0, Wald, Sign-flip and ESF p_values
  inf_obj<-inferenceDataObjectBuilder(test = "sim", dim = 2, n_cov = 2, type = c("w","sf","esf"), component = "nonparametric", f0 = fs.test)
  invisible(capture.output(sol_ref<-smooth.FEM(locations = locations, observations=data, 
                                           covariates = cbind(cov1, cov2),
                                           FEMbasis=FEMbasis,
                                           lambda=lambda,
                                           lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', 
                                           inference.data.object=inf_obj)))
  
  #save(output_CPP,file=paste0(foldername,"/test_7_5.RData"))
  save(sol_ref,file=paste0(foldername,"/test_7_5.RData"))
  
  ### Test 7.6: Inference on f, hypothesis testing and confidence intervals, Wald with new locations  
  mesh_loc = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
  mesh_loc_ref<-refine.mesh.2D(mesh_loc, maximum_area = 0.05)
  
  new_locs <- mesh_loc_ref$nodes[which(mesh_loc_ref$nodesmarkers!=1),]
  new_locs[,1] <- new_locs[,1] + 0.2
  
  inf_obj<-inferenceDataObjectBuilder(test = "sim", interval = "oat", dim = 2, n_cov = 2, type = "w", component = "nonparametric", locations = new_locs)
  invisible(capture.output(sol_ref<-smooth.FEM(locations = locations, observations=data, 
                                           covariates = cbind(cov1, cov2),
                                           FEMbasis=FEMbasis,
                                           lambda.selection.criterion='newton', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', 
                                           inference.data.object=inf_obj)))
  
  #save(output_CPP,file=paste0(foldername,"/test_7_6.RData"))
  save(sol_ref,file=paste0(foldername,"/test_7_6.RData"))

  ### Test 7.7: Inference on both beta and f, hypothesis testing: all implementations p_values 
  inf_obj<-inferenceDataObjectBuilder(test = c("sim", "oat", "sim", "sim", "oat"), dim = 2, n_cov = 2, type = c("w","s","sf","esf","enh-esf"), 
                                      component = c("both", "parametric", "nonparametric", "both", "parametric"), f0 = fs.test, beta0 = c(2,-1))
  
  invisible(capture.output(sol_ref<-smooth.FEM(locations = locations, observations=data, 
                                           covariates = cbind(cov1, cov2),
                                           FEMbasis=FEMbasis, 
                                           lambda=lambda,
                                           lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV', 
                                           inference.data.object=inf_obj)))
  
  #save(output_CPP,file=paste0(foldername,"/test_7_7.RData"))
  save(sol_ref,file=paste0(foldername,"/test_7_7.RData"))
}
