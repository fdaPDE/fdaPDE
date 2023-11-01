if(!dir.exists(test_path("../data/tSR-PDE")))
  dir.create(test_path("../data/tSR-PDE"))

if(!dir.exists(test_path("../data/tSR-PDE/test_13"))){
  dir.create(test_path("../data/tSR-PDE/test_13"))

  options(warn=-1)
  foldername <- test_path("../data/tSR-PDE/test_13/")
  data(horseshoe2D)
  
  mesh = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
  locations = refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
  mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
  
  FEMbasis = create.FEM.basis(mesh)
  
  f<-function(x,y,t)
  {
    K <- (y/0.1*as.double((abs(y)<=0.1 & x>-0.5))+as.double((abs(y)>0.1 | x<=-0.5)))^2
    res=numeric(length =length(x))
    for(i in 1:length(x))
    {
      if(x[i]>=0 && y[i]>0)
        res[i]=cos(t[i])*(0.25*pi+x[i])+(y[i]-0.5)^2
      if(x[i]>=0 && y[i]<=0)
        res[i]=cos(2*t[i])*(-0.25*pi-x[i])+(-y[i]-0.5)^2
      if(x[i]<0 && y[i]>0)
        res[i]=cos(t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
      if(x[i]<0 && y[i]<=0)
        res[i]=cos(2*t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
    }
    res
  }
  
  NumTimeInstants=5
  TimePoints=seq(0,pi,length.out =NumTimeInstants)
  
  space_time_locations = cbind(rep(TimePoints,each=nrow(locations)),rep(locations[,1],NumTimeInstants),rep(locations[,2],NumTimeInstants))
  sol_exact = f(space_time_locations[,2],space_time_locations[,3],space_time_locations[,1])
  
  ndata = length(sol_exact)
  
  # Create covariates
  set.seed(509875)
  cov1 = rnorm(ndata, mean = 1, sd = 2)
  
  # Add error to simulate data
  set.seed(7893475)
  data = sol_exact + 2*cov1 
  data = data + rnorm(ndata, mean = 0, sd =  0.05*diff(range(sol_exact)))
  observations = matrix(data,nrow(locations),NumTimeInstants)
  
  # Set smoothing parameter
  lambdaS = 10^-2
  lambdaT = 10^-2
  
  #### Test 13.1: Without GCV
  invisible(capture.output(sol_ref<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                              observations=observations, 
                              covariates = cov1,
                              FEMbasis=FEMbasis, lambdaS=lambdaS, 
                              lambdaT=lambdaT)))
  
  save(sol_ref, file=paste0(foldername,"/test_13_1.RData"))
 
  #### Test 13.2: exact GCV
  lambdaS = 10^(-1:1)
  lambdaT = 10^(-6:-4)
  invisible(capture.output(sol_ref<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                              observations=observations, 
                              covariates = cov1,
                              FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                              lambda.selection.criterion='grid', 
                              DOF.evaluation='exact', 
                              lambda.selection.lossfunction='GCV')))
  
  save(sol_ref, file=paste0(foldername,"/test_13_2.RData"))
  
  ### Inference
  # Inference obj:
  inf_obj <- inferenceDataObjectBuilder (test='oat', interval='oat',  dim=2, n_cov=1, type=c('w', 's', 'esf', 'enh-esf'), beta0 = 2, component='parametric', n_flip=1000, f_var=T)
  
  #### Test 13.3: overall inference on beta parameters, separable case
  lambdaS = 10^(-1:1)
  lambdaT = 10^(-6:-4)
  invisible(capture.output(sol_ref<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                              observations=observations, 
                              covariates = cov1,
                              FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                              lambda.selection.criterion='grid', DOF.evaluation='exact',
                              lambda.selection.lossfunction='GCV',
                              inference.data.object = inf_obj)))
  
  save(sol_ref, file=paste0(foldername,"/test_13_3.RData"))
  
  #### Test 13.4: overall inference on beta parameters, parabolic case
  lambdaS = 1 #10^(-1:1)
  lambdaT = 1e-6 #10^(-6:-4)
  invisible(capture.output(sol_ref<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                              observations=observations, 
                              covariates = cov1,
                              FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                              lambda.selection.criterion='grid', DOF.evaluation='exact', 
                              lambda.selection.lossfunction='GCV',
                              FLAG_PARABOLIC = TRUE,
                              inference.data.object = inf_obj)))
  
  save(sol_ref, file=paste0(foldername,"/test_13_4.RData"))
}
