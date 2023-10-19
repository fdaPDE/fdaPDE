  options(warn=-1)
  foldername <- test_path("../data/tSR-PDE/test_15/")
  
  set.seed(0)
  
  # Build mesh: Sphere
  data(sphere3Ddata)
  sphere3D<-create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)
  
  FEMbasis <- create.FEM.basis(sphere3D)
  nodesLocations=sphere3D$nodes
  nnodes = nrow(sphere3D$nodes)
  TimeLocations = seq(0,1,length.out = 5)
  Locations = cbind(rep(TimeLocations,each=nnodes),rep(nodesLocations[,1],length(TimeLocations)),rep(nodesLocations[,2],length(TimeLocations)),rep(nodesLocations[,3],length(TimeLocations)))
  
  # Exact test function
  a1 = rnorm(1,mean = 1, sd = 1)
  a2 = rnorm(1,mean = 1, sd = 1)
  a3 = rnorm(1,mean = 1, sd = 1)
  a4 = rnorm(1,mean = 1, sd = 1)
  
  func = function(x)
  {
    a1*sin(2*pi*(x[,1]*x[,2]))+a2*cos(2*pi*x[,2])+a3*cos(2*pi*x[,3])+a4*sin(2*pi*x[,4])
  }
  
  func_evaluation = func(Locations)
  ran=range(func_evaluation)
  
  # Generate locations
  nloc = 1000
  loc=matrix(data=runif(3*nloc, min=-1,max=1),nrow=nloc,ncol=3,byrow=T)
  
  ind=NULL
  for(row in 1:nloc){
    normvec = (loc[row,1]^2+loc[row,2]^2+loc[row,3]^2)
    if(normvec>0.975)   # check points outside the sphere and remove them
      ind = c(ind,row)
  }
  
  loc=loc[-ind,]
  nloc=dim(loc)[1]
  timeloc = seq(0,1,length.out=5)
  loc = cbind(rep(timeloc,each=nloc),rep(loc[,1],length(timeloc)),rep(loc[,2],length(timeloc)),rep(loc[,3],length(timeloc)))
  
  
  # Exact test function - locations different from nodes
  func_evaluation2=func(loc)
  
  cov1=(4*sin(2*pi*Locations[,2])+6*sin((2*pi*Locations[,3])^2))*(1-exp(-Locations[,1]))/3
  cov2=cos(-2*pi*Locations[,4])+2*Locations[,1]*sin(2*pi*Locations[,2])/6
  
  cov1_nonod=(4*sin(2*pi*loc[,2])+6*sin((2*pi*loc[,3])^2))*(1-exp(-loc[,1]))/3
  cov2_nonod=cos(-2*pi*loc[,4])+2*loc[,1]*sin(2*pi*loc[,2])/6
  
  W=cbind(cov1,cov2)
  W2=cbind(cov1_nonod,cov2_nonod)
  
  lambdaS=10^seq(-5.0, -4.0, 0.25)
  lambdaT=10^seq(-1.5, -0.5, 0.25)
  
  lambdaS2=10^seq(-5.5, -4.5, 0.25)
  lambdaT2=10^seq(-1.5, -0.5, 0.25)
  
  lambdaS_par=10^seq(-4.8, -4.4, 0.1)
  lambdaT_par=10^seq(1.4, 1.8, 0.1)
  
  lambdaS_par2=10^seq(-4.4, -4.0, 0.1)
  lambdaT_par2=10^seq(1.4, 1.8, 0.1)
  
  beta_exact= c(0.7,2.0)
  ran = range(func_evaluation)
  data = func_evaluation +rnorm(length(func_evaluation),mean=0,sd=0.05*(ran[2]-ran[1]))
  
  ran = range(func_evaluation2)
  data_noloc = func_evaluation2 +rnorm(length(func_evaluation2),mean=0,sd=0.05*(ran[2]-ran[1]))
  
  ran = range(func_evaluation+ W%*%beta_exact)
  datacov=func_evaluation+ W%*%beta_exact +rnorm(length(func_evaluation),mean=0,sd=0.05*(ran[2]-ran[1]))
  
  data = matrix(data,nnodes,length(TimeLocations))
  data_noloc = matrix(data_noloc,nloc,length(timeloc))
  datacov = matrix(datacov,nnodes,length(TimeLocations))
  ###########################SEPARABLE###########################################
  
  invisible(capture.output(sol_ref <- smooth.FEM.time(observations = data,time_mesh = TimeLocations,
                           FEMbasis = FEMbasis, lambdaS = lambdaS, lambdaT = lambdaT)))
  save(sol_ref, file=paste0(foldername,"/test_15_1.RData"))
  
  invisible(capture.output(sol_ref <- smooth.FEM.time(locations=loc[1:nloc,2:4],time_locations = timeloc,
                                  observations = data_noloc,
                                  time_mesh = timeloc,FEMbasis = FEMbasis, 
                                  lambdaS = lambdaS2, lambdaT = lambdaT2)))
  
  save(sol_ref, file=paste0(foldername,"/test_15_2.RData"))
  
  invisible(capture.output(sol_ref <- smooth.FEM.time(observations = datacov,time_mesh = TimeLocations, covariates = W,
                              FEMbasis = FEMbasis, lambdaS = lambdaS, lambdaT = lambdaT)))
  
  save(sol_ref, file=paste0(foldername,"/test_15_3.RData"))
  
  ##########################################PARABOLIC####################################################
  ### MONOLITIC METHOD
  invisible(capture.output(sol_ref <- smooth.FEM.time(observations = data,time_mesh = TimeLocations,
                           FEMbasis = FEMbasis, lambdaS = lambdaS_par, lambdaT = lambdaT_par, FLAG_PARABOLIC = TRUE)))
  
  save(sol_ref, file=paste0(foldername,"/test_15_4.RData"))
 
  invisible(capture.output(sol_ref <- smooth.FEM.time(locations=loc[1:nloc,2:4],observations = data_noloc,time_mesh = timeloc,
                                  FEMbasis = FEMbasis, lambdaS = lambdaS_par2, lambdaT = lambdaT_par2, FLAG_PARABOLIC = TRUE)))
  
  save(sol_ref, file=paste0(foldername,"/test_15_5.RData"))
  
  invisible(capture.output(sol_ref <- smooth.FEM.time(observations = datacov[,2:length(TimeLocations)],time_mesh = TimeLocations, 
                                                         covariates = W[(1+nnodes):(length(TimeLocations)*nnodes),],
                              FEMbasis = FEMbasis, lambdaS = lambdaS_par, lambdaT = lambdaT_par, FLAG_PARABOLIC = TRUE,
                              IC=func_evaluation[1:nnodes])))
  
  save(sol_ref, file=paste0(foldername,"/test_15_6.RData"))