checkInferenceParameters <- function(R_Inference_Object_Data,checknumber,locations,nodes){
  
  #if no inference is required, construct a dummy inference object, where all parameters are set to nonsense values on purpose
  if(is.null(R_Inference_Object_Data) || R_Inference_Object_Data@definition==0){
    return(
      new("inferenceDataObject", test = as.integer(0), interval = as.integer(0), type = as.integer(0), component = as.integer(0), exact = as.integer(0), enhanced = as.integer(0), dim = as.integer(0), n_cov = as.integer(0),
          locations = matrix(data=0, nrow = 1 ,ncol = 1), locations_indices = as.integer(0), locations_are_nodes = as.integer(0), coeff = matrix(data=0, nrow = 1 ,ncol = 1), beta0 = -1, f0 = function(){}, 
          f0_eval = -1, f_var = as.integer(0), quantile = -1, alpha = 0, n_flip = as.integer(1000), tol_fspai = -1, definition=as.integer(0)))
  }#define a shadow inferenceDataObject in order tell the cpp code not to perform any inferential analysis.
    
  if(!is.null(R_Inference_Object_Data) && is.null(checknumber) && sum(R_Inference_Object_Data@component!=2)!=0){
    warning("Covariates are not defined, inference data are discarded")
    return(
      new("inferenceDataObject", test = as.integer(0), interval = as.integer(0), type = as.integer(0), component = as.integer(0), exact = as.integer(0), enhanced = as.integer(0), dim = as.integer(0), n_cov = as.integer(0),
          locations = matrix(data=0, nrow = 1 ,ncol = 1), locations_indices = as.integer(0), locations_are_nodes = as.integer(0), coeff = matrix(data=0, nrow = 1 ,ncol = 1), beta0 = -1, f0 = function(){}, 
          f0_eval = -1, f_var = as.integer(0), quantile = -1, alpha = 0, n_flip = as.integer(1000), tol_fspai = -1, definition=as.integer(0))) 
  }
  
  #check consistency with covariates dimension (if inference on the linear component is required)
  if(sum(R_Inference_Object_Data@component!=2)!=0 && checknumber!=R_Inference_Object_Data@n_cov){
    stop("Inference data dimension and covariates dimension are not consistent")
  }
  
  #check consistency with locations and evaluate f0 (if inference on the nonparametric component is required)
  if(sum(R_Inference_Object_Data@component!=1)!=0){
    
    #if a vector of indices has been provided
    if(dim(R_Inference_Object_Data@locations)[1]==0){
      
      #check that there aren't indices repetitions and that each index doesn't exceed n_obs, eventually drop them 
      for(j in 1:length(R_Inference_Object_Data@locations_indices)){
        if(R_Inference_Object_Data@locations_indices[j] %in% R_Inference_Object_Data@locations_indices[-j])
          R_Inference_Object_Data@locations_indices = R_Inference_Object_Data@locations_indices[-j]
        else if(R_Inference_Object_Data@locations_indices[j] > dim(locations)[1]){
          warning("Removing a location index exceeding the number of observed locations")
          R_Inference_Object_Data@locations_indices = R_Inference_Object_Data@locations_indices[-j]
        }
      }
      
      #check that the final length of locations indices doesn't exceed the observed ones
      if(length(R_Inference_Object_Data@locations_indices) > dim(locations)[1]){
        warning("Length of 'location_indices' exceeds the number of observed locations, proceeding with the observed locations")
        R_Inference_Object_Data@location_indices = as.integer(seq(1,dim(locations)[1]))
        R_Inference_Object_Data@locations = locations
      }
      else{
        R_Inference_Object_Data@locations = matrix(data = locations[R_Inference_Object_Data@locations_indices,], ncol = R_Inference_Object_Data@dim)
      }
    }
    else{
    # a matrix of locations has been provided  
    if(ncol(R_Inference_Object_Data@locations)==1){
      R_Inference_Object_Data@locations = locations #default case
      R_Inference_Object_Data@locations_indices = as.integer(seq(1,dim(locations)[1]))
      }
    else{
      # only wald implementation can enter here --> in principle new location points are allowed here, so locations_indices is undefined in this case
      # check that the dimension of the problem is consistent
      if(ncol(R_Inference_Object_Data@locations) != ncol(locations))
        stop("Inference data dimension and locations dimension are not consistent")
      # locations_indices is set to a nonsense value
      R_Inference_Object_Data@locations_indices = as.integer(0)
    
      # EVENTUALLY ADD FURTHER CHECKS: are locations inside the domain??
      }
    }
    
    # in wald implementation check that the final number of chosen locations doesn't exceed n-q (to ensure var-cov matrix is invertible)
    if(sum(R_Inference_Object_Data@type==1)!=0){
      if(is.null(checknumber)){
        if(dim(R_Inference_Object_Data@locations)[1] > dim(locations)[1])
          stop("Number of chosen locations is too high, variance-covariance matrix of the estimator is not invertible: 
             decrease the number of locations or choose eigen-sign-flip inferential approach")
      }
      else{
        if(dim(R_Inference_Object_Data@locations)[1] > dim(locations)[1] - checknumber)
          stop("Number of chosen locations is too high, variance-covariance matrix of the estimator is not invertible: 
             decrease the number of locations or choose eigen-sign-flip inferential approach")
      }
    }
    
    locs_are_nodes_before_check = (R_Inference_Object_Data@locations_are_nodes == 1)
    R_Inference_Object_Data@locations_are_nodes = as.integer(1)
    
    # check if the selected locations coincide with nodes
    end_loop = FALSE
    j=1
    while(!end_loop && j <= dim(R_Inference_Object_Data@locations)[1]){
      node_found = FALSE
      l = 1
      while(!node_found && l <= dim(nodes)[1]){
        if(R_Inference_Object_Data@locations[j,1]==nodes[l,1] && R_Inference_Object_Data@locations[j,2]==nodes[l,2]){
          node_found = TRUE
        }
        l = l+1 
      }
      if(!node_found){
        end_loop = TRUE
        R_Inference_Object_Data@locations_are_nodes = as.integer(2)
      }
      j = j+1
    }
    
    if(!locs_are_nodes_before_check && R_Inference_Object_Data@locations_are_nodes!=2 && sum(R_Inference_Object_Data@type==1)!=length(R_Inference_Object_Data@type))
      warning("All the locations selected for inference on the nonparametric component coincide with the nodes. 
              Sign-Flip and Eigen-Sign-Flip tests are performed by area, combining neighborhoods of locations according to the distance induced by the mesh")
    
    if(R_Inference_Object_Data@locations_are_nodes!=1){
      # intervals with sign-flip or eigne-sign-flip implementation are not allowed
      for(k in 1:length(R_Inference_Object_Data@type)){
        if((R_Inference_Object_Data@type[k]==3 || R_Inference_Object_Data@type[k]==4) && R_Inference_Object_Data@interval[k]!=0 && R_Inference_Object_Data@component[k]!=1)
          stop("Sign-Flip and Eigen-Sign-Flip confidence intervals on f are implemented only if the selected locations are a subset of the mesh nodes")
      }
    }
    # finally evaluate f0 at the chosen locations
    f0 <- R_Inference_Object_Data@f0
    dim <- R_Inference_Object_Data@dim
    locs <- R_Inference_Object_Data@locations
    
    if(dim == 2){
      for(i in 1:dim(locs)[1])
        R_Inference_Object_Data@f0_eval <- c(R_Inference_Object_Data@f0_eval, f0(locs[i,1], locs[i,2])) 
    }
    else{
      for(i in 1:dim(locs)[1])
        R_Inference_Object_Data@f0_eval <- c(R_Inference_Object_Data@f0_eval, f0(locs[i,1], locs[i,2], locs[i,3])) 
    }
  }
  else{
    # in case only inference on parametric component is requested, set f0_eval and locations with default values, just for transmission safety
    R_Inference_Object_Data@locations <- matrix(data = 0)
    R_Inference_Object_Data@locations_indices <- as.integer(0)
    R_Inference_Object_Data@f0_eval <- 0
    
  }
  
  return(R_Inference_Object_Data)    #return the original object
}
  
  
  
  
  
  
  
  
  
  
