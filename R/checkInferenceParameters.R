checkInferenceParameters <- function(R_Inference_Object_Data,checknumber){
  
  #################### Consistency Parameter Check #########################
  if(is.null(R_Inference_Object_Data)||R_Inference_Object_Data@definition==0){
    return(
      new("inferenceDataObject", test = as.integer(0), interval =as.integer(0), type = as.integer(0), exact = as.integer(0), dim = as.integer(0), 
          coeff = matrix(data=0, nrow = 1 ,ncol = 1), beta0 = -1, quantile = -1, n_perm = as.integer(1000), definition=as.integer(0)))
  }#define a shadow inferenceDataObject in order tell the cpp code not to perform any Inferential analysis.
    
  if(is.null(checknumber)&!is.null(R_Inference_Object_Data)){
    warning("Covariates are not defined, inference data are discarded")
    return(
      new("inferenceDataObject", test = as.integer(0), interval =as.integer(0), type = as.integer(0), exact = as.integer(0), dim = as.integer(0), 
          coeff = matrix(data=0, nrow = 1 ,ncol = 1), beta0 = -1, quantile = -1, n_perm = as.integer(1000), definition=as.integer(0))) 
  }
  
  #check consistence with covariates dimension.
  if(checknumber!=R_Inference_Object_Data@dim){
    stop("Inference data dimension and covariates dimenion are not consistent")
  }
  
  return(R_Inference_Object_Data)    #return the original object
}
  
  
  
  
  
  
  
  
  
  
