checkInferenceData <- function(R_Inference_Object_Data,checknumber){
  
  #################### Consistency Parameter Check #########################
  if(is.null(R_Inference_Object_Data)|R_Inference_Object_Data@definition==0){
    return(
      new("inferenceDataObject", test = as.integer(0), interval =as.integer(0), type = as.integer(0), exact = as.integer(0), dim = as.integer(0), 
                  coeff = as.integer(-1), beta0 = -1, beta1 = -1, level = -1,definition=as.integer(0)))
  }                                                                             #define a shadow inferenceDataObject in order, so to tell the cpp code not to perform any Inferential analysis
    
  if(is.null(checknumber)&!is.null(R_Inference_Object_Data)){
    warning("Covariates are not defined, inference data are discarded")
    return(
      new("inferenceDataObject", test = as.integer(0), interval =as.integer(0), type = as.integer(0), exact = as.integer(0), dim = as.integer(0), 
                  coeff = as.integer(-1), beta0 = -1, beta1 = -1, level = -1,definition=as.integer(0))) #check consistence with covariates dimension.
  }
  
  if(checknumber!=R_Inference_Object_Data@dim){
    stop("Inference data dimension and covariates dimenion are not consistent")
  }
  
  #Since dimensions are consistent we retrieve the original object in order to take advantage of previously implemented check in inference data constructor
  if(R_Inference_Object_Data@test!=1){                                          #If definition is set to 1, I give for grated that  the object has been created by means of inferenceDataObjectBuilder function. No more checks needed.
  #Test
  if(R_Inference_Object_Data@test==0){
    New_test=NULL
  }else if(R_Inference_Object_Data@test==1){
    New_test="pvalue"
  }else if(R_Inference_Object_Data@test==2){
    New_test="power"
  }else {New_test="other"}
  
  #Interval
  if(R_Inference_Object_Data@interval==0){
    New_interval=NULL
  }else if(R_Inference_Object_Data@interval==1){
    New_interval="one-at-the-time"
  }else if(R_Inference_Object_Data@interval==2){
    New_interval="bonferroni"
  }else {New_interval="other"}
  
  #Type
  if(R_Inference_Object_Data@type==1){
    New_type="wald"
  }else if(R_Inference_Object_Data@type==2){
    New_type="sandwich"
  }else if(R_Inference_Object_Data@type==3){
    New_type="permutational"
  }else{New_type="other"}
 
  #exact
  if(R_Inference_Object_Data@exact==1){
    New_exact="True"
  }else if(R_Inference_Object_Data@exact==2){
    New_exact="False"
  }else {New_exact="other"}
  #None of the others needs to be modified
  #Build the object to be returned, checks carried out by the inferenceDataObjectBuilderr function
  result<-inferenceDataObjectBuilder(test = New_test,
                                     interval = New_interval,
                                     type=New_type,
                                     exact=New_exact,
                                     dim=R_Inference_Object_Data@dim,
                                     coeff=R_Inference_Object_Data@coeff,
                                     beta0=R_Inference_Object_Data@beta0,
                                     beta1=R_Inference_Object_Data@beta1,
                                     level=R_Inference_Object_Data@level)
  return(result)                       #return (if well defined) the object
  }else{
    return(R_Inference_Object_Data)    #return the original object
  }
}
  
  
  
  
  
  
  
  
  
  
