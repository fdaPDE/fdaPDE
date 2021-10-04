#' Class for inference data
#'
#'@slot test A vector of integers taking value 0, 1 or 2; in the first case no test is performed, in the second one-at-the-time tests are performed,
#'in the third a simultaneous test is performed.
#'@slot interval A vector of integers taking value 0, 1, 2 or 3; In the first case no confidence interval is computed, in the second case one-at-the-time confidence intervals are computed, 
#'in the third case simultaneous confidence intervals are computed, in the fourth case Bonferroni confidence intervals are computed.
#'@slot type A vector of integers taking value 1, 2 or 3, corresponding to Wald, speckman or eigen-sign-flip implementation 
#'of the inference analysis.
#'@slot exact An integer taking value 1 or 2. If 1 an exact computation of the test statistics will be performed,
#'whereas if 2 an approximated computation will be carried out.
#'@slot dim Number of covariates taken into account in the linear part of the regression problem.
#'@slot coeff A matrix of numeric coefficients with columns of dimension \code{dim} and each row represents a linear combination of the linear parameters to be tested and/or to be
#' estimated via confidence interval. 
#'@slot beta0 Vector of null hypothesis values for the linear parameters of the model. Used only if \code{test} is not 0.
#'@slot f_var An integer taking value 1 or 2. If 1 local variance estimates for the nonlinear part of the model will be computed,
#'whereas if 2 they will not.
#'@slot quantile Quantile needed for confidence intervals. Used only if interval is not 0.
#'@slot n_flip An integer representing the number of permutations in the case of eigen-sign-flip test.
#'@slot tol_fspai A real number greater than 0 specifying the tolerance for FSPAI algorithm, in case of non-exact inference.
#'@slot definition An integer taking value 0 or 1. If set to 1, the class will be considered as created by the function [inferenceDataObjectBuilder],
#'leading to avoid some of the checks that are performed on inference data within smoothing functions.
#'
#'@description
#' A class that contains all possible information for inference over linear parameters in spatial regression with
#' differential regularization problem. This object can be used as parameter in smoothing function of the fdaPDE library [smooth.FEM].
#' 
#'@details #Warning
#'At least one between test and interval must be nonzero. \code{dim}, \code{coeff} and \code{beta0} need to be coherent.
#'The usage of [inferenceDataObjectBuilder] is recommended for the construction of an object of this class.

inferenceDataObject<-setClass("inferenceDataObject", slots = list(test = "integer", 
                                                                  interval = "integer", 
                                                                  type = "integer", 
                                                                  exact = "integer", 
                                                                  dim = "integer",
                                                                  coeff = "matrix",
                                                                  beta0 = "numeric",
                                                                  f_var = "integer",
                                                                  quantile = "numeric",
                                                                  n_flip = "integer",
                                                                  tol_fspai = "numeric",
                                                                  definition="integer")
                              )

#'Constructor for inferenceDataObject class
#'
#'@param test A list of strings defining the type of test to be performed, with the same length of \code{interval} and \code{type} lists. The default is NULL, and can take values 'one-at-the-time' 
#'or 'simultaneous'. If the value is NULL, no test is performed, and the \code{interval} parameter in the same position of the list needs to be not NULL. If it takes value
#''one-at-the-time', the parameters \code{beta0}, \code{dim}, \code{coeff} are taken into account, and one-at-the-time tests will be performed. If it takes value 
#''simultaneous', a global simultaneous test will be performed.
#'@param interval A list of strings defining the type of confidence intervals to be computed, with the same length of \code{test} and \code{type} lists. The default is NULL,
#' and can take values 'one-at-the-time' 'simultaneous' and 'bonferroni'. If the value is NULL, no interval will be computed, and the \code{test} parameter in the corresponding position needs to be set.
#'  Otherwise one at the time, simultaneous or Bonferroni correction intervals will be computed. If it is not NULL, the parameter \code{level} will be taken into account. Confidence intervals can be 
#'  computed only in Wald or Speckman implementation.
#'@param type A list of strings defining the type of implementation for the inferential analysis, with the same length of \code{test} list and \code{interval} list . The possible values are three:
#''wald'(default), 'speckman' or 'eigen-sign-flip'. If 'eigen-sign-flip' is set, the corresponding \code{interval} position needs to be NULL.
#'@param exact A logical used to decide the method used to estimate the statistics variance.
#'The possible values are: FALSE (default) and TRUE. In the first case an approximate method is used, leading to a lower accuracy, but faster computation.
#'In the second case the evaluation is exact but computationally expensive.
#'@param dim Number of the covariates, defaulted to NULL. (Must be set by the user)
#'@param coeff A matrix, with \code{dim} number of columns, of numeric coefficients, defaulted to NULL. If this parameter is NULL,
#'in the corresponding inferenceDataObject it will be defaulted to an identity matrix. If at least one 'eigen-sign-flip' value is present in \code{type}, needs to be an identity matrix.
#'@param beta0 Vector of real numbers (default NULL). It is used only if the \code{test} parameter is set, and has length the number of rows of matrix \code{coeff}. 
#'If \code{test} is set and \code{beta0} is NULL, will be set to a vector of zeros.
#'@param f_var A logical used to decide whether to estimate the local variance of the nonlinear part of the model.
#'The possible values are: FALSE (default) and TRUE. 
#'@param level Level of significance used to compute quantiles for confidence intervals, defaulted to 0.95. It is taken into account only if \code{interval} is set.
#'@param n_flip Number of flips performed in Eigen-Sign-Flip test, defaulted to 1000. It is taken into account only if at least one position of \code{type} is set to 'eigen-sign-flip'.
#'@param tol_fspai Tolerance for FSPAI algorithm taking value greater than 0, defaulted to 0.05. It is taken into account only if \code{exact} is set to 'False'. The lower is the tolerance, the heavier is the computation.
#'@return The output is a well defined [inferenceDataObject], that can be used as parameter in the [smooth.FEM] function.
#'@description A function that build an [inferenceDataObject]. In the process of construction many checks over the input parameters are carried out so that the output is a well defined object,
#'that can be used as parameter in [smooth.FEM] function. Notice that this constructor ensures well-posedness of the object, but a further check on consistency with smooth.FEM parameters will be carried out inside that function.
#'
#'@usage inferenceDataObjectBuilder<-function(test = NULL, 
#'interval = NULL, 
#'type = 'wald', 
#'exact = FALSE, 
#'dim, 
#'coeff = NULL, 
#'beta0 = NULL, 
#'f_var = FALSE,
#'level = 0.95,
#'n_flip = 1000,
#'tol_fspai = 0.05)
#' @export
#' 
#' 
#' @examples 
#' obj1<-inferenceDataObjectBuilder(test = "simultaneous", interval = NULL, exact = T, dim = 4);
#' obj2<-inferenceDataObjectBuilder(interval = "one-at-the-time", dim = 5, level = 0.99);
#' obj3<-inferenceDataObjectBuilder(test=c('one-at-the-time', 'simultaneous', 'one-at-the-time','none'), interval=c('bonferroni','one-at-the-time','none','simultaneous'),
#'  type=c('wald','speckman','eigen-sign-flip','speckman'),exact=TRUE, dim=2, level=0.99)


inferenceDataObjectBuilder<-function(test = NULL, 
                                interval = NULL, 
                                type = "wald", 
                                exact = F, 
                                dim = NULL, 
                                coeff = NULL, 
                                beta0 = NULL,
                                f_var = F,
                                level = 0.95,
                                n_flip = 1000,
                                tol_fspai = 0.05){
  
  # Preliminary check of parameters input types, translation into numeric representation of default occurrences.
  if(!is.null(test)){
    if(class(test)!="character")
      stop("'test'should be a character: choose one between 'one-at-the-time', 'simultaneous', 'none'")
    if(length(test)==0)
      stop("'test' is zero dimensional, should be a vector taking values among 'ont-at-the-time', 'simultaneous', 'none'")
  }else{test_numeric=as.integer(0)}
  
  if(!is.null(interval)){
    if(class(interval)!="character")
      stop("'interval'should be a character: choose one between 'one-at-the-time' , 'simultaneous', 'bonferroni', 'none'")
    if(length(interval)==0)
      stop("'interval' is zero dimensional, should be a vector taking values among 'one-at-the-time', 'simultaneous', 'bonferroni', 'none'")
  }else{interval_numeric=as.integer(0)}
  
  if(!is.null(type)){
    if(class(type)!="character")
      stop("'type' should be a vector of characters: choose among 'wald', 'speckman' or 'eigen-sign-flip'" )
    if(length(type)==0)
      stop("'type' is zero dimensional, should be a vector taking values among 'wald', 'speckman' or 'eigen-sign-flip'")
  }
  
  if(exact != F){
    if(class(exact)!="logical")
      stop("'exact' should be either TRUE or FALSE ")
    if(length(exact)==0)
      stop("'exact' is zero dimensional, should be either TRUE or FALSE")
    if(exact!=T)
      stop("'exact' should be either TRUE or FALSE")
    exact_numeric=as.integer(1)
  }else{
    exact_numeric=as.integer(2)
  }
  
  if(!is.null(dim)){
    if(class(dim)!="numeric" && class(dim)!="integer")
      stop("'dim' should be an integer or convertible to integer type")
    dim=as.integer(dim)
  }
  
  if(!is.null(coeff)){
    if(class(coeff)[1]!="matrix")
      stop("'coeff' should be of class matrix")
  }
  
  if(!is.null(beta0)){
    if(class(beta0)!="numeric")
      stop("'beta0' should be numeric")
    if(length(beta0)==0)
      stop("'beta0' is zerodimensional")
  }
  
  if(f_var != F){
    if(class(f_var)!="logical")
      stop("'f_var' should be either TRUE or FALSE ")
    if(length(f_var)==0)
      stop("'f_var' is zero dimensional, should be either TRUE or FALSE")
    if(f_var!=T)
      stop("'f_var' should be either TRUE or FALSE")
    f_var_numeric=as.integer(1)
  }else{
    f_var_numeric=as.integer(2)
  }
  
  if(level!=0.95){
    if(class(level)!="numeric")
      stop("'level' should be numeric")
    if(length(level)==0)
      stop("'level' is zerodimensional, should be a positive number between 0 and 1")
  }
  
  if(!is.null(n_flip)){
    if(class(n_flip)!="numeric" && class(n_flip)!="integer")
      stop("'n_flip' should be an integer or convertible to integer type")
    n_flip=as.integer(n_flip)
  }
  
  if(tol_fspai!=0.05){
    if(class(tol_fspai)!="numeric")
      stop("'tol_fspai' should be numeric")
    if(length(tol_fspai)==0)
      stop("'tol_fspai' is zerodimensional, should be a positive number between 0 and 1")
  }
  
  # Check of consistency of parameters. Translation into numeric representation. The checks are repeated for each element of the vectors test, interval and type
  
  if(is.null(test) && is.null(interval))                                        # However, they can be both set
    stop("at least one between 'test' and 'interval' should be not NULL")
  
  if(is.null(dim) || dim <= 0)                                                  # Otherwise the function won't be able to default the coeff neither object nor beta0
    stop("number of covariates is needed")
  else{
    if(is.null(coeff)){                                                         # If it is left as NULL, all the parameters are individually taken into account without any linear combination.
      coeff = diag(1, nrow=dim, ncol=dim)
      }
    else{
      if(dim(coeff)[2]!=dim)
        stop("number of covariates and coefficients do not match")
      for (i in 1:dim(coeff)[1]){
        for(j in 1:dim(coeff)[2]){
        if(class(coeff[i,j])!="numeric")
          stop("'coeff' should be composed by numeric values")
        }
      }
      rm(list = ("i"))
      rm(list = ("j"))
    }
  }
  
  # Consistency check for number of imlementations required (default values assigned in case of NULL)
  n_of_implementations=length(type);
  if(is.null(test)){
    if(length(interval)!=n_of_implementations)
      stop("Intervals vector length is not consistent whith the number of implementation provided")
    test=rep('none',n_of_implementations)
  }else{
    if(is.null(interval)){
      if(length(test)!=n_of_implementations)
        stop("Test vector length is not consistent whith the number of implementation provided")
        interval=rep('none',n_of_implementations)
    }else{
      if(length(interval)!=n_of_implementations)
        stop("Intervals vector length is not consistent whith the number of implementation provided")
      if(length(test)!=n_of_implementations)
        stop("Test vector length is not consistent whith the number of implementation provided")
    }
  }
  
  # Preallocation of space
  test_numeric=rep(0, n_of_implementations)
  interval_numeric=rep(0, n_of_implementations)
  type_numeric=rep(0, n_of_implementations)
  quantile=rep(0, n_of_implementations)
  
  for (index in 1:n_of_implementations){
    
      if(type[index]!="wald" && type[index]!="speckman" && type[index]!="eigen-sign-flip"){
        stop("type should be choosen between 'wald', 'speckman' and 'eigen-sign-flip'")}else{
          if(type[index]=="wald") type_numeric[index]=as.integer(1)
          if(type[index]=="speckman") type_numeric[index]=as.integer(2)
          if(type[index]=="eigen-sign-flip") type_numeric[index]=as.integer(3)
        }
      
      if(test[index]=='none' && interval[index]=='none'){
        stop("At least one between test and interval should be required for each implementation in")
      }
    
      if(test[index]!="one-at-the-time" && test[index]!="simultaneous" &&  test[index]!="none"){
        stop("test should be either 'one-at-the-time', 'simultaneous' or 'none'")}else{
          if(test[index]=="none") test_numeric[index]=as.integer(0)
          if(test[index]=="one-at-the-time") test_numeric[index]=as.integer(1)
          if(test[index]=="simultaneous") test_numeric[index]=as.integer(2)
       }
      
      if(interval[index]!="none" && (type[index]=="wald" || type[index]=="speckman")){
        if(interval[index]!="one-at-the-time" && interval[index]!="simultaneous" && interval[index]!="bonferroni" && interval[index]!="none"){
          stop("interval should be either 'one-at-the-time' 'simultaneous', 'bonferroni' or 'none'")}else{
            if(interval[index]=="none") interval_numeric[index]=as.integer(0)
            if(interval[index]=="one-at-the-time") interval_numeric[index]=as.integer(1)
            if(interval[index]=="simultaneous") interval_numeric[index]=as.integer(2)
            if(interval[index]=="bonferroni") interval_numeric[index]=as.integer(3)
          }
        if(level <= 0 || level > 1)                                                
          stop("level should be a positive value smaller or equal to 1")
      }
      
      if(interval[index]!="none" && type[index]=="eigen-sign-flip"){
        stop("confidence intervals are not implemented in the eigen-sign-flip case")
      }
      
      if(type[index] == "eigen-sign-flip"){
        for(i in 1:dim(coeff)[1]){
          count=0
          for(j in 1:dim(coeff)[2]){
            if(coeff[i,j]!= 0 && coeff[i,j]!=1)
              stop("linear combinations are not allowed in the eigen-sign-flip case")
            count = count + coeff[i,j]
          }
          if(count != 1)
            stop("linear combinations are not allowed in the eigen-sign-flip case")
        }
      }
      
      # Well posedeness check for coeff in simultaneous case;
      if(test[index]=="simultaneous"){
        if(det(coeff %*% t(coeff)) < 0.001){
          stop("coeff is not full rank, variance-covariance matrix of the linear combination not invertible")
        }
      }
        
      if(interval[index]=="simultaneous"){
        if(det(coeff %*% t(coeff)) < 0.001){
          stop("coeff is not full rank, variance-covariance matrix of the linear combination not invertible")
        }
      }
      
      # Build the quantile for Confidence intervals if needed
      if(interval[index]=="none"){
         quantile[index]=0
        }else{
          if(level > 0){
              if(interval[index] == "simultaneous"){ # Simultaneous CI -> Chi-Squared (q) law for statistic
                quantile[index]=qchisq(level, dim(coeff)[1])
              }else{
                if(interval[index]=="one-at-the-time"){  # One at the time CI (each interval has confidence alpha) -> Gaussian law for statistic
                  quantile[index]=qnorm((1-(1-level)/2),0,1)
              }else{ 
                  if(interval[index]== "bonferroni"){# One at the time CI (overall confidence alpha) -> Gaussian law for statistic
                  quantile[index]=qnorm((1-(1-level)/(2*dim(coeff)[1])),0,1)
                  }
              }
            }
          }
      }
      
  } # End of for cycle over the different implementation
  
  if(is.null(beta0))                 # If it left to NULL, beta0 is set to a vector of zeros.
    beta0<-rep(0, dim(coeff)[1])
  else{
    if(length(beta0)!=dim(coeff)[1])
      stop("dimension of 'coeff' and 'beta0' are not consistent")
  }
  
  if(!is.null(n_flip)){
    if(n_flip <= 0)                                                
      stop("number of sign-flips must be a positive value")
  }else {
    n_flip <- as.integer(1000)
  }
  
  if(exact==FALSE){
    if(tol_fspai <= 0 )                                                
      stop("tol_fspai should be a positive value")
  }
  
  
  
  definition=as.integer(1)
  
  # Building the output object, returning it
  result<-new("inferenceDataObject", test = as.integer(test_numeric), interval = as.integer(interval_numeric), type = as.integer(type_numeric), exact = exact_numeric, dim = dim, 
              coeff = coeff, beta0 = beta0, f_var = f_var_numeric, quantile = quantile, n_flip = n_flip, tol_fspai = tol_fspai, definition=definition)
  
  return(result)
}

