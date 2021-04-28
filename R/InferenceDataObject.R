#' Class for inference data
#'
#'@slot test An integer taking value 0, 1 or 2; in the first case no test is performed, in the second one-at-the-time tests are performed,
#'in the third a simultaneous test is performed.
#'@slot interval An integer taking value 0, 1, 2 or 3; In the first case no confidence interval is computed, in the second case one-at-the-time confidence intervals are computed, 
#'in the third case simultaneous confidence intervals are computed, in the fourth case Bonferroni confidence intervals are computed.
#'@slot type An integer taking value 1, 2 or 3, corresponding to Wald, speckman or permutational implementation 
#'of the inference analysis.
#'@slot exact An integer taking value 1 or 2. If 1 an exact computation of the test statistic will be performed,
#'whereas if 2 an approximated computation will be carried out.
#'@slot dim Number of covariates taken into account in the linear part of the regression problem.
#'@slot coeff A matrix of numeric coefficients with columns of dimension \code{dim} and each row represents a linear combination of the linear parameters to be tested and/or to be
#' estimated via confidence interval. 
#'@slot beta0 Vector of null hypothesis values for the linear parameters of the model. Used only if \code{test} is not 0.
#'@slot quantile Quantile needed for confidence intervals. Used only if interval is not 0.
#'@slot n_perm An integer representing the number of permutations in the case of permutational test.
#'@slot definition An integer taking value 0 or 1. If set to 1, the class will be considered as created by the function [inferenceDataObjectBuilder],
#'leading to avoid some of the checks that are performed on inference data within smoothing functions.
#'
#'@description
#' A class that contains all possible information for linear parameters in spatial regression with
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
                                                                  quantile = "numeric",
                                                                  n_perm = "integer",
                                                                  definition="integer")
                              )

#'Constructor for inferenceDataObject class
#'
#'@param test A string defining the type of test to be performed. The default is NULL, and can take values 'one-at-the-time' or 'simultaneous'.
#'If the value is NULL, no test is performed, and the \code{interval} parameter need to be not NULL. If it takes value
#''one-at-the-time', the parameters \code{beta0}, \code{dim}, \code{coeff} are taken into account, and one-at-the-time tests will be performed. If it takes value 
#''simultaneous', a global simultaneous test will be performed.
#'@param interval A string defining the type of confidence intervals to be computed. The default is NULL, and can take value 'one-at-the-time' 'simultaneous' and 'bonferroni'.
#'If the value is NULL, no interval will be computed, and the \code{test} parameter needs to be set. Otherwise one at the time, simultaneous or Bonferroni correction intervals will be computed.
#'If it is not NULL, the parameter \code{level} will be taken into account. Up to now, confidence intervals can be computed only in the Wald implementation.
#'@param type A string defining the type of implementation for the inferential analysis. The possible values are three:
#''wald'(default), 'speckman' or 'permutational', corresponding to the three possible methods developed in Ferraccioli....
#'@param exact A string used to decide the method used to estimate the statistics variance.
#'The possible values are: 'True' and 'False'(default). In the first case the evaluation is exact but computationally very expensive.
#'In the second case an approximate method is used, leading to a lower accuracy, but faster computation.
#'@param dim Number of the covariates, defaulted to NULL.
#'@param coeff A matrix, with \code{dim} number of columns, of numeric coefficients, defaulted to NULL. If this parameter is NULL, it will be defaulted to an identity matrix. 
#'@param beta0 Vector of real numbers (default NULL). It is used only if the \code{test} parameter is set, and has length the number of rows of matrix \code{coeff}. If \code{test} is set and \code{beta0} is NULL,
#'will be set to a vector of zeros.
#'@param level Level of significance, defaulted to 0.05. It is taken into account only if \code{interval} is set.
#'@param n_perm Number of permutations, defaulted to 1000. It is taken into account only if \code{type} is set to "permutational".
#'@return The output is a well defined [inferenceDataObject], that can be used as parameter in the [smooth.FEM] function.
#'@description A function that build an [inferenceDataObject]. In the process of construction many checks over the input parameters are carried out so that the output is a well defined object,
#'that can be used as parameter in [smooth.FEM] function. Notice that this constructor ensures well-posedness of the object, but a further check on consistency with smooth.FEM parameters will be carried out inside that function.
#'
#'@usage inferenceDataObjectBuilder<-function(test = NULL, 
#'interval = NULL, 
#'type = "wald", 
#'exact = "False", 
#'dim = NULL, 
#'coeff = NULL, 
#'beta0 = NULL, 
#'level = 0.05,
#'n_perm = 1000)
#' @export
#' 
#' 
#' @examples 
#' obj1<-inferenceDataObjectBuilder(test = "simultaneous", interval = NULL, dim = 4);
#' obj2<-inferenceDataObjectBuilder(interval = "one-at-the-time", dim = 5, level = 0.01);


inferenceDataObjectBuilder<-function(test = NULL, 
                                interval = NULL, 
                                type = "wald", 
                                exact = "False", 
                                dim = NULL, 
                                coeff = NULL, 
                                beta0 = NULL, 
                                level = 0.05,
                                n_perm = 1000){
  
  # Preliminary check of parameters input types, translation into numeric representation of default occurrences.
  if(!is.null(test)){
    if(class(test)!="character")
      stop("'test'should be a character: choose one between 'one-at-the-time' or 'simultaneous'")
    if(length(test)==0)
      stop("'test' is zero dimensional, should be one between 'ont-at-the-time' or 'simultaneous'")
  }else{test_numeric=as.integer(0)}
  
  if(!is.null(interval)){
    if(class(interval)!="character")
      stop("'interval'should be a character: choose one between 'one-at-the-time' , 'simultaneous' or 'bonferroni'")
    if(length(interval)==0)
      stop("'interval' is zero dimensional, should be one between 'one-at-the-time', 'simultaneous' or 'bonferroni'")
  }else{interval_numeric=as.integer(0)}
  
  if(type!="wald"){
    if(class(type)!="character")
      stop("'type' should be a character: choose one among 'wald', 'speckman' or 'permutational'" )
    if(length(type)==0)
      stop("'type' is zero dimensional, should be one among 'wald', 'speckman' or 'permutational'")
  }
  
  if(exact!="False"){
    if(class(exact)!="character")
      stop("'exact' should be either 'True' or 'False'")
    if(length(exact)==0)
      stop("'exact' is zero dimensional, should be either 'True' or 'False'")
    if(exact!="True")
      stop("'exact' should be either 'True' or 'False'")
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
    if(class(coeff)!="matrix")
      stop("'coeff' should be of class matrix")
  }
  
  if(!is.null(beta0)){
    if(class(beta0)!="numeric")
      stop("'beta0' should be numeric")
    if(length(beta0)==0)
      stop("'beta0' is zerodimensional")
  }
  

  if(level!=0.05){
    if(class(level)!="numeric")
      stop("'level' should be numeric")
    if(length(level)==0)
      stop("'level' is zerodimensional, should be a positive number between 0 and 1")
  }
  
  if(!is.null(n_perm)){
    if(class(n_perm)!="numeric" && class(n_perm)!="integer")
      stop("'n_perm' should be an integer or convertible to integer type")
    n_perm=as.integer(n_perm)
  }
  
  # Check of consistency of parameters. Translation into numeric representation.
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
  
  if(type!="wald" && type!="speckman" && type!="permutational"){
    stop("type should be choosen between 'wald', 'speckman' and 'permutational'")}else{
      if(type=="wald") type_numeric=as.integer(1)
      if(type=="speckman") type_numeric=as.integer(2)
      if(type=="permutational") type_numeric=as.integer(3)
    }
  
  if(!is.null(test)){
    if(test!="one-at-the-time" && test!="simultaneous"){
      stop("test should be either 'one-at-the-time' or 'simultaneous'")}else{
        if(test=="one-at-the-time") test_numeric=as.integer(1)
        if(test=="simultaneous") test_numeric=as.integer(2)
       }
    if(is.null(beta0))                                                          # If it left to NULL, beta0 is set to a vector of zeros.
      beta0<-rep(0, dim(coeff)[1])
    else{
      if(length(beta0)!=dim(coeff)[1])
        stop("dimension of 'coeff' and 'beta0' are not consistent")
    }
  }
  
  if(!is.null(interval) && (type=="wald" || type=="speckman")){
    if(interval!="one-at-the-time" && interval!="simultaneous" && interval!="bonferroni"){
      stop("interval should be either 'one-at-the-time' 'simultaneous' or 'bonferroni'")}else{
        if(interval=="one-at-the-time") interval_numeric=as.integer(1)
        if(interval=="simultaneous") interval_numeric=as.integer(2)
        if(interval=="bonferroni") interval_numeric=as.integer(3)
      }
    
    if(level <= 0 || level > 1)                                                
      stop("level should be a positive value smaller or equal to 1")
  }
  
  if(!is.null(n_perm)){
  if(n_perm <= 0)                                                
    stop("number of permutations must be a positive value")
  }
  
  else {
    n_perm <- as.integer(1000)
  }
  
  
  if(!is.null(interval) && type=="permutational"){
    stop("confidence intervals are not implemented in the permutational case")
  }
  
  if(type == "permutational"){
    for(i in 1:dim(coeff)[1]){
      count=0
      for(j in 1:dim(coeff)[2]){
        if(coeff[i,j]!= 0 && coeff[i,j]!=1)
          stop("linear combinations are not allowed in the permutational case")
        count = count + coeff[i,j]
      }
      if(count != 1)
        stop("linear combinations are not allowed in the permutational case")
    }
  }
  
  if(is.null(beta0)) beta0=0 #won't be used anyway                              # If beta0 is still NULL here, no test is required, and this parameter is not considered. Set to zero in order to compel with the dataInferenceObject class.
  
  # Well posedeness check for coeff in simultaneous case;
  if(!is.null(test)){
  if(test=="simultaneous"){
    if(det(coeff %*% t(coeff)) < 0.001){
      stop("coeff is not full rank, variance-covariance matrix of the linear combination not invertible")
    }
  }
  }
  
  if(!is.null(interval)){
    if(interval=="simultaneous"){
      if(det(coeff %*% t(coeff)) < 0.001){
        stop("coeff is not full rank, variance-covariance matrix of the linear combination not invertible")
      }
    }
  }
  
  definition=as.integer(1)
  
  # Build the quantile for Confidence intervals if needed
  quantile=0
  if(level > 0){
    if(interval_numeric == 1){ # Simultaneous CI -> Chi-Squared (q) law for statistic
      # C++
      #chi_squared distribution(q);
      #quant =std::sqrt(quantile(complement(distribution,alpha)))
      quantile=qchisq(1-level, dim(coeff)[1])
    }
    else{
      if(interval_numeric ==  2 ){  # One at the time CI (each interval has confidence alpha) -> Gaussian law for statistic
        # C++
        #normal distribution(0,1);
        #quant = quantile(complement(distribution,alpha/2));
        quantile=qnorm((1-level/2),0,1)
      }
      else{ 
        if(interval_numeric == 3){# One at the time CI (overall confidence alpha) -> Gaussian law for statistic
        #C++
        #normal distribution(0,1);
        #quant = quantile(complement(distribution,alpha/(2*q)));
        quantile=qnorm((1-level/(2*dim(coeff)[1])),0,1)
        }
      }
    }
  }
  
  # Building the output object, returning it
  result<-new("inferenceDataObject", test = test_numeric, interval =interval_numeric, type = type_numeric, exact = exact_numeric, dim = dim, 
              coeff = coeff, beta0 = beta0, quantile = quantile, n_perm = n_perm, definition=definition)
  
  return(result)
}

