#' Class for inference data
#'
#'@slot test An integer taking value 0, 1 or 2; in the first case no test is performed, in the second a null hypotesis test,
#'in the third a power analysis.
#'@slot interval An integer taking value 0, 1 or 2; In the first case no interval is computed, in the second a one-at-the-time interval, 
#'in the third a boferroni interval.
#'@slot type An integer taking value 1, 2 or 3, corresponding to wald, sandwich or permutational implementation 
#'of the inference analysis.
#'@slot exact An integer taking value 1 or 2. If 1 an exact computation of the test statistic will be performed,
#'whereas if 2 an approximated computation will be carried out.
#'@slot dim Number of covariates taken into account in the linear part of the regression problem-
#'@slot coeff A vector of integer numbers of lengths \code{dim} that can take value 1 or 0. When taking value 1,
#'the corresponding covariate will be teken into account in inference analysis.
#'@slot beta0 Vector of null hypotesis values for the linear parameters of the modeled. Used only if \code{test} is not 0.
#'@slot beta1 Vector of alternative hypotesis values for the linear parameters of the modeled. Used only if \code{test} is not 0.
#'@slot level Level of significance for the confidence intervals. Needs to be between 0 and 1. Is used only if interval is not 0.
#'@slot definition An integer taking value 0 or 1. If set to 1, the class will be considered as created by the function [inferenceDataObjectBuilder],
#'leading to avoid some of the checks that are performed on inference data within smoothing functions.
#'
#'@description
#' A class that contains all possible information for linear parameters in spatial regression with
#' differential regularization problem. This object can be used as parameter in smoothing function of the fdaPDE libray [smooth.FEM].
#' 
#'@details #Warning
#'At least one between test and interval must be nonzero. \code{dim}, \code{coeff} and \code{beta0} need to be coherent.
#'The usage of [inferenceDataObjectBuilder] is reccomended for the construction of an oobject o this class.

inferenceDataObject<-setClass("inferenceDataObject", slots = list(test = "integer", #"character",
                                                                  interval = "integer", #"character",
                                                                  type = "integer", #"character",
                                                                  exact = "integer", #"character",
                                                                  dim = "integer",
                                                                  coeff = "integer",
                                                                  beta0 = "numeric",
                                                                  beta1 = "numeric",
                                                                 level = "numeric",
                                                                 definition="integer")
                              )

#'Constructor for inferenceDataObject class
#'
#'@param test A string defining the type of test to be performed. The default is NULL, and can take value 'p-value' or 'power'.
#'If the value is NULL, no test is performed, and the \code{interval} parameter need to be not NULL. If it takes value
#''p-value', the parameters \code{beta0}, \code{dim}, \code{coeff} are taken into account, and a null hypotesis test will be performed. If it takes value 
#''power', an alternative hypotesis test will be performed, and even the \code{beta1} will be taken  into account.
#'@param interval A string defining the type of confidence interval to be computed. The default is NULL, and can take value 'one-at-the-time' and 'bonferroni'.
#'If the value is NULL, no interval will be computed, and the \code{test} parameter needs to be set. Otherwise a one at the time or bonferroni correction interval will be computed.
#'If it is not NULL, the parameter \code{level} will be taken into account. Up to now, confidence interval can be computed only in the wald implementation.
#'@param type A string defining the type of implementation for the inferential analysis. The possible values are three:
#''wald'(default), 'sandwich' or 'permutational', corresponding to the three possibe methods developed in Ferraccioli....
#'@param exact A string used to decide the method users to estimate the statistics variance.
#'The possible values are: 'True' and 'False'(default). In the first case the evaluation is exact but computationally very expensive.
#'In the second case an apporximate method is used, leading to a lower accuracy, but faster computation.
#'@param dim Number of the covariates, defaulted to NULL.
#'@param coeff A vector, of lengths \code{dim} of integer numbers that can take value 1 or 0, defaulted to NULL. When takes value 1, the corresponding
#'covariate is taken into account for the inferential problem (test or interval). If this parameter is NULL, all the covariates will be taken into account.
#'@param beta0 Vector of real numbers (default NULL). It is used only if the \code{test} parameter is set, and has length the sum of the values in coeff. If \code{test} is set and \code{beta0} is NULL,
#'will be set to a vector of zeros.
#'@param beta1 Vector of real numbers (default NULL). It is used only if the \code{test} parameter is set to 'power'. In thatr case needs to be expressed, and cannot be left NULL.
#'@param level Level of significance, defaulted to 0.05. It is taken into account only if \code{intrval} is set.
#'@param definition A 0-1 integer. If it is set, and the procedure outcomes a valid inferenceDataObjec, some checks will be avoided inside [smooth.FEM].
#'@return The output is a well defined [inferenceDataObject], that can be used as parameter in the [smooth.FEM] function.
#'@description A function that build an [inferenceDataObject]. In the process of construction many checks over the input parameters are carried out so that the output is a well defined object,
#'that can be used as parameter in [smooth.FEM] function. Notice that this constructor ensures well-posedeness of the object, but a further check on consistency with smooth.FEM parameters will be carried out inside that function.
#'
#'@usage inferenceDataObjectBuilder<-function(test = NULL, 
#'interval = NULL, 
#'type = "wald", 
#'exact = "False", 
#'dim = NULL, 
#'coeff = NULL, 
#'beta0 = NULL, 
#'beta1 = NULL, 
#'level = 0.05,
#'definition=1)
#' @export
#' 
#' 
#' @examples 
#' obj1<-inferenceDataObjectBuilder(test = "pvalue", interval = NULL, dim = 4);
#' obj2<-inferenceDataObjectBuilder(interval = "one-at-the-time", dim = 5, level = 0.01);


inferenceDataObjectBuilder<-function(test = NULL, 
                                interval = NULL, 
                                type = "wald", 
                                exact = "False", 
                                dim = NULL, 
                                coeff = NULL, 
                                beta0 = NULL, 
                                beta1 = NULL, 
                                level = 0.05,
                                definition=1){
  
  # Preliminary check of parameters input types, translation into numeric representation of default occurrencies.
  if(!is.null(test)){
    if(class(test)!="character")
      stop("'test'should be a character: choose one between 'pvalue' or 'power'")
    if(length(test)==0)
      stop("'test' is zero dimensional, should be one between 'pvalue' or 'power'")
  }else{test_numeric=as.integer(0)}
  
  if(!is.null(interval)){
    if(class(interval)!="character")
      stop("'interval'should be a character: choose one between 'one-at-the-time' or 'bonferroni'")
    if(length(interval)==0)
      stop("'interval' is zero dimensional, should be one between 'one-at-the-time' or 'bonferroni'")
  }else{interval_numeric=as.integer(0)}
  
  if(type!="wald"){
    if(class(type)!="character")
      stop("'type' should be a character: choose one among 'wald', 'sandwich' or 'permutational'" )
    if(length(type)==0)
      stop("'type' is zero dimensional, should be one among 'wald', 'sandwich' or 'permutational'")
  }
  
  if(exact!="False"){
    if(class(exact)!="character")
      stop("'exact' should be either 'True' or 'False'")
    if(length(exact)==0)
      stop("'exact' is zero dimensional, should be either 'True' or 'False'")
    if(exact!="True")
      stop("'exact' should be either 'True' or 'False'")
    exact_numeric=as.integer(2)
  }else{
    exact_numeric=as.integer(1)
  }
  
  
  if(!is.null(dim)){
    if(length(as.integer(dim))==0)
      stop("'dim' should be an integer or converitble to integer type")
    dim=as.integer(dim)
  }
  
  if(!is.null(coeff)){
    if(length(as.integer(coeff))==0)
      stop("'coeff' should be an integer or converitble to integer type")
  }
  
  if(!is.null(beta0)){
    if(class(beta0)!="numeric")
      stop("'beta0' should be numeric")
    if(length(beta0)==0)
      stop("'beta0' is zerodimensional")
  }
  
  if(!is.null(beta1)){
    if(class(beta1)!="numeric")
      stop("'beta1' should be numeric")
    if(length(beta1)==0)
      stop("'beta1' is zerodimensional")
  }
  
  if(level!=0.05){
    if(class(level)!="numeric")
      stop("'level' should be numeric")
    if(length(level)==0)
      stop("'level' is zerodimensional, should be a positive number between 0 and 1")
  }
  
  # Check of consistency of parameters. Translation into numeric representation.
  if(is.null(test) && is.null(interval))                                        # However, they can be both set
    stop("at least one between 'test' and 'interval' should be not NULL")
  
  if(is.null(dim) || dim <= 0)                                                  # Otherwise the function won't be able to default the coeff neither object nor beta0
    stop("number of covariates is needed")
  else{
    if(is.null(coeff)){                                                         # If it is left as NULL, all the coefficients will be taken into acccount
      coeff = as.integer(rep(1, dim)) 
      count = dim
    }
    else{
      if(length(coeff)!=dim)
        stop("number of covariates and coefficients do not match")
      count = 0
      for (i in 1:dim){
        if(coeff[i]!=0 && coeff[i]!=1)
          stop("'coeff' should be composed by 0-1 values")
        if(coeff[i]==1)
          count = count + 1
      }
      if(count == 0)                                                            # Otherwise object is not needed
        stop("at least one coefficient must be indicated")
      rm(list = ("i"))
    }
  }
  
  if(type!="wald" && type!="sandwich" && type!="permutational"){
    stop("type should be choosen between 'wald', 'sandwich' and 'permutational'")}else{
      if(type=="wald") type_numeric=as.integer(1)
      if(type=="sandwich") type_numeric=as.integer(2)
      if(type=="permutational") type_numeric=as.integer(3)
    }
  
  if(!is.null(test)){
    if(test!="pvalue" && test!="power"){
      stop("test should be either 'pvalue' or 'power'")}else{
        if(test=="pvalue") test_numeric=as.integer(1)
        if(test=="power") test_numeric=as.integer(2)
      }
    if(is.null(beta0))                                                          # If it left to NULL, all the beta0 that are needed w.r.t. coeff will be se to 0 (H0 beta0=0 vs H1 beta0!=0)
      beta0<-rep(0, count)
    else{
      if(length(beta0)!=count)
        stop("dimension of 'coeff' and 'beta0' are not consistent")
    }
    if(test=="power"){
      if(length(beta1)!=count)
        stop("dimension of 'beta1' and 'beta0' do not match")
    }
    rm(list=("count"))
  }
  
  if(!is.null(interval)){
    if(interval!="one-at-the-time" && interval!="bonferroni"){
      stop("interval should be either 'one-at-the-time' or 'bonferroni'")}else{
        if(interval=="one-at-the-time") interval_numeric=as.integer(1)
        if(interval=="bonferroni") interval_numeric=as.integer(2)
      }
    if(level <= 0 || level > 1)                                                
      stop("level should be a positive value smaller or equal to 1")
  }
  if(is.null(beta0)) beta0=0 #won't be used anyway                              # If beta0 is still NULL here, no test is required, and this parameter is not considered. Set to zero in order to compell with the dataInferencebject class.
  if(is.null(beta1)) beta1=0 #won't be used anyway                              # If beta0 is still NULL here, no power is required, and this parameter is not considered. Set to zero in order to compell with the dataInferencebject class.
  definition=as.integer(definition)
  
  # Building the output object, returning it
  result<-new("inferenceDataObject", test = test_numeric, interval =interval_numeric, type = type_numeric, exact = exact_numeric, dim = dim, 
              coeff = coeff, beta0 = beta0, beta1 = beta1, level = level,definition=definition)
  
  return(result)
}

