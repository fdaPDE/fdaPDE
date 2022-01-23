#' Class for inference data
#'
#'@slot test A vector of integers taking value 0, 1 or 2; in the first case no test is performed, in the second one-at-the-time tests are performed,
#'in the third a simultaneous test is performed.
#'@slot interval A vector of integers taking value 0, 1, 2 or 3; In the first case no confidence interval is computed, in the second case one-at-the-time confidence intervals are computed, 
#'in the third case simultaneous confidence intervals are computed, in the fourth case Bonferroni confidence intervals are computed.
#'@slot type A vector of integers taking value 1, 2 or 3, corresponding to Wald, Speckman or Eigen-Sign-Flip implementation 
#'of the inference analysis.
#'@slot exact An integer taking value 1 or 2. If 1 an exact computation of the test statistics will be performed,
#'whereas if 2 an approximated computation will be carried out.
#'@slot enhanced An integer taking value 1 or 2. If 1 an enhanced computation of the beta ESF test statistics will be performed,
#'whereas if 2 the classical beta ESF test will be perfoermed.
#'@slot component A vector of integers taking value 1,2 or 3, indicating whether the inferential analysis should be carried out respectively for the parametric, nonparametric or both the components.  
#'@slot dim Dimension of the problem, it is equal to 2 in the 2D case and equal to 3 in the 2.5D and 3D cases. 
#'@slot n_cov Number of covariates taken into account in the linear part of the regression problem.
#'@slot locations A matrix of numeric coefficients. When nonparametric inference is requested it represents the set of spatial locations for which the inferential analysis should be performed. 
#' The default values is a one-dimensional matrix of value 1 indicating that all the observed location points should be considered. 
#' In the sign-flip and eigen-sign-flip implementations only observed points are allowed. 
#'@slot locations_indices A vector of indices indicating the which spatial points have to be considered among the observed ones for nonparametric inference. If a vector of location indices is provided
#'then the slot 'location' is discarded. 
#'@slot locations_are_nodes An integer taking value 1 or 2; in the first case it indicates that the selected locations to perform inference on f are all coinciding with the nodes; otherwise it takes value 2;  
#'@slot coeff A matrix of numeric coefficients with columns of dimension \code{dim} and each row represents a linear combination of the linear parameters to be tested and/or to be
#' estimated via confidence interval. 
#'@slot beta0 Vector of null hypothesis values for the linear parameters of the model. Used only if \code{test} is not 0.
#'@slot f0 Function representing the expression of the nonparametric component f under the null hypothesis.
#'@slot f0_eval Vector of f0 evaluations at the choosen test locations. It will be eventually set later in checkInferenceParameters, if nonparametric inference is required. 
#'@slot f_var An integer taking value 1 or 2. If 1 local variance estimates for the nonlinear part of the model will be computed,
#'whereas if 2 they will not.
#'@slot quantile Quantile needed for confidence intervals. Used only if interval is not 0.
#'@slot alpha Level vector of Eigen-Sign-Flip confidence intervals. Used only if interval is not 0.
#'@slot n_flip An integer representing the number of permutations in the case of eigen-sign-flip test.
#'@slot tol_fspai A real number greater than 0 specifying the tolerance for FSPAI algorithm, in case of non-exact inference.
#'@slot definition An integer taking value 0 or 1. If set to 1, the class will be considered as created by the function [inferenceDataObjectBuilder],
#'leading to avoid some of the checks that are performed on inference data within smoothing functions.
#'
#'@description
#' A class that contains all possible information for inference over linear parameters and/or nonparametric field in spatial regression with
#' differential regularization problem. This object can be used as parameter in smoothing function of the fdaPDE library [smooth.FEM].
#' 
#'@details #Warning
#'At least one between test and interval must be nonzero. \code{n_cov}, \code{coeff} and \code{beta0} need to be coherent.
#'The usage of [inferenceDataObjectBuilder] is recommended for the construction of an object of this class.

inferenceDataObject<-setClass("inferenceDataObject", slots = list(test = "integer", 
                                                                  interval = "integer", 
                                                                  type = "integer", 
                                                                  component = "integer",
                                                                  exact = "integer",
                                                                  enhanced = "integer",
                                                                  dim = "integer",
                                                                  n_cov = "integer",
                                                                  locations = "matrix",
                                                                  locations_indices = "integer",
                                                                  locations_are_nodes = "integer",
                                                                  coeff = "matrix",
                                                                  beta0 = "numeric",
                                                                  f0 = "function",
                                                                  f0_eval = "numeric",
                                                                  f_var = "integer",
                                                                  quantile = "numeric",
                                                                  alpha = "numeric",
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
#' and can take values 'one-at-the-time', 'simultaneous' (available only for Wald and Speckman cases) and 'bonferroni'. If the value is NULL, no interval will be computed, and the \code{test} parameter in the corresponding position needs to be set.
#'  Otherwise one at the time, simultaneous or Bonferroni correction intervals will be computed. If it is not NULL, the parameter \code{level} will be taken into account.
#'@param type A list of strings defining the type of implementation for the inferential analysis, with the same length of \code{test} list and \code{interval} list . The possible values are three:
#''wald'(default), 'speckman' 'sign-flip' or 'eigen-sign-flip'. If 'eigen-sign-flip' is set, the corresponding \code{interval} position needs to be NULL.
#'@param component A list of strings defining on which model component inference has to be performed. It can take values 'parametric', 'nonparametric' or 'both'. The default is 'parametric'.
#'@param exact A logical used to decide the method used to estimate the statistics variance.
#'The possible values are: FALSE (default) and TRUE. In the first case an approximate method is used, leading to a lower accuracy, but faster computation.
#'In the second case the evaluation is exact but computationally expensive.
#'@param enhanced A logical used to decide wether to perform enhaned on classical ESF tests.
#'The possible values are: FALSE (default) and TRUE. In the first case classical ESF test will be performed.
#'In the second case enhanced ESF procedures will be employed leading to a better Type-I error but slight lower power.
#'@param dim Dimension of the problem, defaulted to NULL. (Must be set by the user)
#'@param n_cov Number of the covariates, defaulted to NULL. (Must be set by the user)
#'@param locations A matrix of the locations of interest when testing the nonparametric component f
#'@param locations_indices A vector of indices indicating the locations to be considered among the observed ones for nonparametric inference. If a vector of indices is provided, then the slot 'locations' is discarded.
#'@param coeff A matrix, with \code{n_cov} number of columns, of numeric coefficients, defaulted to NULL. If this parameter is NULL,
#'in the corresponding inferenceDataObject it will be defaulted to an identity matrix. If at least one 'eigen-sign-flip' value is present in \code{type}, needs to be an identity matrix.
#'@param beta0 Vector of real numbers (default NULL). It is used only if the \code{test} parameter is set, and has length the number of rows of matrix \code{coeff}. 
#'If \code{test} is set and \code{beta0} is NULL, will be set to a vector of zeros.
#'@param f0 A function object representing the expression of the nonparametric component f under the null hypothesis. If NULL, the default is the null function, hence a test on the significance of the nonparametric component is carried out. 
#'@param f_var A logical used to decide whether to estimate the local variance of the nonlinear part of the model.
#'The possible values are: FALSE (default) and TRUE. 
#'@param level A vector containing the level of significance used to compute quantiles for confidence intervals, defaulted to 0.95. It is taken into account only if \code{interval} is set.
#'@param n_flip Number of flips performed in Eigen-Sign-Flip test, defaulted to 1000. It is taken into account only if at least one position of \code{type} is set to 'eigen-sign-flip'.
#'@param tol_fspai Tolerance for FSPAI algorithm taking value greater than 0, defaulted to 0.05. It is taken into account only if \code{exact} is set to 'False'. The lower is the tolerance, the heavier is the computation.
#'@return The output is a well defined [inferenceDataObject], that can be used as parameter in the [smooth.FEM] function.
#'@description A function that build an [inferenceDataObject]. In the process of construction many checks over the input parameters are carried out so that the output is a well defined object,
#'that can be used as parameter in [smooth.FEM] function. Notice that this constructor ensures well-posedness of the object, but a further check on consistency with smooth.FEM parameters will be carried out inside that function.
#'
#'@usage inferenceDataObjectBuilder<-function(test = NULL, 
#'interval = NULL, 
#'type = 'wald', 
#'component = 'parametric',
#'exact = FALSE, 
#'enhanced = FALSE,
#'dim = NULL, 
#'n_cov = NULL,
#'locations = NULL,
#'locations_indices = NULL,
#'coeff = NULL, 
#'beta0 = NULL, 
#'f0 = NULL,
#'f_var = FALSE,
#'level = 0.95,
#'n_flip = 1000,
#'tol_fspai = 0.05)
#' @export
#' 
#' 
#' @examples 
#' obj1<-inferenceDataObjectBuilder(test = "simultaneous", interval = NULL, exact = T, dim = 2, n_cov = 4);
#' obj2<-inferenceDataObjectBuilder(interval = "one-at-the-time", dim = 3, n_cov = 5, level = 0.99);
#' obj3<-inferenceDataObjectBuilder(test=c('one-at-the-time', 'simultaneous', 'one-at-the-time','one-at-the-time'), interval=c('bonferroni','one-at-the-time','none','simultaneous'),
#'  type=c('wald','speckman','eigen-sign-flip','speckman'),exact=TRUE, dim=2, n_cov = 2, level=0.99)
#' obj4<-inferenceDataObjectBuilder(test = "simultaneous", interval = c("one-at-the-time", "simultaneous"), type = "wald", component = c("parametric", "nonparametric"), dim = 2, n_cov = 2, exact = T, locations_indices = 1:10, level = c(0.95, 0.9))


inferenceDataObjectBuilder<-function(test = NULL, 
                                interval = NULL, 
                                type = "wald", 
                                component = "parametric",
                                exact = F,
                                enhanced = F,
                                dim = NULL,
                                n_cov = NULL,
                                locations = NULL,
                                locations_indices = NULL,
                                coeff = NULL, 
                                beta0 = NULL,
                                f0 = NULL,
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
  }
  
  if(!is.null(interval)){
    if(class(interval)!="character")
      stop("'interval'should be a character: choose one between 'one-at-the-time' , 'simultaneous', 'bonferroni', 'none'")
    if(length(interval)==0)
      stop("'interval' is zero dimensional, should be a vector taking values among 'one-at-the-time', 'simultaneous', 'bonferroni', 'none'")
  }
  
  if(length(type)==0)
    stop("'type' is zero dimensional, should be a vector taking values among 'wald', 'speckman', 'sign-flip' or 'eigen-sign-flip'")
 
   if(length(type) > 1 || type!="wald"){
    if(class(type)!="character")
      stop("'type' should be a vector of characters: choose among 'wald', 'speckman', 'sign-flip' or 'eigen-sign-flip'" )
  }
  
  
  if(length(component)==0)
    stop("'component' is zero dimensional, should be a vector of characters taking values among 'parametric', 'nonparametric' or 'both'")
  
  if(length(component) > 1 || component!="parametric"){
    if(class(component)!="character")
      stop("'component' should be a vector of characters with values among 'parametric', 'nonparametric' or 'both'" )
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
  
  if(enhanced != F){
    if(class(enhanced)!="logical")
      stop("'enhanced' should be either TRUE or FALSE ")
    if(length(enhanced)==0)
      stop("'enhanced' is zero dimensional, should be either TRUE or FALSE")
    if(enhanced!=T)
      stop("'enhanced' should be either TRUE or FALSE")
    enhanced_numeric=as.integer(1)
  }else{
    enhanced_numeric=as.integer(2)
  }
  
  if(!is.null(dim)){
    if(class(dim)!="numeric" && class(dim)!="integer")
      stop("'dim' should be an integer or convertible to integer type")
    else if(dim!=2 && dim!=3)
      stop("'dim' should be either 2 or 3")
    
    dim=as.integer(dim)
  }
  
  if(!is.null(n_cov)){
    if(class(n_cov)!="numeric" && class(n_cov)!="integer")
      stop("'dim' should be an integer or convertible to integer type")
    n_cov=as.integer(n_cov)
  }
  
  if(!is.null(locations)){
    if(class(locations)[1]!="matrix")
      stop("'locations' should be of class matrix")
  }
  
  if(!is.null(locations_indices)){
    if(class(locations_indices)!="numeric" && class(locations_indices)!="integer")
      stop("'locations_indices' should be of class numeric or convertible to integer")
    else{
      for(i in 1:length(locations_indices)){
        j <- as.integer(locations_indices[i]) 
        if(locations_indices[i]!=j || locations_indices[i] <= 0)
          stop("'locations_indices' should contain positive integers")
      }
      rm(list = c("i","j"))
    }
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
  
  if(!is.null(f0)){
    if(class(f0)!="function")
      stop("'f0' should be a function")
    if(length(f0)==0)
      stop("'f0' is zerodimensional")
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
  
  
  if(length(level)==0)
    stop("'level' is zerodimensional, should be a vector of positive numbers between 0 and 1")
  if(length(level) > 1 || level!=0.95){
    if(class(level)!="numeric")
      stop("'level' should be numeric")
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
  
  if(sum(component == "parametric")!=length(component)){
    if(is.null(dim))
      stop("dimension of the problem is needed")                                # Otherwise the function won't be able to default f0
  }
  else{
    if(is.null(dim))
      dim = as.integer(0)                                                       # Inference on f is not required, we set dim to zero by default since it won't be required
  }
  
  if(sum(component == "nonparametric")!=length(component)){
    if((is.null(n_cov) || n_cov <= 0))                                          # Otherwise the function won't be able to default the coeff neither object nor beta0
      stop("number of covariates is needed")
    else{
      if(is.null(coeff)){                                                       # If it is left as NULL, all the parameters are individually taken into account without any linear combination.
        coeff = diag(1, nrow=n_cov, ncol=n_cov)
        }
      else{
        if(dim(coeff)[2]!=n_cov)
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
  }else{
    # inference on beta is not required, set n_cov to 0 and coeff to a 1x1 matrix
    n_cov = as.integer(0)
    coeff = matrix(data=0, nrow=1, ncol=1)
  }
  
  if(sum(component == "parametric")==length(component)){
    if(!is.null(locations) || !is.null(locations_indices))
      warning("locations are provided but not used, since inference is requested only on the paramteric component")
    locations_indices = NULL
    locations = matrix(data = 1, nrow = 1, ncol = 1)
  }
  else{
  if(is.null(locations) && is.null(locations_indices)){
    # default case: all observed locations are considered 
    locations <- matrix(data=1, nrow=1, ncol=1)                                 # Just for convenience, the actual default will be set inside checkInferenceParameters
  }
  else if(is.null(locations_indices)){
    if(sum(type == "eigen-sign-flip")==1 || sum(type == "sign-flip")==1)
      stop("For the eigen-sign-flip implementation a vector of location indices is required")
    if(ncol(locations)!=dim)
      stop("number of columns of 'locations' and 'dim' do not match")
  }
  else{
    if(!is.null(locations))
      warning("'locations' are provided but not used, since 'locations_indices' is not NULL")
    # other dimensional checks on 'locations_indices' will be performed later, inside checkInferenceParameters
  }
  }  
  
  # Consistency check for number of implementations required (default values assigned in case of NULL)
  n_of_implementations=max(length(type), length(test), length(interval), length(component));
  
  # Pre-allocation of space
  test_numeric=rep(0, n_of_implementations)
  interval_numeric=rep(0, n_of_implementations)
  component_numeric=rep(0, n_of_implementations)
  type_numeric=rep(0, n_of_implementations)
  quantile=rep(0, n_of_implementations)
  
  if(!is.null(test)){
    if(length(test) < n_of_implementations){
      if(length(test) > 1)
        stop("Test vector length is not consistent with the number of implementations provided")
      else
        test = rep(test, n_of_implementations)
    }
  }else{test = rep("none", n_of_implementations)}
  
  if(!is.null(interval)){
    if(length(interval) < n_of_implementations){
      if(length(interval) > 1)
        stop("Intervals vector length is not consistent with the number of implementations provided")
      else
        interval = rep(interval, n_of_implementations)
    }
  }else{interval = rep("none", n_of_implementations)}
  
  if(!is.null(component)){
    if(length(component) < n_of_implementations){
      if(length(component) > 1)
        stop("Component vector length is not consistent with the number of implementations provided")
      else
        component = rep(component, n_of_implementations)
    }
  }
  
  if(!is.null(type)){
    if(length(type) < n_of_implementations){
      if(length(type) > 1)
        stop("Type vector length is not consistent with the number of implementations provided")
      else
        type = rep(type, n_of_implementations)
    }
  }
  
  # check level consistency
  if(!is.null(level)){
    
    n_of_intervals <- 0
    for(i in 1:n_of_implementations){
      if(interval[i]!="none")
        n_of_intervals <- n_of_intervals + 1
    }
    
    if(length(level)==1){
      level_vec <- rep(level, n_of_implementations)
      alpha <- rep(1-level, n_of_implementations)
    }
    
    else if(length(level)==n_of_intervals){
      level_vec <- numeric(n_of_implementations)
      alpha <- numeric(n_of_implementations)
      for(i in 1:n_of_implementations){
        if(interval[i]!="none")
          level_vec[i] <- level[i]
        else
          level_vec[i] <- 0
        alpha[i] <- 1 - level_vec[i] 
      }
    }
    
    else
      stop("Level vector length is not consistent with the number of confidence interval implementations requested")
    
    level = level_vec
    
    rm(list = ("level_vec"))
    rm(list = ("i"))
    
  }
  
  for (index in 1:n_of_implementations){
    
      if(type[index]!="wald" && type[index]!="speckman" && type[index]!="eigen-sign-flip" && type[index]!="sign-flip"){
        stop("type should be chosen between 'wald', 'speckman' 'sign-flip' and 'eigen-sign-flip'")}else{
          if(type[index]=="wald") type_numeric[index]=as.integer(1)
          if(type[index]=="speckman") type_numeric[index]=as.integer(2)
          if(type[index]=="eigen-sign-flip") type_numeric[index]=as.integer(3)
          if(type[index]=="sign-flip") type_numeric[index]=as.integer(4)
        }
    
    if(component[index]!="parametric" && component[index]!="nonparametric" && component[index]!="both"){
      stop("component should be chosen between 'parametric', 'nonparametric' and 'both'")}else{
        if(component[index]=="parametric") component_numeric[index]=as.integer(1)
        if(component[index]=="nonparametric") component_numeric[index]=as.integer(2)
        if(component[index]=="both") component_numeric[index]=as.integer(3)
      }
      
    if(type[index]=="sign-flip" && component[index]!="nonparametric")
       stop("sign-flip test is implemented only for the nonparametric component")
    
    if(type[index]=="speckman" && component[index]!="parametric")
      stop("speckman test and confidence intervals are implemented only for the parametric component")
    
    if(test[index]=='none' && interval[index]=='none'){
        stop("At least one between test and interval should be required for each implementation")
      }
    
      if(test[index]!="one-at-the-time" && test[index]!="simultaneous" &&  test[index]!="none"){
        stop("test should be either 'one-at-the-time', 'simultaneous' or 'none'")}else{
          if(test[index]=="none") test_numeric[index]=as.integer(0)
          if(test[index]=="one-at-the-time") test_numeric[index]=as.integer(1)
          if(test[index]=="simultaneous") test_numeric[index]=as.integer(2)
       }
      
      if(interval[index]!="none"){
        if(interval[index]!="one-at-the-time" && interval[index]!="simultaneous" && interval[index]!="bonferroni" && interval[index]!="none"){
          stop("interval should be either 'one-at-the-time' 'simultaneous', 'bonferroni' or 'none'")}else{
            if(interval[index]=="none") interval_numeric[index]=as.integer(0)
            if(interval[index]=="one-at-the-time") interval_numeric[index]=as.integer(1)
            if(interval[index]=="simultaneous") interval_numeric[index]=as.integer(2)
            if(interval[index]=="bonferroni") interval_numeric[index]=as.integer(3)
          }
        
        if(level[index] <= 0 || level[index] > 1)                                                
          stop("level should be a positive value smaller or equal to 1")
      }
    
      if(interval[index]!="none" && type[index]!="wald" && component[index]!="parametric")
        stop("confidence intervals are not available with sign-flipping implementation")
      
      if(interval[index]=="simultaneous" && type[index]=="eigen-sign-flip"){
        stop("simultaneous confidence intervals are not implemented in the eigen-sign-flip case")
      }
      
      if(type[index] == "eigen-sign-flip" && component[index]!="nonparametric"){
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
      
      # Well posedness check for coeff in simultaneous case;
      if(test[index]=="simultaneous" && component[index]!="nonparametric"){
        if(det(coeff %*% t(coeff)) < 0.001){
          stop("coeff is not full rank, variance-covariance matrix of the linear combination not invertible")
        }
      }
        
      if(interval[index]=="simultaneous" && component[index]!="nonparametric"){
        if(det(coeff %*% t(coeff)) < 0.001){
          stop("coeff is not full rank, variance-covariance matrix of the linear combination not invertible")
        }
      }
    
      if((test[index]=="one-at-the-time") && component[index] != "parametric"){
        warning("Only simultaneous tests are available for the nonparametric component, proceeding with simultaneous inference")
        test_numeric[index] = as.integer(2)
      }
    
      
      # Build the quantile for confidence intervals if needed 
      if(interval[index]=="none" || component[index]=="nonparametric"){
         quantile[index]=0
        }else{
          if(level[index] > 0){
              if(interval[index] == "simultaneous"){ # Simultaneous CI -> Chi-Squared (q) law for statistic
                quantile[index]=sqrt(qchisq(level[index], dim(coeff)[1]))  
              }else{
                if(interval[index]=="one-at-the-time"){  # One at the time CI (each interval has confidence alpha) -> Gaussian law for statistic
                  quantile[index]=qnorm((1-(1-level[index])/2),0,1)
              }else{ 
                  if(interval[index]== "bonferroni"){# One at the time CI (overall confidence alpha) -> Gaussian law for statistic
                  quantile[index]=qnorm((1-(1-level[index])/(2*dim(coeff)[1])),0,1)
                  }
              }
            }
          }
      }
      
  } # End of for cycle over the different implementation
if(sum(component == "nonparametric")!=length(component)){  
  if(is.null(beta0))                 # If left to NULL, beta0 is set to a vector of zeros.
    beta0<-rep(0, dim(coeff)[1])
  else{
    if(length(beta0)!=dim(coeff)[1])
      stop("dimension of 'coeff' and 'beta0' are not consistent")
  }
}else{
  # inference on beta is not required, just setting it to zero
  beta0<-0
}

  if(sum(component == "parametric")!=length(component)){
  if(is.null(f0)){                  # If left to NULL, f0 is set to a null function.
    if(dim==2){
      f0 <- function(x,y){
        return(0)
      }
    }
    else if(dim==3){
      f0 <- function(x,y,z){
        return(0)
      }
    }
  }
  else{
    args<-formals(f0)
    non_defaulted_args = 0
    for(i in 1:length(args)){
      if(is.name(args[[i]]))
        non_defaulted_args <- non_defaulted_args + 1
    }
    if(non_defaulted_args != dim)
      stop("Number of f0 coordinate arguments is not consistent with the problem dimension 'dim'")
    rm(list = c("non_defaulted_args", "i"))
  }
  }
  else{
    # inference on the nonparametric component is not required, just setting f0 to a zero function
    f0 <- function(){return(0)}
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
  if(!is.null(locations_indices))
    result<-new("inferenceDataObject", test = as.integer(test_numeric), interval = as.integer(interval_numeric), type = as.integer(type_numeric), component = as.integer(component_numeric), exact = exact_numeric, enhanced = enhanced_numeric, dim = dim, n_cov = n_cov,
              locations_indices = as.integer(locations_indices), locations_are_nodes = as.integer(1), coeff = coeff, beta0 = beta0, f0 = f0, f_var = f_var_numeric, quantile = quantile, alpha = alpha, n_flip = n_flip, tol_fspai = tol_fspai, definition=definition)
  else
    result<-new("inferenceDataObject", test = as.integer(test_numeric), interval = as.integer(interval_numeric), type = as.integer(type_numeric), component = as.integer(component_numeric), exact = exact_numeric, enhanced = enhanced_numeric, dim = dim, n_cov = n_cov,
                locations = locations, locations_are_nodes = as.integer(1), coeff = coeff, beta0 = beta0, f0 = f0, f_var = f_var_numeric, quantile = quantile, alpha = alpha, n_flip = n_flip, tol_fspai = tol_fspai, definition=definition)
    
  
  return(result)
}

