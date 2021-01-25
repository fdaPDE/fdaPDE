
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
  
  if(is.null(test) && is.null(interval))
    stop("at least one between 'test' and 'interval' should be not NULL")
  
  if(is.null(dim) || dim <= 0)
    stop("number of covariates is needed")
  else{
    if(is.null(coeff)){
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
      if(count == 0)
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
    if(is.null(beta0))
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
  if(is.null(beta0)) beta0=0 #won't be used anyway
  if(is.null(beta1)) beta1=0 #won't be used anyway
  definition=as.integer(definition)
  
  # attenzione a conversioni 
  result<-new("inferenceDataObject", test = test_numeric, interval =interval_numeric, type = type_numeric, exact = exact_numeric, dim = dim, 
              coeff = coeff, beta0 = beta0, beta1 = beta1, level = level,definition=definition)
  
  return(result)
}

