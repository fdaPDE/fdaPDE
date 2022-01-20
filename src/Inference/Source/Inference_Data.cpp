#include "../Include/Inference_Data.h"
//! Main constructor of the class
/*!
  \param test_Type_ vector parameter used to define the type of tests that are required (if any)
  \param interval_Type_ vector parameter used to define which type of confidence intervals are required (if any)
  \param implementation_Type_ vector parameter used to define which type of implementations (Wald, Speckman, ESF) are used for the tests and intervals computation
  \param component_Type_ vector parameter used to on which component of the model (parametric, nonparametric, both) inference has to be performed
  \param exact_Inference_ parameter for the method used to invert E matrix in Woodbury decomposition for inference
  \param enhanced_Inference_ parameter for the type of ESF test performed
  \param locs_Inference_ matrix that specifies the spatial locations to be considered when performing inference on the nonparametric component
  \param locs_index_Inference_ vector containing the spatial locations indices to be considered among the observed ones for inference on the nonparametric component 
  \param coeff_Inference_ matrix that specifies the linear combinations of the linear parameters to be tested and/or estimated via confidence intervals 
  \param beta_0_ vector for the null hypotesis (if a test on the parametric component is required)
  \param f_0_ vector for the null hypotesis (if a test on the nonparametric component is required)
  \param f_Var_ parameter used to decide whether to compute local f variance or not
  \param inference_Quantile_ vector parameter containing the quantiles to be used for the computation of the confidence intervals (if interval_type is defined)
  \param inference_Alpha_ significance used to compute ESF confidence intervals
  \param n_Flip_ parameter that provides the number of sign-flips to be used for the eigen-sign-flip tests (if they are required)
  \param tol_Fspai_ parameter that provides the tolerance used in the FSPAI algorithm
  \param definition_ parameter used to set definition of the InferenceData object
*/
InferenceData::InferenceData(SEXP test_Type_, SEXP interval_Type_, SEXP implementation_Type_, SEXP exact_Inference_,  SEXP enhanced_Inference_, SEXP coeff_Inference_, SEXP beta_0_, SEXP f_Var_,
			     SEXP inference_Quantile_, SEXP inference_Alpha_, SEXP n_Flip_, SEXP tol_Fspai_, SEXP definition_){
  //test_Type
  UInt size_test_Type=Rf_length(test_Type_);
  test_Type.resize(size_test_Type);
  for(UInt i=0; i<size_test_Type; i++){
    if(INTEGER(test_Type_)[i]==0)
      test_Type[i]="not-defined";
    else if(INTEGER(test_Type_)[i]==1)
      test_Type[i] = "one-at-the-time";
    else if(INTEGER(test_Type_)[i]==2)
      test_Type[i] = "simultaneous";
  }
  
  //interval_Type
  UInt size_interval_Type=Rf_length(interval_Type_);
  interval_Type.resize(size_interval_Type);
  for(UInt i=0; i<size_interval_Type; i++){
    if(INTEGER(interval_Type_)[i]==0)
      interval_Type[i] = "not-defined";
    else if(INTEGER(interval_Type_)[i]==1)
      interval_Type[i] = "one-at-the-time";
    else if(INTEGER(interval_Type_)[i]==2)
      interval_Type[i] = "simultaneous";
    else if(INTEGER(interval_Type_)[i]==3)
      interval_Type[i] = "bonferroni";
  }

  //implementation_Type
  UInt size_implementation_Type=Rf_length(implementation_Type_);
  implementation_Type.resize(size_implementation_Type);
  for(UInt i=0; i<size_implementation_Type; i++){
    if(INTEGER(implementation_Type_)[i]==1)
      implementation_Type[i] = "wald";
    else if(INTEGER(implementation_Type_)[i]==2)
      implementation_Type[i] = "speckman";
    else if(INTEGER(implementation_Type_)[i]==3)
      implementation_Type[i] = "eigen-sign-flip";
    else if(INTEGER(implementation_Type_)[i]==4)
      implementation_Type[i] = "sign-flip";
  }
 
  //exact_Inference
  if(INTEGER(exact_Inference_)[0]==1)
    this->set_exact_inference("exact");
  else
    this->set_exact_inference("non-exact");

  //enhanced_Inference
  if(INTEGER(enhanced_Inference_)[0]==1)
    this->set_enhanced_inference("enhanced");
  else
    this->set_enhanced_inference("classical");

  //coeff_Inference
  UInt n_ = INTEGER(Rf_getAttrib(coeff_Inference_, R_DimSymbol))[0]; // #Rows
  UInt p_ = INTEGER(Rf_getAttrib(coeff_Inference_, R_DimSymbol))[1]; // #Columns
  coeff_Inference.resize(n_, p_);
  for(auto i=0; i<n_; ++i)
    {
      for(auto j=0; j<p_ ; ++j)
	{
	  coeff_Inference(i,j) = REAL(coeff_Inference_)[i+ n_*j];
	}
    }

  //beta_0
  UInt size_beta_0=Rf_length(beta_0_); //We need different sizes for the cases of bad definition
  beta_0.resize(size_beta_0);
  for(UInt i=0;i<size_beta_0;i++){
    beta_0[i]=REAL(beta_0_)[i];
  }

  //f_var
  if(INTEGER(f_Var_)[0]==1)
    this->set_f_Var(true);
  else
    this->set_f_Var(false);
 
  //inference_Quantile
  UInt size_inference_Quantile=Rf_length(inference_Quantile_);
  inference_Quantile.resize(size_inference_Quantile);
  for(UInt i=0;i<size_inference_Quantile;i++){
    inference_Quantile(i)=REAL(inference_Quantile_)[i];
  }

  //inference_Alpha
  UInt size_inference_Alpha=Rf_length(inference_Alpha_);
  inference_Alpha.resize(size_inference_Alpha);
  for(UInt i=0;i<size_inference_Alpha;i++){
    inference_Alpha(i)=REAL(inference_Alpha_)[i];
  }

  //n_flip
  this->set_n_Flip(INTEGER(n_Flip_)[0]);

  //tol_Fspai
  this->set_tol_Fspai(REAL(tol_Fspai_)[0]);

  //definition
  this->set_definition(bool(INTEGER(definition_)[0]));

};

//! Space-only main constructor of the class, with inference for f
InferenceData::InferenceData(SEXP test_Type_, SEXP interval_Type_, SEXP implementation_Type_, SEXP component_Type_,
			     SEXP exact_Inference_, SEXP enhanced_Inference_, SEXP locs_Inference_, SEXP locs_index_Inference_, SEXP coeff_Inference_, SEXP beta_0_, SEXP f_0_, SEXP f_Var_,
			     SEXP inference_Quantile_, SEXP inference_Alpha_, SEXP n_Flip_, SEXP tol_Fspai_, SEXP definition_):InferenceData(test_Type_, interval_Type_, implementation_Type_, exact_Inference_, enhanced_Inference_, coeff_Inference_, beta_0_, f_Var_, inference_Quantile_, inference_Alpha_, n_Flip_, tol_Fspai_, definition_){
 
  //component_Type
  UInt size_component_Type=Rf_length(component_Type_);
  component_Type.resize(size_component_Type);
  for(UInt i=0; i<size_component_Type; i++){
    if(INTEGER(component_Type_)[i]==1)
      component_Type[i] = "parametric";
    else if(INTEGER(component_Type_)[i]==2)
      component_Type[i] = "nonparametric";
    else if(INTEGER(component_Type_)[i]==3)
      component_Type[i] = "both";
  }
  
  //locs_Inference
  UInt m_ = INTEGER(Rf_getAttrib(locs_Inference_, R_DimSymbol))[0]; // #Rows
  UInt d_ = INTEGER(Rf_getAttrib(locs_Inference_, R_DimSymbol))[1]; // #Columns
  locs_Inference.resize(m_, d_);
  for(auto i=0; i<m_; ++i)
    {
      for(auto j=0; j<d_ ; ++j)
	{
	  locs_Inference(i,j) = REAL(locs_Inference_)[i+ m_*j];
	}
    }

  //locs_index_Inference
  UInt size_locs_index=Rf_length(locs_index_Inference_); 
  locs_index_Inference.resize(size_locs_index);
  for(UInt i=0;i<size_locs_index;i++){
    locs_index_Inference[i]=INTEGER(locs_index_Inference_)[i];
  }
  
  //f0_eval
  UInt size_f_0=Rf_length(f_0_); 
  f0_eval.resize(size_f_0);
  for(UInt i=0;i<size_f_0;i++){
    f0_eval[i]=REAL(f_0_)[i];
  }
};

void InferenceData::print_inference_data() const{
  Rprintf("\nInferenceData:\n");
  Rprintf("test_Type:");
  for(UInt i=0; i < test_Type.size(); ++i){
    Rprintf(" %s", test_Type[i].c_str());
  }
  Rprintf("\n");
  Rprintf("interval_Type:");
  for(UInt i=0; i < interval_Type.size(); ++i){
    Rprintf(" %s", interval_Type[i].c_str());
  }
  Rprintf("\n");
  Rprintf("implementation_Type:");
  for(UInt i=0; i < implementation_Type.size(); ++i){
    Rprintf(" %s", implementation_Type[i].c_str());
  }
  Rprintf("\n");
  Rprintf("component_Type:");
  for(UInt i=0; i < component_Type.size(); ++i){
    Rprintf(" %s", component_Type[i].c_str());
  }
  Rprintf("\n");
  
  Rprintf("exact_Inference: %s\n",exact_Inference.c_str());

  Rprintf("enhanced_Inference: %s\n",enhanced_Inference.c_str());
							     
  Rprintf("locs_inference:");
  for(UInt i=0; i < locs_Inference.rows(); ++i){
    for(UInt j=0; j < locs_Inference.cols(); ++j){
      Rprintf(" %f",locs_Inference(i,j));
    }
  }
  Rprintf("\n");

  Rprintf("locs_index_inference: \n");
  for(UInt i=0; i < locs_index_Inference.size(); ++i){
    Rprintf(" %d \n", locs_index_Inference[i]);
  }
  
  
  Rprintf("coeff_inference:");
  for(UInt i=0; i < coeff_Inference.rows(); ++i){
    for(UInt j=0; j < coeff_Inference.cols(); ++j){
      Rprintf(" %f",coeff_Inference(i,j));
    }
  }
  Rprintf("\n");
  Rprintf("beta_0: \n");
  for(UInt i=0; i < beta_0.size(); ++i){
    Rprintf(" %f \n", beta_0(i));
  }

  Rprintf("f0_eval: \n");
  for(UInt i=0; i < f0_eval.size(); ++i){
    Rprintf(" %f \n", f0_eval(i));
  }

  Rprintf("f_var: %d\n",f_Var);
  
  Rprintf("\n");
  Rprintf("inference_Quantile: \n");
  for(UInt i=0; i < inference_Quantile.size(); ++i){
    Rprintf(" %f \n", inference_Quantile(i));
  }
  Rprintf("\n");
  Rprintf("inference_Alpha: \n");
  for(UInt i=0; i < inference_Alpha.size(); ++i){
    Rprintf(" %f \n", inference_Alpha(i));
  }

  Rprintf("n_Flip: %d\n", n_Flip);
  Rprintf("tol_Fspai: %f\n", tol_Fspai);
  Rprintf("definition: %d\n",definition);
};
