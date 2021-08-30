#include "../Include/Inference_Data.h"
//! Main constructor of the class
/*!
 \param test_Type_ vector parameter used to define the type of tests that are required (if any)
 \param interval_Type_ vector parameter used to define which type of confidence intervals are required (if any)
 \param implementation_Type_ vector parameter used to define which type of implementations (Wald, Speckman, ESF) are used for the tests and intervals computation
 \param exact_Inference_ parameter for the method used to invert E matrix in Woodbury decomposition for inference
 \param coeff_Inference_ matrix that specifies the linear combinations of the linear parameters to be tested and/or estimated via confidence intervals 
 \param beta_0_ vector for the null hypotesis (if a test is required)
 \param inference_Quantile_ vector parameter containing the quantiles to be used for the computation of the confidence intervals (if interval_type is defined)
 \param n_flip_ parameter that provides the number of sign-flips to be used for the eigen-sign-flip tests (if they are required)
 \param definition_ parameter used to set definition of the InferenceData object
*/
InferenceData::InferenceData(SEXP test_Type_, SEXP interval_Type_, SEXP implementation_Type_,
			     SEXP exact_Inference_, SEXP coeff_Inference_, SEXP beta_0_,
			     SEXP inference_Quantile_, SEXP n_flip_, SEXP definition_){
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
  }
   //exact_Inference
   if(INTEGER(exact_Inference_)[0]==1)
    this->set_exact_inference("exact");

  else
    this->set_exact_inference("non-exact");

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

 //inference_Quantile
 UInt size_inference_Quantile=Rf_length(inference_Quantile_);
 inference_Quantile.resize(size_inference_Quantile);
 for(UInt i=0;i<size_inference_Quantile;i++){
   inference_Quantile(i)=REAL(inference_Quantile_)[i];
 }

 //definition
 this->set_definition(bool(INTEGER(definition_)[0]));

 //n_perm
 this->set_n_flip(INTEGER(n_flip_)[0]);
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
  Rprintf("exact_Inference: %s\n",exact_Inference.c_str());
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
  Rprintf("\n");
  Rprintf("inference_Quantile: \n");
  for(UInt i=0; i < inference_Quantile.size(); ++i){
    Rprintf(" %f \n", inference_Quantile(i));
  }
  Rprintf("\n");
  Rprintf("n_flip: %d\n", n_flip);
  Rprintf("definition: %d\n",definition);
};
