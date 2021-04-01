#include "../Include/Inference_Data.h"
//! Main constructor of the class
/*!
 \param test_Type_ parameter used to define the type of test is required (if any)
 \param interval_Type_ parameter used to define which type of confidence interval is required (if any)
 \param implementation_Type_ parameter used to define which type of implementation is used for the test and intervals computation
 \param exact_Inference_ parameter for the method used to invert E matrix in Woodbury decomposition for inference
 \param coeff_Inference_ vector that specifies the linear combinations of the linear parameters to be tested and/or estimated via confidence intervals. 
 \param beta_0_ vector for the null hypotesis (if test is defined)
 \param inference_level_ parameter used to set the significance level of the confidence intervals (if interval_type is set)
 \param definition_ parameter used to set definition of the InferenceData object. 
*/
InferenceData::InferenceData(SEXP test_Type_, SEXP interval_Type_, SEXP implementation_Type_,
			     SEXP exact_Inference_, SEXP coeff_Inference_, SEXP beta_0_,
			     SEXP inference_Level_, SEXP n_perm_, SEXP definition_){
  //test_Type
  if(INTEGER(test_Type_)[0]==0)
    this->set_test_type("not-defined");

  else if(INTEGER(test_Type_)[0]==1)
    this->set_test_type("one-at-the-time");

  else if(INTEGER(test_Type_)[0]==2)
    this->set_test_type("simultaneous");
  
  //interval_Type
  if(INTEGER(interval_Type_)[0]==0)
    this->set_interval_type("not-defined");

  else if(INTEGER(interval_Type_)[0]==1)
    this->set_interval_type("one-at-the-time");

  else if(INTEGER(interval_Type_)[0]==2)
    this->set_interval_type("simultaneous");

  else if(INTEGER(interval_Type_)[0]==3)
    this->set_interval_type("bonferroni");

  //implementation_Type
  if(INTEGER(implementation_Type_)[0]==1)
    this->set_implementation_type("wald");

  else if(INTEGER(implementation_Type_)[0]==2)
    this->set_implementation_type("speckman");

  else if(INTEGER(implementation_Type_)[0]==3)
    this->set_implementation_type("permutational");

   //exact_Inference
   if(INTEGER(exact_Inference_)[0]==1)
    this->set_exact_inference(true);

  else
    this->set_exact_inference(false);

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

 //inference_Level
 this->set_inference_level(REAL(inference_Level_)[0]);

 //definition
 this->set_definition(bool(INTEGER(definition_)[0]));

 //n_perm
 this->set_n_perm(INTEGER(n_perm_)[0]);
};

void InferenceData::print_inference_data() const{
  Rprintf("\nInferenceData:\n");
  Rprintf("test_Type: %s\n", test_Type.c_str());
  Rprintf("interval_Type: %s\n", interval_Type.c_str());
  Rprintf("implementation_Type: %s\n",implementation_Type.c_str());
  Rprintf("exact_Inference: %d\n",exact_Inference);
  Rprintf("coeff_inference:");
  for(UInt i=0; i < coeff_Inference.rows(); ++i){
    for(UInt j=0; j < coeff_Inference.cols(); ++j){
      Rprintf(" %f",coeff_Inference(i,j));
  }
  }
  Rprintf("\n");
  Rprintf("beta_0:");
  for(UInt i=0; i < beta_0.size(); ++i){
    Rprintf(" %f", beta_0(i));
  }
  Rprintf("\n");
  Rprintf("inference_Level: %f\n",inference_Level);
  Rprintf("n_perm: %d\n", n_perm);
  Rprintf("definition: %d\n",definition);
};
