#include "InferenceData.h"
//! Main constructor of the class
/*!
 \param test_Type_ parameter used to define the typr of test is rquired (if any)
 \param interval_Type_ parameter used to define which type of confidence interval is required (if any)
 \param implementation_Type_ parameter used to define which type of implementation is used for the test anf interval computation
 \param exact_Inference_ parameter for the method used to invert E matrix in woodbury decomposition for inference
 \param coeff_Inference_ vector tha specifies which covariates are taken into account in the analysis
 \param beta_0_ vector for null hypotesis (if test is defined)
 \param beta_1_ vector for alternative hypotesis (if test is set to power)
 \param inference_level_ parameter used to set the significance level of the confidence interval (if interval_type is set)
 \param definition_ parameter used to set definition
*/
InferenceData::InferenceData(SEXP test_Type_, SEXP interval_Type_, SEXP implementation_Type_,
			     SEXP exact_Inference_, SEXP coeff_Inference_, SEXP beta_0_,
			     SEXP beta_1_, SEXP inference_Level_, SEXP definition_){
  //test_Type
  if(INTEGER(test_Type_)[0]==0)
    this->set_test_type("not-defined");

  else if(INTEGER(test_Type_)[0]==1)
    this->set_test_type("pvalue");

  else if(INTEGER(test_Type_)[0]==2)
    this->set_test_type("power");
  
  //interval_Type
  if(INTEGER(interval_Type_)[0]==0)
    this->set_interval_type("not-defined");

  else if(INTEGER(interval_Type_)[0]==1)
    this->set_interval_type("one-at-the-time");

  else if(INTEGER(interval_Type_)[0]==2)
    this->set_interval_type("bonferroni");

  //implementation_Type
   if(INTEGER(implementation_Type_)[0]==0)
   this->set_implementation_type("not-defined");

   else if(INTEGER(implementation_Type_)[0]==1)
    this->set_implementation_type("wald");

  else if(INTEGER(implementation_Type_)[0]==2)
    this->set_implementation_type("sandwich");

  else if(INTEGER(implementation_Type_)[0]==3)
    this->set_implementation_type("permutational");

   //exact_Inference
 if(INTEGER(exact_Inference_)[0]==1)
    this->set_exact_inference(true);

  else
    this->set_exact_inference(false);

 //coeff_Inference
 UInt size_coeff=INTEGER(Rf_getAttrib(coeff_Inference_, R_DimSymbol))[0];
 coeff_Inference.resize(size);
 for(auto i=0;i<size_coeff;i++){
   if(INTEGER(coeff_Inference_)[i]==0)
     coeff_Inference[i]=false;
   else
     coeff_Inference[i]=true;
 }

 //beta_0
 UInt size_beta_0=INTEGER(Rf_getAttrib(beta_0_, R_DimSymbol))[0]; //We need different sizes for the cases of bad definition
 beta_0.resize(size_beta_0);
 for(auto i=0;i<size;i++){
   beta_0[i]=REAL(beta_0_)[i];
 }

 //beta_1
 UInt size_beta_1=INTEGER(Rf_getAttrib(beta_1_, R_DimSymbol))[0]; //We need different sizes for the cases of bad definition
 beta_1.resize(size_beta_1);
 for(auto i=0;i<size;i++){
   beta_1[i]=REAL(beta_1_)[i];
 }

 //inference_Level
 this->set_inference_level(REAL(inference_Level_)[0]);

 //definition
 this->set_definition(bool(INTEGER(definition_)[0])); 
};

void print_inference_data() const{
  Rprintf("\nInferenceData:\n");
  Rprintf("test_Type: %s\n", test_Type.c_str());
  Rprintf("interval_Type: %s\n", interval_Type.c_str());
  Rprintf("implementation_Type: %s\n",inplementation_Type.c_str());
  Rprintf("exact_Inference: %d\n",exact_Inference);
  Rprintf("coeff_inference:");
  for(auto value : coeff_Inference){
    Rprintf(" %d",value);
  }
  Rprintf("\n");
  Rprintf("beta_0:");
  for(auto value : beta_0){
    Rprintf(" %d",value);
  }
  Rprintf("\n");
  Rprintf("beta_1:");
  for(auto value : beta_1){
    Rprintf(" %d",value);
  }
  Rprintf("\n");
  Rprintf("inference_Level: %f\n",inference_Level);
  Rprintf("definition: %d\n",definition);
}
