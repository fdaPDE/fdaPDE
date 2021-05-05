#include "Inference_Base.h"

template<typename InputHandler>
MatrixXv Inference_Base<InputHandler>::compute_inference_output(void){
  // declare the result Matrix of vectors to be returned
  MatrixXv result;
  
  // get the test_type and interval_type
  std::string test_type = inf_car.getInfData()->get_test_type()[this->pos_impl];
  std::string interval_type = inf_car.getInfData()->get_interval_type()[this->pos_impl];

  // Preallocate space for any case
  UInt p = inf_car.getInfData()->get_coeff_inference().rows();
  result.resize(1, p+1);
  
  // if test_type is not defined, only intervals are required
  if(test_type == "not-defined"){
    result(0).resize(p);
    for(k=0;k<p;k++){
    result(0)(k) = 10e20; // Default value (unfeasible)
    }
    result.rightCols(p) = this->compute_CI();
  }
  // if interval_type is not defined, only test is required
  if(interval_type == "not-defined"){
    result(0) = this->compute_pvalue();
    for(UInt k=0;k<p;k++){
    result(k+1).resize(3);
    result(k+1)(0)=10e20;  // default value (unfeasible)
    result(k+1)(1)=10e20;  // default value (unfeasible)
    result(k+1)(2)=10e20;  // default value (unfeasible)
    } 
  }
  // else, both are required
  else{
    result(0) = this->compute_pvalue();
    result.rightCols(p) = this->compute_CI();
  }
return result;
};

