#include "Solver_Base.h"

template<typename InputHandler>
MatrixXv Solver_Base<InputHandler>::compute_inference_output(void){
  // declare the result Matrix of vectors to be returned
  MatrixXv result;
  
  // get the test_type and interval_type
  std::string test_type = inf_car.getInfData()->get_test_type();
  std::string interval_type = inf_car.getInfData()->get_interval_type();
  
  // if test_type is not defined, only intervals are required
  if(test_type == "not-defined"){
    result = this->compute_CI();
    return result;
  }
  // if interval_type is not defined, only test is required
  if(interval_type == "not-defined"){
    result.resize(1,1);
    result(0) = this->compute_pvalue(); 
    return result;
  }
  // else, both are required
  else{
    UInt q = inf_car.getInfData()->get_coeff_inference().rows();
    result.resize(1, q+1);
    result(0) = this->compute_pvalue();
    result.rightCols(q) = this->compute_CI();
    return result;
  }
};

