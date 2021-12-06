#include "Inference_Base.h"

template<typename InputHandler, typename MatrixType>
MatrixXv Inference_Base<InputHandler, MatrixType>::compute_inference_output(void){
  // declare the result Matrix of vectors to be returned
  MatrixXv result;
  
  // get the test_type and the interval_type 
  std::string test_type = inf_car.getInfData()->get_test_type()[this->pos_impl];
  std::string interval_type = inf_car.getInfData()->get_interval_type()[this->pos_impl];
 
  // Preallocate space for any case
  UInt p = inf_car.getInfData()->get_coeff_inference().rows();
  UInt n_loc = inf_car.getN_loc();
  UInt result_dim = (p > n_loc) ? p : n_loc;   
  result.resize(2, result_dim+1);
  
  // if test_type is not defined, only intervals are required
  if(test_type == "not-defined"){
    // beta pvalues
    result(0).resize(p);
    for(UInt k=0;k<p;k++){
      result(0)(k) = 10e20; // Default value (unfeasible)
    }
    // f pvalue
    result(1,0).resize(1);
    result(1,0)(0) = 10e20;

    result.rightCols(result_dim) = this->compute_CI();
  }

  // if interval_type is not defined, only test is required
  if(interval_type == "not-defined"){
    result.leftCols(1) = this->compute_pvalue();
    for(UInt k=0;k<result_dim;k++){
    for(UInt i=0; i<result.rows(); ++i){
      result(i,k+1).resize(3);
      result(i,k+1)(0)=10e20;  // default value (unfeasible)
      result(i,k+1)(1)=10e20;  // default value (unfeasible)
      result(i,k+1)(2)=10e20;  // default value (unfeasible)
    } 
    }
  }
  // else, both are required
  else{
    result.leftCols(1) = this->compute_pvalue();
    result.rightCols(result_dim) = this->compute_CI();
  }
  return result;
};

template<typename InputHandler, typename MatrixType>
MatrixXv Inference_Base<InputHandler, MatrixType>::compute_pvalue(void){
  // declare the object that will store the p-values
  MatrixXv result;
  result.resize(2,1);
  
  // inference only on beta
  if(this->inf_car.getInfData()->get_component_type()[this->pos_impl] == "parametric"){
  	result(0) = this->compute_beta_pvalue();
        result(1).resize(1);
        result(1)(0) = 10e20; // default value (unfeasible)
  }
  
  // inference only on f
  if(this->inf_car.getInfData()->get_component_type()[this->pos_impl] == "nonparametric"){
	result(1).resize(1);	
	result(1)(0) = this->compute_f_pvalue();
        
	UInt p = inf_car.getInfData()->get_coeff_inference().rows();
   	result(0).resize(p);
    	for(UInt k=0;k<p;k++){
      	  result(0)(k) = 10e20; // Default value (unfeasible)
    	}
  }

  // inference on both beta and f
  if(this->inf_car.getInfData()->get_component_type()[this->pos_impl] == "both"){
  	result(0) = this->compute_beta_pvalue();
        result(1).resize(1);
        result(1)(0) = this->compute_f_pvalue();
  }

  return result; 

};

template<typename InputHandler, typename MatrixType>
MatrixXv Inference_Base<InputHandler, MatrixType>::compute_CI(void){
  // declare the object that will store the intervals
  MatrixXv result;
  UInt p = inf_car.getInfData()->get_coeff_inference().rows();
  UInt n_loc = inf_car.getN_loc();
  UInt res_dim = (p > n_loc) ? p : n_loc; 
  
  result.resize(2,res_dim);
  
  // inference only on beta
  if(this->inf_car.getInfData()->get_component_type()[this->pos_impl] == "parametric"){
  	result.row(0).leftCols(p) = this->compute_beta_CI();
        // for safety
        for(UInt j=0; j < res_dim; ++j){
          	result(1,j).resize(3);
                result(1,j)(0) = 10e20;
                result(1,j)(1) = 10e20;
                result(1,j)(2) = 10e20;
          }
        if(res_dim != p){
          for(UInt j=p; j < res_dim; ++j){
	  	result(0,j).resize(3);
                result(0,j)(0) = 10e20;
                result(0,j)(1) = 10e20;
                result(0,j)(2) = 10e20;
          }
        } 
  }
  
  // inference only on f
  if(this->inf_car.getInfData()->get_component_type()[this->pos_impl] == "nonparametric"){
	result.row(1).leftCols(n_loc) = this->compute_f_CI();
        // for safety
        for(UInt j=0; j < res_dim; ++j){
          	result(0,j).resize(3);
                result(0,j)(0) = 10e20;
                result(0,j)(1) = 10e20;
                result(0,j)(2) = 10e20;
          }
        if(res_dim != n_loc){
          for(UInt j=n_loc; j < res_dim; ++j){
	  	result(1,j).resize(3);
                result(1,j)(0) = 10e20;
                result(1,j)(1) = 10e20;
                result(1,j)(2) = 10e20;
          }
        } 
  }

  // inference on both beta and f 
  if(this->inf_car.getInfData()->get_component_type()[this->pos_impl] == "both"){
  	result.row(0).leftCols(p) = this->compute_beta_CI();
        result.row(1).leftCols(n_loc) = this->compute_f_CI();
        // for safety
        if(res_dim != p){
          for(UInt j=p; j < res_dim; ++j){
		result(0,j).resize(3);
                result(0,j)(0) = 10e20;
                result(0,j)(1) = 10e20;
                result(0,j)(2) = 10e20;
	  }
        } else {
           for(UInt j=n_loc; j < res_dim; ++j){
		result(1,j).resize(3);
                result(1,j)(0) = 10e20;
                result(1,j)(1) = 10e20;
                result(1,j)(2) = 10e20;
           }
        }        
  }

  return result; 
};

template<typename InputHandler, typename MatrixType>
VectorXr Inference_Base<InputHandler, MatrixType>::compute_f_var(void){
  UInt n_obs = inf_car.getN_obs();
  VectorXr result = VectorXr::Zero(n_obs);

  return result;
};
