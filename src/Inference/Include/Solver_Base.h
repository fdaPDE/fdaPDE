#ifndef __SOLVER_BASE_H__
#define __SOLVER_BASE_H__

// HEADERS
#include "../../FdaPDE.h"
#include "../../Regression/Include/Mixed_FE_Regression.h"
#include "../../Regression/Include/Regression_Data.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"
#include "../../Lambda_Optimization/Include/Solution_Builders.h"
#include "Inference_Data.h"
#include "Inference_Carrier.h"
#include "Inverter.h"
#include <memory>

// *** Solver_Base Class ***
//! Hypothesis testing and confidence intervals base class
/*!
  This class provides the basic tools to perform hypothesis testing and/or compute confidence intervals. It contains a reference to an inverter, that manages to compute the invertion of matrixNoCov in an exact or non-exact way; It contains a reference to an Inference_Carrier object that wraps all the information needed to make inference. There is only one public method that calls the proper private methods to compute what is requested by the user.
*/
template<typename InputHandler>
class Solver_Base{
protected:
  std::unique_ptr<Inverse_Base> inverter; 		//!< Pointer to inverter object that computes the inverse of matrixNoCov in exact/non-exact way
  const Inference_Carrier<InputHandler> & inf_car;	//!< Inference carrier that contains all the information needed for inference 
  virtual VectorXr compute_pvalue(void) = 0;		//!< Pure virtual method used to compute the pvalues of the tests 
  virtual MatrixXv compute_CI(void) = 0;		//!< Pure virtual method to compute the confidence intervals
  
public:
  // CONSTUCTOR
  Solver_Base()=delete;	//The default constructor is deleted
  Solver_Base(std::unique_ptr<Inverse_Base> inverter_, const Inference_Carrier<InputHandler> & inf_car_):inverter(std::move(inverter_)), inf_car(inf_car_){}; 
 
  //!< public method that calls the requested functions according to test_type and interval_type
  MatrixXv compute_inference_output (void);

  // DESTRUCTOR
  virtual ~Solver_Base(){};
};


#include "Solver_Base_imp.h"

#endif
