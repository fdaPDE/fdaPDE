#ifndef __INFERENCE_BASE_H__
#define __INFERENCE_BASE_H__

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

// *** Inference_Base Class ***
//! Hypothesis testing and confidence intervals base class
/*!
  This class provides the basic tools to perform hypothesis testing and/or compute confidence intervals. It contains a reference to an inverter, that manages to compute the invertion of matrixNoCov in an exact or non-exact way; It contains a reference to an Inference_Carrier object that wraps all the information needed to make inference. There is only one public method that calls the proper private methods to compute what is requested by the user.
*/
template<typename InputHandler, typename MatrixType>
class Inference_Base{
protected:
  std::shared_ptr<Inverse_Base<MatrixType>> inverter = nullptr;     //!< Pointer to inverter object that computes the inverse of matrixNoCov in exact/non-exact way
  const Inference_Carrier<InputHandler> & inf_car;	//!< Inference carrier that contains all the information needed for inference 
  UInt pos_impl;					//!< Index that gives the position in all the vectors in infecenceData object
  virtual VectorXr compute_pvalue(void) = 0;		//!< Pure virtual method used to compute the pvalues of the tests 
  virtual MatrixXv compute_CI(void) = 0;		//!< Pure virtual method to compute the confidence intervals
  

public:
  // CONSTUCTOR
  Inference_Base()=delete;	//The default constructor is deleted
  Inference_Base(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):inverter(inverter_), inf_car(inf_car_), pos_impl(pos_impl_){}; 
  
  //!< Public method that calls the requested functions according to test_type and interval_type
  MatrixXv compute_inference_output (void);

  //!< Public setter for pos_impl, needed when multiple tests are required
  inline void setpos_impl (UInt pos_impl_){this->pos_impl=pos_impl_;};
  
  //!< Virtual public method that computes exact GCV, implemented only for Wald 
  inline virtual Real compute_GCV_from_inference(void) const {return 0;};

  // DESTRUCTOR
  virtual ~Inference_Base(){};
};


#include "Inference_Base_imp.h"

#endif
