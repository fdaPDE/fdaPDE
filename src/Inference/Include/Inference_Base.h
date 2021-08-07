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
//! Hypothesis testing and confidence intervals abstract base class
/*!
  This template class provides the basic tools to perform hypothesis testing and/or compute confidence intervals. It contains a shared pointer to an inverter, that manages the invertion of matrixNoCov in an exact or non-exact way; It contains a reference to an Inference_Carrier object that wraps all the information needed to make inference; It contains an integer pos_impl that indicates the index position of the current inferential procedure to be carried out. It is needed to take the correct information from the vector parameters in the Inference_Carrier object. There are two pure virtual protected methods for the computation of p-values and confidence intervals; there is a main public method that calls the proper functions according to the current test and interval types; then there is a public setter for the index position that is needed when multiple inferential procedures are required and a virtual method that will be actually implemented only in the derived Wald class used to compute exact GCV. 
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
