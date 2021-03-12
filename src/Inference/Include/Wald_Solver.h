#ifndef __WALD_SOLVER_H__
#define __WALD_SOLVER_H__

// HEADERS
#include "../../FdaPDE.h"
#include "../../Regression/Include/Mixed_FE_Regression.h"
#include "../../Regression/Include/Regression_Data.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"
#include "../../Lambda_Optimization/Include/Solution_Builders.h"
#include "Inference_Data.h"
#include "Inference_Carrier.h"
#include "Inverter.h"

// *** Wald_Solver Class ***
//! Hypothesis testing and confidence intervals using Wald implementation
/*!
  This class performes hypothesis testing and/or computes confidence intervals using a Wald-type approach. It contains a reference to an inverter, that manages to compute the invertion of matrixNoCov in an exact or non-exact way; It contains a reference to an Inference_Carrier object that wraps all the information needed to make inference. There is only one public method that calls the proper private methods to compute what is requested by the user.
*/
template<typename InputHandler>
class Wald_Solver{
private:
  Inverse_Base & inverter; 				//!< Inverter object that computes the inverse of matrixNoCov in exact/non-exact way
  const Inference_Carrier<InputHandler> & inf_car;	//!< Inference carrier that contains all the information needed for inference 
  MatrixXr S;						//!< Smoothing matrix 
  MatrixXr S_t;   					//!< Transpose of the smoothing matrix
  bool is_S_computed = false;				//!< Boolean that tells whether S has been computed or not
  MatrixXr V;						//!< Variance-Covariance matrix of the beta parameters
  bool is_V_computed = false;				//!< Boolean that tells whether V has been computed or not
  void compute_S(void);					//!< Method used to compute S
  void compute_V(void);					//!< Method used to compute V
  VectorXr compute_pvalue(void);			//!< Method used to compute the pvalues of the tests 
  MatrixXv compute_CI(void);				//!< Method to compute the confidence intervals

public:
  // CONSTUCTOR
  Wald_Solver()=delete;	//The default constructor is deleted
  Wald_Solver(Inverse_Base & inverter_, const Inference_Carrier<InputHandler> & inf_car_):inverter(inverter_), inf_car(inf_car_){}; 
  
  // GETTERS
  inline const MatrixXr * getSp (void) const {return &this->S;}      //!< Getter of Sp \return Sp
  inline const MatrixXr * getS_tp (void) const {return &this->S_t;}  //!< Getter of Sp_tp \return Sp_tp
  inline const MatrixXr * getVp (void) const {return &this->V;}      //!< Getter of Vp \ return Vp

  //!< public method that calls the requested functions according to test_type and interval_type
  MatrixXv compute_inference_output (void);
  void print_for_debug(void) const;
};


#include "Wald_Solver_imp.h"

#endif
