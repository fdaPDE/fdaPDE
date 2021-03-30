#ifndef __SPECKMAN_SOLVER_H__
#define __SPECKMAN_SOLVER_H__

// HEADERS
#include "../../FdaPDE.h"
#include "../../Regression/Include/Mixed_FE_Regression.h"
#include "../../Regression/Include/Regression_Data.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"
#include "../../Lambda_Optimization/Include/Solution_Builders.h"
#include "Inference_Data.h"
#include "Inference_Carrier.h"
#include "Inverter.h"

// *** Speckman_Solver Class ***
//! Hypothesis testing and confidence intervals using Speckman implementation
/*!
  This class performes hypothesis testing and/or computes confidence intervals using a Speckman approach. It contains a reference to an inverter, that manages to compute the invertion of matrixNoCov in an exact or non-exact way; It contains a reference to an Inference_Carrier object that wraps all the information needed to make inference. There is only one public method that calls the proper private methods to compute what is requested by the user.
*/
template<typename InputHandler>
class Speckman_Solver{
private:
  Inverse_Base & inverter; 				//!< Inverter object that computes the inverse of matrixNoCov in exact/non-exact way
  const Inference_Carrier<InputHandler> & inf_car;	//!< Inference carrier that contains all the information needed for inference 
  MatrixXr B;						//!< Matrix Psi*(Psi^t * Psi + lambda*R)^-1*Psi^t 
  MatrixXr Lambda;   					//!< I - B
  bool is_Lambda_computed = false;			//!< Boolean that tells whether Lambda has been computed or not
  MatrixXr V;						//!< Variance-Covariance matrix of the beta parameters
  bool is_V_computed = false;				//!< Boolean that tells whether V has been computed or not
  void compute_Lambda(void);				//!< Method used to compute Lambda
  void compute_V(void);					//!< Method used to compute V
  VectorXr compute_pvalue(void);			//!< Method used to compute the pvalues of the tests 
  MatrixXv compute_CI(void);				//!< Method to compute the confidence intervals
  
public:
  // CONSTUCTOR
  Speckman_Solver()=delete;	//The default constructor is deleted
  Speckman_Solver(Inverse_Base & inverter_, const Inference_Carrier<InputHandler> & inf_car_):inverter(inverter_), inf_car(inf_car_){}; 
  
  // GETTERS
  inline const MatrixXr * getLambdap (void) const {return &this->Lambda;}      	//!< Getter of Lambdap \return Lambdap
  inline const MatrixXr * getBp (void) const {return &this->B;}  		//!< Getter of Bp \return Bp
  inline const MatrixXr * getVp (void) const {return &this->V;}     	 	//!< Getter of Vp \ return Vp
  
  //!< public method that calls the requested functions according to test_type and interval_type
  MatrixXv compute_inference_output (void);
  void print_for_debug(void) const;
};


#include "Speckman_Solver_imp.h"

#endif
