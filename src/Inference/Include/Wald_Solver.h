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
#include "Solver_Base.h"
#include <memory>

// *** Wald_Solver Class ***
//! Hypothesis testing and confidence intervals using Wald implementation
/*!
  This class performes hypothesis testing and/or computes confidence intervals using a Wald-type approach. It contains a reference to an inverter, that manages to compute the invertion of matrixNoCov in an exact or non-exact way; It contains a reference to an Inference_Carrier object that wraps all the information needed to make inference. There is only one public method that calls the proper private methods to compute what is requested by the user.
*/
template<typename InputHandler>
class Wald_Solver:public Solver_Base<InputHandler>{
private:
  MatrixXr S;						//!< Smoothing matrix 
  MatrixXr S_t;   					//!< Transpose of the smoothing matrix
  bool is_S_computed = false;				//!< Boolean that tells whether S has been computed or not
  MatrixXr V;						//!< Variance-Covariance matrix of the beta parameters
  bool is_V_computed = false;				//!< Boolean that tells whether V has been computed or not
  void compute_S(void);					//!< Method used to compute S
  void compute_V(void);					//!< Method used to compute V
  VectorXr compute_pvalue(void) override;		//!< Method used to compute the pvalues of the tests 
  MatrixXv compute_CI(void) override;			//!< Method to compute the confidence intervals
  Real compute_sigma_hat_sq(void) const;                //!< Method to compute the estimator of the variance of the residuals 
  
public:
  // CONSTUCTOR
  Wald_Solver()=delete;	//The default constructor is deleted
  Wald_Solver(std::unique_ptr<Inverse_Base> inverter_, const Inference_Carrier<InputHandler> & inf_car_):Solver_Base<InputHandler>(std::move(inverter_), inf_car_){}; 
  Wald_Solver(Wald_Solver & rhs) = delete; //The default copy constructor is deleted
  Wald_Solver(Wald_Solver && rhs):inverter(std::move(rhs.inverter)), inf_car(rhs.inf_car), S(std::move(rhs.S)), S_t(std::move(rhs.S_t)), is_S_computed(rhs.is_S_computed), V(std::move(rhs.V)), is_V_computed(rhs.is_V_computed){}; //Definition of the move constructor
  Wald_Solver & operator=(Wald_Solver && rhs) = delete; //The move assignment operator is deleted
  
  
  // GETTERS
  inline const MatrixXr * getSp (void) const {return &this->S;}      //!< Getter of Sp \return Sp
  inline const MatrixXr * getS_tp (void) const {return &this->S_t;}  //!< Getter of Sp_tp \return Sp_tp
  inline const MatrixXr * getVp (void) const {return &this->V;}      //!< Getter of Vp \ return Vp
  
  void print_for_debug(void) const;
};


#include "Wald_Solver_imp.h"

#endif
