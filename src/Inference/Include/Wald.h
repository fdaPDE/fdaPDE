#ifndef __WALD_H__
#define __WALD_H__

// HEADERS
#include "../../FdaPDE.h"
#include "../../Regression/Include/Mixed_FE_Regression.h"
#include "../../Regression/Include/Regression_Data.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"
#include "../../Lambda_Optimization/Include/Solution_Builders.h"
#include "Inference_Data.h"
#include "Inference_Carrier.h"
#include "Inverter.h"
#include "Inference_Base.h"
#include <memory>

// *** Wald Class ***
//! Hypothesis testing and confidence intervals using Wald implementation
/*!
  This class performes hypothesis testing and/or computes confidence intervals using a Wald-type approach. It contains a reference to an inverter, that manages to compute the invertion of matrixNoCov in an exact or non-exact way; It contains a reference to an Inference_Carrier object that wraps all the information needed to make inference. There is only one public method that calls the proper private methods to compute what is requested by the user.
*/
template<typename InputHandler>
class Wald:public Inference_Base<InputHandler>{
private:
  MatrixXr S;						//!< Smoothing matrix 
  MatrixXr S_t;   					//!< Transpose of the smoothing matrix
  Real tr_S=0; 						//!< Trace of smoothing matrix, needed for the variance-covariance matrix (V) and eventually GCV computation
  Real sigma_hat_sq; 					//!< Estimator for the variance of the residuals (SSres/(n_obs-(q+tr_S)))
  bool is_S_computed = false;				//!< Boolean that tells whether S has been computed or not
  MatrixXr V;						//!< Variance-Covariance matrix of the beta parameters
  bool is_V_computed = false;				//!< Boolean that tells whether V has been computed or not
  void compute_S(void);					//!< Method used to compute S
  void compute_V(void);					//!< Method used to compute V
  VectorXr compute_pvalue(void) override;		//!< Method used to compute the pvalues of the tests 
  MatrixXv compute_CI(void) override;			//!< Method to compute the confidence intervals
  void compute_sigma_hat_sq(void);                      //!< Method to compute the estimator of the variance of the residuals 
  
public:
  // CONSTUCTOR
  Wald()=delete;	//The default constructor is deleted
  Wald(std::unique_ptr<Inverse_Base> inverter_, const Inference_Carrier<InputHandler> & inf_car_):Inference_Base<InputHandler>(std::move(inverter_), inf_car_){}; 
  Wald(Wald & rhs) = delete; //The default copy constructor is deleted
  inline Wald(Wald && rhs):S(std::move(rhs.S)), S_t(std::move(rhs.S_t)), is_S_computed(rhs.is_S_computed), V(std::move(rhs.V)), is_V_computed(rhs.is_V_computed){this->inverter=std::move(rhs.inverter); this->inf_car=rhs.inf_car;}; //Definition of the move constructor

  Wald & operator=(Wald && rhs) = delete; //The move assignment operator is deleted

  Real compute_GCV_from_inference(void) const override; //!< Needed to compute exact GCV in case Wald test is required and GCV exact is not provided by lambda optimization (Run after S computation)
  
  
  // GETTERS
  inline const MatrixXr * getSp (void) const {return &this->S;}      //!< Getter of Sp \return Sp
  inline const MatrixXr * getS_tp (void) const {return &this->S_t;}  //!< Getter of Sp_tp \return Sp_tp
  inline const MatrixXr * getVp (void) const {return &this->V;}      //!< Getter of Vp \ return Vp
  
  void print_for_debug(void) const;
};


#include "Wald_imp.h"

#endif
