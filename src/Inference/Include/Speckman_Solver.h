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
#include "Solver_Base.h"
#include <memory>

// *** Speckman_Solver Class ***
//! Hypothesis testing and confidence intervals using Speckman implementation
/*!
  This class performes hypothesis testing and/or computes confidence intervals using a Speckman approach. It contains a reference to an inverter, that manages to compute the invertion of matrixNoCov in an exact or non-exact way; It contains a reference to an Inference_Carrier object that wraps all the information needed to make inference. There is only one public method that calls the proper private methods to compute what is requested by the user.
*/
template<typename InputHandler>
class Speckman_Solver:public Solver_Base<InputHandler>{
private:
  MatrixXr B;						//!< Matrix Psi*(Psi^t * Psi + lambda*R)^-1*Psi^t 
  MatrixXr Lambda2;   					//!< (I - B)^2
  bool is_Lambda2_computed = false;			//!< Boolean that tells whether Lambda^2 has been computed or not
  MatrixXr V;						//!< Variance-Covariance matrix of the beta parameters
  bool is_V_computed = false;				//!< Boolean that tells whether WLW has been computed or not
  Eigen::PartialPivLU<MatrixXr> WLW_dec; 		//!< Decomposition of [W^t * Lambda^2 * W] 
  bool is_WLW_computed=false; 				//!< Boolean that tells whether Lambda has been computed or not
  void compute_Lambda2(void);				//!< Method used to compute Lambda^2
  void compute_V(void);					//!< Method used to compute V
  void compute_WLW_dec(void); 				//!< Method that computes the decomposition for WLW
  VectorXr compute_beta_hat(void);               	//!< Method used to compute beta estimates for the Speckman test
  VectorXr compute_pvalue(void) override;		//!< Method used to compute the pvalues of the tests 
  MatrixXv compute_CI(void) override;			//!< Method to compute the confidence intervals
  
public:
  // CONSTUCTOR
  Speckman_Solver()=delete;	//The default constructor is deleted
  Speckman_Solver(std::unique_ptr<Inverse_Base> inverter_, const Inference_Carrier<InputHandler> & inf_car_):Solver_Base<InputHandler>(std::move(inverter_), inf_car_){}; 
  
  // GETTERS
  inline const MatrixXr * getLambda2p (void) const {return &this->Lambda2;}     //!< Getter of Lambda2p \return Lambda2p
  inline const MatrixXr * getBp (void) const {return &this->B;}  		//!< Getter of Bp \return Bp
  inline const MatrixXr * getVp (void) const {return &this->V;}     	 	//!< Getter of Vp \ return Vp
  
  void print_for_debug(void) const;
};


#include "Speckman_Solver_imp.h"

#endif
