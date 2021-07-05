#ifndef __SPECKMAN_H__
#define __SPECKMAN_H__

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

// *** Speckman Class ***
//! Hypothesis testing and confidence intervals using Speckman implementation
/*!
  This class performes hypothesis testing and/or computes confidence intervals using a Speckman approach. It contains a reference to an inverter, that manages to compute the invertion of matrixNoCov in an exact or non-exact way; It contains a reference to an Inference_Carrier object that wraps all the information needed to make inference. There is only one public method that calls the proper private methods to compute what is requested by the user.
*/
template<typename InputHandler, typename MatrixType>
class Speckman_Base:public Inference_Base<InputHandler, MatrixType>{
protected: 
  MatrixXr Lambda2;   					//!< (I - Psi*(Psi^t * Psi + lambda*R)^-1*Psi^t)^2
  bool is_Lambda2_computed = false;			//!< Boolean that tells whether Lambda^2 has been computed or not
  MatrixXr V;						//!< Variance-Covariance matrix of the beta parameters
  bool is_V_computed = false;				//!< Boolean that tells whether WLW has been computed or not
  Eigen::PartialPivLU<MatrixXr> WLW_dec; 		//!< Decomposition of [W^t * Lambda^2 * W] 
  bool is_WLW_computed=false; 				//!< Boolean that tells whether Lambda has been computed or not
  virtual void compute_Lambda2(void) = 0;		//!< Method used to compute Lambda^2, either in an exact or non-exact way
  void compute_V(void);					//!< Method used to compute V
  void compute_WLW_dec(void); 				//!< Method that computes the decomposition for WLW
  VectorXr compute_beta_hat(void);               	//!< Method used to compute beta estimates for the Speckman test
  VectorXr compute_pvalue(void) override;		//!< Method used to compute the pvalues of the tests 
  MatrixXv compute_CI(void) override;			//!< Method to compute the confidence intervals
  
public:
  // CONSTUCTOR
  Speckman_Base()=delete;	//The default constructor is deleted
  Speckman_Base(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Inference_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; 
  
  // DESTRUCTOR
  virtual ~ Speckman_Base() {};

  // GETTERS
  inline const MatrixXr * getLambda2p (void) const {return &this->Lambda2;}     //!< Getter of Lambda2p \return Lambda2p
  inline const MatrixXr * getVp (void) const {return &this->V;}     	 	//!< Getter of Vp \ return Vp
  
  void print_for_debug(void) const;
};


template<typename InputHandler, typename MatrixType>
class Speckman_Exact:public Speckman_Base<InputHandler, MatrixType>{
private: 
  void compute_Lambda2(void) override;
public:
  // CONSTUCTOR
  Speckman_Exact()=delete;	//The default constructor is deleted
  Speckman_Exact(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Speckman_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; 

};

template<typename InputHandler, typename MatrixType>
class Speckman_Non_Exact:public Speckman_Base<InputHandler, MatrixType>{
private: 
  void compute_Lambda2(void) override;
public:
  // CONSTUCTOR
  Speckman_Non_Exact()=delete;	//The default constructor is deleted
  Speckman_Non_Exact(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Speckman_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; 

};


#include "Speckman_imp.h"

#endif
