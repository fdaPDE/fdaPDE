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

// *** Wald_Base Class ***
//! Hypothesis testing and confidence intervals using Wald implementation
/*!
  This template class is an abstract base class to perform hypothesis testing and/or compute confidence intervals using a Wald-type approach. Beyond all the objects and methods inherited from the abstract base inference class, it stores the smoothing matrix S, its trace, the estimator of the residuals variance, the variance-covariance matrix of the beta parameters V, alongside with some convenient boolean objects. It overrides the methods that specify how to compute the p-values and the confidence intervals, according to the Wald apporach. It has a pure virtual method for the computation of the smoothing matrix S, since it relies on the inversion of MatrixNoCov in an exact or non-exact way. Moreover it also overrides the method for the computation of the exact GCV, since it can be computed in a straight-forward way after having computed the smoothing matrix S and its trace. 
*/
template<typename InputHandler, typename MatrixType>
class Wald_Base:public Inference_Base<InputHandler, MatrixType>{
protected:
  MatrixXr S;						//!< Smoothing matrix 
  MatrixXr Partial_S;                                   //!< (Psi^t*Q*Psi + lambda*P)^-1 * Psi^t computed only if local f variance is required 
  Real tr_S=0; 						//!< Trace of smoothing matrix, needed for the variance-covariance matrix (V) and eventually GCV computation
  Real sigma_hat_sq; 					//!< Estimator for the variance of the residuals (SSres/(n_obs-(q+tr_S)))
  bool is_sigma_hat_sq_computed = false;                //!< Boolean that tells whether sigma_hat_sq has been computed or not
  bool is_S_computed = false;				//!< Boolean that tells whether S has been computed or not
  MatrixXr V;						//!< Variance-Covariance matrix of the beta parameters
  bool is_V_computed = false;				//!< Boolean that tells whether V has been computed or not
  virtual void compute_S(void) = 0;			//!< Pure virtual method used to compute S, either in an exact or non-exact way 
  void compute_V(void);					//!< Method used to compute V
  MatrixXv V_f;                                         //!< Variance-Covariance matrix of f_hat estimator
  bool is_V_f_computed = false;                         //!< Boolean that tells whether V_f has been computed or not
  void compute_V_f(void);                               //!< Method used to compute V_f
  void compute_sigma_hat_sq(void);                      //!< Method to compute the estimator of the variance of the residuals 

  // methods that compute pvalues and/or CI on beta and on f respectively
  VectorXr compute_beta_pvalue(void) override;
  Real compute_f_pvalue(void) override;
  MatrixXv compute_beta_CI(void) override;
  MatrixXv compute_f_CI(void) override;
  
public:
  // CONSTUCTOR
  Wald_Base()=delete;	//The default constructor is deleted
  Wald_Base(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Inference_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; //Main constructor of the class
  
  virtual ~ Wald_Base(){};
  
  Real compute_GCV_from_inference(void) const override; //!< Needed to compute exact GCV in case Wald test is required and GCV exact is not provided by lambda optimization (Run after S computation)
  VectorXr compute_f_var(void) override; //!< Needed to compute local f variance if required
  
  // GETTERS
  inline const MatrixXr * getSp (void) const {return &this->S;}      //!< Getter of Sp \return Sp
  inline const MatrixXr * getVp (void) const {return &this->V;}      //!< Getter of Vp \ return Vp
};


// *** Wald_Exact Class ***
//! Hypothesis testing and confidence intervals using Wald implementation in an exact way 
/*!
   This template class derives from the Wald_Base class and it overrides the method that manages the computation of the smoothing matrix S, relying on an exact inversion of the MatrixNoCov. 
*/
template<typename InputHandler, typename MatrixType>
class Wald_Exact:public Wald_Base<InputHandler, MatrixType>{
private: 
  void compute_S(void) override;
public:
  // CONSTUCTOR
  Wald_Exact()=delete;	//The default constructor is deleted
  Wald_Exact(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Wald_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; 
};

// *** Wald_Non_Exact Class ***
//! Hypothesis testing and confidence intervals using Wald implementation in a non-exact way
/*!
   This template class derives from the Wald_Base class and it overrides the method that manages the computation of the smoothing matrix S, relying on an approximated inversion of the MatrixNoCov.
*/
template<typename InputHandler, typename MatrixType>
class Wald_Non_Exact:public Wald_Base<InputHandler, MatrixType>{
private: 
  void compute_S(void) override;
public:
  // CONSTUCTOR
  Wald_Non_Exact()=delete;	//The default constructor is deleted
  Wald_Non_Exact(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Wald_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; 
};


#include "Wald_imp.h"

#endif
