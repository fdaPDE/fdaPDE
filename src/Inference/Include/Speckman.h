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

// *** Speckman_Base Class ***
//! Hypothesis testing and confidence intervals using Speckman implementation
/*!
This template class is an abstract base class to perform hypothesis testing and/or compute confidence intervals using a Speckman correction approach. Beyond all the objects and methods inherited from the abstract base inference class, it stores the matrix Lambda squared, whose type is given by the template parameter MatrixType which can be either a dense or a sparse matrix depending on the inversion exactness of the MatrixNoCov; it stores the variance-covariance matrix V of the beta parameters, the vector of estimated beta parameters via Speckaman estimator, the LU decomposition of W^T *Lambda^2 * W, alongside with some covenient boolean objects. It overrides the methods that specify how to compute the p-values and the confidence intervals, according to the Speckman apporach. The methods related to the nonparametric component are overriden, but they are not actually implemented. It has a pure virtual method for the computation of Lambda squared, since it relies on the inversion of MatrixNoCov in an exact or non-exact way. Moreover it has also a method for the computation of the estimators beta_hat required by the Speckman inferential approach. 
\tparam InputHandler the type of regression problem needed to determine the MixedFERegressionBase object type in Inference_Carrier<InputHandler>
\tparam MatrixType the type of matrix (MatrixXr or SpMat) used to store diffferent objects related to the smoother Lambda. SpMat type is related to approximated inference computation.
*/
template<typename InputHandler, typename MatrixType>
class Speckman_Base:public Inference_Base<InputHandler, MatrixType>{
protected: 
  MatrixType Lambda2;   				//!< (I - Psi*(Psi^t * Psi + lambda*R)^-1*Psi^t)^2
  bool is_Lambda2_computed = false;			//!< Boolean that tells whether Lambda^2 has been computed or not
  MatrixXr V;						//!< Variance-Covariance matrix of the beta parameters
  bool is_V_computed = false;				//!< Boolean that tells whether V has been computed or not
  VectorXr beta_hat;                                    //!< Vector of estimated beta parameters via Speckman estimator
  bool is_beta_hat_computed = false;                    //!< Boolean that tells whether beta_hat has been computed or not
  Eigen::PartialPivLU<MatrixXr> WLW_dec; 		//!< Decomposition of [W^t * Lambda^2 * W] 
  bool is_WLW_computed=false; 				//!< Boolean that tells whether WLW decomposition has been computed or not
  virtual void compute_Lambda2(void) = 0;		//!< Pure virtual method used to compute Lambda^2, either in an exact or non-exact way
  void compute_V(void);					//!< Method used to compute V
  void compute_WLW_dec(void); 				//!< Method that computes the decomposition for WLW
  void compute_beta_hat(void);               	        //!< Method used to compute beta estimates for the Speckman test
  
  // methods that compute pvalues and/or CI on beta and on f respectively
  VectorXr compute_beta_pvalue(void) override;
  Real compute_f_pvalue(void) override;
  MatrixXv compute_beta_CI(void) override;
  MatrixXv compute_f_CI(void) override;
  
public:
  // CONSTUCTOR
  Speckman_Base()=delete;	//The default constructor is deleted
  Speckman_Base(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Inference_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; //Main constructor of the class
  
  // DESTRUCTOR
  virtual ~ Speckman_Base() {};

  // GETTERS
  inline const MatrixXr * getLambda2p (void) const {return &this->Lambda2;}     //!< Getter of Lambda2p \return Lambda2p
  inline const MatrixXr * getVp (void) const {return &this->V;}     	 	//!< Getter of Vp \ return Vp
  inline const VectorXr * getBeta_hatp (void) const {return &this->beta_hat;}   //!< Getter of beta_hatp \ return beta_hatp
};

// *** Speckman_Exact Class ***
//! Hypothesis testing and confidence intervals using Speckman implementation in an exact way 
/*!
   This template class derives from the Speckman_Base class and it overrides the method that manages the computation of the matrix Lambda2, relying on an exact inversion of the MatrixNoCov. 
*/
template<typename InputHandler, typename MatrixType>
class Speckman_Exact:public Speckman_Base<InputHandler, MatrixType>{
private: 
  void compute_Lambda2(void) override;
public:
  // CONSTUCTOR
  Speckman_Exact()=delete;	//The default constructor is deleted
  Speckman_Exact(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Speckman_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; 

};

// *** Speckman_Non_Exact Class ***
//! Hypothesis testing and confidence intervals using Speckman implementation in a non-exact way 
/*!
   This template class derives from the Speckman_Base class and it overrides the method that manages the computation of the matrix Lambda2, relying on an approximated inversion of the MatrixNoCov. 
*/
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
