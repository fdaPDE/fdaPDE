#ifndef __EIGEN_SIGN_FLIP_H__
#define __EIGEN_SIGN_FLIP_H__

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

// *** Eigen_Sign_Flip_Base Class ***
//! Hypothesis testing using eigen-sign-flip implementation
/*!
  This template class is an abstract base class to perform hypothesis testing using an eigen-sign-flip approach. Beyond all the objects and methods inherited from the abstract base inference class, it stores the matrix Lambda, whose type is given by the template parameter MatrixType which can be either a dense or a sparse matrix depending on the inversion exactness of the MatrixNoCov; it stores the partial residuals under the null hypothesis and a boolean indicating if the matrix Lambda has been computed. It overrides the methods that specify how to compute the p-values and the confidence intervals, according to the eigen-sign-flip apporach, but the latter is not actually implemented. It has a pure virtual method for the computation of Lambda, since it relies on the inversion of MatrixNoCov in an exact or non-exact way.
*/
template<typename InputHandler, typename MatrixType>
class Eigen_Sign_Flip_Base:public Inference_Base<InputHandler, MatrixType>{
protected:
  MatrixXr Partial_res_H0; 				//!< Contains: z - W^t * beta_0
  MatrixType Lambda;   					//!< I - Psi*(Psi^t * Psi + lambda*R)^-1*Psi^t
  bool is_Lambda_computed = false;			//!< Boolean that tells whether Lambda has been computed or not
  VectorXr Speckman_aux_ranges;                         //!< Speckman auxiliary CI ranges needed for CI method initialization
  bool is_speckman_aux_computed = false;                //!< Boolean that tells whether Speckman auxiliary ranges have been computed or not
  virtual void compute_Lambda(void) = 0;		//!< Pure virtual method used to compute Lambda, either in an exact or non-exact way
  void Compute_speckman_aux(void);                      //!< Auxiliary function for CI that computes the speckman ranges
  
  // methods that compute pvalues and/or CI on beta and on f respectively
  VectorXr compute_beta_pvalue(void) override;
  Real compute_f_pvalue(void) override;
  MatrixXv compute_beta_CI(void) override;
  MatrixXv compute_f_CI(void) override;
  
  Real compute_CI_aux_pvalue(const VectorXr &, const MatrixXr &, const  MatrixXr &) const;  //!< Computes the p value for a generic beta proposed in the research algorithm for CI 
  
public:
  // CONSTUCTOR
  Eigen_Sign_Flip_Base()=delete;	//The default constructor is deleted
  Eigen_Sign_Flip_Base(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Inference_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; //Main constructor of the class

  // DESTRUCTOR
  virtual ~ Eigen_Sign_Flip_Base() {};
    
  // GETTERS
  inline const MatrixXr * getLambdap (void) const {return &this->Lambda;}      	                        //!< Getter of Lambdap \return Lambdap
  inline const MatrixXr * getPartial_res_H0p (void) const {return &this->Partial_res_H0;}  		//!< Getter of Partial_res_H0p \return Partial_res_H0p
};

// *** Eigen_Sign_Flip_Exact Class ***
//! Hypothesis testing using Eigen-Sign-Flip implementation in an exact way 
/*!
  This template class derives from the Eigen_Sign_Flip_Base class and it overrides the method that manages the computation of the matrix Lambda, relying on an exact inversion of the MatrixNoCov. 
*/
template<typename InputHandler, typename MatrixType>
class Eigen_Sign_Flip_Exact:public Eigen_Sign_Flip_Base<InputHandler, MatrixType>{
private: 
  void compute_Lambda(void) override;
public:
  // CONSTUCTOR
  Eigen_Sign_Flip_Exact()=delete;	//The default constructor is deleted
  Eigen_Sign_Flip_Exact(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Eigen_Sign_Flip_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; 
};

// *** Eigen_Sign_Flip_Non_Exact Class ***
//! Hypothesis testing using Eigen-Sign-Flip implementation in a non-exact way 
/*!
  This template class derives from the Eigen_Sign_Flip_Base class and it overrides the method that manages the computation of the matrix Lambda, relying on an approximated inversion of the MatrixNoCov. 
*/
template<typename InputHandler, typename MatrixType>
class Eigen_Sign_Flip_Non_Exact:public Eigen_Sign_Flip_Base<InputHandler, MatrixType>{
private: 
  void compute_Lambda(void) override;
public:
  // CONSTUCTOR
  Eigen_Sign_Flip_Non_Exact()=delete;	//The default constructor is deleted
  Eigen_Sign_Flip_Non_Exact(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Eigen_Sign_Flip_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; 
};

//We define the naive operator to compare two vectors, namely v > u if and only if 
//each component of v (in absolute value) is strictly greater than 
//the corresponding component of u (in absolute value)
inline bool operator > (VectorXr v, VectorXr u){
  UInt q=v.size();
  if(u.size()!=q){
    Rprintf("Errore: dimensioni non combaciano");
    return false;
  }
  for (UInt i=0; i< q; i++){
    if(fabs(v(i))<=fabs(u(i))){
      return false;
    }
  }
  return true;
};
#include "Eigen_Sign_Flip_imp.h"

#endif
