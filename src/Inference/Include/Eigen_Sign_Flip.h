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

// *** Eigen_Sign_Flip Class ***
//! Hypothesis testing using implementation eigen sign-flip
/*!
  This class performes hypothesis testing using a eigen sign-flip approach. It contains a reference to an inverter, that manages to compute the invertion of matrixNoCov in an exact or non-exact way; It contains a reference to an Inference_Carrier object that wraps all the information needed to make inference. There is only one public method that calls the proper private methods to compute what is requested by the user.
*/
template<typename InputHandler, typename MatrixType>
class Eigen_Sign_Flip_Base:public Inference_Base<InputHandler, MatrixType>{
protected:
  MatrixXr Partial_res_H0; 				//!< Contains: z - W^t * beta_0
  MatrixXr Lambda;   					//!< I - Psi*(Psi^t * Psi + lambda*R)^-1*Psi^t
  bool is_Lambda_computed = false;			//!< Boolean that tells whether Lambda has been computed or not
  virtual void compute_Lambda(void) = 0;		//!< Method used to compute Lambda, either in an exact or non-exact way
  VectorXr compute_pvalue(void) override;		//!< Method used to compute the pvalues of the tests 
  MatrixXv compute_CI(void) override;			//!< Method to compute the confidence intervals (not implemented yet)
  
public:
  // CONSTUCTOR
  Eigen_Sign_Flip_Base()=delete;	//The default constructor is deleted
  Eigen_Sign_Flip_Base(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Inference_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; 

  // DESTRUCTOR
  virtual ~Eigen_Sign_Flip_Base() {};
    
  // GETTERS
  inline const MatrixXr * getLambdap (void) const {return &this->Lambda;}      	                        //!< Getter of Lambdap \return Lambdap
  inline const MatrixXr * getPartial_res_H0p (void) const {return &this->Partial_res_H0;}  		//!< Getter of Partial_res_H0p \return Partial_res_H0p
  
  void print_for_debug(void) const;
};

template<typename InputHandler, typename MatrixType>
class Eigen_Sign_Flip_Exact:public Eigen_Sign_Flip_Base<InputHandler, MatrixType>{
private: 
  void compute_Lambda(void) override;
public:
  // CONSTUCTOR
  Eigen_Sign_Flip_Exact()=delete;	//The default constructor is deleted
  Eigen_Sign_Flip_Exact(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Eigen_Sign_Flip_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; 
};

template<typename InputHandler, typename MatrixType>
class Eigen_Sign_Flip_Non_Exact:public Eigen_Sign_Flip_Base<InputHandler, MatrixType>{
private: 
  void compute_Lambda(void) override;
public:
  // CONSTUCTOR
  Eigen_Sign_Flip_Non_Exact()=delete;	//The default constructor is deleted
  Eigen_Sign_Flip_Non_Exact(std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_):Eigen_Sign_Flip_Base<InputHandler, MatrixType>(inverter_, inf_car_, pos_impl_){}; 
};

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
