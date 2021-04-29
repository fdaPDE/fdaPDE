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
template<typename InputHandler>
class Eigen_Sign_Flip:public Inference_Base<InputHandler>{
private:
  MatrixXr Partial_res_H0; 				//!< Contains: z - W^t * beta_0
  MatrixXr B;						//!< Matrix Psi*(Psi^t * Psi + lambda*R)^-1*Psi^t 
  MatrixXr Lambda;   					//!< I - B
  bool is_Lambda_computed = false;			//!< Boolean that tells whether Lambda has been computed or not
  void compute_Lambda(void);				//!< Method used to compute Lambda
  VectorXr compute_pvalue(void) override;		//!< Method used to compute the pvalues of the tests 
  MatrixXv compute_CI(void) override;			//!< Method to compute the confidence intervals (not implemented yet)
  
public:
  // CONSTUCTOR
  Eigen_Sign_Flip()=delete;	//The default constructor is deleted
  Eigen_Sign_Flip(std::shared_ptr<Inverse_Base> inverter_, const Inference_Carrier<InputHandler> & inf_car_):Inference_Base<InputHandler>(inverter_, inf_car_){}; 
  Eigen_Sign_Flip(Eigen_Sign_Flip & rhs) = delete; //The default copy constructor is deleted
  inline Eigen_Sign_Flip(Eigen_Sign_Flip && rhs): Partial_res_H0(std::move(rhs.Partial_res_H0)), B(std::move(rhs.B)), Lambda(std::move(rhs.Lambda)), is_Lambda_computed(rhs.is_Lambda_computed){this->inverter=std::move(rhs.inverter); this->inf_car=rhs.inf_car;}; //Definition of the move constructor
  Eigen_Sign_Flip & operator=(Eigen_Sign_Flip && rhs) = delete; //The move assignment operator is deleted
  
  
  // GETTERS
  inline const MatrixXr * getLambdap (void) const {return &this->Lambda;}      	//!< Getter of Lambdap \return Lambdap
  inline const MatrixXr * getBp (void) const {return &this->B;}  		//!< Getter of Bp \return Bp
  
  void print_for_debug(void) const;
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
