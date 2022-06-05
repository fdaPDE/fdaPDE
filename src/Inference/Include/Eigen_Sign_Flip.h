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
//! Hypothesis testing and confidence intervals using eigen-sign-flip implementation
/*!
  This template class is an abstract base class to perform hypothesis testing and computing confidence intervals using an eigen-sign-flip approach. Beyond all the objects and methods inherited from the abstract base inference class, it stores the matrix Lambda, whose type is given by the template parameter MatrixType which can be either a dense or a sparse matrix depending on the inversion exactness of the MatrixNoCov; it stores the partial residuals under the null hypothesis and a boolean indicating if the matrix Lambda has been computed. It overrides the methods that specify how to compute the p-values and the confidence intervals, according to the eigen-sign-flip apporach. It has a pure virtual method for the computation of Lambda, since it relies on the inversion of MatrixNoCov in an exact or non-exact way (not implemented in this version).
  \tparam InputHandler the type of regression problem needed to determine the MixedFERegressionBase object type in Inference_Carrier<InputHandler>
  \tparam MatrixType the type of matrix (MatrixXr or SpMat) used to store diffferent objects related to the smoothers S and Lambda. SpMat type is related to approximated inference computation (not implemented in this version).
*/
template<typename InputHandler, typename MatrixType>
class Eigen_Sign_Flip_Base:public Inference_Base<InputHandler, MatrixType>{
protected:
  MatrixXr Partial_res_H0; 				//!< Contains: z - W^t * beta_0
  VectorXr Partial_f_res_H0;                            //!< Contains: Q_loc*(z_loc - f_0)
  MatrixType Lambda;   					//!< I - Psi*(Psi^t * Psi + lambda*R)^-1*Psi^t
  bool is_Lambda_computed = false;			//!< Boolean that tells whether Lambda has been computed or not
  VectorXr Speckman_aux_ranges;                         //!< Speckman auxiliary CI ranges needed for CI method initialization (for beta)
  VectorXr Wald_aux_ranges;                             //!< Wald auxiliary CI ranges needed for CI method initialization (for f) 
  bool is_speckman_aux_computed = false;                //!< Boolean that tells whether Speckman auxiliary ranges have been computed or not
  bool is_wald_aux_computed = false;                    //!< Boolean that tells whether Wald auxiliary ranges have been computed or not
  virtual void compute_Lambda(void) = 0;		//!< Pure virtual method used to compute Lambda, either in an exact or non-exact way
  void Compute_speckman_aux(void);                      //!< Auxiliary function for beta CI that computes the speckman ranges
  void Compute_wald_aux(void);                          //!< Auxiliary function for f CI that computes the wald ranges
  
  // methods that compute pvalues and/or CI on beta and on f respectively
  VectorXr compute_beta_pvalue(void) override;
  Real compute_f_pvalue(void) override;
  MatrixXv compute_beta_CI(void) override;
  MatrixXv compute_f_CI(void) override;
  
  Real compute_CI_aux_beta_pvalue(const VectorXr &, const MatrixXr &, const VectorXr &, const  MatrixXr &) const;  //!< Computes the unilateral p value for a generic beta proposed in the research algorithm for CI 
  Real compute_CI_aux_f_pvalue(const VectorXr &, const UInt current_index) const;  //!< Computes the unilateral p value for a generic value of f in a given point proposed in the research algorithm for CI  
    
  
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


//We define the naive operator to compare unilaterally two vectors, namely v > u if and only if 
//each component of v is strictly greater than 
//the corresponding component of u
inline bool is_Unilaterally_Greater (VectorXr v, VectorXr u){
  UInt q=v.size();
  if(u.size()!=q){
    Rprintf("Error: in Eigen-Sign-Flip procedure two vectors of different length have been compared");
    return false;
  }
  for (UInt i=0; i< q; i++){
    if(v(i)<=u(i)){
      return false;
    }
  }
  return true;
};

//We define the naive operator to compare unilaterally two vectors, namely v < u if and only if 
//each component of v is strictly greater than 
//the corresponding component of u
inline bool is_Unilaterally_Smaller (VectorXr v, VectorXr u){
  UInt q=v.size();
  if(u.size()!=q){
    Rprintf("Error: in Eigen-Sign-Flip procedure two vectors of different length have been compared");
    return false;
  }
  for (UInt i=0; i< q; i++){
    if(v(i)>=u(i)){
      return false;
    }
  }
  return true;
};

// minimum between VectorXr (needed in one-at-the-time tests)
inline VectorXr min(const VectorXr & v, const VectorXr & u){
  VectorXr result;
  result.resize(v.size());
  for(UInt i=0; i<v.size(); ++i){
    result(i)=std::min(v(i),u(i));
  }
  return result;
}

#include "Eigen_Sign_Flip_imp.h"
#endif
