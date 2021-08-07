#ifndef __INFERENCE_INVERTER_H__
#define __INFERENCE_INVERTER_H__

// HEADERS
#include "../../FdaPDE.h"
#include "Inference_Carrier.h"

// *** Inverse_Base Class ***
//! Class for the inversion of sparse matrices in inference framework
/*!
  \tparam MatrixType the type of the inverse of MatrixNoCov, it will be either MatrixXr or SpMat 
  This template class provides the base for the classes that are used to compute the exact or non exact inverse of sparse matrices needed for inferential work
  The class is pure virtual, as one needs to specify the policy of inversion by using one of its derivation (inverse_Exact or inverse_Approximated)
*/
template<typename MatrixType>
class Inverse_Base {
protected:
  MatrixType E_inv; 			//!< Matrix that stores the inverse when it has been computed
  bool inverse_computed = false;

public:
  // Constructor
  Inverse_Base() = default;		//!< Default constructor
		
  // Getters
  inline const MatrixType * getInv(void){if(inverse_computed==false){Compute_Inv();} return &E_inv;}; 	//!< Getter for the Inverse Matrix [Need to call Compute_inv before calling this]
  virtual inline bool get_status_inverse (void) const {return true;} 					//!< Virtual method that will be overriden by the non exact derived class to get the FSPAI status
  virtual void Compute_Inv (void) = 0; 						     //!< Pure virtual function for the computation of the inverse matrix, takes as input MatrixNoCov decomposition
  
  // Destructor
  virtual ~Inverse_Base() {};

};

// *** Inverse_Exact Class ***
//! Class for the exact inversion of sparse matrices in inference framework
/*!
  This class is used to compute the exact inverse of the sparse MatrixNoCov, needed to recover the exact inverse of the system matrix
*/
class Inverse_Exact : public Inverse_Base<MatrixXr> {
private:
  const SpMat * Ep;			        //!< Const pointer to the MatrixNoCov
  const Eigen::SparseLU<SpMat> * E_decp; 	//!< Const pointer to the (already computed) decomposition of MatrixNoCov
		
public:
  // Constructor
  Inverse_Exact()=delete; 										//!< Default constructor deleted
  Inverse_Exact(const SpMat * Ep_, const Eigen::SparseLU<SpMat> * E_decp_): Ep(Ep_),E_decp(E_decp_){}; 	//!< Main constructor

  void Compute_Inv(void) override;                                                                      //!< Function for the exact computation of the inverse matrix
}; 

// *** Inverse_Non_Exact Class ***
//! Class for the approximate inversion of sparse matrices in inference framework
/*!
  \tparam InputHandler the regression data type of the problem 
  This class is used to compute the sparse approximate inverse of the sparse matrix [Psi^T Psi + lambda*R] via FSPAI algorithm
*/
template<typename InputHandler>
class Inverse_Non_Exact : public Inverse_Base<SpMat> {
private:
  const Inference_Carrier<InputHandler> & inf_car; 	  //!< Reference to inference carrier
  SpMat E_tilde;					  //!< [Psi^T Psi + lambda*R] matrix
  SpMat R0_inv_tilde;                                     //!< Approximated inverse of the mass matrix in space as sparse matrix
  bool status_R0_inv_tilde = false;                       //!< Boolean that states whether the FSPAI computation of R0_inv went well
  bool status_E_tilde_inv = false;                        //!< Boolean that states whether the FSPAI computation of E_tilde_inv (E_inv) went well
     
  void pre_Inverse(void);                                 //!< Method that computes matrices R0_tilde and E_tilde
		
public:
  // Constructor
  Inverse_Non_Exact()=delete; 									//!< Default constructor deleted
  Inverse_Non_Exact(const Inference_Carrier<InputHandler> & inf_car_): inf_car(inf_car_){}; 	//!< Main constructor

  // Getter         
  inline bool get_status_inverse (void) const override {return status_E_tilde_inv;}; 		//!< Getter for the status of the inverse after FSPAI computation \return status_E_tilde_inv

  void Compute_Inv(void) override; 								//!< Function for the non-exact computation of the inverse matrix
};

#include "Inverter_imp.h"


#endif 

