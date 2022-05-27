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
  The class is pure virtual, as one needs to specify the policy of inversion by using one of its derivation (inverse_Exact or inverse_Non_Exact)
  In this CRAN version only Exact inversion is actually implemented
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
  inline const MatrixType * getInv(void){if(inverse_computed==false){Compute_Inv();} return &E_inv;}; 	//!< Getter for the Inverse Matrix [Needs to call Compute_Inv() before calling this]
  virtual inline bool get_status_inverse (void) const {return true;} 					//!< Virtual method that will be overriden by the Non_Exact derived class to get the FSPAI status
  virtual void Compute_Inv (void) = 0; 						     			//!< Pure virtual function for the computation of the inverse matrix
  
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

#endif 

