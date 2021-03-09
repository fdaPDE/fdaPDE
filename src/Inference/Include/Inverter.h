#ifndef __INFERENCE_INVERTER_H__
#define __INFERENCE_INVERTER_H__

// HEADERS
#include "../../FdaPDE.h"
#include "Inference_Carrier.h"

// *** inverse_Base Class ***
//! Class for the inversion of sparse matrices in inference framework
/*!
 This class provides the base for the classes that are used to compute the exact or non exact inverse of sparse matrices needed for inferential work
 The class is pure virtual, as one needs to specify the policy of inversion by using one of its derivation (inverse_Exact or inverse_Approximated)
*/
class Inverse_Base {
	protected:
		MatrixXr E_inv; //!< Matrix that stores the inverse when it has been computed
		bool inverse_computed = false;

	public:
		// Getter
		inline const MatrixXr * getInv(const Eigen::SparseLU<SpMat> * E_decp, const SpMat * Ep){if(inverse_computed==false){Compute_Inv(E_decp,Ep);} return &E_inv;}; //!< Getter for the Inverse Matrix [Need to call Compute_inv before calling this]
		virtual void Compute_Inv (const Eigen::SparseLU<SpMat> * E_decp, const SpMat * Ep) = 0; //!< Virtual function for the computation of the inverse matrix, takes as input matrixNoCov decomposition /param E_decp and matrixNoCov /param Ep

};

// *** inverse_Exact Class ***
//! Class for the exact inversion of sparse matrices in inference framework
/*!
 This class provides is used to compute the exact inverse of the spase noCovMatrix, needed to recover the exact inverse of the system matrix
*/
class Inverse_Exact : public Inverse_Base {
	public:
		void Compute_Inv(const Eigen::SparseLU<SpMat> * E_decp, const SpMat * Ep) override; //!< Function for the exact computation of the inverse matrix, takes as imput matrixNoCov decomposition /param E_decp and matrixNoCov /param Ep
};


#endif 
