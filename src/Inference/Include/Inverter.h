#ifndef __INFERENCE_INVERTER_H__
#define __INFERENCE_INVERTER_H__

// HEADERS
#include "../../FdaPDE.h"
#include "Inference_Carrier.h"
#include "Inference_Carrrier_imp.h"

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
		inline const MatrixXr * getInv(void) const {if(inverse_computed==false){Compute_Inv();} return &E_inv;}; //!< Getter for the Inverse Matrix [Need to call Compute_inv before calling this]
		virtual void Compute_Inv (const Inference_Carrier & inf_car) = 0; //!< Virtual function for the computation of the inverse matrix, takes as imput inference Carrier /param inf_car

}

// *** inverse_Exact Class ***
//! Class for the exact inversion of sparse matrices in inference framework
/*!
 This class provides is used to compute the exact inverse of the spase noCovMatrix, needed to recover the exact inverse of the system matrix
*/
class Inverese_Exact : public inverse_Base {
	public:
		void Compute_Inv (const Inference_Carrier & inf_car) override; //!< Function for the exact computation of the inverse matrix, takes as imput inference Carrier /param inf_car
}

#endif 
