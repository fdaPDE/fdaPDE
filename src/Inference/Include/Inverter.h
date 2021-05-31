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
		MatrixXr E_inv; 			//!< Matrix that stores the inverse when it has been computed
		bool inverse_computed = false;

	public:
		// Constructor
		Inverse_Base() = default;		//!< Default constructor
		
		// Getter
		inline const MatrixXr * getInv(void){if(inverse_computed==false){Compute_Inv();} return &E_inv;}; 	//!< Getter for the Inverse Matrix [Need to call Compute_inv before calling this]
		virtual void Compute_Inv (void) = 0; 						//!< Virtual function for the computation of the inverse matrix, takes as input matrixNoCov decomposition
		void print_for_debug(void) const;
		virtual ~Inverse_Base() {};

};

// *** inverse_Exact Class ***
//! Class for the exact inversion of sparse matrices in inference framework
/*!
 This class provided is used to compute the exact inverse of the sparse noCovMatrix, needed to recover the exact inverse of the system matrix
*/
class Inverse_Exact : public Inverse_Base {
	private:
		const SpMat * Ep;			//!< Const pointer to the MatrixNoCov
		const Eigen::SparseLU<SpMat> * E_decp; 	//!< Const pointer to the (already computed) decomposition of MatrixNoCov
		
	public:
		// Constructor
		Inverse_Exact()=delete; 										//!< Default constructor deleted
		Inverse_Exact(const SpMat * Ep_, const Eigen::SparseLU<SpMat> * E_decp_): Ep(Ep_),E_decp(E_decp_){}; 	//!< Main constructor

		void Compute_Inv(void) override;									//!< Function for the exact computation of the inverse matrix

// *** inverse_Non_Exact Class ***
//! Class for the approximate inversion of sparse matrices in inference framework
/*!
 This class provided computes the sparse approximate inverse of the sparse matrix [Psi^T Psi + lambda*R] via FSPAI algorithm
*/
template<typename InputHandler>
class Inverse_Non_Exact : public Inverse_Base {
	private:
		const Inference_Carrier<InputHandler> & inf_car; 	//!< Refernce to inference carrier
		//const SpMat E_tilde;					//!< [Psi^T Psi + lambda*R] matrix
		// ++ eventual parameters TO BE MODIFIED
		
	public:
		// Constructor
		Inverse_Non_Exact()=delete; 									//!< Default constructor deleted
		Inverse_Non_Exact(const Inference_Carrier<InputHandler> & inf_car_): inf_car(inf_car_){}; 	//!< Main constructor

		void Compute_Inv(void) override; 								//!< Function for the non-exact computation of the inverse matrix
};


#endif 

