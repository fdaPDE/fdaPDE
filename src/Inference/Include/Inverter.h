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
		const SpMat * Ep;		//!< Sparse matrix that needs to be inverted
		MatrixXr E_inv; 		//!< Matrix that stores the inverse when it has been computed
		bool inverse_computed = false;

	public:
		// Constructors
		Inverse_Base()=delete;
		Inverse_Base(const SpMat * Ep_):Ep(Ep_){};

		// Getter
		inline const MatrixXr * getInv(void){if(inverse_computed==false){Compute_Inv();} return &E_inv;}; 	//!< Getter for the Inverse Matrix [Call Compute_inv before calling this]

		virtual void Compute_Inv (void) = 0; 									//!< Virtual function for the computation of the inverse matrix
		void print_for_debug(void) const;
		virtual ~Inverse_Base() {};

};

// *** inverse_Exact Class ***
//! Class for the exact inversion of sparse matrices in inference framework
/*!
 This class provides is used to compute the exact inverse of the spase noCovMatrix
*/
class Inverse_Exact : public Inverse_Base {
	private:
		const Eigen::SparseLU<SpMat> * E_decp;
	public:
		// Constructor
		Inverse_Exact()=delete;
		Inverse_Exact(const SpMat * Ep_, const Eigen::SparseLU<SpMat> * E_decp_):Inverse_Base(Ep_),E_decp(E_decp_){}; 

		void Compute_Inv(void) override; 										//!< Function for the exact computation of the inverse matrix
};

// *** inverse_Non_Exact Class ***
//! Class for the approximate inversion of sparse matrices in inference framework
/*!
 This class provides is used to compute approximate inverse of the spase noCovMatrix
*/
class Inverse_Non_Exact : public Inverse_Base {
	public:
		//Constructor
		Inverse_Non_Exact()=delete;
		Inverse_Non_Exact(const SpMat * Ep_):Inverse_Base(Ep_){};

		void Compute_Inv(void) override; //!< Function for the non-exact computation of the inverse matrix, using BICGSTAB algorithm
};


#endif 
