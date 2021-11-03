#ifndef __INFERENCE_DATA_H__
#define __INFERENCE_DATA_H__

#include "../../FdaPDE.h"
#include <string>
#include <vector>
//!  Class for inference data
/*!
 * This class collects all the data needed for inferential analysis over the linear estimated parameter of the model,
 * it is constructed using the input defined in R.
*/
class InferenceData
{
	private:
		// Type of analysis required
  		std::vector<std::string> test_Type		= {"not-defined"};	        //!< Values: not-defined [default], one-at-the-time, simultaneous
  		std::vector<std::string> interval_Type  	= {"not-defined"};	        //!< Values: not-defined [default], one-at-the-time, simultaneous, bonferroni
  		std::vector<std::string> implementation_Type  	= {"wald"}; 	                //!< Values: wald [default], speckman, eigen-sign-flip
  		std::string exact_Inference			= "non-exact";		        //!< Values: non-exact [default], exact 
		// Parameters needed
                MatrixXr coeff_Inference;	        					//!< Matrix of coefficients for the linear combinations of parameters 
                VectorXr beta_0;             	              					//!< Values for the null hypostesis, if test_Type != not-defined
                bool f_Var                                      = false;             	        //!< Defines whether the local f variance has to be computed or not
  		VectorXr inference_Quantile;		 					//!< Quantiles needed for confidence intervals if interval_Type != not-defined
		Real inference_Alpha                            = 0.05; 			//!< Significance used in ESF confidence interval computation;
  		bool definition					= false;			//!< Defines whether the inference analysis needs to be carried out or not
                long int n_Flip 				= 1000; 			//!< Number of sign-flips if eigen-sign-flip tests are required
		Real tol_Fspai 					= 0.05; 			//!< Tolerance given in input to the FSPAI algorithm

	public:
  	//Constructors
		//! Default constructor
  		InferenceData() = default;

  		InferenceData(SEXP test_Type_, SEXP interval_Type_, SEXP implementation_Type_,
				SEXP exact_Inference_, SEXP coeff_Inference_, SEXP beta_0_, SEXP f_Var_,
			        SEXP inference_Quantile_, SEXP inference_Alpha_, SEXP n_Flip_, SEXP tol_Fspai_, SEXP definition_);
        //Setters
  		inline void set_test_type(const std::vector<std::string> & test_Type_){test_Type = test_Type_;};		 //!< Setter for test_Type \param test_Type_ new test_Type
  		inline void set_interval_type(const std::vector<std::string> & interval_Type_){interval_Type = interval_Type_;}; //!< Setter for interval_Type \param interval_Type_ new interval_Type
  		inline void set_implementation_type(const std::vector<std::string> & implementation_Type_){implementation_Type = implementation_Type_;}; //!< Setter for implementation_Type \param implementation_Type_ new implementation_Type	
  		inline void set_exact_inference(const std::string && exact_Inference_){exact_Inference = exact_Inference_;};	//!< Setter for exact_Inference \param exact_Inference_ new exact_Inference
  		inline void set_coeff_inference(const MatrixXr & coeff_inf){coeff_Inference = coeff_inf;};		        //!< Setter for coeff_Inference \param coeff_Inference_ new coeff_Inference
  		inline void set_beta_0(const VectorXr & beta_0_){beta_0 = beta_0_;};					        //!< Setter for beta0 \param beta0_ new beta0
                inline void set_f_Var(const bool & f_Var_){f_Var = f_Var_;};				                        //!< Setter for f_Var \param f_Var_ new f_Var
  		inline void set_inference_quantile(const VectorXr & inference_Quantile_){inference_Quantile = inference_Quantile_;};//!< Setter for inference_Quantile \param inference_Quantile_ new inference_Quantile
		inline void set_inference_Alpha(Real inference_Alpha_){inference_Alpha=inference_Alpha_;}; 			//!< Setter for inference_Alpha \param inference_Alpha_ new inference_Alpha
  		inline void set_definition(const bool & definition_){definition = definition_;};				//!< Setter for definition \param definition_ new definition
		inline void set_n_Flip(long int n_Flip_){n_Flip=n_Flip_;}; 							//!< Setter for n_Flip \param n_Flip_ new n_Flip
		inline void set_tol_Fspai(Real tol_Fspai_){tol_Fspai=tol_Fspai_;}; 						//!< Setter for tol_Fspai \param tol_Fspai_ new tol_Fspai

  	//Getters
  		inline std::vector<std::string> get_test_type() const{return this->test_Type;}; 				//!< Getter for test_Type \return test_Type
  		inline std::vector<std::string> get_interval_type() const{return this->interval_Type;};			        //!< Getter for interval_Type \return interval_Type
  		inline std::vector<std::string> get_implementation_type() const{return this->implementation_Type;};		//!< Getter for implementation_Type \return implementation_Tupe
 		inline std::string get_exact_inference() const{return this->exact_Inference;};			                //!< Getter for exact_Inference \return exact_Inference
 		inline MatrixXr get_coeff_inference() const{return this->coeff_Inference;};		                        //!< Getter for coeff_Inference \return coeff_Inference
  		inline VectorXr get_beta_0() const{return this->beta_0;};				                        //!< Getter for beta0 \return beta0
                inline bool get_f_var() const{return this->f_Var;};					                        //!< Getter for f_var \return f_var
  		inline VectorXr get_inference_quantile() const{return this->inference_Quantile;};			        //!< Getter for inference_Quantile \return inference_Quantile
		inline Real get_inference_Alpha() const{return this->inference_Alpha;}; 			                //!< Getter for inference_Alpha \return inference_Alpha
  		inline bool get_definition() const{return this->definition;};					                //!< Getter for definition \return definition
		inline long int get_n_Flip() const{return this->n_Flip;}; 					                //!< Getter for n_Flip \return n_Flip
		inline Real get_tol_Fspai() const{return this->tol_Fspai;}; 					                //!< Getter for tol_Fspai \return tol_Fspai

	//For debugging
  		void print_inference_data() const;
    
};



#endif /* __INFERENCE_DATA_H__ */
