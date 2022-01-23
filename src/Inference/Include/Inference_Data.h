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
  		std::vector<std::string> implementation_Type  	= {"wald"}; 	                //!< Values: wald [default], speckman, sign-flip, eigen-sign-flip
		std::vector<std::string> component_Type		= {"parametric"};		//!< Values: parametric [default], nonparametric, both
  		std::string exact_Inference			= "non-exact";		        //!< Values: non-exact [default], exact 
		std::string enhanced_Inference			= "classical";		        //!< Values: classical [default], enhanced 
		// Parameters needed
		MatrixXr locs_Inference;							//!< Matrix of spatial locations to be considered for nonparametric inference
		std::vector<UInt> locs_index_Inference;						//!< Vector of location indices to be considered for nonparametric inference
                bool locations_are_nodes                        = false;                        //!< Boolean defining whether the selected locations are a subset of mesh nodes
                MatrixXr coeff_Inference;	        					//!< Matrix of coefficients for the linear combinations of parameters 
                VectorXr beta_0;             	              					//!< Values for the null hypostesis, if test_Type != not-defined
		VectorXr f0_eval;								//!< Evaluations of the field f under the null hypothesis in the selected locations 
                bool f_Var                                      = false;             	        //!< Defines whether the local f variance has to be computed or not
  		VectorXr inference_Quantile;		 					//!< Quantiles needed for confidence intervals if interval_Type != not-defined
		VectorXr inference_Alpha                        = VectorXr::Constant(1,0.05); 	//!< Significance used in ESF confidence interval computation;
  		bool definition					= false;			//!< Defines whether the inference analysis needs to be carried out or not
                long int n_Flip 				= 1000; 			//!< Number of sign-flips if eigen-sign-flip tests are required
		Real tol_Fspai 					= 0.05; 			//!< Tolerance given in input to the FSPAI algorithm

	public:
  	//Constructors
		//! Default constructor
  		InferenceData() = default;
                //! Space constructor (with inference for f)
  		InferenceData(SEXP test_Type_, SEXP interval_Type_, SEXP implementation_Type_, SEXP component_Type_,
				SEXP exact_Inference_, SEXP enhanced_Inference_,SEXP locs_Inference_, SEXP locs_index_Inference_, SEXP locs_are_nodes, SEXP coeff_Inference_, SEXP beta_0_, SEXP f0_eval_, SEXP f_Var_,
			        SEXP inference_Quantile_, SEXP inference_Alpha_, SEXP n_Flip_, SEXP tol_Fspai_, SEXP definition_);
                //! Space-time constructor (without inference for f --> not implemented yet)
                InferenceData(SEXP test_Type_, SEXP interval_Type_, SEXP implementation_Type_, SEXP exact_Inference_, SEXP enhanced_Inference_, SEXP coeff_Inference_, SEXP beta_0_, SEXP f_Var_,
			        SEXP inference_Quantile_, SEXP inference_Alpha_, SEXP n_Flip_, SEXP tol_Fspai_, SEXP definition_);
                
                
        //Setters
  		inline void set_test_type(const std::vector<std::string> & test_Type_){test_Type = test_Type_;};		 //!< Setter for test_Type \param test_Type_ new test_Type
  		inline void set_interval_type(const std::vector<std::string> & interval_Type_){interval_Type = interval_Type_;}; //!< Setter for interval_Type \param interval_Type_ new interval_Type
  		inline void set_implementation_type(const std::vector<std::string> & implementation_Type_){implementation_Type = implementation_Type_;}; //!< Setter for implementation_Type \param implementation_Type_ new implementation_Type	
		inline void set_component_type(const std::vector<std::string> & component_Type_){component_Type = component_Type_;}; //!< Setter for component_Type \param component_Type_ new component_Type
  		inline void set_exact_inference(const std::string && exact_Inference_){exact_Inference = exact_Inference_;};	//!< Setter for exact_Inference \param exact_Inference_ new exact_Inference
		inline void set_enhanced_inference(const std::string && enhanced_Inference_){enhanced_Inference = enhanced_Inference_;}; //!< Setter for enhanced_Inference \param enhanced_Inference_ new enhanced_Inference
                inline void set_locs_inference(const MatrixXr & locs_inf){locs_Inference = locs_inf;};		        	//!< Setter for locs_Inference \param locs_inf new locs_Inference
		inline void set_locs_index_inference(const std::vector<UInt> & locs_ind_inf){locs_index_Inference = locs_ind_inf;}; //!< Setter for locs_index_Inference \param locs_ind_inf new locs_index_Inference
                inline void set_locs_are_nodes_inference(const bool & locs_are_nodes){locations_are_nodes = locs_are_nodes;};   //!< Setter for locations_are_nodes \param new locs_are_nodes
  		inline void set_coeff_inference(const MatrixXr & coeff_inf){coeff_Inference = coeff_inf;};		        //!< Setter for coeff_Inference \param coeff_inf new coeff_Inference
  		inline void set_beta_0(const VectorXr & beta_0_){beta_0 = beta_0_;};					        //!< Setter for beta0 \param beta0_ new beta0
		inline void set_f_0(const VectorXr & f_0_){f0_eval = f_0_;};					        	//!< Setter for f0_eval \param f_0_ new f0_eval
                inline void set_f_Var(const bool & f_Var_){f_Var = f_Var_;};				                        //!< Setter for f_Var \param f_Var_ new f_Var
  		inline void set_inference_quantile(const VectorXr & inference_Quantile_){inference_Quantile = inference_Quantile_;};//!< Setter for inference_Quantile \param inference_Quantile_ new inference_Quantile
		inline void set_inference_Alpha(const VectorXr & inference_Alpha_){inference_Alpha=inference_Alpha_;}; 		//!< Setter for inference_Alpha \param inference_Alpha_ new inference_Alpha
  		inline void set_definition(const bool & definition_){definition = definition_;};				//!< Setter for definition \param definition_ new definition
		inline void set_n_Flip(long int n_Flip_){n_Flip=n_Flip_;}; 							//!< Setter for n_Flip \param n_Flip_ new n_Flip
		inline void set_tol_Fspai(Real tol_Fspai_){tol_Fspai=tol_Fspai_;}; 						//!< Setter for tol_Fspai \param tol_Fspai_ new tol_Fspai

  	//Getters
  		inline std::vector<std::string> get_test_type() const{return this->test_Type;}; 				//!< Getter for test_Type \return test_Type
  		inline std::vector<std::string> get_interval_type() const{return this->interval_Type;};			        //!< Getter for interval_Type \return interval_Type
  		inline std::vector<std::string> get_implementation_type() const{return this->implementation_Type;};		//!< Getter for implementation_Type \return implementation_Type
		inline std::vector<std::string> get_component_type() const{return this->component_Type;};			//!< Getter for component_Type \return component_Type
 		inline std::string get_exact_inference() const{return this->exact_Inference;};			                //!< Getter for exact_Inference \return exact_Inference
		inline std::string get_enhanced_inference() const{return this->enhanced_Inference;};		                //!< Getter for enhanced_Inference \return enhanced_Inference
		inline MatrixXr get_locs_inference() const{return this->locs_Inference;};		                        //!< Getter for locs_Inference \return locs_Inference
		inline std::vector<UInt> get_locs_index_inference() const{return this->locs_index_Inference;};		        //!< Getter for locs_index_Inference \return locs_index_Inference
                inline bool get_locs_are_nodes_inference() const{return this->locations_are_nodes;};                            //!< Getter for locations_are_nodes \return locations_are_nodes
 		inline MatrixXr get_coeff_inference() const{return this->coeff_Inference;};		                        //!< Getter for coeff_Inference \return coeff_Inference
  		inline VectorXr get_beta_0() const{return this->beta_0;};				                        //!< Getter for beta0 \return beta0
		inline VectorXr get_f_0() const{return this->f0_eval;};				                        	//!< Getter for f0_eval \return f0_eval
                inline bool get_f_var() const{return this->f_Var;};					                        //!< Getter for f_var \return f_var
  		inline VectorXr get_inference_quantile() const{return this->inference_Quantile;};			        //!< Getter for inference_Quantile \return inference_Quantile
		inline VectorXr get_inference_alpha() const{return this->inference_Alpha;}; 			                //!< Getter for inference_Alpha \return inference_Alpha
  		inline bool get_definition() const{return this->definition;};					                //!< Getter for definition \return definition
		inline long int get_n_Flip() const{return this->n_Flip;}; 					                //!< Getter for n_Flip \return n_Flip
		inline Real get_tol_Fspai() const{return this->tol_Fspai;}; 					                //!< Getter for tol_Fspai \return tol_Fspai

	//For debugging
  		void print_inference_data() const;
    
};



#endif /* __INFERENCE_DATA_H__ */
