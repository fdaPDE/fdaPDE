#ifndef __INFERENCE_DATA_H___
#define __INFERENCE_DATA_H__

#include "../../FdaPDE.h"
#include <string>
//!  Class for inference data
/*!
 * This class collects all the data needed for inferential analysis over the linear estimated parameter of the model,
 * is constructed using the input defined in R.
*/
class InferenceData
{
	private:
		// Type of analysis required
  		std::string test_Type            	= "not-defined";			//!< Values: not-defined [default], one-at-the-time, simultaneous
  		std::string interval_Type		= "not-defined";			//!< Values: not-defined [default], one-at-the-time, simultaneous, bonferroni
  		std::string implementation_Type		= "wald"; 				//!< Values: wald [default], sandwich, permutational
  		bool exact_Inference			= false;				//!< If true the inverse of E matrix of woodbury decomposition of the system is computed exactly
		// Parameters needed
                MatrixXr coeff_Inference;	        					//!< If position j is true, the j-th covariate is taken into account in the inferential analysis
                VectorXr beta_0;             	              					//!< Values for the null hypostesis, if test_Type != not-defined
  		Real inference_Level			= 0; 					//!< Significance level for confidence intervals if interval_Type != not-defined
  		bool definition				= false;				//!< Defines whether the inference analysis needs to be carried out or not

	public:
  	//Constructors
		//! Default constructor
  		InferenceData() = default;

  		InferenceData(SEXP test_Type_, SEXP interval_Type_, SEXP implementation_Type_,
				SEXP exact_Inference_, SEXP coeff_Inference_, SEXP beta_0_, 
			        SEXP inference_Level_, SEXP definition_);
        //Setters
  		inline void set_test_type(const std::string && test_Type_){test_Type = test_Type_;};				//!< Setter for test_Type \param test_Type_ new test_Type
  		inline void set_interval_type(const std::string && interval_Type_){interval_Type = interval_Type_;};		//!< Setter for interval_Type \param interval_Type_ new interval_Type
  		inline void set_implementation_type(const std::string && implementation_Type_){implementation_Type = implementation_Type_;};	//!< Setter for implementation_Type \param implementation_Type_ new implementation_Type	
  		inline void set_exact_inference(const bool & exact_Inference_){exact_Inference = exact_Inference_;};		//!< Setter for exact_Inference \param exact_Inference_ new exact_Inference
  		inline void set_coeff_inference(const MatrixXr & coeff_inf){coeff_Inference = coeff_inf;};		        //!< Setter for coeff_Inference \param coeff_Inference_ new coeff_Inference
  		inline void set_beta_0(const VectorXr & beta_0_){beta_0 = beta_0_;};					        //!< Setter for beta0 \param beta0_ new beta0
  		inline void set_inference_level(const Real & inference_Level_){inference_Level = inference_Level_;};		//!< Setter for inference_level \param inference_level_ new inference_level
  		inline void set_definition(const bool & definition_){definition = definition_;};				//!< Setter for definition \param definition_ new definition

  //Getters
  inline std::string get_test_type() const{return this->test_Type;}; 				//!< Getter for test_Type \return test_Type
  inline std::string get_interval_type() const{return this->interval_Type;};			//!< Getter for interval_Type \return interval_Type
  inline std::string get_implementation_type() const{return this->implementation_Type;};	//!< Getter for implementation_Type \return implementation_Tupe
  inline bool get_exact_inference() const{return this->exact_Inference;};			//!< Getter for exact_Inference \return exact_Inference
  inline MatrixXr get_coeff_inference() const{return this->coeff_Inference;};		        //!< Getter for coeff_Inference \return coeff_Inference
  inline VectorXr get_beta_0() const{return this->beta_0;};				        //!< Getter for beta0 \return beta0
  inline Real get_inference_level() const{return this->inference_Level;};			//!< Getter for inference_Level \return inference_Level
  inline bool get_definition() const{return this->definition;};					//!< Getter for definition \return definition

  //For debugging
  void print_inference_data() const;
    
};



#endif /* __INFERENCE_DATA_H__ */
