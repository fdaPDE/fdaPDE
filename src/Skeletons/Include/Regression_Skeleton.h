#ifndef __REGRESSION_SKELETON_H__ 
#define __REGRESSION_SKELETON_H__

#include "../../FdaPDE.h"
#include "../../Lambda_Optimization/Include/Carrier.h"
#include "../../Lambda_Optimization/Include/Grid_Evaluator.h"
#include "../../Lambda_Optimization/Include/Lambda_Optimizer.h"
#include "../../Lambda_Optimization/Include/Newton.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"
#include "../../Lambda_Optimization/Include/Optimization_Methods_Factory.h"
#include "../../Lambda_Optimization/Include/Solution_Builders.h"
#include "../../Inference/Include/Inference_Data.h"
#include "../../Inference/Include/Inference_Carrier.h"
#include "../../Inference/Include/Inverter.h"
#include "../../Inference/Include/Wald.h"
#include "../../Inference/Include/Speckman.h"
#include "../../Inference/Include/Eigen_Sign_Flip.h"
#include "../../Inference/Include/Inference_Factory.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../Regression/Include/Mixed_FE_Regression.h"
#include "../../Global_Utilities/Include/FSPAI_Wrapper.h"
#include <memory>

template<typename CarrierType>
std::pair<MatrixXr, output_Data> optimizer_method_selection(CarrierType & carrier);
template<typename EvaluationType, typename CarrierType>
std::pair<MatrixXr, output_Data> optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier);
template<typename InputHandler>
void inference_wrapper(const OptimizationData & opt_data, output_Data & output, const Inference_Carrier<InputHandler> & inf_car, MatrixXv & inference_output);
template<typename InputHandler>
void lambda_inference_selection(const OptimizationData & optimizationData, const output_Data & output, const InferenceData & inferenceData, MixedFERegression<InputHandler> & regression, Real & lambda_inference);

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
SEXP regression_skeleton(InputHandler & regressionData, OptimizationData & optimizationData, InferenceData & inferenceData, SEXP Rmesh)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, regressionData.getSearch());	// Create the mesh
	MixedFERegression<InputHandler> regression(regressionData, optimizationData, mesh.num_nodes()); // Define the mixed object

	regression.preapply(mesh); // preliminary apply (preapply) to store all problem matrices

	std::pair<MatrixXr, output_Data> solution_bricks;	// Prepare solution to be filled
	MatrixXv inference_Output; 				// Matrix that will store the output from inference, i.e. p-values and/or intervals
        Real lambda_inference = 0; 				// Will store the value of the optimal lambda

	// Build the Carrier according to problem type
	if(regression.isSV())
	{
		if(regressionData.getNumberOfRegions()>0)
		{
			//Rprintf("Areal-forced\n");
			Carrier<InputHandler,Forced,Areal>
				carrier = CarrierBuilder<InputHandler>::build_forced_areal_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler, Forced,Areal>>(carrier);
                        
              		lambda_inference_selection(optimizationData, solution_bricks.second, inferenceData, regression, lambda_inference); // Set lambda for inference
		}
		else
		{
			//Rprintf("Pointwise-forced\n");
			Carrier<InputHandler,Forced>
				carrier = CarrierBuilder<InputHandler>::build_forced_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler,Forced>>(carrier);

                        lambda_inference_selection(optimizationData, solution_bricks.second, inferenceData, regression, lambda_inference); // Set lambda for inference
		}
	}
	else
	{
		if(regressionData.getNumberOfRegions()>0)
		{
			//Rprintf("Areal\n");
			Carrier<InputHandler,Areal>
				carrier = CarrierBuilder<InputHandler>::build_areal_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler,Areal>>(carrier);

                        lambda_inference_selection(optimizationData, solution_bricks.second, inferenceData, regression, lambda_inference); // Set lambda for inference
		}
		else
		{
			//Rprintf("Pointwise\n");
			Carrier<InputHandler>
				carrier = CarrierBuilder<InputHandler>::build_plain_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler>>(carrier);

                        lambda_inference_selection(optimizationData, solution_bricks.second, inferenceData, regression, lambda_inference); // Set lambda for inference	
		}

	}
        //Inference
	if(inferenceData.get_definition()==true){ 
		//only if inference is actually required
		Inference_Carrier<InputHandler> inf_car(&regressionData, &regression, &solution_bricks.second, &inferenceData, lambda_inference); //Carrier for inference Data
		inference_wrapper(optimizationData, solution_bricks.second, inf_car, inference_Output);    
       }
 	return Solution_Builders::build_solution_plain_regression<InputHandler, ORDER, mydim, ndim>(solution_bricks.first,solution_bricks.second,mesh,regressionData,inference_Output,inferenceData);
}

//! Function to select the right optimization method
/*
 \tparam CarrierType the type of Carrier to be employed
 \param carrier the Carrier used for the methods
 \return the solution to pass to the Solution_Builders
*/
template<typename CarrierType>
std::pair<MatrixXr, output_Data> optimizer_method_selection(CarrierType & carrier)
{
	// Build the optimizer
	const OptimizationData * optr = carrier.get_opt_data();
	if(optr->get_loss_function() == "GCV" && optr->get_DOF_evaluation() == "exact")
	{
		//Rprintf("GCV exact\n");
		GCV_Exact<CarrierType, 1> optim(carrier);
		return optimizer_strategy_selection<GCV_Exact<CarrierType, 1>, CarrierType>(optim, carrier);
	}
	else if(optr->get_loss_function() == "GCV" && (optr->get_DOF_evaluation() == "stochastic" || optr->get_DOF_evaluation() == "not_required"))
	{
		//Rprintf("GCV stochastic\n");
		GCV_Stochastic<CarrierType, 1> optim(carrier, true);
		return optimizer_strategy_selection<GCV_Stochastic<CarrierType, 1>, CarrierType>(optim, carrier);
	}
	else // if(optr->get_loss_function() == "unused" && optr->get_DOF_evaluation() == "not_required")
	{
		//Rprintf("Pure evaluation\n");
		GCV_Stochastic<CarrierType, 1> optim(carrier, false);

		timer Time_partial; // Of the sole optimization
		Time_partial.start();
		// Rprintf("WARNING: start taking time\n");

		// Get the solution
		output_Data output;
		output.z_hat.resize(carrier.get_psip()->rows(),carrier.get_opt_data()->get_size_S());
		output.lambda_vec = carrier.get_opt_data()->get_lambda_S();
		MatrixXr solution;
 		MatrixXv betas;
		betas.resize(carrier.get_opt_data()->get_size_S(),1);

		for(UInt j=0; j<carrier.get_opt_data()->get_size_S(); j++)
		{
			if(j==0)
			{
				MatrixXr sol = carrier.apply(carrier.get_opt_data()->get_lambda_S()[j]);
				solution.resize(sol.rows(),carrier.get_opt_data()->get_size_S());
				solution.col(j) = sol;
			}
			else
			{
				solution.col(j) = carrier.apply(carrier.get_opt_data()->get_lambda_S()[j]);
			}
			optim.combine_output_prediction(solution.topRows(solution.rows()/2).col(j),output,j);
			if(carrier.get_model()->getBeta().cols()>0 & carrier.get_model()->getBeta().rows()>0)
				betas.coeffRef(j,0)=carrier.get_model()->getBeta().coeffRef(0,0);
		}

		// Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		output.time_partial = T.tv_sec + 1e-9*T.tv_nsec;

                // postponed after apply in order to have betas computed
        output.betas = betas;

        return {solution, output};
	}
}

//! Function to apply the optimization strategy, grid or Newton
/*
\\tparam EvaluationType optimization type to be used
 \tparam CarrierType the type of Carrier to be employed
 \param optim EvaluationType containing the class related to the function to be optimized,together with the method (exact or stochastic)
 \param carrier the Carrier used for the methods
 \return the solution to pass to the Solution_Builders
*/
template<typename EvaluationType, typename CarrierType>
std::pair<MatrixXr, output_Data> optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier)
{
	// Build wrapper and newton method
	Function_Wrapper<Real, Real, Real, Real, EvaluationType> Fun(optim);
	typedef Function_Wrapper<Real, Real, Real, Real, EvaluationType> FunWr;

	const OptimizationData * optr = carrier.get_opt_data();
	if(optr->get_criterion() == "grid")
	{
		timer Time_partial; // Of the sole optimization
		Time_partial.start();
		// Rprintf("WARNING: start taking time\n");

		// this will be used when grid will be correctly implemented, also for return elements

		Eval_GCV<EvaluationType> eval(Fun, optr->get_lambda_S());
		output_Data output = eval.Get_optimization_vectorial();

		// Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		// Get the solution
		MatrixXr solution = carrier.apply(output.lambda_sol);

		output.time_partial = T.tv_sec + 1e-9*T.tv_nsec;

                //postponed after apply in order to have betas computed
                output.betas = carrier.get_model()->getBeta();

                return {solution, output};

	}
	else // 'not_required' optimization can't enter here!! [checked in R code]
	{
		std::unique_ptr<Opt_methods<Real,Real,EvaluationType>> optim_p =
			Opt_method_factory<Real, Real, EvaluationType>::create_Opt_method(optr->get_criterion(), Fun);

                // Compute optimal lambda
		Checker ch;
		std::vector<Real> lambda_v_;
		std::vector<Real> GCV_v_;
		Real lambda = optr->get_initial_lambda_S();

		if(lambda<=0)
		{
			lambda = -1.0;
		}

		timer Time_partial; // Of the sole optimization
		Time_partial.start();
		// Rprintf("WARNING: start taking time\n");

		std::pair<Real, UInt> lambda_couple = optim_p->compute(lambda, optr->get_stopping_criterion_tol(), 40, ch, GCV_v_, lambda_v_);

		//Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		// Get the solution
		// to compute f and g hat
		MatrixXr solution = carrier.apply(lambda_couple.first);

		// postponed after apply in order to have betas computed
		// now the last values in GCV_exact are the correct ones, related to the last iteration
		output_Data  output = Fun.get_output(lambda_couple, T, GCV_v_, lambda_v_, ch.which());
		// the copy is necessary for the bulders outside

		return {solution, output};
	}
}

//! Function to select the right inference method
/*
  \tparam InputHandler the type of regression problem
  \param opt_data the object containing optimization data
  \param output the object containing the solution of the optimization problem 
  \param inf_car the object wrapping all the objects needed to make inference
  \param inference_output the object to be filled with inference output 
  \return void
*/
template<typename InputHandler>
void inference_wrapper(const OptimizationData & opt_data, output_Data & output, const Inference_Carrier<InputHandler> & inf_car, MatrixXv & inference_output)
{
  UInt n_implementations = inf_car.getInfData()->get_implementation_type().size();
  UInt p = inf_car.getInfData()->get_coeff_inference().rows();

  inference_output.resize(n_implementations, p+1);

  if(inf_car.getInfData()->get_exact_inference() == "exact"){
    // Select the right policy for inversion of MatrixNoCov
    std::shared_ptr<Inverse_Base<MatrixXr>> inference_Inverter = std::make_shared<Inverse_Exact>(inf_car.getEp(), inf_car.getE_decp());

    for(UInt i=0; i<n_implementations; ++i){
      // Factory instantiation for solver: using factory provided in Inference_Factory.h
      std::shared_ptr<Inference_Base<InputHandler,MatrixXr>> inference_Solver = Inference_Factory<InputHandler,MatrixXr>::create_inference_method(inf_car.getInfData()->get_implementation_type()[i], inference_Inverter, inf_car, i); // Selects the right implementation and solves the inferential problems		
		
      inference_output.row(i) = inference_Solver->compute_inference_output();

      if(inf_car.getInfData()->get_implementation_type()[i]=="wald" && opt_data.get_loss_function()=="unused" && opt_data.get_size_S()==1){
	output.GCV_opt=inference_Solver->compute_GCV_from_inference(); // Computing GCV if Wald has being called is an almost zero-cost function, since tr(S) hase been already computed
      }
    }
  }
  else{
    // Select the right policy for inversion of MatrixNoCov
    std::shared_ptr<Inverse_Base<SpMat>> inference_Inverter = std::make_shared<Inverse_Non_Exact<InputHandler>>(inf_car);

    for(UInt i=0; i<n_implementations; ++i){
      // Factory instantiation for solver: using factory provided in Inference_Factory.h
      std::shared_ptr<Inference_Base<InputHandler,SpMat>> inference_Solver = Inference_Factory<InputHandler,SpMat>::create_inference_method(inf_car.getInfData()->get_implementation_type()[i], inference_Inverter, inf_car, i); // Selects the right implementation and solves the inferential problems		
		
      inference_output.row(i) = inference_Solver->compute_inference_output();


    }
  }

  return;
  
}

//! Function that sets the correct lambda needed for inferential operations
/*
  \tparam InputHandler the type of regression problem
  \param optimization_data the object containing optimization data
  \param output_Data the object containing the solution of the optimization problem
  \param inferenceData the object containing the data needed for for inference
  \param regression the object containing the model of the problem
  \param lambda_inference the lambda that will be used to compute the optimal model and the right inferential solutions
  \return void
*/
template<typename InputHandler>
void lambda_inference_selection (const OptimizationData & optimizationData, const output_Data & output, const InferenceData & inferenceData, MixedFERegression<InputHandler> & regression, Real & lambda_inference){
	if(inferenceData.get_definition()==true && optimizationData.get_loss_function()!="unused"){
		lambda_inference = output.lambda_sol;
		if(optimizationData.get_last_lS_used() != lambda_inference){
			regression.build_regression_inference(lambda_inference);
			// for debug only 
			Rprintf("I'm computing again the matrices in Mixed_FERegression\n");
			}
	}else{ 		// supposing we have only one lambda, otherwise inference gets discarded in smoothing.R
		if(inferenceData.get_definition()==true){
			lambda_inference = optimizationData.get_last_lS_used();
			}
		}
	return; 
}

#endif
