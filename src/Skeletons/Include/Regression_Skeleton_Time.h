#ifndef __REGRESSION_SKELETON_TIME_H__
#define __REGRESSION_SKELETON_TIME_H__

#include "../../FdaPDE.h"
#include "../../Lambda_Optimization/Include/Carrier.h"
#include "../../Lambda_Optimization/Include/Grid_Evaluator.h"
#include "../../Lambda_Optimization/Include/Lambda_Optimizer.h"
#include "../../Lambda_Optimization/Include/Newton.h"
#include "../../Lambda_Optimization/Include/Optimization_Methods_Factory.h"
#include "../../Lambda_Optimization/Include/Solution_Builders.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../Regression/Include/Mixed_FE_Regression.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"
#include "../../Global_Utilities/Include/Lambda.h"
#include "../../Inference/Include/Inference_Data.h"
#include "../../Inference/Include/Inference_Carrier.h"
#include "../../Inference/Include/Inverter.h"
#include "../../Inference/Include/Wald.h"
#include "../../Inference/Include/Speckman.h"
#include "../../Inference/Include/Eigen_Sign_Flip.h"
#include "../../Inference/Include/Inference_Factory.h"
#include <memory>

template<typename CarrierType>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Temporal, CarrierType>::value>, t_type>::value,
			std::pair<MatrixXr, output_Data<2>> >::type optimizer_method_selection(CarrierType & carrier);

template<typename EvaluationType, typename CarrierType, UInt size>
typename std::enable_if<size==2, std::pair<MatrixXr, output_Data<2>>>::type
optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier);
template<class GCV_type, typename CarrierType>
std::pair<MatrixXr, output_Data<2>> parabolic_routine(CarrierType & carrier);
template<typename InputHandler>
void lambda_inference_selection(const OptimizationData & optimizationData, const  output_Data<2> & output, const InferenceData & inferenceData, MixedFERegression<InputHandler> & regression, Real & lambda_inference_S, Real & lambda_inference_T);
template<typename InputHandler>
void inference_wrapper_time(const OptimizationData & opt_data, const output_Data<2> & output, const Inference_Carrier<InputHandler> & inf_car, MatrixXv & inference_output);
template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void compute_nonparametric_inference_matrices_time(const MeshHandler<ORDER, mydim, ndim>  & mesh_, const std::vector<Real> meshTime_, const InputHandler & regressionData_, InferenceData & inferenceData_, Inference_Carrier<InputHandler> & inf_car_)

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
SEXP regression_skeleton_time(InputHandler & regressionData, OptimizationData & optimizationData, InferenceData & inferenceData, SEXP Rmesh, SEXP Rmesh_time)
{
  MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, regressionData.getSearch());//! load the mesh
  UInt n_time = Rf_length(Rmesh_time);
  std::vector<Real> mesh_time(n_time);
  for(UInt i=0; i<n_time; ++i)
    {
      mesh_time[i] = REAL(Rmesh_time)[i];
    }
	
  MixedFERegression<InputHandler> regression(mesh_time, regressionData, optimizationData, mesh.num_nodes());//! load data in a C++ object
	
  regression.preapply(mesh); // preliminary apply (preapply) to store all problem matrices

  std::pair<MatrixXr, output_Data<2>> solution_bricks;	// Prepare solution to be filled
  MatrixXv inference_Output; 				// Matrix that will store the output from inference, i.e. p-values and/or intervals
  Real lambda_inference_S = 0; 				// Will store the value of the optimal lambda_S
  Real lambda_inference_T = 0; 				// Will store the value of the optimal lambda_T
	
  // Build the Carrier according to problem type
  if(regression.isSV())
    {
      if(regressionData.getNumberOfRegions()>0)
	{
	  //Rprintf("Temporal, Areal-forced\n");
	  Carrier<InputHandler,Temporal,Forced,Areal>
	    carrier = CarrierBuilder<InputHandler>::build_temporal_forced_areal_carrier(regressionData, regression, optimizationData);
	  solution_bricks = optimizer_method_selection<Carrier<InputHandler,Temporal,Forced,Areal>>(carrier);

	  lambda_inference_selection(optimizationData, solution_bricks.second, inferenceData, regression, lambda_inference_S, lambda_inference_T); // Set lambdas for inference
	}
      else
	{
	  //Rprintf("Temporal, Pointwise-forced\n");
	  Carrier<InputHandler,Temporal,Forced>
	    carrier = CarrierBuilder<InputHandler>::build_temporal_forced_carrier(regressionData, regression, optimizationData);
	  solution_bricks = optimizer_method_selection<Carrier<InputHandler,Temporal,Forced>>(carrier);

	  lambda_inference_selection(optimizationData, solution_bricks.second, inferenceData, regression, lambda_inference_S, lambda_inference_T); // Set lambdas for inference
	}
    }
  else
    {
      if(regressionData.getNumberOfRegions()>0)
	{
	  //Rprintf("Temporal, Areal\n");
	  Carrier<InputHandler,Temporal,Areal>
	    carrier = CarrierBuilder<InputHandler>::build_temporal_areal_carrier(regressionData, regression, optimizationData);
	  solution_bricks = optimizer_method_selection<Carrier<InputHandler,Temporal,Areal>>(carrier);

	  lambda_inference_selection(optimizationData, solution_bricks.second, inferenceData, regression, lambda_inference_S, lambda_inference_T); // Set lambdas for inference
	}
      else
	{
	  //Rprintf("Temporal, Pointwise\n");
	  Carrier<InputHandler,Temporal>
	    carrier = CarrierBuilder<InputHandler>::build_temporal_plain_carrier(regressionData, regression, optimizationData);
	  solution_bricks = optimizer_method_selection<Carrier<InputHandler,Temporal>>(carrier);

	  lambda_inference_selection(optimizationData, solution_bricks.second, inferenceData, regression, lambda_inference_S, lambda_inference_T); // Set lambdas for inference
	}
    }

  //!Inference 
  if(inferenceData.get_definition()==true){ 

    //!Only if inference is actually required
    Inference_Carrier<InputHandler> inf_car(&regressionData, &regression, &solution_bricks.second,  &inferenceData, lambda_inference_S, lambda_inference_T); //Carrier for inference
    
    //get the component on which inference is required
    const std::vector<std::string> inf_component = inferenceData.get_component_type(); 

    //if nonparametric inference is required
    if(std::find(inf_component.begin(), inf_component.end(), "nonparametric") != inf_component.end() || 
       std::find(inf_component.begin(), inf_component.end(), "both") != inf_component.end()){
      // set the solution of the system inside the inference carrier
      inf_car.setSolutionp(&(solution_bricks.first));
      // compute other local matrices according to the implementation
      compute_nonparametric_inference_matrices_time<InputHandler, ORDER, mydim, ndim>(mesh, mesh_time, regressionData, inferenceData, inf_car);
    }

    
    inference_wrapper_time(optimizationData, solution_bricks.second, inf_car, inference_Output);    
  }
	
  //Se si riesce, usare plain anche nel caso temporal
  return Solution_Builders::build_solution_temporal_regression<InputHandler, ORDER, mydim, ndim>(solution_bricks.first, solution_bricks.second, mesh, regressionData, regression, inference_Output,inferenceData);
}



//! Function to select the right optimization method
/*
  \tparam CarrierType the type of Carrier to be employed
  \param carrier the Carrier used for the methods
  \return the solution to pass to the Solution_Builders
*/
template<typename CarrierType>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Temporal, CarrierType>::value>, t_type>::value,
			std::pair<MatrixXr, output_Data<2>> >::type optimizer_method_selection(CarrierType & carrier)
{
  // Build the optimizer
  const OptimizationData * optr = carrier.get_opt_data();
	
  if(optr->get_loss_function() == "GCV" && optr->get_DOF_evaluation() == "exact")
    {
      //Rprintf("GCV exact\n");
      if(carrier.get_flagParabolic())	
	return parabolic_routine<GCV_Exact<CarrierType, 1>>(carrier);
      else
	{
	  GCV_Exact<CarrierType, 2> optim(carrier);
	  return optimizer_strategy_selection<GCV_Exact<CarrierType, 2>, CarrierType, 2>(optim, carrier);
	}

    }
  else if(optr->get_loss_function() == "GCV" && (optr->get_DOF_evaluation() == "stochastic" || optr->get_DOF_evaluation() == "not_required"))
    {		
      //Rprintf("GCV stochastic\n");
      if(carrier.get_flagParabolic())
	return parabolic_routine<GCV_Stochastic<CarrierType, 1>>(carrier);
      else
	{
	  GCV_Stochastic<CarrierType, 2> optim(carrier, true);
	  return optimizer_strategy_selection<GCV_Stochastic<CarrierType, 2>, CarrierType, 2>(optim, carrier);
	}
    }
  else // if(optr->get_loss_function() == "unused" && optr->get_DOF_evaluation() == "not_required")
    {
      //Rprintf("Pure evaluation\n");
      GCV_Stochastic<CarrierType, 2> optim(carrier, false);

      timer Time_partial; // Of the sole optimization
      Time_partial.start();
      // Rprintf("WARNING: start taking time\n");

      // Get the solution
      output_Data<2> output;
      output.size_S = carrier.get_opt_data()->get_size_S();
      output.size_T = carrier.get_opt_data()->get_size_T();
      UInt lambdas_count = carrier.get_opt_data()->get_size_S()*carrier.get_opt_data()->get_size_T();
      output.z_hat.resize(carrier.get_psip()->rows(), lambdas_count);
      output.lambda_vec.reserve(lambdas_count);
      output.lambda_vec.clear();
      output.GCV_evals.resize(lambdas_count, -1);
      output.dof.resize(lambdas_count, -1);
      MatrixXr solution;
      MatrixXv betas;
      betas.resize(lambdas_count, 1);
      for(UInt j=0; j<carrier.get_opt_data()->get_size_T(); j++)
	{
	  for(UInt i=0; i<carrier.get_opt_data()->get_size_S(); i++)
	    {
	      Real lambdaS = carrier.get_opt_data()->get_lambda_S()[i];
	      Real lambdaT = carrier.get_opt_data()->get_lambda_T()[j];
	      lambda::type<2> lambda = lambda::make_pair(lambdaS, lambdaT);
	      output.lambda_vec.push_back(lambda);
	      UInt couple_index = j*carrier.get_opt_data()->get_size_S()+i;
	      if(i==0 && j==0)
		{
		  MatrixXr sol = carrier.apply(lambda); //apply_iterative
		  solution.resize(sol.rows(),lambdas_count);
		  solution.col(couple_index) = sol;
		}
	      else
		{
		  solution.col(couple_index) = carrier.apply(lambda); //apply_iterative
		}
	      optim.combine_output_prediction(solution.topRows(solution.rows()/2).col(couple_index),output,couple_index);
	      if(carrier.get_model()->getBeta().cols()>0 && carrier.get_model()->getBeta().rows()>0)
		betas.coeffRef(couple_index,0) = carrier.get_model()->getBeta().coeffRef(0,0);
	    }
	}

      // Rprintf("WARNING: partial time after the optimization method\n");
      timespec T = Time_partial.stop();

      output.time_partial = T.tv_sec + 1e-9*T.tv_nsec;

      // postponed after apply in order to have betas computed
      output.betas = betas;
		
      return {solution, output};
    }
}

template<class GCV_type, typename CarrierType>
std::pair<MatrixXr, output_Data<2>> parabolic_routine(CarrierType & carrier)
{
  timer Time_partial;
  Time_partial.start();

  std::pair<MatrixXr, output_Data<1>> opt_res;
  std::pair<MatrixXr, output_Data<2>> res;
	
  const OptimizationData * optr = carrier.get_opt_data();

  std::vector<Real> lambdaS = optr->get_lambda_S();
  std::vector<Real> lambdaT = optr->get_lambda_T();
  if(optr->get_criterion() != "grid"){	// set initial data for Newton
    lambdaS[0] = optr->get_initial_lambda_S();
    lambdaT[0] = optr->get_initial_lambda_T();
  }
  res.second.lambda_vec.reserve(res.second.size_S*res.second.size_T);
  res.second.GCV_evals.reserve(res.second.size_S*res.second.size_T);
  res.second.lambda_vec.clear();
  res.second.GCV_evals.clear();
				
  for(UInt j=0; j<optr->get_size_T(); j++)
    {			
      GCV_type optim(carrier, lambdaT[j]);
      std::pair<MatrixXr, output_Data<1>> par_res = 
	optimizer_strategy_selection<GCV_type, CarrierType, 1>(optim, carrier);

      UInt dim = par_res.second.lambda_vec.size();

      for(UInt i=0; i<dim; i++)
	res.second.lambda_vec.push_back(
					lambda::make_pair(par_res.second.lambda_vec[i], lambdaT[j]));
		
      if(par_res.second.GCV_opt < opt_res.second.GCV_opt || j == 0)
	{
	  opt_res = par_res;
	  res.second.lambda_sol =
	    lambda::make_pair(par_res.second.lambda_sol, lambdaT[j]);
	  res.second.lambda_pos = j*optr->get_size_S()+par_res.second.lambda_pos;
	  res.second.termination = par_res.second.termination;
	}
		
      res.second.rmse.insert(res.second.rmse.end(), par_res.second.rmse.begin(), par_res.second.rmse.end());
      res.second.dof.insert(res.second.dof.end(), par_res.second.dof.begin(), par_res.second.dof.end());
      res.second.GCV_evals.insert(res.second.GCV_evals.end(), par_res.second.GCV_evals.begin(), par_res.second.GCV_evals.end());
		
      res.second.n_it += par_res.second.n_it;
    }

  if(optr->get_criterion() == "grid")
    res.second.size_S = optr->get_size_S();
  else
    res.second.size_S = res.second.lambda_vec.size();  // size_S is the number of lamdaS visited by Newton method
  res.second.size_T = optr->get_size_T();
  res.first = opt_res.first;
  res.second.content = opt_res.second.content;
  timespec T = Time_partial.stop();
  res.second.time_partial = T.tv_sec + 1e-9*T.tv_nsec;
  res.second.z_hat = opt_res.second.z_hat;
  res.second.sigma_hat_sq = opt_res.second.sigma_hat_sq;
  res.second.betas = opt_res.second.betas;
  res.second.GCV_opt = opt_res.second.GCV_opt;
	
  return res;
}

//! Function to apply the optimization strategy, grid or Newton
/*
  \\tparam EvaluationType optimization type to be used
  \tparam CarrierType the type of Carrier to be employed
  \param optim EvaluationType containing the class related to the function to be optimized,together with the method (exact or stochastic)
  \param carrier the Carrier used for the methods
  \return the solution to pass to the Solution_Builders
*/
template<typename EvaluationType, typename CarrierType, UInt size>
typename std::enable_if<size==2, std::pair<MatrixXr, output_Data<2>>>::type
optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier)
{
  // Build wrapper and newton method
  Function_Wrapper<lambda::type<2>, Real, lambda::type<2>, MatrixXr, EvaluationType> Fun(optim);
  typedef Function_Wrapper<lambda::type<2>, Real, lambda::type<2>, MatrixXr, EvaluationType> FunWr;
	
  const OptimizationData * optr = carrier.get_opt_data();
  if(optr->get_criterion() == "grid")
    {
      timer Time_partial; // Of the sole optimization
      Time_partial.start();
      // Rprintf("WARNING: start taking time\n");

      // this will be used when grid will be correctly implemented, also for return elements
		
      UInt lambdas_count = carrier.get_opt_data()->get_size_S()*carrier.get_opt_data()->get_size_T();
      std::vector<lambda::type<2>> lambda_vec;
      lambda_vec.reserve(lambdas_count);
      for(UInt j=0; j<carrier.get_opt_data()->get_size_T(); j++)
	{
	  Real lambdaT = carrier.get_opt_data()->get_lambda_T()[j];
	  for(UInt i=0; i<carrier.get_opt_data()->get_size_S(); i++)
	    {
	      Real lambdaS = carrier.get_opt_data()->get_lambda_S()[i];
	      lambda::type<2> lambda = lambda::make_pair(lambdaS, lambdaT);
	      lambda_vec.push_back(lambda);
	    }
	}
		
      Eval_GCV<lambda::type<2>, MatrixXr, EvaluationType> eval(Fun, lambda_vec);
      output_Data<2> output = eval.Get_optimization_vectorial();

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
      std::unique_ptr<Opt_methods<lambda::type<2>,MatrixXr,EvaluationType>> optim_p =
	Opt_method_factory<lambda::type<2>,MatrixXr,EvaluationType>::create_Opt_method(optr->get_criterion(), Fun);			

      // Choose initial lambdaS, lambdaT with grid
      Real lambdaS = optr->get_initial_lambda_S();
      Real lambdaT = optr->get_initial_lambda_T();

      lambda::type<2> lambda_init = lambda::make_pair(lambdaS, lambdaT);
      std::vector<lambda::type<2>> lambda_vec;
      if(lambda_init(0)>0 && lambda_init(1)>0){
	lambda_vec.reserve(5);
	lambda_vec.push_back(lambda_init);
      }
      else
	lambda_vec.reserve(4);
      lambda_vec.push_back(lambda::make_pair(exp(-3), exp(-9)));
      lambda_vec.push_back(lambda::make_pair(exp(-5), exp(-7)));
      lambda_vec.push_back(lambda::make_pair(exp(-7), exp(-5)));
      lambda_vec.push_back(lambda::make_pair(exp(-9), exp(-3)));

		
      UInt dim = lambda_vec.size();
      lambda::type<2> lambda_min;
      Real GCV_min = -1.0;

      for (UInt i=0; i<dim; i++)
	{
	  Rprintf("Pre-Newton grid: evaluating %d/%d\n", i+1, dim);
	  Real evaluation = Fun.evaluate_f(lambda_vec[i]); //only scalar functions;

	  if (evaluation<GCV_min || i==0)
	    {
	      GCV_min = evaluation;
	      lambda_min = lambda_vec[i];
	    }
	}

      // Compute optimal lambda
      Checker ch;
      std::vector<lambda::type<2>> lambda_v_;
      std::vector<Real> GCV_v_;

      timer Time_partial; // Of the sole optimization
      Time_partial.start();
      // Rprintf("WARNING: start taking time\n");

      std::pair<lambda::type<2>, UInt> lambda_couple = optim_p->compute(lambda_min, optr->get_stopping_criterion_tol(), 40, ch, GCV_v_, lambda_v_);

      //Rprintf("WARNING: partial time after the optimization method\n");
      timespec T = Time_partial.stop();

      // Get the solution
      // to compute f and g hat
      MatrixXr solution = carrier.apply(lambda_couple.first); //apply_iterative (potrebbe esistere nel GCV un apply_to_b iterative?)

      // postponed after apply in order to have betas computed
      // now the last values in GCV_exact are the correct ones, related to the last iteration
      output_Data<2> output = Fun.get_output(lambda_couple, T, GCV_v_, lambda_v_, ch.which());
      // the copy is necessary for the bulders outside

      return {solution, output};
    }
}



//! Function to select the right inference method
/*
  \tparam InputHandler the type of regression problem
  \param opt_data the object containing optimization data
  \param output the output_Data object containing the solution to be returned
  \param inf_car the object wrapping all the objects needed to make inference
  \param inference_output the object to be filled with inference output 
  \return void
*/
template<typename InputHandler>
void inference_wrapper_time(const OptimizationData & opt_data, output_Data<2> & output, const Inference_Carrier<InputHandler> & inf_car, MatrixXv & inference_output)
{
  UInt n_implementations = inf_car.getInfData()->get_implementation_type().size();
  UInt p = inf_car.getInfData()->get_coeff_inference().rows();
  UInt n_loc = inf_car.getN_loc();

  UInt out_dim = (p > n_loc) ? p : n_loc; 

  // since only inference on beta is implemented for ST, the rows corresponding to f inference will be empty to be coherent with the size of inference_Solver's output
  inference_output.resize(2*n_implementations+1, out_dim+1);

  // Select the right policy for inversion of MatrixNoCov
  std::shared_ptr<Inverse_Base<MatrixXr>> inference_Inverter = std::make_shared<Inverse_Exact>(inf_car.getEp(), inf_car.getE_decp());

  for(UInt i=0; i<n_implementations; ++i){
    // Factory instantiation for solver: using factory provided in Inference_Factory.h
    std::shared_ptr<Inference_Base<InputHandler,MatrixXr>> inference_Solver = Inference_Factory<InputHandler,MatrixXr>::create_inference_method(inf_car.getInfData()->get_implementation_type()[i], inference_Inverter, inf_car, i); // Selects the right implementation and solves the inferential problems
    inference_output.middleRows(2*i,2) = inference_Solver->compute_inference_output();

    if(inf_car.getInfData()->get_implementation_type()[i]=="wald" && opt_data.get_loss_function()=="unused" && opt_data.get_size_S()==1 && opt_data.get_size_T()==1){
      output.GCV_opt=inference_Solver->compute_GCV_from_inference(); // Computing GCV if Wald has being called is an almost zero-cost function, since tr(S) hase been already computed
    }
  }
    
  // Check if local f variance has to be computed
  if(inf_car.getInfData()->get_f_var()){
    std::shared_ptr<Inference_Base<InputHandler,MatrixXr>> inference_Solver = Inference_Factory<InputHandler,MatrixXr>::create_inference_method("wald", inference_Inverter, inf_car, n_implementations);
    inference_output(2*n_implementations,0) = inference_Solver->compute_f_var();
  }

  return;
  
}

//! Function that sets the correct lambda needed for inferential operations
/*
  \tparam InputHandler the type of regression problem
  \param optimization_data the object containing optimization data
  \param output the output_Data object containing the solution of the optimization problem
  \param inferenceData the object containing the data needed for inference
  \param regression the object containing the model of the problem
  \param lambda_inference_S the spatial smoothing parameter that will be used to compute the optimal model and the right inferential solutions
  \param lambda_inference_T the temporal smoothing parameter that will be used to compute the optimal model and the right inferential solutions
  \return void
*/
template<typename InputHandler>
void lambda_inference_selection (const OptimizationData & optimizationData, const output_Data<2> & output, const InferenceData & inferenceData, MixedFERegression<InputHandler> & regression, Real & lambda_inference_S, Real & lambda_inference_T){
  if(inferenceData.get_definition()==true && optimizationData.get_loss_function()!="unused"){
    lambda_inference_S = output.lambda_sol(0);
    lambda_inference_T = output.lambda_sol(1);
    if(optimizationData.get_last_lS_used() != lambda_inference_S || optimizationData.get_last_lT_used() != lambda_inference_T){
      regression.build_regression_inference(lambda_inference_S, lambda_inference_T);
    }
  }else{ 		// supposing we have only one lambda when GCV is unused, otherwise inference gets discarded in smoothing.R
    if(inferenceData.get_definition()==true){
      lambda_inference_S = optimizationData.get_last_lS_used();
      lambda_inference_T = optimizationData.get_last_lT_used();
    }
  }
  return; 
}

//! Function that evaluates the spatial basis functions in a set of new location points, needed for inference on f
/*
  \tparam InputHandler the type of regression problem
  \tparam ORDER the order of the mesh 
  \tparam mydim specifies if the mesh lie in R^2 or R^3
  \tparam ndim specifies if the local dimension is 2 or 3
  \param mesh_ the mesh of the problem
  \param regressionData_ the object containing regression informations 
  \param inferenceData_ the object containing the data needed for inference
  \param inf_car_ the inference carrier object to be modified 
  \return void
*/
template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void compute_nonparametric_inference_matrices_time(const MeshHandler<ORDER, mydim, ndim>  & mesh_, const std::vector<Real> meshTime_, const InputHandler & regressionData_, InferenceData & inferenceData_, Inference_Carrier<InputHandler> & inf_car_){
  // if a matrix of locations has been provided, compute spatial Psi_loc by directly evaluating the spatial basis functions in the provided points
  
  // define the psi matrix
  SpMat psi_temp;
  
  if((inferenceData_.get_locs_index_inference())[0] == -1){
    // first fetch the dimensions
    UInt nnodes = mesh_.num_nodes();
    UInt nlocations = (inferenceData_.get_locs_inference()).rows();
    psi_temp.resize(nlocations, nnodes);
		
    constexpr UInt EL_NNODES = how_many_nodes(ORDER,mydim);
    Eigen::Matrix<Real,EL_NNODES,1> coefficients;    // Dummy for point evaluation
    Real evaluator;                // Dummy for evaluation storage

    for(UInt i=0; i<nlocations;i++)
      { // Update psi looping on all locations
	// [[GM missing a defaulted else, raising a WARNING!]]
	VectorXr coords = (inferenceData_.get_locs_inference()).row(i);
	Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.findLocation(Point<ndim>(i, coords));

	// Search the element containing the point
	if(tri_activated.getId() == Identifier::NVAL)
	  { // If not found
	    Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
	  }
	else
	  {
	    for(UInt node=0; node<EL_NNODES ; ++node)
	      {// Loop on all the nodes of the found element and update the related entries of Psi
		// Define vector of all zeros but "node" component (necessary for function evaluate_point)
		coefficients = Eigen::Matrix<Real,EL_NNODES,1>::Zero();
		coefficients(node) = 1; //Activates only current base-node
		// Evaluate psi in the node
		evaluator = tri_activated.evaluate_point(Point<ndim>(i, coords), coefficients);
		// Insert the value in the column given by the GLOBAL indexing of the evaluated NODE
		psi_temp.insert(i, tri_activated[node].getId()) = evaluator;
	      }
	  }
      } // End of for loop
        
    psi_temp.makeCompressed();   
           
    // update the inference carrier
    inf_car_.setPsi_loc(psi);
    inf_car_.setN_loc(psi.rows());                 	
  }
  else{
    // the locations are chosen among the observed ones, hence psi can be extracted from Psi
    const std::vector<UInt> row_indices = inferenceData_.get_locs_index_inference();

    // vector that converts global indices into local indices, common for both Psi_loc and W_loc
    VectorXi rel_rows = VectorXi::Constant(inf_car_.getPsip()->rows(), -1);
    for(UInt i=0; i < row_indices.size(); ++i){
      rel_rows(row_indices[i]) = i; 
    } 
	
    UInt nnodes = mesh_.num_nodes();
    UInt nlocations = row_indices.size();

    psi_temp.resize(nlocations, nnodes);

    if(nlocations == inf_car_.getPsip()->rows()){
      psi_temp = *(inf_car_.getPsip());
    }
    else{

      // vector storing the non zero elements of Psi to be inserted in psi --> they are at most Psi.nonZeros()
      std::vector<coeff> coefficients;
      coefficients.reserve(inf_car_.getPsip()->nonZeros());

      // loop over the nonzero elements
      for (UInt k = 0; k < inf_car_.getPsip()->outerSize(); ++k){
	for (SpMat::InnerIterator it(*(inf_car_.getPsip()),k); it; ++it)
	  {
	    if(std::find(row_indices.begin(), row_indices.end(), it.row()) != row_indices.end()){
	      coefficients.push_back(coeff(rel_rows(it.row()), it.col(), it.value()));
	    }
	  }
      }

      psi_temp.setFromTriplets(coefficients.begin(), coefficients.end());
      psi_temp.makeCompressed();
    }
    
  }
  // Now we compute the temporal Phi_loc and then we compute the kroenecker product
  VectorXr time_locs_inf = inf_car.getInfData()->get_time_locs_inf();
    
  Spline<SPLINE_DEGREE, ORDER_DERIVATIVE> spline(meshTime_);
  UInt M = spline.num_knots()-SPLINE_DEGREE-1;
  UInt m = time_locs_inf.size();
    
  phi_.resize(m, M);
  Real value;

  for(UInt i=0; i<m; ++i)
    {
      for(UInt j=0; j<M; ++j)
	{
	  value = spline.BasisFunction(j, time_locs_inf[i]);
	  if(value!=0)
	    {
	      phi_.coeffRef(i,j) = value;
	    }
	}
    }
  phi_.makeCompressed();
    
  // set the spatial Psi_loc into the inference carrier
  inf_car_.setPsi_loc(kroneckerProduct(phi, psi_temp));
  inf_car_.setN_loc(psi.rows());
}

#endif
