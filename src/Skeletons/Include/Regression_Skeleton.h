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
#include <algorithm>
#include <set>

template<typename CarrierType>
std::pair<MatrixXr, output_Data> optimizer_method_selection(CarrierType & carrier);
template<typename EvaluationType, typename CarrierType>
std::pair<MatrixXr, output_Data> optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier);
template<typename InputHandler>
void inference_wrapper_space(const OptimizationData & opt_data, output_Data & output, const Inference_Carrier<InputHandler> & inf_car, MatrixXv & inference_output);
template<typename InputHandler>
void lambda_inference_selection(const OptimizationData & optimizationData, const output_Data & output, const InferenceData & inferenceData, MixedFERegression<InputHandler> & regression, Real & lambda_inference);
template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void compute_nonparametric_inference_matrices(const MeshHandler<ORDER, mydim, ndim>  & mesh, const InputHandler & regressionData, const InferenceData & inferenceData, Inference_Carrier<InputHandler> & inf_car);

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

	  //get the component on which inference is required
	  const std::vector<std::string> inf_component = inferenceData.get_component_type(); 

	  //if nonparametric inference is required
	  if(std::find(inf_component.begin(), inf_component.end(), "nonparametric") != inf_component.end() || 
	     std::find(inf_component.begin(), inf_component.end(), "both") != inf_component.end()){
	    // set the solution of the system inside the inference carrier
	    inf_car.setSolutionp(&(solution_bricks.first));
	    // compute other local matrices according to the implementation
	    compute_nonparametric_inference_matrices<InputHandler, ORDER, mydim, ndim>(mesh, regressionData, inferenceData, inf_car);
	  }

	  inference_wrapper_space(optimizationData, solution_bricks.second, inf_car, inference_Output);    
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
void inference_wrapper_space(const OptimizationData & opt_data, output_Data & output, const Inference_Carrier<InputHandler> & inf_car, MatrixXv & inference_output)
{
  // Get the number of implementations
  UInt n_implementations = inf_car.getInfData()->get_implementation_type().size();

  UInt p = inf_car.getInfData()->get_coeff_inference().rows();
  UInt n_loc = inf_car.getN_loc();

  UInt out_dim = (p > n_loc) ? p : n_loc; 

  inference_output.resize(2*n_implementations+1, out_dim+1);

  if(inf_car.getInfData()->get_exact_inference() == "exact"){
    // Select the right policy for inversion of MatrixNoCov
    std::shared_ptr<Inverse_Base<MatrixXr>> inference_Inverter = std::make_shared<Inverse_Exact>(inf_car.getEp(), inf_car.getE_decp());

    for(UInt i=0; i<n_implementations; ++i){
      // Factory instantiation for solver: using factory provided in Inference_Factory.h
      std::shared_ptr<Inference_Base<InputHandler,MatrixXr>> inference_Solver = Inference_Factory<InputHandler,MatrixXr>::create_inference_method(inf_car.getInfData()->get_implementation_type()[i], inference_Inverter, inf_car, i); // Selects the right implementation and solves the inferential problems		
		
      inference_output.middleRows(2*i,2) = inference_Solver->compute_inference_output();

      if(inf_car.getInfData()->get_implementation_type()[i]=="wald" && opt_data.get_loss_function()=="unused" && opt_data.get_size_S()==1){
	output.GCV_opt=inference_Solver->compute_GCV_from_inference(); // Computing GCV if Wald has being called is an almost zero-cost function, since tr(S) hase been already computed
      }
    }
    
    // Check if local f variance has to be computed
    if(inf_car.getInfData()->get_f_var()){
      std::shared_ptr<Inference_Base<InputHandler,MatrixXr>> inference_Solver = Inference_Factory<InputHandler,MatrixXr>::create_inference_method("wald", inference_Inverter, inf_car, n_implementations);
      inference_output(2*n_implementations,0) = inference_Solver->compute_f_var();
    }
  }
  else{
    // Select the right policy for inversion of MatrixNoCov
    std::shared_ptr<Inverse_Base<SpMat>> inference_Inverter = std::make_shared<Inverse_Non_Exact<InputHandler>>(inf_car);

    for(UInt i=0; i<n_implementations; ++i){
      // Factory instantiation for solver: using factory provided in Inference_Factory.h
      std::shared_ptr<Inference_Base<InputHandler,SpMat>> inference_Solver = Inference_Factory<InputHandler,SpMat>::create_inference_method(inf_car.getInfData()->get_implementation_type()[i], inference_Inverter, inf_car, i); // Selects the right implementation and solves the inferential problems		
		
      inference_output.middleRows(2*i,2) = inference_Solver->compute_inference_output();


    }
    
    // Check if local f variance has to be computed
    if(inf_car.getInfData()->get_f_var()){
      std::shared_ptr<Inference_Base<InputHandler,SpMat>> inference_Solver = Inference_Factory<InputHandler,SpMat>::create_inference_method("wald", inference_Inverter, inf_car, n_implementations);
      inference_output(2*n_implementations,0) = inference_Solver->compute_f_var();
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
			//Rprintf("I'm computing again the matrices in Mixed_FERegression\n");
			}
	}else{ 		// supposing we have only one lambda, otherwise inference gets discarded in smoothing.R
		if(inferenceData.get_definition()==true){
			lambda_inference = optimizationData.get_last_lS_used();
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
  \param inferenceData_ the object containing the data needed for for inference
  \param inf_car_ the inference carrier object to be modified 
  \return void
*/
template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void compute_nonparametric_inference_matrices(const MeshHandler<ORDER, mydim, ndim>  & mesh_, const InputHandler & regressionData_, const InferenceData & inferenceData_, Inference_Carrier<InputHandler> & inf_car_){
  // if a matrix of locations has been provided, compute Psi_loc by directly evaluating the spatial basis functions in the provided points
  // only wald implementation can enter here, no other additional matrices are needed
  if((inferenceData_.get_locs_index_inference())[0] == -1){
    // define the psi matrix
    SpMat psi;
    // first fetch the dimensions
    UInt nnodes = mesh_.num_nodes();
    UInt nlocations = (inferenceData_.get_locs_inference()).rows();
    psi.resize(nlocations, nnodes);
		
    constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;
    Element<Nodes, mydim, ndim> tri_activated;	// Dummy for element search
    Eigen::Matrix<Real,Nodes,1> coefficients;	// Dummy for point evaluation
    Real evaluator;				// Dummy for evaluation storage

    for(UInt i=0; i<nlocations;++i)
      { // Update psi looping on all locations
	// [[GM missing a defaulted else, raising a WARNING!]]
	VectorXr coords = (inferenceData_.get_locs_inference()).row(i);
	if(regressionData_.getSearch() == 1)
	  { // Use Naive search
	    tri_activated = mesh_.findLocationNaive(Point<ndim>(i, coords));
	  }
	else if(regressionData_.getSearch() == 2)
	  { // Use Tree search (default)
	    tri_activated = mesh_.findLocationTree(Point<ndim>(i, coords));
	  }

	// Search the element containing the point
	if(tri_activated.getId() == Identifier::NVAL)
	  { // If not found
	    Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform inference\n", i+1);
	  }
	else
	  { // tri_activated.getId() found, it's action might be felt a priori by all the psi of the element, one for each node
				
	    for(UInt node=0; node<Nodes ; ++node)
	      {// Loop on all the nodes of the found element and update the related entries of Psi
		// Define vector of all zeros but "node" component (necessary for function evaluate_point)
		coefficients = Eigen::Matrix<Real,Nodes,1>::Zero();
		coefficients(node) = 1; //Activates only current base-node
					// Evaluate psi in the node
		evaluator = tri_activated.evaluate_point(Point<ndim>(i, coords), coefficients);
		// Insert the value in the column given by the GLOBAL indexing of the evaluated NODE
		psi.insert(i, tri_activated[node].getId()) = evaluator;
	      }
	  }
      } // End of for loop
		
    psi.makeCompressed();  
           
    // update the inference carrier
    inf_car_.setPsi_loc(psi);
    inf_car_.setN_loc(psi.rows());                 	
  }
  else{
    // the locations are chosen among the observed ones, hence psi can be extracted from Psi
    SpMat psi; 
    const std::vector<UInt> row_indices = inferenceData_.get_locs_index_inference();

    // vector that converts global indices into local indices, common for both Psi_loc and W_loc
    VectorXi rel_rows = VectorXi::Constant(inf_car_.getPsip()->rows(), -1);
    for(UInt i=0; i < row_indices.size(); ++i){
      rel_rows(row_indices[i]) = i; 
    } 
	
    UInt nnodes = mesh_.num_nodes();
    UInt nlocations = row_indices.size();

    psi.resize(nlocations, nnodes);

    if(nlocations == inf_car_.getPsip()->rows()){
      psi = *(inf_car_.getPsip());
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

      psi.setFromTriplets(coefficients.begin(), coefficients.end());
      psi.makeCompressed();
    }

    // set it into the inference carrier
    inf_car_.setPsi_loc(psi);
    inf_car_.setN_loc(psi.rows());
         
    

    // in the sign-flip and eigen-sign-flip cases, additional matrices have to be computed
    const std::vector<std::string> implementation_type = inferenceData_.get_implementation_type();

    if(std::find(implementation_type.begin(), implementation_type.end(), "sign-flip") != implementation_type.end() ||
       std::find(implementation_type.begin(), implementation_type.end(), "eigen-sign-flip") != implementation_type.end()){
      // reduced vector of observations		
      VectorXr z_loc; 
      z_loc.resize(row_indices.size());
		
      for(UInt i=0; i < inf_car_.getZp()->size(); ++i){
	if(std::find(row_indices.begin(), row_indices.end(), i) != row_indices.end())
	  z_loc(rel_rows(i)) = inf_car_.getZp()->coeff(i);
      }
		
      inf_car_.setZ_loc(z_loc);
                		
      // reduced design matrix, only if there are covariates
      if(inf_car_.getRegData()->getCovariates()->rows()!=0) {
        MatrixXr W_loc; 
        W_loc.resize(row_indices.size(), inf_car_.getWp()->cols());

        for(UInt i=0; i < inf_car_.getWp()->cols(); ++i){
	  for(UInt j=0; j < inf_car_.getWp()->rows(); ++j){
	    if(std::find(row_indices.begin(), row_indices.end(), j) != row_indices.end())
	      W_loc(rel_rows(j), i) = inf_car_.getWp()->coeff(j,i);
	  }
        }
		
        inf_car_.setW_loc(W_loc);
      }
      // if the selected locations coincide with the nodes compute the matrix that groups closer locations, according to the distance induced by the mesh elements
      if(inferenceData_.get_locs_are_nodes_inference()){
	MatrixXr Group_locs = MatrixXr::Zero(mesh_.num_nodes(), inferenceData_.get_locs_inference().rows());
	// get the selected location indices
	std::vector<UInt> locations_index = inferenceData_.get_locs_index_inference(); 
	// declare the vector that contains, for each location-node, the corresponding neighbors' indices
	std::vector<std::set<UInt>> NearestIndices; 
	NearestIndices.resize(mesh_.num_nodes()); 
       
	for(UInt k=0; k<mesh_.num_nodes(); ++k){
	  // prepare the object to insert
	  std::set<UInt> neighbors;
 
	  // loop on the mesh elements 
	  for(auto i=0; i < mesh_.num_elements(); ++i){
	    auto elem = mesh_.getElement(i);
	    // check if the current point is inside the current element
	    if(elem.isPointInside(mesh_.getPoint(k))){
	      // loop on all the points in the current element and insert them into the set of neighbors
	      for(auto it = elem.begin(); it != elem.end(); ++it){
		neighbors.insert(it->id());
	      }
	    }
	  }
	  // insert the set of neighbors in the final vector
	  NearestIndices[k] = neighbors;
	}

	// compute the group matrix
	for(UInt i=0; i < locations_index.size(); ++i){
	  for(UInt j : NearestIndices[locations_index[i]]){
            if(rel_rows(j)>=0) // only if it belongs to the selected locations
	    	Group_locs(locations_index[i], rel_rows(j)) = 1;
          }
	}
    
	// set it into inference carrier 
	inf_car_.setGroup_loc(Group_locs);
       
      }
    }
  }
	
}

#endif
