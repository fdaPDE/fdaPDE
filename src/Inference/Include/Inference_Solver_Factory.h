#ifndef __INFERENCE_SOLVER_FACTORY_H__
#define __INFERENCE_SOLVER_FACTORY_H__

// HEADERS
#include "Solver_Base.h"
#include "Wald_Solver.h"
#include "Speckman_Solver.h"
#include "Eigen_Sign_Flip_Solver.h"
#include "../../Global_Utilities/Include/Make_Unique.h"
#include <memory>

//! A Factory class: A class for the choice of implementation for the computation of inferential objects.
/* \tparam InputHandler RegressionData of the problem
 */
template<typename InputHandler>
class Inference_Solver_Factory
{
 public:
  //! A method that takes as parameter a string and builds a pointer to the right implementation object
  /*!
    \param implementation_type type of implementation required
    \param inverter_ class demanded to the computation of MatrixNoCov inverse
    \param inf_car_ inference carrier object
    \return std::unique_ptr to the chosen solver
  */
  static std::unique_ptr<Solver_Base<InputHandler>> create_inference_solver_method(const std::string & implementation_type_, const std::string & exact_inference_, Inference_Carrier<InputHandler> inf_car_)
  {
    if(implementation_type_=="wald")
      return make_unique<Wald_Solver<InputHandler>>(exact_inference_, inf_car_);
    if(implementation_type_=="speckman")
      return make_unique<Speckman_Solver<InputHandler>>(exact_inference_, inf_car_);
    if(implementation_type_=="permutational")
      return make_unique<Eigen_Sign_Flip_Solver<InputHandler>>(exact_inference_, inf_car_);
    else // deafult Wald
      {
	Rprintf("Method not found, using Newton_fd");
	return make_unique<Wald_Solver<InputHandler>>(exact_inference_, inf_car_);
      }
  }
};

#endif
