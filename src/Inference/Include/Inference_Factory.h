#ifndef __INFERENCE_FACTORY_H__
#define __INFERENCE_FACTORY_H__

// HEADERS
#include "Inference_Base.h"
#include "Wald.h"
#include "Speckman.h"
#include "Eigen_Sign_Flip.h"
#include "../../Global_Utilities/Include/Make_Unique.h"
#include <memory>

//! A Factory class: A class for the choice of implementation for the computation of inferential objects.
/* \tparam InputHandler RegressionData of the problem
 */
template<typename InputHandler>
class Inference_Factory
{
public:
  //! A method that takes as parameter a string and builds a pointer to the right implementation object
  /*!
    \param implementation_type type of implementation required
    \param inverter_ class demanded to the computation of MatrixNoCov inverse
    \param inf_car_ inference carrier object
    \return std::unique_ptr to the chosen solver
  */
  static std::unique_ptr<Inference_Base<InputHandler>> create_inference_method(const std::string & implementation_type_, std::unique_ptr<Inverse_Base> inverter_, const Inference_Carrier<InputHandler> & inf_car_)
  {
    if(implementation_type_=="wald")
      return make_unique<Wald<InputHandler>>(std::move(inverter_), inf_car_);
    if(implementation_type_=="speckman")
      return make_unique<Speckman<InputHandler>>(std::move(inverter_), inf_car_);
    if(implementation_type_=="permutational")
      return make_unique<Eigen_Sign_Flip<InputHandler>>(std::move(inverter_), inf_car_);
    else // deafult Wald
      {
	Rprintf("Method not found, using Newton_fd");
	return make_unique<Wald<InputHandler>>(std::move(inverter_), inf_car_);
      }
  }
};

#endif
