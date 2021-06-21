#ifndef __INVERTER_FACTORY_H__
#define __INVERTER_FACTORY_H__

// HEADERS
#include "Inverter.h"
#include "Inference_Carrier.h"
#include "../../Global_Utilities/Include/Make_Unique.h"
#include "../../Global_Utilities/Include/Make_Shared.h"
#include <memory>

//! A Factory class: A class for the choice of the method used to invert MatrixNoCov in inferenctial problems.
template<typename InputHandler>
class Inverter_Factory
{
public:
  //! A method that takes as parameter a string and builds a pointer to object required
  /*!
    \param policy a string code to decide which pointer to create
    \param inf_car inference carrier needed to extract the matrices to be inverted
    \return a pointer to the selected policy
  */
  static std::shared_ptr<Inverse_Base> create_inverter_method(const Inference_Carrier<InputHandler> & inf_car)
  {
    const std::string policy = inf_car.getInfData()->get_exact_inference();
    if(policy=="exact")
      return fdaPDE::make_shared<Inverse_Exact>(inf_car.getEp(), inf_car.getE_decp());

    if(policy!="non-exact"){
      Rprintf("Method not found, using non-exact");
    }

    return fdaPDE::make_shared<Inverse_Non_Exact<InputHandler>>(inf_car);
  }
};

#endif
