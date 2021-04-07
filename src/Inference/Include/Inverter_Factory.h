#ifndef __INVERTER_FACTORY_H__
#define __INVERTER_FACTORY_H__

// HEADERS
#include "Inverter.h"
#include "../../Global_Utilities/Include/Make_Unique.h"
#include <memory>

//! A Factory class: A class for the choice of the method used to invert MatrixNoCov in inferenctial problems.

class Inverter_Factory
{
public:
  //! A method that takes as parameter a string and builds a pointer to object required
  /*!
    \param policy a string code to decide which pointer to create
    \return a pointer to the selected policy
  */
  static std::unique_ptr<Inverse_Base> create_inverter_method(const std::string & policy)
  {
    if(policy=="exact")
      return make_unique<Inverse_Exact>();
    // if(policy=="non-exact-1")
    //        return make_unique<Inverse_Non_Exact_1>();
    else // default is fd
      {
	Rprintf("Method not found, using non-exact-1");
	return make_unique<Inverse_Exact>(); //// OCCHIO DA SOSTITUIRE
      }
  }
};

#endif
