#ifndef __INVERTER_FACTORY_H__
#define __INVERTER_FACTORY_H__

// HEADERS
#include "../../FdaPDE.h"
#include "Inverter.h"
#include "../../Global_Utilities/Include/Make_Unique.h"
#include "../../Global_Utilities/Include/Make_Shared.h"
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

  static std::shared_ptr<Inverse_Base> create_inverter_method(const std::string & policy, const SpMat * Ep, const Eigen::SparseLU<SpMat> * E_decp)
  {
    if(policy=="exact")
      return make_shared<Inverse_Exact>(Ep,E_decp);
    if(policy=="non-exact")
      return make_shared<Inverse_Non_Exact>(Ep);
    else // default is fd
      {
	Rprintf("Method not found, using non-exact");
	return make_shared<Inverse_Non_Exact>(Ep);
      }
  }
};

#endif
