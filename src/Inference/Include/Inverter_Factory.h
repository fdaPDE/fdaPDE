#ifndef __INVERTER_FACTORY_H__
#define __INVERTER_FACTORY_H__

// HEADERS
#include "Inverter.h"
#include "Inference_Carrier.h"
#include "../../Global_Utilities/Include/Make_Unique.h"
#include "../../Global_Utilities/Include/Make_Shared.h"
#include <memory>

//! A Factory class: A class for the choice of the method used to invert MatrixNoCov in inferenctial problems.
template<typename InputHandler, MatrixType>
class Inverter_Factory
{
private: 
  static std::map<std::string,std::shared_ptr<Inverter_Base<MatrixType>>>& get_Factory_Store(void)
  {
    static std::map<std::string,std::shared_ptr<Inverter_Base<MatrixType>>> factory_Store; // initialize the static map
    return factory_Store;
  }
public:
  //! A method that takes as parameter a string and builds a pointer to object required
  /*!
    \param policy a string code to decide which pointer to create
    \param inf_car inference carrier needed to extract the matrices to be inverted
    \return a pointer to the selected policy
  */
  static std::shared_ptr<Inverse_Base> create_inverter_method(const Inference_Carrier<InputHandler> & inf_car)
  {
    std::map<std::string,std::shared_ptr<Inverter_Base<MatrixType>>> factory_Store=get_Factory_Store(); // Get the static factory
    
    const std::string policy = inf_car.getInfData()->get_exact_inference();
    if(policy=="exact"){
      auto It = factory_Store.find("exact");
      if(It==factory_Store.end()){
	factory_Store.insert(std::make_pair<std::string, std::shared_ptr<Inverter_Base<MatrixType>>>("exact", fdaPDE::make_shared<Inverse_Exact<MatrixType>>(inf_car.getEp(), inf_car.getE_decp())));
      }
      return factory_Store["exact"];
    }

    if(policy=="non-exact"){
      auto It = factory_Store.find("non-exact");
      if(It==factory_Store.end()){
	factory_Store.insert(std::make_pair<std::string, std::shared_ptr<Inverter_Base<MatrixType>>>("non-exact", fdaPDE::make_shared<Inverse_Non_Exact<InputHandler, MatrixType>>(inf_car)));
      }
      return factory_Store["non-exact"];
    }

    else{
      Rprintf("Method not found, using non-exact");
      if(It==factory_Store.end()){
	factory_Store.insert(std::make_pair<std::string, std::shared_ptr<Inverter_Base<MatrixType>>>("non-exact", fdaPDE::make_shared<Inverse_Non_Exact<InputHandler, MatrixType>>(inf_car)));
      }
      return factory_Store["non-exact"];
    }

    
  }
};

#endif

