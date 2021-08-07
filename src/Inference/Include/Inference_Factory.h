#ifndef __INFERENCE_FACTORY_H__
#define __INFERENCE_FACTORY_H__

// HEADERS
#include "Inference_Base.h"
#include "Wald.h"
#include "Speckman.h"
#include "Eigen_Sign_Flip.h"
#include <memory>
#include <map>
#include <type_traits>

//! A Factory class: A class for the choice of implementation for the computation of inferential objects.
/* \tparam InputHandler RegressionData of the problem
 * \tparam MatrixType the parameter that characterises the inverter object
 */
template<typename InputHandler, typename MatrixType>
class Inference_Factory
{
private:
  //!< Getter of the Inference_Factory
  static std::map<std::string,std::shared_ptr<Inference_Base<InputHandler, MatrixType>>>& get_Factory_Store(void)
  {
    static std::map<std::string,std::shared_ptr<Inference_Base<InputHandler, MatrixType>>> factory_Store; // initialize the static map
    return factory_Store;
  }

public:
  //! A method that takes as parameter a string and builds a pointer to the right implementation object
  /*!
    \param implementation_type_ type of implementation required
    \param inverter_ class demanded to the computation of MatrixNoCov inverse
    \param inf_car_ inference carrier object
    \param pos_impl_ the position index of the current inferential procedure
    \return std::shared_ptr to the chosen solver
  */
  static std::shared_ptr<Inference_Base<InputHandler, MatrixType>> create_inference_method(const std::string & implementation_type_, std::shared_ptr<Inverse_Base<MatrixType>> inverter_, const Inference_Carrier<InputHandler> & inf_car_, UInt pos_impl_)
  {
    std::map<std::string,std::shared_ptr<Inference_Base<InputHandler, MatrixType>>> factory_Store=get_Factory_Store(); // Get the static factory
    
    if(implementation_type_=="wald"){
      // check if it is exact or non-exact according to MatrixType
      if(std::is_same<MatrixType, SpMat>::value==false){
        // look if the object is already present in the factory
        auto It = factory_Store.find("wald_exact");
        // if not, insert the new object
	if(It==factory_Store.end()){
	  factory_Store.insert(std::make_pair<std::string, std::shared_ptr<Inference_Base<InputHandler, MatrixType>>>("wald_exact", std::make_shared<Wald_Exact<InputHandler, MatrixType>>(inverter_, inf_car_, pos_impl_)));
	}else{
          // if it is already in the factory, just update the pos_impl
	  It->second->setpos_impl(pos_impl_);
	}
	return factory_Store["wald_exact"];
      }else{
        // look if the object is already present in the factory
	auto It = factory_Store.find("wald_non_exact");
        // if not, insert the new object
	if(It==factory_Store.end()){
	  factory_Store.insert(std::make_pair<std::string, std::shared_ptr<Inference_Base<InputHandler, MatrixType>>>("wald_non_exact", std::make_shared<Wald_Non_Exact<InputHandler, MatrixType>>(inverter_, inf_car_, pos_impl_)));
	}else{
          // if it is already in the factory, just update the pos_impl
	  It->second->setpos_impl(pos_impl_);
	}
	return factory_Store["wald_non_exact"];
      }
    }

    if(implementation_type_=="speckman"){
      // check if it is exact or non-exact according to MatrixType
      if(std::is_same<MatrixType, SpMat>::value==false){
        // look if the object is already present in the factory
	auto It = factory_Store.find("speckman_exact");
        // if not, insert the new object
	if(It==factory_Store.end()){
	  factory_Store.insert(std::make_pair<std::string, std::shared_ptr<Inference_Base<InputHandler, MatrixType>>>("speckman_exact", std::make_shared<Speckman_Exact<InputHandler, MatrixType>>(inverter_, inf_car_, pos_impl_)));
	}else{
          // if it is already in the factory, just update the pos_impl
	  It->second->setpos_impl(pos_impl_);
	}
	return factory_Store["speckman_exact"];
      }else{
        // look if the object is already present in the factory
	auto It = factory_Store.find("speckman_non_exact");
        // if not, insert the new object
	if(It==factory_Store.end()){
	  factory_Store.insert(std::make_pair<std::string, std::shared_ptr<Inference_Base<InputHandler, MatrixType>>>("speckman_non_exact", std::make_shared<Speckman_Non_Exact<InputHandler, MatrixType>>(inverter_, inf_car_, pos_impl_)));
	}else{
          // if it is already in the factory, just update the pos_impl
	  It->second->setpos_impl(pos_impl_);
	}
	return factory_Store["speckman_non_exact"];
      }
    }

    if(implementation_type_=="eigen-sign-flip"){
      // check if it is exact or non-exact according to MatrixType
      if(std::is_same<MatrixType, SpMat>::value==false){
        // look if the object is already present in the factory
	auto It = factory_Store.find("eigen-sign-flip_exact");
        // if not, insert the new object
	if(It==factory_Store.end()){
	  factory_Store.insert(std::make_pair<std::string, std::shared_ptr<Inference_Base<InputHandler, MatrixType>>>("eigen-sign-flip_exact", std::make_shared<Eigen_Sign_Flip_Exact<InputHandler, MatrixType>>(inverter_, inf_car_, pos_impl_)));
	}else{
          // if it is already in the factory, just update the pos_impl
	  It->second->setpos_impl(pos_impl_);
	}
	return factory_Store["eigen-sign-flip_exact"];
      }else{
        // look if the object is already present in the factory
	auto It = factory_Store.find("eigen-sign-flip_non_exact");
        // if not, insert the new object
	if(It==factory_Store.end()){
	  factory_Store.insert(std::make_pair<std::string, std::shared_ptr<Inference_Base<InputHandler, MatrixType>>>("eigen-sign-flip_non_exact", std::make_shared<Eigen_Sign_Flip_Non_Exact<InputHandler, MatrixType>>(inverter_, inf_car_, pos_impl_)));
	}else{
          // if it is already in the factory, just update the pos_impl
	  It->second->setpos_impl(pos_impl_);
	}
	return factory_Store["eigen-sign-flip_non_exact"]; 
      }
    }
    else // default Wald
      {
	if(std::is_same<MatrixType, SpMat>::value==false){
    
	  Rprintf("Implementation not found, using wald exact");
	  auto It = factory_Store.find("wald_exact");
	  if(It==factory_Store.end()){
	    factory_Store.insert(std::make_pair<std::string, std::shared_ptr<Inference_Base<InputHandler, MatrixType>>>("wald_exact", std::make_shared<Wald_Exact<InputHandler, MatrixType>>(inverter_, inf_car_, pos_impl_)));
	  }else{
	    It->second->setpos_impl(pos_impl_);
	  }
	  return factory_Store["wald_exact"];
	}else{
	  Rprintf("Implementation not found, using wald non exact");
	  auto It = factory_Store.find("wald_non_exact");
	  if(It==factory_Store.end()){
	    factory_Store.insert(std::make_pair<std::string, std::shared_ptr<Inference_Base<InputHandler, MatrixType>>>("wald_non_exact", std::make_shared<Wald_Non_Exact<InputHandler, MatrixType>>(inverter_, inf_car_, pos_impl_)));
	  }else{
	    It->second->setpos_impl(pos_impl_);
	  }
	  return factory_Store["wald_non_exact"];
	}
      }
  }
};

#endif
