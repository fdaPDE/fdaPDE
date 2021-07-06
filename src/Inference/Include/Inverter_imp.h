#include "Inverter.h"
#include "../../Global_Utilities/Include/FSPAI_Wrapper.h"

template<typename InputHandler>
void Inverse_Non_Exact<InputHandler>::pre_Inverse(void){
  SpMat R0 = *(inf_car.getR0p());
  R0.makeCompressed();
  bool status = FSPAI_Wrapper(R0, this->R0_inv_tilde);
  if(!status){
    return; // Add a boolean as a member of Inverse_Non_Exact for well-posedness
  }
  this->R0_inv_tilde.makeCompressed();
  this->E_tilde = (*inf_car.getPsi_tp())*(*inf_car.getPsip()) + inf_car.getlambda()*(inf_car.getR1p()->transpose())*(this->R0_inv_tilde)*(*inf_car.getR1p());
  this->E_tilde.makeCompressed();
};


template<typename InputHandler>
void Inverse_Non_Exact<InputHandler>::Compute_Inv(void){
  if(!this->inverse_computed){
  this->pre_Inverse();
  bool status = FSPAI_Wrapper(this->E_tilde, this->E_inv);
  if(!status){
    return; // Add a boolean as a member of Inverse_Non_Exact for well-posedness
  }
  }
  
  this->inverse_computed = true;
  this->E_inv.makeCompressed();

  return;
};
