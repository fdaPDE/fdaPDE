#include "Inverter.h"
#include "../../Global_Utilities/Include/FSPAI_Wrapper.h"

template<typename InputHandler>
void Inverse_Non_Exact<InputHandler>::pre_Inverse(void){
  SpMat R0 = *(inf_car.getR0p());
  R0.makeCompressed();

  status_R0 = FSPAI_Wrapper(R0, this->R0_inv_tilde);
  if(!status_R0_inv_tilde){
    return;
  }
  this->R0_inv_tilde.makeCompressed();

  this->E_tilde = (*inf_car.getPsi_tp())*(*inf_car.getPsip()) + inf_car.getlambda()*(inf_car.getR1p()->transpose())*(this->R0_inv_tilde)*(*inf_car.getR1p());
  this->E_tilde.makeCompressed();

  return;
};


template<typename InputHandler>
void Inverse_Non_Exact<InputHandler>::Compute_Inv(void){
  if(!this->inverse_computed){
    this->pre_Inverse();

    if(!status_R0_inv_tilde){
      status_E_tilde_inv=false;
      return;
    }
  
    status_E_tilde_inv = FSPAI_Wrapper(this->E_tilde, this->E_inv);
  
    if(!status_E_tilde_inv){
      return;
    }
    this->inverse_computed = true;
  }
  
  this->E_inv.makeCompressed();

  return;
};
