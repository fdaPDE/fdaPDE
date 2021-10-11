#include "Inverter.h"
#include "../../Global_Utilities/Include/FSPAI_Wrapper.h"

template<typename InputHandler>
void Inverse_Non_Exact<InputHandler>::pre_Inverse(void){
  // get the mass matrix from the inference carrier
  SpMat R0 = *(inf_car.getR0p());
  R0.makeCompressed();

  // call the FSPAI to get a sparse approximation of R0_inv_tilde
  status_R0_inv_tilde = FSPAI_Wrapper(R0, this->R0_inv_tilde, inf_car.getInfData()->get_tol_Fspai());
  if(!status_R0_inv_tilde){
    return;
  }
  this->R0_inv_tilde.makeCompressed();
 
  UInt n_nodes = inf_car.getN_nodes();
  
  // compute the sparse matrix E_tilde
  if(!(inf_car.getRegData()->isSpaceTime())){ // Space only penalization
    this->E_tilde = (inf_car.getEp()->block(0,0, n_nodes, n_nodes)) + inf_car.getlambda_S()*(inf_car.getR1p()->transpose())*(this->R0_inv_tilde)*(*inf_car.getR1p());
  }else{
    if(inf_car.getRegData()->getFlagParabolic()){// Parabolic Case
      this->E_tilde = (inf_car.getEp()->block(0,0, n_nodes, n_nodes)) + inf_car.getlambda_S()*(*inf_car.getR1p() + inf_car.getlambda_T()* (*inf_car.getLR0kp())).transpose()*(this->R0_inv_tilde)*(*inf_car.getR1p() + inf_car.getlambda_T()* (*inf_car.getLR0kp()));
    }else{
      this->E_tilde = (inf_car.getEp()->block(0,0, n_nodes, n_nodes)) + inf_car.getlambda_S()*(inf_car.getR1p()->transpose())*(this->R0_inv_tilde)*(*inf_car.getR1p());
    }
  }
  this->E_tilde.makeCompressed();

  return;
};


template<typename InputHandler>
void Inverse_Non_Exact<InputHandler>::Compute_Inv(void){
  if(!this->inverse_computed){
    // compute R0_tilde_inv and E_tilde
    this->pre_Inverse();

    if(!status_R0_inv_tilde){
      status_E_tilde_inv=false;
      return;
    }
    // call the FSPAI to get a sparse approximation of E_inv
    status_E_tilde_inv = FSPAI_Wrapper(this->E_tilde, this->E_inv, inf_car.getInfData()->get_tol_Fspai());
  
    if(!status_E_tilde_inv){
      return;
    }
    this->inverse_computed = true;
  }
  
  this->E_inv.makeCompressed();

  return;
};
