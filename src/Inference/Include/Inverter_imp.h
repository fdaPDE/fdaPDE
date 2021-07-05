#include "Inverter.h"
#include "../../Global_Utilities/Include/FSPAI_Wrapper.h"

template<typename MatrixType>
void Inverse_Exact<MatrixType>::Compute_Inv(void){
  if(!inverse_computed){
    this->E_inv=this->E_decp->solve(MatrixXr::Identity(Ep->rows(),Ep->cols())); //Solve directly the system for an identity matrix
    this->inverse_computed=true;
  }
  
  return;
};

template<typename InputHandler, typename MatrixType>
void Inverse_Non_Exact<InputHandler, MatrixType>::pre_Inverse(void){
  SpMat R0 = *(inf_car.getR0p());
  bool status = FSPAI_Wrapper(R0, this->R0_inv_tilde);
  if(!status){
    return; // Add a boolean as a member of Inverse_Non_Exact for well-posedness
  }
  this->E_tilde = (*inf_car.getPsi_tp())*(*inf_car.getPsip()) + inf_car.getlambda()*(inf_car.getR1p()->transpose())*(this->R0_inv_tilde)*(*inf_car.getR1p());
};


template<typename InputHandler, typename MatrixType>
void Inverse_Non_Exact<InputHandler, MatrixType>::Compute_Inv(void){

  bool status = FSPAI_Wrapper(this->E_tilde, this->E_inv);
  if(!status){
    return; // Add a boolean as a member of Inverse_Non_Exact for well-posedness
  }
 
  //this->E_inv.makecompressed();

  return;
};
