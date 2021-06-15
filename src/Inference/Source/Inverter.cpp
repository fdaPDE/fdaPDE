#include "../Include/Inverter.h"
#include "../../Global_Utilities/Include/FSPAI_Wrapper.h"


void Inverse_Base::print_for_debug(void) const {
  
  Rprintf( "Inverse computed: %d \n ", inverse_computed);
  if(inverse_computed){
    Rprintf( "Inverse E_inv (only some samples): \n");
    for (UInt i=0; i<10; i++){
      Rprintf( "E_inv( %d, %d):  %f \n", 10*i, 20*i, E_inv(10*i,20*i));
    } 
  }
  return;
};

void Inverse_Exact::Compute_Inv(void){
  if(!inverse_computed){
    this->E_inv=this->E_decp->solve(MatrixXr::Identity(Ep->rows(),Ep->cols())); //Solve directly the system for an identity matrix
    this->inverse_computed=true;
  }
  
  return;
};


void Inverse_Non_Exact::pre_Inverse(void){
  SpMat R0 = *(inf_car.getR0p());
  bool status = FSPAI_wrapper(R0, this->R0_inv_tilde);
  if(!status){
    return; // Add a boolean as a member of Inverse_Non_Exact for well-posedness
  }
  this->E_tilde = (*inf_car.getPsi_tp())*(*inf_car.getPsip()) + inf_car.getlambda()*(inf_car.getR1p()->transpose())*(this->R0_inv_tilde)*(*inf_car.getR1p());
};
