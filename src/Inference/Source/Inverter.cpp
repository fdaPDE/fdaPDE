#include "../Include/Inverter.h"

void Inverse_Exact::Compute_Inv(void){
  if(!this->inverse_computed){
    this->E_inv=this->E_decp->solve(MatrixXr::Identity(Ep->rows(),Ep->cols())); //Solve directly the system for an identity matrix
    this->inverse_computed=true;
  }
  
  return;
};
