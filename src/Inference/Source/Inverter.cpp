#include "../Include/Inverter.h"


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

void Inverse_Exact::Compute_Inv(const Eigen::SparseLU<SpMat> * E_decp, const SpMat * Ep){
  if(!inverse_computed){
  E_inv=E_decp->solve(MatrixXr::Identity(Ep->rows(),Ep->cols())); //Solve directly the system for an identity matrix
  inverse_computed=true;
  }
  return;
};

