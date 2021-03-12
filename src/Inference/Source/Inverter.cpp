#include "../Include/Inverter.h"

void Inverse_Base::print_for_debug(void) const {

  std::cout << "Inverse computed: " << inverse_computed << std::endl;
  if(inverse_computed){
    std::cout << "The inverse is (only some samples): \n" << std::endl;
    for (UInt i=0; i<10; i++){
      std::cout << "E_inv(" << 10*i << "," << 20*i << "): " << E_inv(10*i,20*i) << std::endl;
    }
  }
  return;
}

void Inverse_Exact::Compute_Inv(const Eigen::SparseLU<SpMat> * E_decp, const SpMat * Ep){
  E_inv=E_decp->solve(MatrixXr::Identity(Ep->rows(),Ep->cols())); //Solve directly the system for an identity matrix
  return;
}

