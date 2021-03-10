#include "../Include/Inverter.h"

void Inverse_Exact::Compute_Inv(const Eigen::SparseLU<SpMat> * E_decp, const SpMat * Ep){
  E_inv=E_decp->solve(MatrixXr::Identity(Ep->rows(),Ep->cols())); //Solve directly the system for an identity matrix
  return;
}
