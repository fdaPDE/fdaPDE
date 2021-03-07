#include "Inverter.h"

void Inverse_Exact::Compute_inv(const inference_Carrier & inf_car){
  E_inv=inf_car.getE_decp()->solve(MatrixXr::identity(inf_car.getEp()->rows(),(inf_car.getEp()->cols()));
  return;
}
