#include "Eigen_Sign_Flip_Solver.h"
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include <random>

using boost::math::chi_squared;
using boost::math::normal;
using namespace boost::math; 

template<typename InputHandler> 
void Eigen_Sign_Flip_Solver<InputHandler>::compute_Lambda(void){
  // call the inverter to compute the inverse of the sparse matrix E of the Woodbury decomposition
  inverter.Compute_Inv(inf_car.getE_decp(), inf_car.getEp());
  // extract the inverse of E
  const MatrixXr * E_inv = inverter.getInv(inf_car.getE_decp(), inf_car.getEp());
  
  UInt n_obs = inf_car.getN_obs();
  UInt n_nodes = inf_car.getN_nodes();
  B.resize(n_obs, n_obs);
  const SpMat * Psi = inf_car.getPsip();
  const SpMat * Psi_t = inf_car.getPsi_tp();
  UInt p = inf_car.getp(); 
  
  B = (*Psi)*((*E_inv).block(0,0, n_nodes, n_nodes)*(*Psi_t));
  
  Lambda = (MatrixXr::Identity(B.rows(),B.cols()) - B);
  is_Lambda_computed = true;
  // For Debug Only
  Rprintf("Lambda computed: %d \n", is_Lambda_computed); 
  
  if(is_Lambda_computed==true){
    Rprintf("Matrix Lambda is (only some samples): \n");
    for (UInt i=0; i<10; i++){
      Rprintf( "Lambda( %d, %d):  %f \n", 10*i, 20*i, Lambda(10*i,20*i));
    }
  }
  return; 
};

template<typename InputHandler> 
VectorXr Eigen_Sign_Flip_Solver<InputHandler>::compute_pvalue(void){
  // declare the vector that will store the p-values
  VectorXr result;
  result.resize(1);
  
  // one-at-the-time tests
  unsigned long int n_perm=inf_car.getInfData()->get_n_perm();
  
  // get the value of the parameters under the null hypothesis
  VectorXr beta_0 = inf_car.getInfData()->get_beta_0(); 
  
  // compute Lambda
  if(!is_Lambda_computed){
    compute_Lambda();
  }
  
  // compute eigenvectors and eigenvalues of Lambda
  Eigen::SelfAdjointEigenSolver<MatrixXr> Lambda_dec(Lambda);
  
  // extract covariates matrices
  const MatrixXr * W = inf_car.getWp();
  const MatrixXr W_t = W->transpose();
  
  
  // compute the partial residuals
  Partial_res_H0 = *(inf_car.getZp()) - (*W) * (beta_0);
  
  // compute the vectors needed for the statistic
  MatrixXr TildeX = W_t * Lambda_dec.eigenvectors()*Lambda_dec.eigenvalues().asDiagonal();   	// W^t * V * D
  VectorXr Tilder = Lambda_dec.eigenvectors().transpose()*Partial_res_H0;   			// V^t * partial_res_H0
  
  // Observed statistic
  VectorXr stat=TildeX*Tilder;
  VectorXr stat_perm=stat;
  std::default_random_engine eng;
  std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)
  Real count=0;
  VectorXr Pi; // Sign flip matrix preallocation
  Pi.resize(TildeX.cols());
  VectorXr Tilder_perm=Tilder;
  
  
  
  for(unsigned long int i=0;i<n_perm;i++){
    for(unsigned long int j=0;j<TildeX.cols();j++){
      UInt flip=2*distr(eng)-1;
      Tilder_perm(j)=Tilder(j)*flip;
    }
    stat_perm=TildeX*Tilder_perm; // Flipped statistic
    // If(fabs(stat_perm) > fabs(stat)){ ++count; } ///// >??? fabs????
  }
  
  Real pval = count/n_perm;
  result(0) = pval;
  
  return result;
  
};

template<typename InputHandler>
MatrixXv Eigen_Sign_Flip_Solver<InputHandler>::compute_inference_output(void){
  // declare the result Matrix of vectors to be returned
  MatrixXv result;
  
  // prepare the space for the output (just one-at-the-time p_value case)
  result.resize(1,1);
  result(0) = compute_pvalue(); 
  return result;
  
};


template<typename InputHandler>
void Eigen_Sign_Flip_Solver<InputHandler>::print_for_debug(void) const {
  int aaaa=1;
  return;
};
