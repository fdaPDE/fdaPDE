#include "Eigen_Sign_Flip.h"
#include <cmath>
#include <random>

template<typename InputHandler, typename MatrixType> 
VectorXr Eigen_Sign_Flip_Base<InputHandler, MatrixType>::compute_pvalue(void){
  
  // extract matrix C  
  MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
  UInt p = C.rows(); 
  
  // declare the vector that will store the p-values
  VectorXr result;
  
  // get the number of permutations
  unsigned long int n_perm=this->inf_car.getInfData()->get_n_perm();
  
  // get the value of the parameters under the null hypothesis
  VectorXr beta_0 = this->inf_car.getInfData()->get_beta_0(); 
  
  // compute Lambda
  if(!is_Lambda_computed){
    compute_Lambda();
  }
  
  // compute eigenvectors and eigenvalues of Lambda
  Eigen::SelfAdjointEigenSolver<MatrixXr> Lambda_dec(Lambda);
  
  // extract covariates matrices
  const MatrixXr * W = this->inf_car.getWp();
  
  // simultaneous test
  if(this->inf_car.getInfData()->get_test_type()[this->pos_impl] == "simultaneous"){
    // Store beta_hat
    VectorXr beta_hat = (*(this->inf_car.getBeta_hatp()))(0);
    
    // Build auxiliary matrix for residuals computation
    MatrixXr C_out = MatrixXr::Zero(C.cols(), C.cols());
    
    for(UInt j = 0; j < C.cols(); ++j){
      UInt count = 0; 
      for(UInt i = 0; i < p; ++i){
	if(C(i,j) == 1)
	  count++;
      }
      if(count == 0)
	C_out(j,j) = 1; 
    }
    
    // compute the partial residuals
    Partial_res_H0 = *(this->inf_car.getZp()) - (*W) * (C.transpose()) * (beta_0) - (*W) * (C_out) * beta_hat;
    
    // compute the vectors needed for the statistic
    MatrixXr TildeX = (C * W->transpose()) * Lambda_dec.eigenvectors()*Lambda_dec.eigenvalues().asDiagonal();   	// W^t * V * D
    VectorXr Tilder = Lambda_dec.eigenvectors().transpose()*Partial_res_H0;   			        		// V^t * partial_res_H0
    
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
      if(stat_perm > stat){ ++count; } 
    }
    
    
    Real pval = count/n_perm;
    
    result.resize(p); // Allocate more space so that R receives a well defined object (different implementations may require higher number of pvalues)
    result(0) = pval;
    for(UInt k=1;k<p;k++){
      result(k)=10e20;
    }
  }
  else{
    
    // one-at-the-time tests

    // Store beta_hat
    VectorXr beta_hat = (*(this->inf_car.getBeta_hatp()))(0);
    
    Partial_res_H0.resize(Lambda.cols(), p);
    for(UInt i=0; i<p; ++i){
      // Build auxiliary vector for residuals computation
      VectorXr other_covariates = MatrixXr::Ones(beta_hat.size(),1)-C.row(i);

      // compute the partial residuals
      Partial_res_H0.col(i) = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (beta_0); // (z-W*beta_hat(non in test)-W*beta_0(in test))
    }
    // compute the vectors needed for the statistic
    MatrixXr TildeX = (C * W->transpose()) * Lambda_dec.eigenvectors()*Lambda_dec.eigenvalues().asDiagonal();   	// W^t * V * D
    MatrixXr Tilder = Lambda_dec.eigenvectors().transpose()*Partial_res_H0;   			        		// V^t * partial_res_H0
    
    // Observed statistic
    MatrixXr stat=TildeX*Tilder;
    MatrixXr stat_perm=stat;
    std::default_random_engine eng;
    std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)
    VectorXr count = VectorXr::Zero(p);
    
    MatrixXr Tilder_perm=Tilder;
    
    for(unsigned long int i=0;i<n_perm;i++){
      for(unsigned long int j=0;j<TildeX.cols();j++){
	UInt flip=2*distr(eng)-1;
	Tilder_perm.row(j)=Tilder.row(j)*flip;
      }
      stat_perm=TildeX*Tilder_perm; // Flipped statistic
      
      for(UInt k=0; k<p; ++k){
	if(fabs(stat_perm(k,k)) > fabs(stat(k,k))){
	  ++count(k);
	} 
      } 
    }
    
    VectorXr pval = count/n_perm;
    
    result.resize(p);
    result = pval;
  } 
  return result;
  
};

template<typename InputHandler, typename MatrixType>
MatrixXv Eigen_Sign_Flip_Base<InputHandler, MatrixType>::compute_CI(void){
  MatrixXv result;

  // Not implemented yet //

  return result;
  
};

template<typename InputHandler, typename MatrixType> 
void Eigen_Sign_Flip_Exact<InputHandler, MatrixType>::compute_Lambda(void){
  this->inverter->Compute_Inv();
  // extract the inverse of E
  const MatrixType * E_inv = this->inverter->getInv();
  
  UInt n_obs = this->inf_car.getN_obs();
  UInt n_nodes = this->inf_car.getN_nodes();
  const SpMat * Psi = this->inf_car.getPsip();
  const SpMat * Psi_t = this->inf_car.getPsi_tp();
  UInt q = this->inf_car.getq(); 
  
  this->Lambda.resize(n_obs,n_obs);
  this->Lambda = (MatrixXr::Identity(n_obs,n_obs) - (*Psi)*((*E_inv).block(0,0, n_nodes, n_nodes)*(*Psi_t)));
  this->is_Lambda_computed = true;
  
  return; 
};

template<typename InputHandler, typename MatrixType> 
void Eigen_Sign_Flip_Non_Exact<InputHandler, MatrixType>::compute_Lambda(void){
  this->inverter->Compute_Inv();
  // extract the inverse of E
  const MatrixType * E_tilde_inv = this->inverter->getInv();
  
  UInt n_obs = this->inf_car.getN_obs();
  UInt n_nodes = this->inf_car.getN_nodes();
  const SpMat * Psi = this->inf_car.getPsip();
  const SpMat * Psi_t = this->inf_car.getPsi_tp();
  UInt q = this->inf_car.getq(); 
  
  this->Lambda.resize(n_obs,n_obs);
  this->Lambda = (MatrixXr::Identity(n_obs,n_obs) - (*Psi)*((*E_tilde_inv)*(*Psi_t)));
  this->is_Lambda_computed = true;
  
  return; 
};

template<typename InputHandler, typename MatrixType>
void Eigen_Sign_Flip_Base<InputHandler, MatrixType>::print_for_debug(void) const {
  return;
};


