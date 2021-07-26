#include "Speckman.h"
#include <cmath>
#include <type_traits>

template<typename InputHandler, typename MatrixType> 
void Speckman_Base<InputHandler, MatrixType>::compute_V(){

  if(!is_Lambda2_computed){
    this->compute_Lambda2();
  }
  if(!is_WLW_computed){
    compute_WLW_dec();
  }
  
  // get the residuals needed
  VectorXr eps_hat = (*(this->inf_car.getZp())) - (*(this->inf_car.getZ_hatp()));
  // build squared residuals
  VectorXr Res2=eps_hat.array()*eps_hat.array();
  
  // resize the variance-covariance matrix
  UInt q = this->inf_car.getq();
  V.resize(q,q);
  
  const MatrixXr * W = this->inf_car.getWp();
  const MatrixXr W_t = W->transpose();
  
  MatrixXr diag = Res2.asDiagonal();
  
  V = (WLW_dec).solve((W_t)*Lambda2*Res2.asDiagonal()*Lambda2*(*W)*(WLW_dec).solve(MatrixXr::Identity(q,q)));
  is_V_computed = true;
  
  return;
};

template<typename InputHandler, typename MatrixType>
void Speckman_Base<InputHandler, MatrixType>::compute_WLW_dec(void){
  if(!is_Lambda2_computed){
    this->compute_Lambda2();
  }
  
  // compute the decomposition of W^T*Lambda^2*W
  const MatrixXr * W = this->inf_car.getWp();
  
  WLW_dec.compute(W->transpose()*Lambda2*(*W));
  is_WLW_computed=true;
};

template<typename InputHandler, typename MatrixType> 
VectorXr Speckman_Base<InputHandler, MatrixType>::compute_beta_hat(void){
  if(!is_WLW_computed){
    compute_WLW_dec();
  }
  const MatrixXr * W = this->inf_car.getWp();
  
  VectorXr beta = WLW_dec.solve(W->transpose()*Lambda2*(*(this->inf_car.getZp())));
  
  return beta; 
  
};

template<typename InputHandler, typename MatrixType> 
VectorXr Speckman_Base<InputHandler, MatrixType>::compute_pvalue(void){
  // declare the vector that will store the p-values
  VectorXr result;

  // compute the variance-covariance matrix if needed
  if(!is_Lambda2_computed){
    this->compute_Lambda2();
    if(!is_Lambda2_computed){     // Failed computation of E_tilde_inv/E_inv, returning, unfeasible p_values
      Rprintf("error: failed FSPAI inversion in p_values computation, discarding inference");
      MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
      result.resize(C.rows());
    
      for(UInt k=0;k<C.rows();k++){
	result(k)==10e20;
      }
      return result;
    }
  }
  if(!is_V_computed){
    compute_V();
  }
  
  // simultaneous test
  if(this->inf_car.getInfData()->get_test_type()[this->pos_impl] == "simultaneous"){
    // get the matrix of coefficients for the test
    MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
    // get the value of the parameters under the null hypothesis
    VectorXr beta_0 = this->inf_car.getInfData()->get_beta_0();
    // get the estimates of the parameters
    VectorXr beta_hat = compute_beta_hat();
    // compute the difference
    VectorXr diff = C*beta_hat - beta_0; 
    
    MatrixXr Sigma = C*V*C.transpose();
    
    // compute the LU factorization of Sigma
    Eigen::PartialPivLU<MatrixXr> Sigma_dec;
    Sigma_dec.compute(Sigma);
    
    // compute the test statistic
    Real stat = diff.adjoint() * Sigma_dec.solve(diff);
    
    result.resize(C.rows()); // Allocate more space so that R receives a well defined object (different implementations may require higher number of pvalues)
    result(0) = stat;
    
    for(UInt k=1;k<C.rows();k++){
      result(k)==10e20;
    }
    
    return result;
  }
  
  // one-at-the-time tests
  else{
    // get the matrix of coefficients for the tests
    MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
    Real p = C.rows();
    result.resize(p);
    
    // get the value of the parameters under the null hypothesis
    VectorXr beta_0 = this->inf_car.getInfData()->get_beta_0();
    // get the estimates of the parameters
    VectorXr beta_hat = compute_beta_hat(); 
    
    // for each row of C matrix
    for(UInt i=0; i<p; ++i){
      VectorXr col = C.row(i);
      Real difference = col.adjoint()*beta_hat - beta_0(i);
      Real sigma = col.adjoint()*V*col;
      
      // compute the test statistic
      Real stat = difference/std::sqrt(sigma);
      
      result(i) = stat;	
    }
    
    return result;
  }
  
};

template<typename InputHandler, typename MatrixType> 
MatrixXv Speckman_Base<InputHandler, MatrixType>::compute_CI(void){
  
  // compute Lambda2 and V if needed
  if(!is_Lambda2_computed){
    compute_Lambda2();
    if(!is_Lambda2_computed){
      Rprintf("error: failed FSPAI inversion in confidence intervals computation, discarding inference");
      MatrixXv result;
      MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
      for(UInt i=0; i<C.rows(); ++i){
	result(i).resize(3);

	//Central element
	result(i)(1)=10e20;

	// compute the limits of the interval
	result(i)(0) = 10e20;
	result(i)(2) = 10e20; 	
      }
      return result;
    }
  }
  if(!is_V_computed){
    compute_V();
  }
  // get the matrix of coefficients
  MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
  
  // get the estimates of the parameters
  VectorXr beta_hat = compute_beta_hat();
  
  // declare the matrix that will store the p-values
  UInt p=C.rows();
  MatrixXv result;
  result.resize(1,p);

  // Extract the quantile needed for the computation of upper and lower bounds
  Real quant = this->inf_car.getInfData()->get_inference_quantile()[this->pos_impl];
  
  // for each row of C matrix
  for(UInt i=0; i<p; ++i){
    result(i).resize(3);
    VectorXr col = C.row(i);
    
    //Central element
    result(i)(1)=col.adjoint()*beta_hat;
    
    // compute the standard deviation of the linear combination and half range of the interval
    Real sd_comb = std::sqrt(col.adjoint()*V*col);
    Real half_range=sd_comb*quant;
    
    // compute the limits of the interval
    result(i)(0) = result(i)(1) - half_range; 
    result(i)(2) = result(i)(1) + half_range; 	
  }
  
  return result;
};

template<typename InputHandler, typename MatrixType> 
void Speckman_Exact<InputHandler, MatrixType>::compute_Lambda2(void){
  this->inverter->Compute_Inv();
  // extract the inverse of E
  const MatrixType * E_inv = this->inverter->getInv();
  
  UInt n_obs = this->inf_car.getN_obs();
  UInt n_nodes = this->inf_car.getN_nodes();
  
  const SpMat * Psi = this->inf_car.getPsip();
  const SpMat * Psi_t = this->inf_car.getPsi_tp(); 
  
  this->Lambda2.resize(n_obs, n_obs);
  this->Lambda2 = (MatrixXr::Identity(n_obs,n_obs) - (*Psi)*((*E_inv).block(0,0, n_nodes, n_nodes)*(*Psi_t)))*(MatrixXr::Identity(n_obs,n_obs) - (*Psi)*((*E_inv).block(0,0, n_nodes, n_nodes)*(*Psi_t)));
  this->is_Lambda2_computed = true;
  
  return; 
};

template<typename InputHandler, typename MatrixType> 
void Speckman_Non_Exact<InputHandler, MatrixType>::compute_Lambda2(void){
  this->inverter->Compute_Inv();

  if(this->inverter->get_status_inverse()==false){
    this->is_Lambda2_computed = false;
    return;
  } 

  // extract the inverse of E_tilde
  const MatrixType * E_tilde_inv = this->inverter->getInv();
  
  UInt n_obs = this->inf_car.getN_obs();
  UInt n_nodes = this->inf_car.getN_nodes();
  
  const SpMat * Psi = this->inf_car.getPsip();
  const SpMat * Psi_t = this->inf_car.getPsi_tp(); 
  
  this->Lambda2.resize(n_obs, n_obs);
  SpMat Identity(n_obs, n_obs);
  Identity.setIdentity();
  this->Lambda2 = (Identity - (*Psi)*((*E_tilde_inv)*(*Psi_t)))*(Identity - (*Psi)*((*E_tilde_inv)*(*Psi_t)));
  //if(std::is_same<MatrixType, SpMat>::value==true){
  //  this->Lambda2.makeCompressed();
  //}
  this->is_Lambda2_computed = true;
  
  return; 
};
