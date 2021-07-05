#include "Wald.h"
#include <cmath>

template<typename InputHandler, MatrixType> 
void Wald_Base<InputHandler, MatrixType>::compute_sigma_hat_sq(void){
  if(is_S_computed==true){
    VectorXr eps_hat = (*(this->inf_car.getZp())) - (*(this->inf_car.getZ_hatp()));
    Real SS_res = eps_hat.squaredNorm();
  
    UInt n = this->inf_car.getN_obs();
    UInt q = this->inf_car.getq();
    tr_S = this->S.trace();
    sigma_hat_sq = SS_res/(n - (q + tr_S));
  }
  
  
  return; 
};

template<typename InputHandler, MatrixType> 
void Wald_Base<InputHandler, MatrixType>::compute_V(){
  // resize the variance-covariance matrix
  UInt q = this->inf_car.getq();
  V.resize(q,q);
  
  const MatrixXr * W = this->inf_car.getWp();
  const Eigen::PartialPivLU<MatrixXr> * WtW_decp = this->inf_car.getWtW_decp();
  
  this->compute_sigma_hat_sq();
  V = this->sigma_hat_sq*((*WtW_decp).solve(MatrixXr::Identity(q,q)) + (*WtW_decp).solve(W->transpose()*S*S.transpose()*(*W)*(*WtW_decp).solve(MatrixXr::Identity(q,q))));
  is_V_computed = true;
  
  return;
};

template<typename InputHandler, MatrixType> 
VectorXr Wald_Base<InputHandler, MatrixType>::compute_pvalue(void){
  // declare the vector that will store the p-values
  VectorXr result;
  
  // simultaneous test
  if(this->inf_car.getInfData()->get_test_type()[this->pos_impl] == "simultaneous"){
    // get the matrix of coefficients for the test
    MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
    // get the value of the parameters under the null hypothesis
    VectorXr beta_0 = this->inf_car.getInfData()->get_beta_0();
    // get the estimates of the parameters
    VectorXr beta_hat = (*(this->inf_car.getBeta_hatp()))(0);
    // compute the difference
    VectorXr diff = C*beta_hat - beta_0; 
    
    // compute the variance-covariance matrix if needed
    if(!is_S_computed){
      compute_S();
    }
    if(!is_V_computed){
      compute_V();
    }
    
    MatrixXr Sigma = C*V*C.transpose();
    // compute the LU factorization of Sigma
    Eigen::PartialPivLU<MatrixXr> Sigma_dec;
    Sigma_dec.compute(Sigma);
    
    // compute the test statistic
    Real stat = diff.adjoint() * Sigma_dec.solve(diff);
    
    //chi_squared distribution(C.rows());
    
    // compute the p-value
    //Real pval = cdf(complement(distribution, stat));
    
    result.resize(C.rows()); // Allocate more space so that R receives a well defined object (different implementations may require higher number of pvalues)
    //result(0) = pval;
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
    VectorXr beta_hat = (*(this->inf_car.getBeta_hatp()))(0);
    
    // compute S and V if needed
    if(!is_S_computed){
      compute_S();
    }
    if(!is_V_computed){
      compute_V();
    }
      
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

template<typename InputHandler, MatrixType> 
MatrixXv Wald_Base<InputHandler, MatrixType>::compute_CI(void){
  
  // get the matrix of coefficients
  MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
  
  // get the estimates of the parameters
  VectorXr beta_hat = (*(this->inf_car.getBeta_hatp()))(0);
  
  // declare the matrix that will store the p-values
  UInt p=C.rows();
  MatrixXv result;
  result.resize(1,p);

  // Extract the quantile needed for computing the upper and lower bounds
  Real quant = this->inf_car.getInfData()->get_inference_quantile()[this->pos_impl];
  
  // compute S and V if needed
  if(!is_S_computed){
    compute_S();
  }
  if(!is_V_computed){
    compute_V();
  }
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

template<typename InputHandler, MatrixType>
Real Wald_Base<InputHandler, MatrixType>::compute_GCV_from_inference(void) const {
  UInt n_obs =this->inf_car.getN_obs();
  UInt q = this->inf_car.getq();
  if(this->is_S_computed==true){
    return sigma_hat_sq * n_obs /(n_obs - q - tr_S);
  } else{
    return -1; // S has not been computed, returning default value
  }
};

template<typename InputHandler, MatrixType> 
void Wald_Exact<InputHandler, MatrixType>::compute_S(void){
  // compute the inverse of the system matrix M by reconstructing the Woodbury decomposition
  this->inverter->Compute_Inv();

  MatrixXr M_inv;
  M_inv.resize(this->inverter->getInv()->rows(), this->inverter->getInv()->cols());

  const MatrixXr * E_inv = this->inverter->getInv();
  const MatrixXr * U = this->inf_car.getUp();
  const MatrixXr * V = this->inf_car.getVp();
  const Eigen::PartialPivLU<MatrixXr> * G_decp = this->inf_car.getG_decp();
  
  M_inv = *E_inv - (*E_inv)*(*U)*((*G_decp).solve((*V)*(*E_inv)));
  
  UInt n_obs = this->inf_car.getN_obs();
  UInt n_nodes = this->inf_car.getN_nodes();
  S.resize(n_obs, n_obs);
  const SpMat * Psi = this->inf_car.getPsip();
  const SpMat * Psi_t = this->inf_car.getPsi_tp();
  UInt q = this->inf_car.getq(); 
  MatrixXr Q = MatrixXr::Identity(q, q) - *(this->inf_car.getHp()); 
  
  S = (*Psi)*M_inv.block(0,0, n_nodes, n_nodes)*((*Psi_t)*Q);
  
  is_S_computed = true;
  
  return; 
};

template<typename InputHandler, MatrixType> 
void Wald_Non_Exact<InputHandler, MatrixType>::compute_S(void){
  // compute the inverse of the system matrix M by reconstructing the Woodbury decomposition
  this->inverter->Compute_Inv();

  MatrixXr M_tilde_inv;
  M_tilde_inv.resize(this->inverter->getInv()->rows(), this->inverter->getInv()->cols());

  UInt n_obs = this->inf_car.getN_obs();
  UInt n_nodes = this->inf_car.getN_nodes();  
  UInt q = this->inf_car.getq(); 

  const MatrixType * E_tilde_inv = this->inverter->getInv();
  const MatrixXr * U_tilde = *(this->inf_car.getUp()).topRows(n_nodes);
  const MatrixXr * V_tilde = *(this->inf_car.getVp()).leftCols(n_nodes);
  const MatrixXr C_tilde = -*(this->inf_car.getWtW_decp()).solve(MatrixXr::Identity(q, q));
  const MatrixXr G_tilde = C_tilde + (*V_tilde)*(*E_tilde_inv)*(*U_tilde);
  const Eigen::PartialPivLU<MatrixXr> G_tilde_decp; 
  G_tilde_decp.compute(G_tilde);
  
  M_tilde_inv = *E_tilde_inv - (*E_tilde_inv)*(*U_tilde)*((*G_tilde_decp).solve((*V_tilde)*(*E_tilde_inv)));
  
  S.resize(n_obs, n_obs);
  const SpMat * Psi = this->inf_car.getPsip();
  const SpMat * Psi_t = this->inf_car.getPsi_tp();
  MatrixXr Q = MatrixXr::Identity(q, q) - *(this->inf_car.getHp()); 
  
  S = (*Psi)*M_tilde_inv*((*Psi_t)*Q);
  
  is_S_computed = true;
  
  return; 
};


template<typename InputHandler, MatrixType>
void Wald_Base<InputHandler, MatrixType>::print_for_debug(void) const {
  
  Rprintf("S computed: %d \n", is_S_computed); 
  
  if(is_S_computed==true){
    Rprintf("Matrix Smoothing S is (only some samples): \n");
    for (UInt i=0; i<10; i++){
      Rprintf( "S( %d, %d):  %f \n", 10*i, 20*i, S(10*i,20*i));
    }
  }
  
  Rprintf("V computed: %d \n" , is_V_computed); 
  
  if(is_V_computed==true){
    Rprintf( "Matrix variance V is: \n");
    for(UInt i=0; i < V.rows(); ++i){
      for(UInt j=0; j < V.cols(); ++j){
	Rprintf(" %f",V(i,j));
      }
    }
    
  } 
  
  return;
};


