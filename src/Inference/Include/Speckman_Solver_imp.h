#include "Speckman_Solver.h"
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <cmath>

using boost::math::chi_squared;
using boost::math::normal;
using namespace boost::math; 

template<typename InputHandler> 
void Speckman_Solver<InputHandler>::compute_Lambda2(void){
  // call the inverter to compute the inverse of the sparse matrix E of the Woodbury decomposition
  this->inverter.Compute_Inv(this->inf_car.getE_decp(), this->inf_car.getEp());
  // extract the inverse of E
  const MatrixXr * E_inv = this->inverter.getInv(this->inf_car.getE_decp(), this->inf_car.getEp());
  
  UInt n_obs = this->inf_car.getN_obs();
  UInt n_nodes = this->inf_car.getN_nodes();
  B.resize(n_obs, n_obs);
  const SpMat * Psi = this->inf_car.getPsip();
  const SpMat * Psi_t = this->inf_car.getPsi_tp();
  UInt p = this->inf_car.getp(); 
  
  B = (*Psi)*((*E_inv).block(0,0, n_nodes, n_nodes)*(*Psi_t));
  
  Lambda2 = (MatrixXr::Identity(B.rows(),B.cols()) - B)*(MatrixXr::Identity(B.rows(),B.cols()) - B);
  is_Lambda2_computed = true;
  
  return; 
};

template<typename InputHandler> 
void Speckman_Solver<InputHandler>::compute_V(){
  // get the residuals needed
  VectorXr eps_hat = (*(this->inf_car.getZp())) - (*(this->inf_car.getZ_hatp()));
  // build squared residuals
  VectorXr Res2=eps_hat.array()*eps_hat.array();
  
  // resize the variance-covariance matrix
  UInt p = this->inf_car.getp();
  V.resize(p,p);
  
  const MatrixXr * W = this->inf_car.getWp();
  const MatrixXr W_t = W->transpose();
  
  MatrixXr diag = Res2.asDiagonal();
  
  V = (WLW_dec).solve((W_t)*Lambda2*Res2.asDiagonal()*Lambda2*(*W)*(WLW_dec).solve(MatrixXr::Identity(p,p)));
  is_V_computed = true;
  
  return;
};

template<typename InputHandler>
void Speckman_Solver<InputHandler>::compute_WLW_dec(void){
  if(!is_Lambda2_computed){
    compute_Lambda2();
  }
  
  // compute the decomposition of W^T*Lambda^2*W
  const MatrixXr * W = this->inf_car.getWp();
  const MatrixXr W_t = W->transpose();
  
  WLW_dec.compute(W_t*Lambda2*(*W));
  is_WLW_computed=true;
};

template<typename InputHandler> 
VectorXr Speckman_Solver<InputHandler>::compute_beta_hat(void){
  if(!is_WLW_computed){
    compute_WLW_dec();
  }
  const MatrixXr * W = this->inf_car.getWp();
  const MatrixXr W_t = W->transpose();
  
  VectorXr beta = WLW_dec.solve(W_t*Lambda2*(*(this->inf_car.getZp())));
  
  return beta; 
  
};

template<typename InputHandler> 
VectorXr Speckman_Solver<InputHandler>::compute_pvalue(void){
  // declare the vector that will store the p-values
  VectorXr result;
  
  // simultaneous test
  if(this->inf_car.getInfData()->get_test_type() == "simultaneous"){
    // get the matrix of coefficients for the test
    MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
    MatrixXr C_t = C.transpose();
    // get the value of the parameters under the null hypothesis
    VectorXr beta_0 = this->inf_car.getInfData()->get_beta_0();
    // get the estimates of the parameters
    VectorXr beta_hat = compute_beta_hat();
    // compute the difference
    VectorXr diff = C*beta_hat - beta_0; 
    
    // compute the variance-covariance matrix if needed
    if(!is_Lambda2_computed){
      compute_Lambda2();
    }
    if(!is_V_computed){
      compute_V();
    }
    
    MatrixXr Sigma = C*V*C_t;
    // compute the LU factorization of Sigma
    Eigen::PartialPivLU<MatrixXr> Sigma_dec;
    Sigma_dec.compute(Sigma);
    
    // compute the test statistic
    Real stat = diff.adjoint() * Sigma_dec.solve(diff);
    
    chi_squared distribution(C.rows());
    
    // compute the p-value
    Real pval = cdf(complement(distribution, stat));
    
    result.resize(1);
    result(0) = pval;
    return result;
  }
  
  // one-at-the-time tests
  else{
    // get the matrix of coefficients for the tests
    MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
    Real q = C.rows();
    result.resize(q);
    
    // get the value of the parameters under the null hypothesis
    VectorXr beta_0 = this->inf_car.getInfData()->get_beta_0();
    // get the estimates of the parameters
    VectorXr beta_hat = compute_beta_hat(); 
    
    // compute Lambda2 and V if needed
    if(!is_Lambda2_computed){
      compute_Lambda2();
    }
    if(!is_V_computed){
      compute_V();
    }
    
    
    /* // FOR DEBUG TEMPORARY
       this->print_for_debug();
       inverter.print_for_debug(); */
    
    // for each row of C matrix
    for(UInt i=0; i<q; ++i){
      VectorXr col = C.row(i);
      Real difference = col.adjoint()*beta_hat - beta_0(i);
      Real sigma = col.adjoint()*V*col;
      // compute the test statistic
      Real stat = difference/std::sqrt(sigma);
      normal distribution(0,1);
      // compute the pvalue
      Real pval = 2*cdf(complement(distribution, fabs(stat)));
      result(i) = pval; 	
    }
    
    return result;
  }
  
};

template<typename InputHandler> 
MatrixXv Speckman_Solver<InputHandler>::compute_CI(void){
  
  // get the matrix of coefficients
  MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
  
  // get the estimates of the parameters
  VectorXr beta_hat = compute_beta_hat();
  
  // declare the matrix that will store the p-values
  Real alpha=this->inf_car.getInfData()->get_inference_level();
  Real quant=0;
  UInt q=C.rows();
  MatrixXv result;
  result.resize(1,q);
  
  
  // simultaneous confidence interval (overall confidence aplha)
  if(this->inf_car.getInfData()->get_interval_type() == "simultaneous"){
    chi_squared distribution(q);
    quant =std::sqrt(quantile(complement(distribution,alpha)));
  }
  else{
    // one-at-the-time confidence intervals (each interval has confidence alpha)
    if(this->inf_car.getInfData()->get_interval_type() == "one-at-the-time"){
      normal distribution(0,1);
      quant = quantile(complement(distribution,alpha/2));
    }
    // Bonferroni confidence intervals (overall confidence approximately alpha)
    else{
      normal distribution(0,1);
      quant = quantile(complement(distribution,alpha/(2*q)));
    }
  }
  
  // compute Lambda2 and V if needed
  if(!is_Lambda2_computed){
    compute_Lambda2();
  }
  if(!is_V_computed){
    compute_V();
  }
  // for each row of C matrix
  for(UInt i=0; i<q; ++i){
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


template<typename InputHandler>
void Speckman_Solver<InputHandler>::print_for_debug(void) const {
  int aaaa=1;
  return;
};
