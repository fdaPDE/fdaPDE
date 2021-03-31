#include "Speckman_Solver.h"
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <cmath>

using boost::math::chi_squared;
using boost::math::normal;
using namespace boost::math; 

template<typename InputHandler> 
void Speckman_Solver<InputHandler>::compute_Lambda(void){
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
  
  Lambda = MatrixXr::Identity(B.rows(),B.cols()) - B;
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
void Speckman_Solver<InputHandler>::compute_V(){
  // get the residuals needed
  VectorXr eps_hat = (*inf_car.getZp()) - (*inf_car.getZ_hatp());
  // build squared residuals
  VectorXr Res2=eps_hat.array()*eps_hat.array();
  
    Rprintf("Vector Res2 is (only some samples): \n");
    for (UInt i=0; i<10; i++){
      Rprintf( "Res2( %d):  %f \n", 10*i, Res2(10*i));
    }
    
    // resize the variance-covariance matrix
    UInt p = inf_car.getp();
    V.resize(p,p);
    
    const MatrixXr * W = inf_car.getWp();
    const MatrixXr W_t = W->transpose();
    const Eigen::PartialPivLU<MatrixXr> * WtW_decp = inf_car.getWtW_decp();
    
    MatrixXr diag = Res2.asDiagonal();
    
    Rprintf("Vector diag is (only some samples): \n");
    for (UInt i=0; i<10; i++){
      Rprintf( "diag( %d, %d):  %f \n", 10*i, 10*i, diag(10*i, 10*i));
    }
    Rprintf("diag( %d, %d):  %f \n", 10, 20, diag(10, 20));
    
    
  
    V = (*WtW_decp).solve((W_t)*Lambda*Lambda*Res2.asDiagonal()*Lambda*Lambda*(*W)*(*WtW_decp).solve(MatrixXr::Identity(p,p)));
    is_V_computed = true;
    //For debug only
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

template<typename InputHandler> 
VectorXr Speckman_Solver<InputHandler>::compute_beta_hat(void){
 if(!is_Lambda_computed){
   compute_Lambda();
 }
 // compute the decomposition of W^T*Lambda^2*W
 const MatrixXr * W = inf_car.getWp();
 const MatrixXr W_t = W->transpose();
 
 Eigen::PartialPivLU<MatrixXr> WLW_decp; 
 WLW_decp.compute(W_t*Lambda*Lambda*(*W));
 
 VectorXr beta = WLW_decp.solve(W_t*Lambda*Lambda*(*inf_car.getZp()));
 
 return beta; 
 
};

template<typename InputHandler> 
VectorXr Speckman_Solver<InputHandler>::compute_pvalue(void){
  // declare the vector that will store the p-values
  VectorXr result;
  
  // simultaneous test
  if(inf_car.getInfData()->get_test_type() == "simultaneous"){
    // get the matrix of coefficients for the test
    MatrixXr C = inf_car.getInfData()->get_coeff_inference();
    MatrixXr C_t = C.transpose();
    // get the value of the parameters under the null hypothesis
    VectorXr beta_0 = inf_car.getInfData()->get_beta_0();
    // get the estimates of the parameters
    VectorXr beta_hat = compute_beta_hat();
    // compute the difference
    VectorXr diff = C*beta_hat - beta_0; 
    
    // compute the variance-covariance matrix if needed
    if(!is_Lambda_computed){
      compute_Lambda();
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
    MatrixXr C = inf_car.getInfData()->get_coeff_inference();
    Real q = C.rows();
    result.resize(q);
    
    // get the value of the parameters under the null hypothesis
    VectorXr beta_0 = inf_car.getInfData()->get_beta_0();
    // get the estimates of the parameters
    VectorXr beta_hat = compute_beta_hat(); 
    
    Rprintf("Vector beta_hat: \n");
    for (UInt i=0; i<beta_hat.size(); i++){
      Rprintf( "beta_hat( %d):  %f \n", i, beta_hat(i));
    }
    
    // compute Lambda and V if needed
    if(!is_Lambda_computed){
      compute_Lambda();
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
  MatrixXr C = inf_car.getInfData()->get_coeff_inference();
  
  // get the estimates of the parameters
  VectorXr beta_hat = compute_beta_hat();
  
  // declare the matrix that will store the p-values
  Real alpha=inf_car.getInfData()->get_inference_level();
  Real quant=0;
  UInt q=C.rows();
  MatrixXv result;
  result.resize(1,q);
  
  
  // simultaneous confidence interval (overall confidence aplha)
  if(inf_car.getInfData()->get_interval_type() == "simultaneous"){
    chi_squared distribution(q);
    quant =std::sqrt(quantile(complement(distribution,alpha)));
  }
  else{
    // one-at-the-time confidence intervals (each interval has confidence alpha)
    if(inf_car.getInfData()->get_interval_type() == "one-at-the-time"){
      normal distribution(0,1);
      quant = quantile(complement(distribution,alpha/2));
    }
    // Bonferroni confidence intervals (overall confidence approximately alpha)
    else{
      normal distribution(0,1);
      quant = quantile(complement(distribution,alpha/(2*q)));
    }
  }
  
  // compute Lambda and V if needed
  if(!is_Lambda_computed){
    compute_Lambda();
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
MatrixXv Speckman_Solver<InputHandler>::compute_inference_output(void){
  // declare the result Matrix of vectors to be returned
  MatrixXv result;
  
  // get the test_type and interval_type
  std::string test_type = inf_car.getInfData()->get_test_type();
  std::string interval_type = inf_car.getInfData()->get_interval_type();
  
  // if test_type is not defined, only intervals are required
  if(test_type == "not-defined"){
    result = compute_CI();
    return result;
  }
  // if interval_type is not defined, only test is required
  if(interval_type == "not-defined"){
    result.resize(1,1);
    result(0) = compute_pvalue(); 
    return result;
  }
  // else, both are required
  else{
    UInt q = inf_car.getInfData()->get_coeff_inference().rows();
    result.resize(1, q+1);
    result(0) = compute_pvalue();
    result.rightCols(q) = compute_CI();
    return result;
  }
};


template<typename InputHandler>
void Speckman_Solver<InputHandler>::print_for_debug(void) const {
  int aaaa=1;
  return;
};
