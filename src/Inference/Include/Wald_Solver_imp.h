#include "Wald_Solver.h"
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <cmath>

using boost::math::chi_squared;
using boost::math::normal;
using boost::math; 

void template<InputHandler> Wald_Solver<InputHandler>::compute_S(void){
  // call the inverter to compute the inverse of the sparse matrix E of the Woodbury decomposition
  inverter.compute_inverse();
  // compute the inverse of the system matrix M by reconstructing the Woodbury decomposition
  MatrixXr M_inv;
  M_inv.resize(inverter.getE_inv()->rows(), inverter.getE_inv()->cols());
  const MatrixXr * E_inv = inverter.getE_inv();
  const MatrixXr * U = inf_car.getUp();
  const MatrixXr * V = inf_car.getVp();
  const Eigen::PartialPivLU<MatrixXr> * G_decp = inf_car.getG_decp();

  M_inv = *E_inv - (*E_inv)*(*U)*((*G_decp).solve((*V)*(*E_inv)));

  Uint n_obs = inf_car.getN_obs();
  Uint n_nodes = inf_car.getN_nodes();
  S.resize(n_obs, n_obs);
  const Spmat * Psi = inf_car.getPsip();
  const Spmat * Psi_t = inf_car.getPsi_tp();
  Uint p = inf_car.getp(); 
  MatrixXr Q = MatrixXr::Identity(p, p) - *(inf_car.getHp()); 

  S = (*Psi)*M_inv.system.block<n_nodes, n_nodes>(0,0)*(-(*Psi_t)*Q);

  S_t.resize(n_obs, n_obs);
  S_t = this->S.transpose();
  is_S_computed = true;
  return; 
};

void template<InputHandler> Wald_Solver<InputHandler>::compute_V(){
  // get the variance of the residuals estimators from the inference carrier
  const Real var_res = inf_car.getVar_res();
  // resize the variance-covariance matrix
  Uint p = inf_car.getp();
  V.resize(p,p);

  const MatrixXr * W = inf_car.getWp();
  const MatrixXr W_t = W->transpose();
  const Eigen::PartialPivLU<MatrixXr> * WtW_decp = inf_car.getWtW_decp();
  
  V = var_res*((*WtW_decp).solve(MatrixXr::Identity(p,p)) + (*WtW_decp).solve(W_t*S*S_t*(*W)*(*WtW_decp).solve(MatrixXr::Identity(p,p))));
  is_V_computed = true;

  return;
};

VectorXr template<InputHandler> Wald_Solver<InputHandler>::compute_pvalue(void){
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
    VectorXr beta_hat = *(inf_car.getBeta_hatp());
    // compute the difference
    VectorXr diff = C*beta_hat - beta_0; 
    MatrixXr diff_t = diff.transpose();

    // compute the variance-covariance matrix if needed
    if(!is_S_computed){
      compute_S();
    }
    if(!is_V_computed){
      compute_V();
    }

    MatrixXr Sigma = C*V*C_t;
    // compute the LU factorization of Sigma
    Eigen::PartialPivLU<MatrixXr> Sigma_dec;
    Sigma_dec.compute(Sigma);

    // compute the test statistic
    Real stat = diff_t * Sigma_dec.solve(diff);

    chi_squared distribution(C.rows());

    // compute the p-value
    Real pval = cdf(complement(distribution), stat);

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
    VectorXr beta_hat = *(inf_car.getBeta_hatp());

    // compute S and V if needed
      if(!is_S_computed){
        compute_S();
      }
      if(!is_V_computed){
        compute_V();
      }

    // for each row of C matrix
    for(Uint i=0; i<q; ++i){
      MatrixXr row = C.row(i);
      Real diff = row*beta_hat - beta_0(i);
      Real sigma = row*V*col;
      // compute the test statistic
      Real stat = diff/sqrt(sigma);
      normal distribution(0,1);
      // compute the pvalue
      Real pval = 2*cdf(complement(distr), fabs(stat));
      result(i) = pval; 	
    }

    return result;
  }

};

MatrixXv template<InputHandler> Wald_Solver<InputHandler>::compute_inference_output(void){
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
    Uint q = inf_car.getInfData()->get_coeff_inference().rows();
    result.resize(1, q+1);
    result(0) = compute_pvalue();
    result.rightCols(q) = compute_CI();
    return result;
  }
  
  
};
