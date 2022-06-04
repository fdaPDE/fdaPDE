#include "Wald.h"
#include <cmath>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>

template<typename InputHandler, typename MatrixType> 
void Wald_Base<InputHandler, MatrixType>::compute_sigma_hat_sq(void){
  
  VectorXr eps_hat = (*(this->inf_car.getZp())) - (this->inf_car.getZ_hat());
  Real SS_res = eps_hat.squaredNorm();
  UInt n = this->inf_car.getN_obs();
  
  if(this->inf_car.getRegData()->getCovariates()->rows()!=0){ //case with covariates
    //check if S has been computed 
		   if(is_S_computed==false){
		     this->compute_S();
		   }

    UInt q = this->inf_car.getq();
    tr_S = this->S.trace();
    sigma_hat_sq = SS_res/(n - (q + tr_S));
  }

  else{//no covariates case
    //check if B has been computed 
		   if(is_B_computed==false){
		     this->compute_B();
		   }

    Real tr_B = this->B.trace();
    sigma_hat_sq = SS_res/(n - tr_B);
  }
  
  is_sigma_hat_sq_computed = true;
  
  return; 
};

template<typename InputHandler, typename MatrixType> 
void Wald_Base<InputHandler, MatrixType>::compute_V(void){
  //check if S has been computed 
		 if(is_S_computed==false){
		   this->compute_S();
		 }
  
  // resize the variance-covariance matrix
  UInt q = this->inf_car.getq();
  V.resize(q,q);
  
  const MatrixXr * W = this->inf_car.getWp();
  const Eigen::PartialPivLU<MatrixXr> * WtW_decp = this->inf_car.getWtW_decp();
  
  if(is_sigma_hat_sq_computed==false){
    this->compute_sigma_hat_sq();
  }
  V = this->sigma_hat_sq*((*WtW_decp).solve(MatrixXr::Identity(q,q)) + (*WtW_decp).solve(W->transpose()*S*S.transpose()*(*W)*(*WtW_decp).solve(MatrixXr::Identity(q,q))));
  is_V_computed = true;
  
  return;
};

template<typename InputHandler, typename MatrixType> 
void Wald_Base<InputHandler, MatrixType>::compute_V_f(void){
  // make sure that B and Partial_S has been computed
  if(this->inf_car.getRegData()->getCovariates()->rows()==0){ // case with no covariates
    if(!is_B_computed){
      compute_B();
    }
    if(!is_B_computed){
      // something went wrong with FSPAI inversion 
	this->is_V_f_computed = false;
      return; 
    }
  }
  else{ // case with covariates
    if(!is_S_computed){
      compute_S();
    }
    if(!is_S_computed){
      // something went wrong with FSPAI inversion 
	this->is_V_f_computed = false;
      return; 
    }
  }
  // make sure sigma_hat_sq has been computed
  if(is_sigma_hat_sq_computed==false){
    compute_sigma_hat_sq();
  }

  // compute variance-covariance matrix of f_hat
  this->V_f = this->sigma_hat_sq * (this->Partial_S) * (this->Partial_S.transpose());
  this->is_V_f_computed = true;

  return;
   
};

template<typename InputHandler, typename MatrixType> 
VectorXr Wald_Base<InputHandler, MatrixType>::compute_beta_pvalue(void){

  VectorXr result;

  // compute the variance-covariance matrix if needed
  if(!is_S_computed){
    compute_S();
    if(!is_S_computed){     // Failed computation of E_tilde_inv/E_inv, returning, unfeasible p_values
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
    VectorXr beta_hat = (*(this->inf_car.getBeta_hatp()));
    // compute the difference
    VectorXr diff = C*beta_hat - beta_0; 
    
    MatrixXr Sigma = C*V*C.transpose();
    // compute the LU factorization of Sigma
    Eigen::PartialPivLU<MatrixXr> Sigma_dec;
    Sigma_dec.compute(Sigma);
    
    // compute the test statistic
    Real stat = diff.adjoint() * Sigma_dec.solve(diff);

    // Allocate more space so that R receives a well defined object (different implementations may require higher number of pvalues)
    result.resize(C.rows());
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
    VectorXr beta_hat = (*(this->inf_car.getBeta_hatp()));
      
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
Real Wald_Base<InputHandler, MatrixType>::compute_f_pvalue(void){
  
  Real result;
  // check that variance-covariance matrix of f_hat has been computed
  if(!is_V_f_computed){
    this->compute_V_f();
  }

  // check that FSPAI inversion went well (if non-exact inference was required)
  if(!is_V_f_computed){
    Rprintf("error: failed FSPAI inversion in p_values computation, discarding inference"); 
    result = 10e20;
    return result; 
  }
  
  // get the estimator f_hat 
       UInt nnodes = this->inf_car.getN_nodes();
  const VectorXr f_hat = this->inf_car.getSolutionp()->topRows(nnodes);

  // get f0 
       const VectorXr f_0 = this->inf_car.getInfData()->get_f_0();

  // get Psi_loc
  const SpMat Psi_loc = this->inf_car.getPsi_loc();
  
  // compute local estimator 
       VectorXr f_loc_hat = Psi_loc * f_hat; 

  // derive the variance-covariance matrix of f_loc_hat
  MatrixXr V_f_loc = Psi_loc * this->V_f * Psi_loc.transpose();

  // compute eigenvalue decomposition of V_f_loc
  Eigen::SelfAdjointEigenSolver<MatrixXr> V_f_loc_eig(V_f_loc);
  MatrixXr eig_values = V_f_loc_eig.eigenvalues().asDiagonal();

  UInt k = 0; 
  UInt max_it = eig_values.cols();
    
  Real threshold = 0.01; 
  bool stop = false;

  while(!stop){
    if(k >= max_it || eig_values(k,k) > threshold)
      stop = true;
    ++k;
  }

  // check that at least one eigenvalue > threshold is found (this may happen with FSPAI implementation) 
       if(max_it - k + 1 == 0){
	 Rprintf("error: cannot invert variance-covariance matrix in Wald-type inference for f, returning NA");
	 result = 10e20;
	 return result;  
       }

  MatrixXr eig_values_red = eig_values.bottomRightCorner(max_it - k + 1, max_it - k + 1);
  MatrixXr eig_vector_red = V_f_loc_eig.eigenvectors().rightCols(max_it - k + 1);

  // Rank-r pseudo inverse
  MatrixXr eig_inv = eig_values_red; 
  eig_inv.diagonal() = eig_values_red.diagonal().array().inverse();
  MatrixXr V_f_loc_red_inv = eig_vector_red * eig_inv * eig_vector_red.transpose();
  
  // compute the test statistics
  Real result_red = (f_loc_hat - f_0).transpose() * V_f_loc_red_inv * (f_loc_hat - f_0);
  
  // generate the distribution
  boost::math::chi_squared dist(max_it - k + 1);
  
  // compute the pvalue
  result = boost::math::cdf(complement(dist,result_red));
  
  return result; 	
};

template<typename InputHandler, typename MatrixType> 
MatrixXv Wald_Base<InputHandler, MatrixType>::compute_beta_CI(void){
  
  // compute S and V if needed
  if(!is_S_computed){
    compute_S();
    if(!is_S_computed){   // Failed inversion of E_tilde_inv/E_inv, returning unfeasible CI
      Rprintf("error: failed FSPAI inversion in confidence intervals computation, discarding inference");
      MatrixXv result;
      MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
      result.resize(1,C.rows());
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
  VectorXr beta_hat = (*(this->inf_car.getBeta_hatp()));
  
  // declare the matrix that will store the confidence intervals
  UInt p=C.rows();
  MatrixXv result;
  result.resize(1,p);

  // Extract the quantile needed for computing the upper and lower bounds
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
MatrixXv Wald_Base<InputHandler, MatrixType>::compute_f_CI(void){
  // declare the object that will store the confidence intervals
  MatrixXv result;
  UInt n_loc = this->inf_car.getN_loc();
  result.resize(1, n_loc);
  
  // make sure that variance-covariance matrix of f_hat has been computed
  if(!is_V_f_computed){
    this->compute_V_f();
  }

  // check that FSPAI inversion went well (if non-exact inference was required)
  if(!is_V_f_computed){
    Rprintf("error: failed FSPAI inversion in p_values computation, discarding inference"); 
    for(UInt i=0; i < n_loc; ++i){
      result(i).resize(3); 
      result(i)(1) = 10e20; 
      result(i)(2) = 10e20;
      result(i)(3) = 10e20;
    }
    return result; 
  }
  
  // get the estimator f_hat 
       UInt nnodes = this->inf_car.getN_nodes();
  const VectorXr f_hat = this->inf_car.getSolutionp()->topRows(nnodes);

  // get Psi_loc
  const SpMat Psi_loc = this->inf_car.getPsi_loc();
  
  // compute local estimator 
       VectorXr f_loc_hat = Psi_loc * f_hat; 

  // derive the variance-covariance matrix of f_loc_hat
  MatrixXr V_f_loc = Psi_loc * this->V_f * Psi_loc.transpose();

  // compute the quantile 
       Real alpha = this->inf_car.getInfData()->get_inference_alpha()[this->pos_impl];
  Real quant;
  
  //ont-at-the-time confidence intervals is the only possible implementation here
  // generate the distribution
  boost::math::normal dist; // default is a standard normal
  quant = boost::math::quantile(complement(dist, alpha/2));
  

  for(UInt i=0; i<n_loc; ++i){
    result(i).resize(3);
    // central element
    result(i)(1)=f_loc_hat(i);
    
    // compute the half range of the interval
    Real sd = std::sqrt(V_f_loc(i,i));
    Real half_range=sd*quant;
    
    // compute the limits of the interval
    result(i)(0) = result(i)(1) - half_range; 
    result(i)(2) = result(i)(1) + half_range; 	
  }
  
  return result;

};

template<typename InputHandler, typename MatrixType>
Real Wald_Base<InputHandler, MatrixType>::compute_GCV_from_inference(void) const {
  UInt n_obs =this->inf_car.getN_obs();
  UInt q = this->inf_car.getq();
  if(this->is_S_computed==true){
    return sigma_hat_sq * n_obs /(n_obs - q - tr_S);
  } else{
    return -1; // S has not been computed, returning default value
  }
};


template<typename InputHandler, typename MatrixType>
VectorXr Wald_Base<InputHandler, MatrixType>::compute_f_var(void){
  UInt n_obs = this->inf_car.getN_obs();
  VectorXr result; 
  result.resize(n_obs);

  if(is_S_computed==false){
    this->compute_S();
  }
  
  if(is_sigma_hat_sq_computed==false){
    compute_sigma_hat_sq();
  }

  const SpMat * Psi = this->inf_car.getPsip();
  const SpMat * Psi_t = this->inf_car.getPsi_tp();

  for(long int i=0; i<n_obs; ++i){
    result(i) = sigma_hat_sq*((*Psi).row(i))*Partial_S*Partial_S.transpose()*((*Psi_t).col(i));
  }
  
  return result; 
};


template<typename InputHandler, typename MatrixType> 
void Wald_Exact<InputHandler, MatrixType>::compute_S(void){
  // compute the inverse of the system matrix M by reconstructing the Woodbury decomposition
  this->inverter->Compute_Inv();

  MatrixXr M_inv;
  M_inv.resize(this->inverter->getInv()->rows(), this->inverter->getInv()->cols());

  const MatrixType * E_inv = this->inverter->getInv();
  const VectorXr * A = this->inf_car.getAp();
  const MatrixXr * U = this->inf_car.getUp();
  const MatrixXr * V = this->inf_car.getVp();
  const Eigen::PartialPivLU<MatrixXr> * G_decp = this->inf_car.getG_decp();
  
  M_inv = *E_inv - (*E_inv)*(*U)*((*G_decp).solve((*V)*(*E_inv)));
  
  UInt n_obs = this->inf_car.getN_obs();
  UInt n_nodes = this->inf_car.getN_nodes();
  this->S.resize(n_obs, n_obs);
  const SpMat * Psi = this->inf_car.getPsip();
  const SpMat * Psi_t = this->inf_car.getPsi_tp();
  MatrixXr Q = MatrixXr::Identity(n_obs, n_obs) - *(this->inf_car.getHp()); 

  if(this->inf_car.getInfData()->get_f_var() || (this->inf_car.getInfData()->get_component_type())[this->pos_impl]!="parametric"){
    this->Partial_S.resize(n_nodes, n_obs);
    
    if(this->inf_car.getRegData()->getNumberOfRegions()>0){
      this->Partial_S = M_inv.block(0,0, n_nodes, n_nodes)*(*Psi_t)*(A->asDiagonal())*Q;
    }else{
      this->Partial_S = M_inv.block(0,0, n_nodes, n_nodes)*(*Psi_t)*Q;
    }
  }
  else{
    this->Partial_S.resize(1,1);
    this->Partial_S(0) = 0;
  }

  if(this->inf_car.getRegData()->getNumberOfRegions()>0){
    this->S = (*Psi)*M_inv.block(0,0, n_nodes, n_nodes)*((*Psi_t)*(A->asDiagonal())*Q);
  }else{
    this->S = (*Psi)*M_inv.block(0,0, n_nodes, n_nodes)*((*Psi_t)*Q);
  }
  
  this->is_S_computed = true;
  
  return; 
};

template<typename InputHandler, typename MatrixType> 
void Wald_Exact<InputHandler, MatrixType>::compute_B(void){
  this->inverter->Compute_Inv();
  // extract the inverse of E
  const MatrixType * E_inv = this->inverter->getInv();
  
  UInt n_obs = this->inf_car.getN_obs();
  UInt n_nodes = this->inf_car.getN_nodes();
  const SpMat * Psi = this->inf_car.getPsip();
  const SpMat * Psi_t = this->inf_car.getPsi_tp();
  MatrixXr Q = MatrixXr::Identity(n_obs, n_obs) - *(this->inf_car.getHp()); 
  const VectorXr * A = this->inf_car.getAp();
  
  this->B.resize(n_obs,n_obs);
  this->Partial_S.resize(n_nodes, n_obs);
  
  if(this->inf_car.getRegData()->getNumberOfRegions()>0){
    this->B = (*Psi)*((*E_inv).block(0,0, n_nodes, n_nodes)*(*Psi_t)*(A->asDiagonal()));
    this->Partial_S = (*E_inv).block(0,0, n_nodes, n_nodes)*(*Psi_t)*(A->asDiagonal())*Q;
  }else{
    this->B = (*Psi)*((*E_inv).block(0,0, n_nodes, n_nodes)*(*Psi_t));
    this->Partial_S = (*E_inv).block(0,0, n_nodes, n_nodes)*(*Psi_t)*Q;
  }

  this->is_B_computed = true;

  return; 
};
