#include "Eigen_Sign_Flip.h"
#include <cmath>
#include <random>
#include <type_traits>
#include <vector>

template<typename InputHandler, typename MatrixType>
void Eigen_Sign_Flip_Base<InputHandler, MatrixType>::Compute_speckman_aux(void){
    
  //check if Lambda has been computed
  if(!is_Lambda_computed){
    this->compute_Lambda();
  }

  // extract W and W^T
  const MatrixXr * W = this->inf_car.getWp();
  const MatrixXr W_t = W->transpose();
  
  // Decomposition of [W^t * Lambda^2 * W] 
  Eigen::PartialPivLU<MatrixXr> WLW_dec; 
  WLW_dec.compute((W_t)*((this->Lambda)*(this->Lambda))*(*W));
  
  // get the residuals needed
  VectorXr eps_hat = (*(this->inf_car.getZp())) - (this->inf_car.getZ_hat());
  // build squared residuals
  VectorXr Res2=eps_hat.array()*eps_hat.array();
  
  // resize the variance-covariance matrix
  UInt q = this->inf_car.getq();
  MatrixXr V;
  V.resize(q,q);
  
  
  MatrixXr diag = Res2.asDiagonal();
  
  V = (WLW_dec).solve((W_t)*(Lambda*Lambda)*Res2.asDiagonal()*(Lambda*Lambda)*(*W)*(WLW_dec).solve(MatrixXr::Identity(q,q))); // V = [(W*Lambda2*W)^-1 * Res2 * (W*Lambda2*W)^-1]

  // Extract the quantile needed for the computation of upper and lower bounds
  Real quant = this->inf_car.getInfData()->get_inference_quantile()[this->pos_impl];

  // extract matrix C 
  MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
  UInt p = C.rows(); 
  
  Speckman_aux_ranges.resize(p);
 
  // for each row of C matrix
  for(UInt i=0; i<p; ++i){
    VectorXr col = C.row(i);
    
    // compute the standard deviation of the linear combination and half range of the interval
    Real sd_comb = std::sqrt(col.adjoint()*V*col);
    Real half_range=sd_comb*quant;
    
    // save the half range
    Speckman_aux_ranges(i)=half_range;  	
  }

  this->is_speckman_aux_computed = true;

  return;
}

template<typename InputHandler, typename MatrixType> 
Real Eigen_Sign_Flip_Base<InputHandler, MatrixType>::compute_CI_aux_pvalue(const VectorXr & partial_res_H0_CI, const MatrixXr & TildeX, const  MatrixXr & Tilder_star) const {
  
  // extract matrix C 
  // (in the eigen-sign-flip case we cannot have linear combinations, but we can have at most one 1 for each column of C) 
  UInt p = TildeX.rows(); 

  // declare the vector that will store the p-values
  Real result;
    
  // compute the vectors needed for the statistic 
  MatrixXr Tilder = Tilder_star * partial_res_H0_CI;

  // Observed statistic
  MatrixXr stat=TildeX*Tilder;
  MatrixXr stat_perm=stat;

  // Random sign-flips
  std::default_random_engine eng;
  std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)
  Real count = 0;
    
  MatrixXr Tilder_perm=Tilder;
    
  // get the number of flips
  unsigned long int n_flip=this->inf_car.getInfData()->get_n_Flip();

  for(unsigned long int i=0;i<n_flip;i++){
    for(unsigned long int j=0;j<TildeX.cols();j++){
      UInt flip=2*distr(eng)-1;
      Tilder_perm.row(j)=Tilder.row(j)*flip;
    }
    stat_perm=TildeX*Tilder_perm; // Flipped statistic
      
    for(UInt k=0; k<p; ++k){
      if(fabs(stat_perm(k,k)) > fabs(stat(k,k))){
	++count;
      } 
    } 
  }
    
  Real pval = count/n_flip;

  result = pval;

  return result;
  
};

template<typename InputHandler, typename MatrixType> 
VectorXr Eigen_Sign_Flip_Base<InputHandler, MatrixType>::compute_pvalue(void){
  
  // extract matrix C 
  // (in the eigen-sign-flip case we cannot have linear combinations, but we can have at most one 1 for each column of C) 
  MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
  UInt p = C.rows(); 
  
  // declare the vector that will store the p-values
  VectorXr result;
  
  // get the number of flips
  unsigned long int n_flip=this->inf_car.getInfData()->get_n_Flip();
  
  // get the value of the parameters under the null hypothesis
  VectorXr beta_0 = this->inf_car.getInfData()->get_beta_0(); 
  
  // compute Lambda if necessary
  if(!is_Lambda_computed){
    this->compute_Lambda();
    if(!is_Lambda_computed){
      Rprintf("error: failed FSPAI inversion in p_values computation, discarding inference");
      result.resize(p);
    
      for(UInt k=0;k<p;k++){
	result(k)==10e20;
      }
      return result;
    }
  }
  
  // compute eigenvectors and eigenvalues of Lambda
  Eigen::SelfAdjointEigenSolver<MatrixType> Lambda_dec(Lambda);
  
  // extract covariates matrices
  const MatrixXr * W = this->inf_car.getWp();
  
  // simultaneous test
  if(this->inf_car.getInfData()->get_test_type()[this->pos_impl] == "simultaneous"){
    // Store beta_hat
    VectorXr beta_hat = (*(this->inf_car.getBeta_hatp()));
    
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
    
    //Random sign-flips
    std::default_random_engine eng;
    std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)
    Real count=0;
    VectorXr Pi; // Sign flip matrix preallocation
    Pi.resize(TildeX.cols());
    VectorXr Tilder_perm=Tilder;
    
    for(unsigned long int i=0;i<n_flip;i++){
      for(unsigned long int j=0;j<TildeX.cols();j++){
	UInt flip=2*distr(eng)-1;
	Tilder_perm(j)=Tilder(j)*flip;
      }
      stat_perm=TildeX*Tilder_perm; // Flipped statistic
      if(stat_perm > stat){ ++count; } //Here we use the custom-operator defined in Eigen_Sign_Flip.h 
    }
    
    
    Real pval = count/n_flip;
    
    result.resize(p); // Allocate more space so that R receives a well defined object (different implementations may require higher number of pvalues)
    result(0) = pval;
    for(UInt k=1;k<p;k++){
      result(k)=10e20;
    }
  }
  else{
    
    // one-at-the-time tests

    // Store beta_hat
    VectorXr beta_hat = (*(this->inf_car.getBeta_hatp()));
    
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

    // Random sign-flips
    std::default_random_engine eng;
    std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)VectorXr other_covariates = MatrixXr::Ones(beta_hat.size(),1)-C.row(i);
    VectorXr count = VectorXr::Zero(p);
    
    MatrixXr Tilder_perm=Tilder;
    
    for(unsigned long int i=0;i<n_flip;i++){
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
    
    VectorXr pval = count/n_flip;
    
    result.resize(p);
    result = pval;
  } 
  return result;
  
};

template<typename InputHandler, typename MatrixType>
MatrixXv Eigen_Sign_Flip_Base<InputHandler, MatrixType>::compute_CI(void){

  // compute Lambda
  if(!is_Lambda_computed){
    compute_Lambda();
    if(!is_Lambda_computed){
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
  
  // Store beta_hat
  VectorXr beta_hat = (*(this->inf_car.getBeta_hatp()));

  // get the matrix of coefficients
  MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
  UInt p=C.rows();

  // compute eigenvectors and eigenvalues of Lambda
  Eigen::SelfAdjointEigenSolver<MatrixType> Lambda_dec(Lambda);
  
  // extract covariates matrices
  const MatrixXr * W = this->inf_car.getWp();

  // declare the matrix that will store the intervals
  MatrixXv result;
  result.resize(1,p);

  // compute the initial ranges from speckman's CI (initial guess for )
  if(!is_speckman_aux_computed){
    this->Compute_speckman_aux();
  }
 
  // this vector will store the tolerance for each interval upper/lower limit
  VectorXr ESF_bisection_tolerances = 0.1*Speckman_aux_ranges; // 0.05 of the speckman CI as maximum tolerance

  // define storage structures for bisection algorithms
  VectorXr UU; // Upper limit for Upper bound
  UU.resize(p);
  VectorXr UL; // Lower limit for Upper bound
  UL.resize(p);
  VectorXr LU; // Upper limit for Lower bound
  LU.resize(p);
  VectorXr LL; // Lower limit for Lower bound
  LL.resize(p);


 
  // speckman intervals initialization
  for(UInt i=0; i<p; ++i){
    result(i).resize(3);
    Real half_range = Speckman_aux_ranges(i);
    
    // compute the limits of the interval
    result(i)(0) = beta_hat(i) - half_range;
    LU(i)=result(i)(0)+0.5*half_range;
    LL(i)=result(i)(0)-0.5*half_range;
    result(i)(2) = beta_hat(i) + half_range;
    UU(i)=result(i)(2)+0.5*half_range;
    UL(i)=result(i)(2)-0.5*half_range; 	
  }

  // define booleans used to unserstand which CI need to be computed forward on
  std::vector<bool> converged_up(p,false);
  std::vector<bool> converged_low(p,false);
  bool all_betas_converged=false;

  // matrix that stores p_values of the bounds at actual step
  MatrixXr local_p_values;
  local_p_values.resize(4,p);
  
  // compute the vectors needed for the statistic
  MatrixXr TildeX = (C * W->transpose()) * Lambda_dec.eigenvectors()*Lambda_dec.eigenvalues().asDiagonal();   	// W^t * V * D
  MatrixXr Tilder_star = Lambda_dec.eigenvectors().transpose();   			        		// V^t

  VectorXr Partial_res_H0_CI;
  Partial_res_H0_CI.resize(Lambda.cols());

  // fill the p_values matrix
  for (UInt i=0; i<p; i++){
  
    // Build auxiliary vector for residuals computation
    VectorXr other_covariates = MatrixXr::Ones(beta_hat.size(),1)-C.row(i);

    MatrixXr TildeX_loc= TildeX.row(i);

    // compute the partial residuals and p value
    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (UU(i)); // (z-W*beta_hat(non in test)-W*UU[i](in test))
    local_p_values(0,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

    // compute the partial residuals and p value
    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (UL(i)); // (z-W*beta_hat(non in test)-W*UL[i](in test))
    local_p_values(1,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

    // compute the partial residuals and p value
    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (LU(i)); // (z-W*beta_hat(non in test)-W*LU[i](in test))
    local_p_values(2,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

    // compute the partial residuals and p value
    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (LL(i)); // (z-W*beta_hat(non in test)-W*LL[i](in test))
    local_p_values(3,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

  }

  // extract the CI significance (1-confidence)
  Real alpha=0;

  if(this->inf_car.getInfData()->get_interval_type()[this->pos_impl]=="one-at-the-time"){
    alpha=0.5*(this->inf_car.getInfData()->get_inference_alpha());
  }else{
    alpha=0.5/p * (this->inf_car.getInfData()->get_inference_alpha());
  }
    
  UInt Max_Iter=20;
  UInt Count_Iter=0;
  while(!all_betas_converged & Count_Iter<Max_Iter){
  
    // Compute all p_values (only those needed)
    for (UInt i=0; i<p; i++){

      VectorXr other_covariates = MatrixXr::Ones(beta_hat.size(),1)-C.row(i);

      MatrixXr TildeX_loc= TildeX.row(i);
  
      if(!converged_up[i]){
	if(local_p_values(0,i)>alpha){ // Upper-Upper bound excessively tight

	  UU(i)=UU(i)+0.5*(UU(i)-UL(i));
  
	  // compute the partial residuals
	  Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (UU(i)); // (z-W*beta_hat(non in test)-W*UU[i](in test))
	  local_p_values(0,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);
  
	}else{
 
	  if(local_p_values(1,i)<alpha){ // Upper-Lower bound excessively tight
	    UL(i)=beta_hat(i)+0.5*(UL(i)-beta_hat(i));
  
	    // compute the partial residuals
	    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (UL(i)); // (z-W*beta_hat(non in test)-W*UL[i](in test))
	    local_p_values(1,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

	  }else{//both the Upper bounds are well defined

	    if(UU(i)-UL(i)<ESF_bisection_tolerances(i)){

	      converged_up[i]=true;

	    }else{

	      Real proposal=0.5*(UU(i)+UL(i));
   
	      // compute the partial residuals
	      Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (proposal); // (z-W*beta_hat(non in test)-W*proposal)
	      Real prop_p_val=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

	      if(prop_p_val<=alpha){UU(i)=proposal; local_p_values(0,i)=prop_p_val;}else{UL(i)=proposal;local_p_values(1,i)=prop_p_val;}
	    }
	  }
	}
      }


      if(!converged_low[i]){
	if(local_p_values(2,i)<alpha){ // Lower-Upper bound excessively tight

	  LU(i)=beta_hat(i)-0.5*(beta_hat(i)-LU(i));
  
	  // compute the partial residuals
	  Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (LU(i)); // (z-W*beta_hat(non in test)-W*LU[i](in test))
	  local_p_values(2,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);
  
	}else{
 
	  if(local_p_values(3,i)>alpha){ // Lower-Lower bound excessively tight
	    LL(i)=beta_hat(i)-0.5*(UL(i)-LL(i));
  
	    // compute the partial residuals
	    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (LL(i)); // (z-W*beta_hat(non in test)-W*LL[i](in test))
	    local_p_values(3,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

	  }else{//both the Upper bounds are well defined

	    if(LU(i)-LL(i)<ESF_bisection_tolerances(i)){

	      converged_low[i]=true;

	    }else{

	      Real proposal=0.5*(LU(i)+LL(i));
   
	      // compute the partial residuals
	      Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (proposal); // (z-W*beta_hat(non in test)-W*proposal)
	      Real prop_p_val=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

	      if(prop_p_val<=alpha){LL(i)=proposal; local_p_values(3,i)=prop_p_val;}else{LU(i)=proposal;local_p_values(2,i)=prop_p_val;}
	    }
	  }
	}
      }
    }
    all_betas_converged =true;
    for(UInt j=0; j<p; j++){

      if(!converged_up[j] || !converged_low[j]){
	all_betas_converged=false;
      }
    }

    Count_Iter++;

  }
   
  // for each row of C matrix
  for(UInt i=0; i<p; ++i){
    result(i).resize(3);
    
    if(Count_Iter < Max_Iter){ // No discrepancy between beta_hat(i) and ESF, bisection converged
      // Central element
      result(i)(1)=beta_hat(i);
      
      // Limits of the interval
      result(i)(0) = 0.5*(LU(i)+LL(i));
      result(i)(2) = 0.5*(UU(i)+UL(i)); 
    }else{ // Not converged in time, give a warning in R
      // Central element
      result(i)(1)=10e20;
      
      // Limits of the interval
      result(i)(0) = 10e20;
      result(i)(2) = 10e20; 
    }
  }
  
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
  
  if(this->inverter->get_status_inverse()==false){
    this->is_Lambda_computed=false;
    return;
  }
  
  // extract the inverse of E
  const MatrixType * E_tilde_inv = this->inverter->getInv();
  
  UInt n_obs = this->inf_car.getN_obs();
  UInt n_nodes = this->inf_car.getN_nodes();
  const SpMat * Psi = this->inf_car.getPsip();
  const SpMat * Psi_t = this->inf_car.getPsi_tp();
  UInt q = this->inf_car.getq(); 
  
  this->Lambda.resize(n_obs,n_obs);
  SpMat Identity(n_obs, n_obs);
  Identity.setIdentity();
  this->Lambda = (Identity - (*Psi)*((*E_tilde_inv)*(*Psi_t)));
  this->is_Lambda_computed = true;
  
  return; 
};


