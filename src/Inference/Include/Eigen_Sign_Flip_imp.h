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
Real Eigen_Sign_Flip_Base<InputHandler, MatrixType>::compute_CI_aux_pvalue(const VectorXr & partial_res_H0_CI, const MatrixXr & TildeX, const VectorXr & Tilder_hat, const  MatrixXr & Tilder_star, const UInt Direction) const {
  
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

  UInt n_obs = this->inf_car.getN_obs();

  // Estimate the standard error
  VectorXr eps_hat = (*(this->inf_car.getZp())) - (this->inf_car.getZ_hat());
  Real SS_res = eps_hat.squaredNorm();
  Real Sigma_hat = std::sqrt(SS_res/(n_obs-1));

  Real threshold = 5*Sigma_hat; // This threshold is used to determine how many components will not be flipped: we drop those that show large alpha_hat w.r.t. the expected standar error
  UInt N_Eig_Out=0;
    

  // Random sign-flips
  std::random_device rd; 
  std::default_random_engine eng{rd()};
  std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)
  Real count_Up = 0;   // Counter for the number of flipped statistics that are larger the observed statistic
  Real count_Down = 0; // Counter for the number of flipped statistics that are smaller the observed statistic
    
  MatrixXr Tilder_perm=Tilder;
    
  // get the number of flips
  unsigned long int n_flip=this->inf_car.getInfData()->get_n_Flip();

  for(unsigned long int i=0;i<n_flip;i++){
    N_Eig_Out=0;
    for(unsigned long int j=0;j<TildeX.cols();j++){
      UInt flip;
      if((this->inf_car.getInfData()->get_enhanced_inference()=="enhanced") & (N_Eig_Out<n_obs/2) & (fabs(Tilder_hat(j))>threshold)){
	flip=1;
	++N_Eig_Out;
      }else{
	flip=2*distr(eng)-1;
      }
      Tilder_perm(j)=Tilder(j)*flip;
    }
    stat_perm=TildeX*Tilder_perm; // Flipped statistic
    if(is_Unilaterally_Greater(stat_perm,stat)){ ++count_Up;}else{ //Here we use the custom-operator defined in Eigen_Sign_Flip.h
      if(is_Unilaterally_Smaller(stat_perm,stat)){ ++count_Down;} //Here we use the custom-operator defined in Eigen_Sign_Flip.h 
    }
  }
    
  Real pval_up = count_Up/n_flip;     // This is the unilateral p_value in the case of H0 b=b_0 vs H1 b>b_0
  Real pval_down = count_Down/n_flip; // This is the unilateral p_value in the case of H0 b=b_0 vs H1 b<b_0

  if(Direction==1){
    result = pval_up;
  }else{
    result = pval_down;
  }

  return result;
  
};

template<typename InputHandler, typename MatrixType> 
VectorXr Eigen_Sign_Flip_Base<InputHandler, MatrixType>::compute_beta_pvalue(void){
  
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
        
    UInt n_obs = this->inf_car.getN_obs();

    // Seclect eigenvalues that will not be flipped basing on the estimated bias carried
    VectorXr Tilder_hat = Lambda_dec.eigenvectors().transpose()* (*(this->inf_car.getZp()) - (*W)* beta_hat); // This vector represents Tilder using only beta_hat, needed for bias estimation

    // Estimate the standard error
    VectorXr eps_hat = (*(this->inf_car.getZp())) - (this->inf_car.getZ_hat());
    Real SS_res = eps_hat.squaredNorm();
    Real Sigma_hat = std::sqrt(SS_res/(n_obs-1));

    Real threshold = 5*Sigma_hat; // This threshold is used to determine how many components will not be flipped: we drop those that show large alpha_hat w.r.t. the expected standar error
    UInt N_Eig_Out=0;
    
    // Observed statistic
    VectorXr stat=TildeX*Tilder;
    VectorXr stat_perm=stat;
    
    //Random sign-flips
    std::random_device rd; 
    std::default_random_engine eng{rd()};
    std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)
    Real count_Up=0;
    Real count_Down=0;

    VectorXr Tilder_perm=Tilder;
    
    for(unsigned long int i=0;i<n_flip;i++){
      N_Eig_Out=0;
      for(unsigned long int j=0;j<TildeX.cols();j++){
	UInt flip;
        if((this->inf_car.getInfData()->get_enhanced_inference()=="enhanced") & (N_Eig_Out<n_obs/2) & (fabs(Tilder_hat(j))>threshold)){
	  flip=1;
          ++N_Eig_Out;
        }else{
	  flip=2*distr(eng)-1;
        }
	Tilder_perm(j)=Tilder(j)*flip;
      }
      stat_perm=TildeX*Tilder_perm; // Flipped statistic
      if(is_Unilaterally_Greater(stat_perm,stat)){ ++count_Up;}else{ //Here we use the custom-operator defined in Eigen_Sign_Flip.h
	if(is_Unilaterally_Smaller(stat_perm,stat)){ ++count_Down;} //Here we use the custom-operator defined in Eigen_Sign_Flip.h 
      }
    }
    
    Real pval_Up = count_Up/n_flip;
    Real pval_Down = count_Down/n_flip;

    result.resize(p); // Allocate more space so that R receives a well defined object (different implementations may require higher number of pvalues)
    result(0) = 2*std::min(pval_Up,pval_Down); // Obtain the blateral p_value starting from the unilateral
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
    
    UInt n_obs = this->inf_car.getN_obs();

    // Seclect eigenvalues that will not be flipped basing on the estimated bias carried
    VectorXr Tilder_hat = Lambda_dec.eigenvectors().transpose()* (*(this->inf_car.getZp()) - (*W)* beta_hat); // This vector represents Tilder using only beta_hat, needed for bias estimation

    // Estimate the standard error
    VectorXr eps_hat = (*(this->inf_car.getZp())) - (this->inf_car.getZ_hat());
    Real SS_res = eps_hat.squaredNorm();
    Real Sigma_hat = std::sqrt(SS_res/(n_obs-1));

    Real threshold = 5*Sigma_hat; // This threshold is used to determine how many components will not be flipped: we drop those that show large alpha_hat w.r.t. the expected standar error
    UInt N_Eig_Out=0;

    // Observed statistic
    MatrixXr stat=TildeX*Tilder;
    MatrixXr stat_perm=stat;

    // Random sign-flips
    std::random_device rd; 
    std::default_random_engine eng{rd()};
    std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)VectorXr other_covariates = MatrixXr::Ones(beta_hat.size(),1)-C.row(i);
    VectorXr count_Up = VectorXr::Zero(p);
    VectorXr count_Down = VectorXr::Zero(p);
    
    MatrixXr Tilder_perm=Tilder;
    
    for(unsigned long int i=0;i<n_flip;i++){
      N_Eig_Out=0;
      for(unsigned long int j=0;j<TildeX.cols();j++){
	UInt flip;
	if((this->inf_car.getInfData()->get_enhanced_inference()=="enhanced") & (N_Eig_Out<n_obs/2) & (fabs(Tilder_hat(j))>threshold)){
	  flip=1 ;
	}else{
	  flip=2*distr(eng)-1;
	}
	Tilder_perm.row(j)=Tilder.row(j)*flip;
      }
      stat_perm=TildeX*Tilder_perm; // Flipped statistic
      
      for(UInt k=0; k<p; ++k){
	if(stat_perm(k,k) > stat(k,k)){
	  ++count_Up(k);
	}else{
	  ++count_Down(k);
	}
      } 
    }
    
    VectorXr pval_Up = count_Up/n_flip;
    VectorXr pval_Down = count_Down/n_flip;
    
    result.resize(p);
    result = 2*min(pval_Up,pval_Down); // Obtain the blateral p_value starting from the unilateral
  } 
  return result;
  
};

template<typename InputHandler, typename MatrixType>
Real Eigen_Sign_Flip_Base<InputHandler, MatrixType>::compute_f_pvalue(void){
  // declare the result
  Real p_value; 

  // get all the necessary matrices from the inf_car
  const MatrixXr W_loc = this->inf_car.getW_loc();
  const SpMat Psi_loc = this->inf_car.getPsi_loc();
  const VectorXr Z_loc = this->inf_car.getZ_loc(); 
  const UInt n_loc = this->inf_car.getN_loc(); 

  const VectorXr f_0 = this->inf_car.getInfData()->get_f_0(); 
  unsigned long int n_flip=this->inf_car.getInfData()->get_n_Flip();

  // compute Q_loc
  MatrixXr Q_loc; 
  Q_loc.resize(n_loc, n_loc);
  
  if(this->inf_car.getRegData()->getCovariates()->rows()==0){ // no covariates case
    Q_loc = MatrixXr::Identity(n_loc, n_loc);
  }
  else{ // covariates case
    MatrixXr WtW_loc = W_loc.transpose() * W_loc; 
    Q_loc = MatrixXr::Identity(n_loc, n_loc) - W_loc * WtW_loc.ldlt().solve(W_loc.transpose());
  }

  // compute the residuals under H0
  this->Partial_f_res_H0 = Q_loc*(Z_loc - f_0); 

  // Matrix that groups close location points, needed only when locations are nodes
  MatrixXr Group_res = MatrixXr::Constant(this->inf_car.getInfData()->get_locs_inference().rows(), this->inf_car.getInfData()->get_locs_inference().rows(), 0);

  if(this->inf_car.getInfData()->get_locs_are_nodes_inference()){ 
    // fix the number of residuals to combine 
    UInt k = 5; 

    // get the selected locations
    MatrixXr locations = this->inf_car.getInfData()->get_locs_inference();
    std::vector<UInt> locations_index = this->inf_car.getInfData()->get_locs_index_inference(); 

    std::vector<std::vector<UInt>> NearestIndex; 
    NearestIndex.resize(locations.rows()); 

    for(UInt i=0; i<locations.rows(); ++i){

      Real x0 = locations(i,0);
      Real y0 = locations(i,1);

      NearestIndex[i].resize(k);
      NearestIndex[i] = this->compute_k_closest_points_2D(k, x0, y0, locations, locations_index); 
    }


    // vector that converts global indices into local indices
    VectorXi rel_rows = VectorXi::Constant(this->inf_car.getPsip()->rows(), -1);
    for(UInt i=0; i < locations_index.size(); ++i){
      rel_rows(locations_index[i]) = i; 
    } 

    for(UInt a=0; a < locations.rows(); ++a){
      for(UInt b=0; b < k; ++b)
        Group_res(a, rel_rows(NearestIndex[a][b])) = 1;
    }
   
  }
  
  // eigen-sign-flip implementation
  if(this->inf_car.getInfData()->get_implementation_type()[this->pos_impl] == "eigen-sign-flip"){
    // compute Q_loc decomposition
    Eigen::SelfAdjointEigenSolver<MatrixXr> Q_dec(Q_loc);
    
    MatrixXr V = Q_dec.eigenvectors();
    
    if(this->inf_car.getRegData()->getCovariates()->rows()!=0){
      UInt q = this->inf_car.getq();
      V = Q_dec.eigenvectors().rightCols(n_loc - q);
    }
   
    // update residuals to be flipped
    VectorXr V_Partial_f_res_H0 = V.transpose() * this->Partial_f_res_H0; 
    
    // observed statistics
    VectorXr T; 
 
    if(this->inf_car.getInfData()->get_locs_are_nodes_inference()){
      T = Group_res * V * V_Partial_f_res_H0; 
    }
    else{
      T = Psi_loc.transpose() * V * V_Partial_f_res_H0;
    }
    
    VectorXr T_perm = T;

    // observed final statistic (for simultaneous test) 
    Real T_comb;
    
    VectorXr T_sq = T.array() * T.array();
    T_comb = T_sq.sum();
    Real T_comb_perm = T_comb; 

    // random sign-flips
    std::random_device rd; 
    std::default_random_engine eng{rd()};
    std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)
    Real count=0;
    VectorXr res_perm = V_Partial_f_res_H0; 

    for(unsigned long int i=0; i < n_flip; ++i){
      for(unsigned long int j=0; j < V_Partial_f_res_H0.size(); ++j){
	UInt flip=2*distr(eng)-1;
	res_perm(j)=V_Partial_f_res_H0(j)*flip;
      }
      if(this->inf_car.getInfData()->get_locs_are_nodes_inference()){
        T_perm = Group_res * V * res_perm; 
      }
      else{
        T_perm = Psi_loc.transpose() * V * res_perm;
      }
      // flipped statistics
      VectorXr T_perm_sq = T_perm.array() * T_perm.array();
      T_comb_perm = T_perm_sq.sum();
      
      if(T_comb_perm >= T_comb){ ++count;} 
    }
    p_value = count/n_flip;
    
  }
  else{ // sign-flip implementation
    // observed statistics
    VectorXr T; 
    if(this->inf_car.getInfData()->get_locs_are_nodes_inference()){
      T = Group_res * this->Partial_f_res_H0; 
    }
    else{
      T = Psi_loc.transpose() * this->Partial_f_res_H0;
    }
    
    VectorXr T_perm = T;

    // observed final statistic (for simultaneous test) 
    Real T_comb;
    
    VectorXr T_sq = T.array() * T.array();
    T_comb = T_sq.sum();
    Real T_comb_perm = T_comb; 
    
    // random sign-flips
    std::random_device rd; 
    std::default_random_engine eng{rd()};
    std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)
    Real count=0;
    VectorXr res_perm = this->Partial_f_res_H0; 

    for(unsigned long int i=0; i < n_flip; ++i){
      for(unsigned long int j=0; j < n_loc; ++j){
	UInt flip=2*distr(eng)-1;
	res_perm(j)=this->Partial_f_res_H0(j)*flip;
      }
      if(this->inf_car.getInfData()->get_locs_are_nodes_inference()){
        T_perm = Group_res * res_perm; 
      }
      else{
        T_perm = Psi_loc.transpose() * res_perm; 
      }
      // flipped statistic
      VectorXr T_perm_sq = T_perm.array() * T_perm.array();
      T_comb_perm = T_perm_sq.sum();
      
      if(T_comb_perm >= T_comb){ ++count;} 
    }
    p_value = count/n_flip; 
  } 

  return p_value; 
};


template<typename InputHandler, typename MatrixType>
MatrixXv Eigen_Sign_Flip_Base<InputHandler, MatrixType>::compute_f_CI(void){
  // Eigen-Sign-Flip CI are not implemented for f
  // this function won't be called 
  MatrixXv null_mat; 
  null_mat.resize(1,1);
  null_mat(0) = VectorXr::Constant(3,0);
  return null_mat; 
};

template<typename InputHandler, typename MatrixType>
MatrixXv Eigen_Sign_Flip_Base<InputHandler, MatrixType>::compute_beta_CI(void){

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
  VectorXr ESF_bisection_tolerances = 0.2*Speckman_aux_ranges; // 0.1 of the speckman CI as maximum tolerance

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

  // Seclect eigenvalues that will not be flipped basing on the estimated bias carried
  VectorXr Tilder_hat = Lambda_dec.eigenvectors().transpose()* (*(this->inf_car.getZp()) - (*W)* beta_hat); // This vector represents Tilder using only beta_hat, needed for bias estimation

  VectorXr Partial_res_H0_CI;
  Partial_res_H0_CI.resize(Lambda.cols());

  // fill the p_values matrix
  for (UInt i=0; i<p; i++){
  
    // Build auxiliary vector for residuals computation
    VectorXr other_covariates = MatrixXr::Ones(beta_hat.size(),1)-C.row(i);

    MatrixXr TildeX_loc= TildeX.row(i);

    // compute the partial residuals and p value
    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (UU(i)); // (z-W*beta_hat(non in test)-W*UU[i](in test))
    local_p_values(0,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star, 2);

    // compute the partial residuals and p value
    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (UL(i)); // (z-W*beta_hat(non in test)-W*UL[i](in test))
    local_p_values(1,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star, 2);

    // compute the partial residuals and p value
    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (LU(i)); // (z-W*beta_hat(non in test)-W*LU[i](in test))
    local_p_values(2,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star, 1);

    // compute the partial residuals and p value
    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (LL(i)); // (z-W*beta_hat(non in test)-W*LL[i](in test))
    local_p_values(3,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star, 1);

  }

  // extract the CI significance (1-confidence)
  Real alpha=0;

  if(this->inf_car.getInfData()->get_interval_type()[this->pos_impl]=="one-at-the-time"){
    alpha=0.5*(this->inf_car.getInfData()->get_inference_alpha()(this->pos_impl));
  }else{
    alpha=0.5/p * (this->inf_car.getInfData()->get_inference_alpha()(this->pos_impl));
  }
    
  UInt Max_Iter=50;
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
	  local_p_values(0,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star, 2);
  
	}else{
 
	  if(local_p_values(1,i)<alpha){ // Upper-Lower bound excessively tight
	    UL(i)=beta_hat(i)+0.5*(UL(i)-beta_hat(i));
  
	    // compute the partial residuals
	    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (UL(i)); // (z-W*beta_hat(non in test)-W*UL[i](in test))
	    local_p_values(1,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star, 2);

	  }else{//both the Upper bounds are well defined

	    if(UU(i)-UL(i)<ESF_bisection_tolerances(i)){

	      converged_up[i]=true;

	    }else{

	      Real proposal=0.5*(UU(i)+UL(i));
   
	      // compute the partial residuals
	      Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (proposal); // (z-W*beta_hat(non in test)-W*proposal)
	      Real prop_p_val=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star, 2);

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
	  local_p_values(2,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star, 1);
  
	}else{
 
	  if(local_p_values(3,i)>alpha){ // Lower-Lower bound excessively tight
	    LL(i)=LL(i)-0.5*(LU(i)-LL(i));
  
	    // compute the partial residuals
	    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (LL(i)); // (z-W*beta_hat(non in test)-W*LL[i](in test))
	    local_p_values(3,i)=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star, 1);

	  }else{//both the Upper bounds are well defined

	    if(LU(i)-LL(i)<ESF_bisection_tolerances(i)){

	      converged_low[i]=true;

	    }else{

	      Real proposal=0.5*(LU(i)+LL(i));
   
	      // compute the partial residuals
	      Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (other_covariates) * (beta_hat) - (*W) * (C.row(i)) * (proposal); // (z-W*beta_hat(non in test)-W*proposal)
	      Real prop_p_val=compute_CI_aux_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star, 1);

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


template<typename InputHandler, typename MatrixType> 
const std::vector<UInt> Eigen_Sign_Flip_Base<InputHandler, MatrixType>::compute_k_closest_points_2D(const UInt k, const Real x0, const Real y0, const MatrixXr & locations, const std::vector<UInt> & locations_index) const{

 // prepare the object to be returned
 std::vector<UInt> result;
 result.reserve(k);
 
 // compute the matrix containing the distances in the first column, and the point index in the second one
 std::vector<std::pair<Real, UInt>> distances; 
 distances.resize(locations.rows());

 for(UInt l=0; l<locations.rows(); ++l){
   Real d = (locations(l,0) - x0)*(locations(l,0) - x0) + (locations(l,1) - y0)*(locations(l,1) - y0);
   d = std::sqrt(d); 
  
   distances[l] = std::make_pair(d, locations_index[l]); 
 }

 // sort the rows according to increasing distance (by default it sorts according to the first element)
 std::sort(distances.begin(), distances.end());
 
 for(UInt i=0; i<k; ++i){
   result.push_back(distances[i].second);
 }

 return result;

};


