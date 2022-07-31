#include "Eigen_Sign_Flip.h"
#include "Inference_Factory.h"
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
void Eigen_Sign_Flip_Base<InputHandler, MatrixType>::Compute_wald_aux(void){

  // create or get an already existing Wald inference object 
  std::shared_ptr<Inference_Base<InputHandler,MatrixType>> wald_obj = Inference_Factory<InputHandler,MatrixType>::create_inference_method("wald", this->inverter, this->inf_car, this->pos_impl);
  MatrixXv wald_output = wald_obj->compute_inference_output(); 
  // extract the confidence intervals 
  UInt n_loc = this->inf_car.getN_loc(); 
  UInt p = this->inf_car.getInfData()->get_coeff_inference().rows();
  UInt result_dim = (p > n_loc) ? p : n_loc; 
  // this matrix contains the overall CI, also the ones for beta if they were originally required as retrieved by pos_impl from inferenceData 
  MatrixXv wald_CI = wald_output.rightCols(result_dim);
  // extract the confidence intervals for f
  MatrixXv wald_f_CI = (wald_CI.row(1)).leftCols(n_loc);
  
  Wald_aux_ranges.resize(n_loc);
  // for each location point 
  for(UInt i=0; i<n_loc; ++i){
    // get the half range 
    Real half_range = wald_f_CI(i)(2) - wald_f_CI(i)(1);
    // save the half range
    Wald_aux_ranges(i) = half_range; 
  }
  
  this->is_wald_aux_computed = true;

  return;
}

template<typename InputHandler, typename MatrixType> 
Real Eigen_Sign_Flip_Base<InputHandler, MatrixType>::compute_CI_aux_beta_pvalue(const VectorXr & partial_res_H0_CI, const MatrixXr & TildeX, const VectorXr & Tilder_hat, const  MatrixXr & Tilder_star) const {
  
  // declare the vector that will store the p-values
  Real result;
    
  // compute the vectors needed for the statistic 
  VectorXr Tilder = Tilder_star * partial_res_H0_CI;

  // Initialize observed statistic and sign_flipped statistic
  MatrixXr stat_temp = TildeX*Tilder;
  Real stat=stat_temp(0);
  Real stat_flip=stat;

  UInt n_obs = this->inf_car.getN_obs();

  // Estimate the standard error
  VectorXr eps_hat = (*(this->inf_car.getZp())) - (this->inf_car.getZ_hat());
  Real SS_res = eps_hat.squaredNorm();
  Real Sigma_hat = std::sqrt(SS_res/(n_obs-1));

  Real threshold = 10*Sigma_hat; // This threshold is used to determine how many components will not be flipped: we drop those that show large alpha_hat w.r.t. the expected standar error
  UInt N_Eig_Out=0; // Initialize number of biased components that will be fixed in Enhanced ESF p_value computation
    

  // Random sign-flips
  std::random_device rd; 
  std::default_random_engine eng{rd()};
  std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)
  Real count_Up = 0;   // Counter for the number of flipped statistics that are larger the observed statistic
  Real count_Down = 0; // Counter for the number of flipped statistics that are smaller the observed statistic
    
  VectorXr Tilder_perm=Tilder;
    
  // get the number of flips
  unsigned long int n_flip=this->inf_car.getInfData()->get_n_Flip();

  for(unsigned long int i=0;i<n_flip;i++){
    N_Eig_Out=0;
    for(unsigned long int j=0;j<TildeX.cols();j++){
      UInt flip;
      if((this->inf_car.getInfData()->get_enhanced_inference()[this->pos_impl]==true) & (N_Eig_Out<n_obs/2) & (fabs(Tilder_hat(j))>threshold)){ // Enhanced ESF test has been required and component biased
	flip=1; // Fix the biased component
	++N_Eig_Out;
      }else{
	flip=2*distr(eng)-1;
      }
      Tilder_perm(j)=Tilder(j)*flip;
    }
    MatrixXr stat_flip_temp = TildeX*Tilder_perm; 
    stat_flip= stat_flip_temp(0);// Flipped statistic
    if(stat_flip > stat){ ++count_Up;}else{ 
      if(stat_flip < stat){ ++count_Down;}  
    }
  }
    
  Real pval_Up = count_Up/n_flip;     
  Real pval_Down = count_Down/n_flip; 

  result = std::min(pval_Up, pval_Down); // Selecting the correct unilateral p_value 

  return result;
  
};

template<typename InputHandler, typename MatrixType> 
Real Eigen_Sign_Flip_Base<InputHandler, MatrixType>::compute_CI_aux_f_pvalue(const VectorXr & partial_res_H0_CI, const UInt current_index) const {
  // declare the result
  Real result; 

  // get all the necessary matrices from the inf_car
  const MatrixXr W_loc = this->inf_car.getW_loc();
  const SpMat Psi_loc = this->inf_car.getPsi_loc();
  const VectorXr Z_loc = this->inf_car.getZ_loc(); 
  const UInt n_loc = this->inf_car.getN_loc(); 
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

  // Matrix that groups close location points, needed only when locations are nodes
  const MatrixXr Group_res = this->inf_car.getGroup_loc();

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
    const VectorXr V_partial_res_H0_CI = V.transpose() * partial_res_H0_CI; 

    // observed statistics
    VectorXr T; 
 
    if(this->inf_car.getInfData()->get_locs_are_nodes_inference()){
      T = Group_res * V * V_partial_res_H0_CI; 
    }
    else{
      T = Psi_loc.transpose() * V * V_partial_res_H0_CI;
    }
    
    VectorXr T_perm = T;

    // Random sign-flips
    std::random_device rd; 
    std::default_random_engine eng{rd()};
    std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)
    Real count_Up = 0;   // Counter for the number of flipped statistics that are larger the observed statistic
    Real count_Down = 0; // Counter for the number of flipped statistics that are smaller the observed statistic
    VectorXr res_perm = V_partial_res_H0_CI; 
 
    for(unsigned long int i=0; i < n_flip; ++i){
      for(unsigned long int j=0; j < V_partial_res_H0_CI.size(); ++j){
	UInt flip=2*distr(eng)-1;
	res_perm(j)=V_partial_res_H0_CI(j)*flip;
      }
      if(this->inf_car.getInfData()->get_locs_are_nodes_inference()){
        T_perm = Group_res * V * res_perm; 
      }
      else{
        T_perm = Psi_loc.transpose() * V * res_perm;
      }
      // flipped statistics (one-at-the-time tests on the nodes)
      if(T_perm(current_index) > T(current_index)){++count_Up;}else{
        if(T_perm(current_index) < T(current_index)){++count_Down;}
      }     
    }
    Real pval_Up = count_Up/n_flip;
    Real pval_Down = count_Down/n_flip;

    result = std::min(pval_Up, pval_Down);

  }
  else{ //sign-flip implementation
      
    // observed statistics
    VectorXr T; 
 
    if(this->inf_car.getInfData()->get_locs_are_nodes_inference()){
      T = Group_res * partial_res_H0_CI; 
    }
    else{
      T = Psi_loc.transpose() * partial_res_H0_CI;
    }
    
    VectorXr T_perm = T;

    // Random sign-flips
    std::random_device rd; 
    std::default_random_engine eng{rd()};
    std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)
    Real count_Up = 0;   // Counter for the number of flipped statistics that are larger the observed statistic
    Real count_Down = 0; // Counter for the number of flipped statistics that are smaller the observed statistic
    VectorXr res_perm = partial_res_H0_CI; 
 
    for(unsigned long int i=0; i < n_flip; ++i){
      for(unsigned long int j=0; j < partial_res_H0_CI.size(); ++j){
	UInt flip=2*distr(eng)-1;
	res_perm(j)=partial_res_H0_CI(j)*flip;
      }
      if(this->inf_car.getInfData()->get_locs_are_nodes_inference()){
        T_perm = Group_res * res_perm; 
      }
      else{
        T_perm = Psi_loc.transpose() * res_perm;
      }
      // flipped statistics (one-at-the-time tests on the nodes)
      if(T_perm(current_index) > T(current_index)){++count_Up;}else{
        if(T_perm(current_index) < T(current_index)){++count_Down;}
      }     
    }
    Real pval_Up = count_Up/n_flip;
    Real pval_Down = count_Down/n_flip;

    result = std::min(pval_Up, pval_Down);
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

  // Store beta_hat
  VectorXr beta_hat = (*(this->inf_car.getBeta_hatp()));
  VectorXr beta_hat_mod=beta_hat;
  
  // simultaneous test
  if(this->inf_car.getInfData()->get_test_type()[this->pos_impl] == "simultaneous"){    
    // extract the current betas in test
    for(UInt i=0; i<p; i++){
      for(UInt j=0; j<C.cols(); j++){
        if(C(i,j)>0){beta_hat_mod[j]=beta_0[j];}
      }
    }
    
    // compute the partial residuals
    Partial_res_H0 = *(this->inf_car.getZp()) - (*W) * beta_hat_mod;
    
    // compute the vectors needed for the statistic
    MatrixXr TildeX = (C * W->transpose()) * Lambda_dec.eigenvectors()*Lambda_dec.eigenvalues().asDiagonal();   	// W^t * V * D
    VectorXr Tilder = Lambda_dec.eigenvectors().transpose()*Partial_res_H0;   			        		// V^t * partial_res_H0
        
    UInt n_obs = this->inf_car.getN_obs();

    // Prepare vectors for enhanced-ESF if requested
    VectorXr Tilder_hat = Lambda_dec.eigenvectors().transpose()* (*(this->inf_car.getZp()) - (*W)* beta_hat); // This vector represents Tilder using only beta_hat, needed for bias estimation

    // Estimate the standard error
    VectorXr eps_hat = (*(this->inf_car.getZp())) - (this->inf_car.getZ_hat());
    Real SS_res = eps_hat.squaredNorm();
    Real Sigma_hat = std::sqrt(SS_res/(n_obs-1));

    Real threshold = 10*Sigma_hat; // This threshold is used to determine how many components will not be flipped: we drop those that show large alpha_hat w.r.t. the expected standar error
    UInt N_Eig_Out=0; // It will store the number of biased components that will be kept fixed if enhanced-ESF is required
    
    // Initialize observed statistic and sign-flipped statistic
    VectorXr stat=TildeX*Tilder;
    VectorXr stat_flip=stat;
    
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
        if((this->inf_car.getInfData()->get_enhanced_inference()[this->pos_impl]==true) & (N_Eig_Out<n_obs/2) & (fabs(Tilder_hat(j))>threshold)){ // If enhanced-ESF is required and component is biased
	  flip=1; // Fix the biased component
          ++N_Eig_Out;
        }else{
	  flip=2*distr(eng)-1;
        }
	Tilder_perm(j)=Tilder(j)*flip;
      }
      stat_flip=TildeX*Tilder_perm; // Flipped statistic
      if(is_Unilaterally_Greater(stat_flip,stat)){ ++count_Up;}else{ //Here we use the custom-operator defined in Eigen_Sign_Flip.h
	if(is_Unilaterally_Smaller(stat_flip,stat)){ ++count_Down;} //Here we use the custom-operator defined in Eigen_Sign_Flip.h 
      }
    }
    
    Real pval_Up = count_Up/n_flip;
    Real pval_Down = count_Down/n_flip;

    result.resize(p); // Allocate more space so that R receives a well defined object (different implementations may require higher number of pvalues)
    result(0) = 2*std::min(pval_Up,pval_Down); // Obtain the bilateral p_value starting from the unilateral
    for(UInt k=1;k<p;k++){
      result(k)=10e20;
    }
  }
  else{
    
    // one-at-the-time tests    
    Partial_res_H0.resize(Lambda.cols(), p);
    for(UInt i=0; i<p; ++i){
      // Extract the current beta in test
      beta_hat_mod = beta_hat;

      for(UInt j=0; j<C.cols(); j++){
        if(C(i,j)>0){beta_hat_mod[j]=beta_0[j];}
      }
      // compute the partial residuals
      Partial_res_H0.col(i) = *(this->inf_car.getZp()) - (*W) * beta_hat_mod; // (z-W*beta_hat(non in test)-W*beta_0(in test))
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

    Real threshold = 10*Sigma_hat; // This threshold is used to determine how many components will not be flipped: we drop those that show large alpha_hat w.r.t. the expected standar error
    UInt N_Eig_Out=0; // It will store the number of biased components that will be kept fixed if enhanced-ESF is required

    // Observed statistic
    MatrixXr stat=TildeX*Tilder;
    MatrixXr stat_flip=stat;

    // Random sign-flips
    std::random_device rd; 
    std::default_random_engine eng{rd()};
    std::uniform_int_distribution<> distr{0,1}; // Bernoulli(1/2)
    VectorXr count_Up = VectorXr::Zero(p);
    VectorXr count_Down = VectorXr::Zero(p);
    
    MatrixXr Tilder_perm=Tilder;
    
    for(unsigned long int i=0;i<n_flip;i++){
      N_Eig_Out=0;
      for(unsigned long int j=0;j<TildeX.cols();j++){
	UInt flip;
	if((this->inf_car.getInfData()->get_enhanced_inference()[this->pos_impl]==true) & (N_Eig_Out<n_obs/2) & (fabs(Tilder_hat(j))>threshold)){ // If enhanced-ESF is required and component is biased
	  flip=1; // Fix the biased component
          ++N_Eig_Out;
	}else{
	  flip=2*distr(eng)-1;
	}
	Tilder_perm.row(j)=Tilder.row(j)*flip;
      }
      stat_flip=TildeX*Tilder_perm; // Flipped statistic
      
      for(UInt k=0; k<p; ++k){
	if(stat_flip(k,k) > stat(k,k)){
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
  MatrixXr Group_res = this->inf_car.getGroup_loc(); 
  
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
  // This function will be called only when the chosen locations are a subset of the mesh nodes  
  
  // compute/get all the needed objects from the inference carrier 
  const MatrixXr W_loc = this->inf_car.getW_loc();
  const SpMat Psi_loc = this->inf_car.getPsi_loc();
  const VectorXr Z_loc = this->inf_car.getZ_loc(); 
  const UInt n_loc = this->inf_car.getN_loc(); 
  const std::vector<UInt> sub_locations_index = this->inf_car.getInfData()->get_locs_index_inference(); 
  
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

  // Matrix that groups close location points, needed only when locations are nodes
  MatrixXr Group_res = this->inf_car.getGroup_loc();

  // get the estimator for f in the selected locations 
  UInt nnodes = this->inf_car.getN_nodes();
  const VectorXr f_hat = this->inf_car.getSolutionp()->topRows(nnodes);
  const VectorXr f_hat_loc = Psi_loc * f_hat; 

  // 2) declare the matrix that will store the confidence intervals 
  MatrixXv result; 
  result.resize(1, n_loc);

  // 3) initialization steps
  // compute initial ranges from Wald confidence intervals (initial guess) 
  if(!is_wald_aux_computed){
    this->Compute_wald_aux(); 
  }
  
  // this vector will store the tolerance for each interval upper/lower limit
  VectorXr bisection_tolerances = 0.2*Wald_aux_ranges; 

  // define storage structures for bisection algorithms
  VectorXr UU; // Upper limit for Upper bound
  UU.resize(n_loc);
  VectorXr UL; // Lower limit for Upper bound
  UL.resize(n_loc);
  VectorXr LU; // Upper limit for Lower bound
  LU.resize(n_loc);
  VectorXr LL; // Lower limit for Lower bound
  LL.resize(n_loc);

  // wald intervals initialization
  for(UInt i=0; i<n_loc; ++i){
    result(i).resize(3);
    Real half_range = Wald_aux_ranges(i);
    
    // compute the limits of the interval
    result(i)(0) = f_hat_loc(i) - half_range;
    LU(i)=result(i)(0)+0.5*half_range;
    LL(i)=result(i)(0)-0.5*half_range;
    result(i)(2) = f_hat_loc(i) + half_range;
    UU(i)=result(i)(2)+0.5*half_range;
    UL(i)=result(i)(2)-0.5*half_range; 	
  }

  // define booleans used to unserstand which CI need to be computed forward on
  std::vector<bool> converged_up(n_loc,false);
  std::vector<bool> converged_low(n_loc,false);
  bool all_f_converged=false;

  // matrix that stores p_values of the bounds at actual step
  MatrixXr local_p_values;
  local_p_values.resize(4,n_loc);

  // compute the vectors needed for the statistic
  
  VectorXr Partial_res_H0_CI;
  Partial_res_H0_CI.resize(n_loc);

  // fill the p_values matrix
  for (UInt i=0; i<n_loc; i++){

    VectorXr f_hat_loc_UU = f_hat_loc;
    VectorXr f_hat_loc_UL = f_hat_loc;
    VectorXr f_hat_loc_LU = f_hat_loc;
    VectorXr f_hat_loc_LL = f_hat_loc; 
    
    for(UInt j=0; j < n_loc; ++j){
      if(Group_res(sub_locations_index[i],j) == 1){
	f_hat_loc_UU(j) = UU(i);  
	f_hat_loc_UL(j) = UL(i); 
	f_hat_loc_LU(j) = LU(i); 
	f_hat_loc_LL(j) = LL(i); 
      }
    }
    
    // compute the partial residuals and p value
    Partial_res_H0_CI = Q_loc*(Z_loc - f_hat_loc_UU); 
    local_p_values(0,i)=compute_CI_aux_f_pvalue(Partial_res_H0_CI, sub_locations_index[i]);

    // compute the partial residuals and p value
    Partial_res_H0_CI = Q_loc*(Z_loc - f_hat_loc_UL);
    local_p_values(1,i)=compute_CI_aux_f_pvalue(Partial_res_H0_CI, sub_locations_index[i]);

    // compute the partial residuals and p value
    Partial_res_H0_CI = Q_loc*(Z_loc - f_hat_loc_LU);
    local_p_values(2,i)=compute_CI_aux_f_pvalue(Partial_res_H0_CI, sub_locations_index[i]);

    // compute the partial residuals and p value
    Partial_res_H0_CI = Q_loc*(Z_loc - f_hat_loc_LL);
    local_p_values(3,i)=compute_CI_aux_f_pvalue(Partial_res_H0_CI, sub_locations_index[i]);


  }

  // extract the CI significance (1-confidence)
  Real alpha=0;
  // only one-at-the-time CI are available for f 
  alpha=0.5*(this->inf_car.getInfData()->get_inference_alpha()(this->pos_impl));
 
  UInt Max_Iter=50;
  UInt Count_Iter=0;
  while((!all_f_converged) && (Count_Iter<Max_Iter)){
  
    // Compute all p_values (only those needed)
    for (UInt i=0; i<n_loc; i++){

      VectorXr f_hat_loc_UU = f_hat_loc;
      VectorXr f_hat_loc_UL = f_hat_loc;
      VectorXr f_hat_loc_LU = f_hat_loc;
      VectorXr f_hat_loc_LL = f_hat_loc; 
  
      if(!converged_up[i]){
	if(local_p_values(0,i)>alpha){ // Upper-Upper bound excessively tight

	  UU(i)=UU(i)+0.5*(UU(i)-UL(i));
	  for(UInt j=0; j < n_loc; ++j){
	    if(Group_res(sub_locations_index[i],j) == 1){
	      f_hat_loc_UU(j) = UU(i);  
	    }
	  }
          // compute the partial residuals
	  Partial_res_H0_CI = Q_loc*(Z_loc - f_hat_loc_UU); 
	  local_p_values(0,i)=compute_CI_aux_f_pvalue(Partial_res_H0_CI, sub_locations_index[i]);
  
	}else{
 
	  if(local_p_values(1,i)<alpha){ // Upper-Lower bound excessively tight
	    UL(i)=f_hat_loc(i)+0.5*(UL(i)-f_hat_loc(i));
	    for(UInt j=0; j < n_loc; ++j){
	      if(Group_res(sub_locations_index[i],j) == 1){
		f_hat_loc_UL(j) = UL(i);  
	      }
	    } 
  
	    // compute the partial residuals
	    Partial_res_H0_CI = Q_loc*(Z_loc - f_hat_loc_UL); 
	    local_p_values(1,i)=compute_CI_aux_f_pvalue(Partial_res_H0_CI, sub_locations_index[i]);

	  }else{//both the Upper bounds are well defined

	    if(UU(i)-UL(i)<bisection_tolerances(i)){

	      converged_up[i]=true;

	    }else{

	      Real proposal=0.5*(UU(i)+UL(i));
              VectorXr f_hat_loc_proposal = f_hat_loc; 

              for(UInt j=0; j < n_loc; ++j){
        	if(Group_res(sub_locations_index[i],j) == 1){
           	  f_hat_loc_proposal(j) = proposal;
       		}
	      }
 
	      // compute the partial residuals
	      Partial_res_H0_CI = Q_loc*(Z_loc - f_hat_loc_proposal);
	      Real prop_p_val=compute_CI_aux_f_pvalue(Partial_res_H0_CI, sub_locations_index[i]);


	      if(prop_p_val<=alpha){UU(i)=proposal; local_p_values(0,i)=prop_p_val;}else{UL(i)=proposal;local_p_values(1,i)=prop_p_val;}
	    }
	  }
	}
      }


      if(!converged_low[i]){
	if(local_p_values(2,i)<alpha){ // Lower-Upper bound excessively tight

	  LU(i)=f_hat_loc(i)-0.5*(f_hat_loc(i)-LU(i));
          for(UInt j=0; j < n_loc; ++j){
	    if(Group_res(sub_locations_index[i],j) == 1){
	      f_hat_loc_LU(j) = LU(i);  
	    }
	  }
  
	  // compute the partial residuals
	  Partial_res_H0_CI = Q_loc*(Z_loc - f_hat_loc_LU);
	  local_p_values(2,i)=compute_CI_aux_f_pvalue(Partial_res_H0_CI, sub_locations_index[i]);
  
	}else{
 
	  if(local_p_values(3,i)>alpha){ // Lower-Lower bound excessively tight
	    LL(i)=LL(i)-0.5*(LU(i)-LL(i));
            for(UInt j=0; j < n_loc; ++j){
	      if(Group_res(sub_locations_index[i],j) == 1){
		f_hat_loc_LL(j) = LL(i);  
	      }
	    }
  
	    // compute the partial residuals
	    Partial_res_H0_CI = Q_loc*(Z_loc - f_hat_loc_LL);
	    local_p_values(3,i)=compute_CI_aux_f_pvalue(Partial_res_H0_CI, sub_locations_index[i]);

	  }else{//both the Upper bounds are well defined

	    if(LU(i)-LL(i)<bisection_tolerances(i)){

	      converged_low[i]=true;

	    }else{

	      Real proposal=0.5*(LU(i)+LL(i));
              VectorXr f_hat_loc_proposal = f_hat_loc; 

              for(UInt j=0; j < n_loc; ++j){
        	if(Group_res(sub_locations_index[i],j) == 1){
           	  f_hat_loc_proposal(j) = proposal;
       		}
	      }
 
   
	      // compute the partial residuals
	      Partial_res_H0_CI = Q_loc*(Z_loc - f_hat_loc_proposal);
	      Real prop_p_val=compute_CI_aux_f_pvalue(Partial_res_H0_CI, sub_locations_index[i]);

	      if(prop_p_val<=alpha){LL(i)=proposal; local_p_values(3,i)=prop_p_val;}else{LU(i)=proposal;local_p_values(2,i)=prop_p_val;}
	    }
	  }
	}
      }
    }
    all_f_converged =true;
    for(UInt j=0; j<n_loc; j++){

      if(!converged_up[j] || !converged_low[j]){
	all_f_converged=false;
      }
    }

    Count_Iter++;

  }

  // give a warning in R if at least one interval did not converge
  if(!all_f_converged)
    Rprintf("warning: some sign-flip confidence intervals for the nonparametric component did not converge");
   
  // for each point
  for(UInt i=0; i<n_loc; ++i){
    result(i).resize(3);

    if(!converged_up[i] || !converged_low[i]){
      // Central element
      result(i)(1)=10e20;
      // Limits of the interval
      result(i)(0) = 10e20;
      result(i)(2) = 10e20; 
    }
    
    else{ // No discrepancy between f_hat_loc(i) and ESF, bisection converged
      // Central element
      result(i)(1)=f_hat_loc(i);
      // Limits of the interval
      result(i)(0) = 0.5*(LU(i)+LL(i));
      result(i)(2) = 0.5*(UU(i)+UL(i)); 
    }
  }

  return result;  


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
  VectorXr beta_hat_mod=beta_hat;

  // get the matrix of coefficients
  MatrixXr C = this->inf_car.getInfData()->get_coeff_inference();
  UInt p=C.rows();
  std::vector<UInt> beta_in_test; // In this vector are stored the respective positions of the beta we are testing in the actual test
  beta_in_test.resize(p);
  for(UInt i=0; i<p; i++){
    for(UInt j=0; j<C.cols(); j++){
      if(C(i,j)>0){beta_in_test[i]=j;}
    }
  }

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
    result(i)(0) = beta_hat(beta_in_test[i]) - half_range;
    LU(i)=result(i)(0)+0.5*half_range;
    LL(i)=result(i)(0)-0.5*half_range;
    result(i)(2) = beta_hat(beta_in_test[i]) + half_range;
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
  MatrixXr TildeX = (W->transpose()) * Lambda_dec.eigenvectors()*Lambda_dec.eigenvalues().asDiagonal();   	// W^t * V * D
  MatrixXr Tilder_star = Lambda_dec.eigenvectors().transpose();   			        		// V^t

  // Seclect eigenvalues that will not be flipped basing on the estimated bias carried
  VectorXr Tilder_hat = Lambda_dec.eigenvectors().transpose()* (*(this->inf_car.getZp()) - (*W)* beta_hat); // This vector represents Tilder using only beta_hat, needed for bias estimation

  VectorXr Partial_res_H0_CI;
  Partial_res_H0_CI.resize(Lambda.cols());

  // fill the p_values matrix
  for (UInt i=0; i<p; i++){
    MatrixXr TildeX_loc= TildeX.row(beta_in_test[i]);
    beta_hat_mod = beta_hat;

    // compute the partial residuals and p value
    beta_hat_mod(beta_in_test[i])=UU(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*UU[i](in test))
    local_p_values(0,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star);

    // compute the partial residuals and p value
    beta_hat_mod(beta_in_test[i])=UL(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*UL[i](in test))
    local_p_values(1,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star);

    // compute the partial residuals and p value
    beta_hat_mod(beta_in_test[i])=LU(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*LU[i](in test))
    local_p_values(2,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star);

    // compute the partial residuals and p value
    beta_hat_mod(beta_in_test[i])=LL(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*LL[i](in test))
    local_p_values(3,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star);

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
  while((!all_betas_converged) & (Count_Iter<Max_Iter)){
  
    // Compute all p_values (only those needed)
    for (UInt i=0; i<p; i++){

      MatrixXr TildeX_loc= TildeX.row(beta_in_test[i]);
      beta_hat_mod = beta_hat;
  
      if(!converged_up[i]){
	if(local_p_values(0,i)>alpha){ // Upper-Upper bound excessively tight

	  UU(i)=UU(i)+0.5*(UU(i)-UL(i));
  
	  // compute the partial residuals
          beta_hat_mod(beta_in_test[i])=UU(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
	  Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*UU[i](in test))
	  local_p_values(0,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star);
  
	}else{
 
	  if(local_p_values(1,i)<alpha){ // Upper-Lower bound excessively tight
	    UL(i)=beta_hat(beta_in_test[i])+0.5*(UL(i)-beta_hat(beta_in_test[i]));
  
	    // compute the partial residuals
            beta_hat_mod(beta_in_test[i])=UL(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
	    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*UL[i](in test))
	    local_p_values(1,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star);

	  }else{//both the Upper bounds are well defined

	    if(UU(i)-UL(i)<ESF_bisection_tolerances(i)){

	      converged_up[i]=true;

	    }else{

	      Real proposal=0.5*(UU(i)+UL(i));
   
	      // compute the partial residuals
              beta_hat_mod(beta_in_test[i])=proposal; // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
	      Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*proposal)
	      Real prop_p_val=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star);

	      if(prop_p_val<=alpha){UU(i)=proposal; local_p_values(0,i)=prop_p_val;}else{UL(i)=proposal;local_p_values(1,i)=prop_p_val;}
	    }
	  }
	}
      }


      if(!converged_low[i]){
	if(local_p_values(2,i)<alpha){ // Lower-Upper bound excessively tight

	  LU(i)=beta_hat(beta_in_test[i])-0.5*(beta_hat(beta_in_test[i])-LU(i));
  
	  // compute the partial residuals
          beta_hat_mod(beta_in_test[i])=LU(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
	  Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*LU[i](in test))
	  local_p_values(2,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star);
  
	}else{
 
	  if(local_p_values(3,i)>alpha){ // Lower-Lower bound excessively tight
	    LL(i)=LL(i)-0.5*(LU(i)-LL(i));
  
	    // compute the partial residuals
            beta_hat_mod(beta_in_test[i])=LL(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
	    Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*LL[i](in test))
	    local_p_values(3,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star);

	  }else{//both the Upper bounds are well defined

	    if(LU(i)-LL(i)<ESF_bisection_tolerances(i)){

	      converged_low[i]=true;

	    }else{

	      Real proposal=0.5*(LU(i)+LL(i));
   
	      // compute the partial residuals
              beta_hat_mod(beta_in_test[i])=proposal; // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
	      Partial_res_H0_CI = *(this->inf_car.getZp()) - (*W) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*proposal)
	      Real prop_p_val=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_hat, Tilder_star);

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
      result(i)(1)=beta_hat(beta_in_test[i]);
      
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
  // extract areal matrix (VectorXr)
  const VectorXr * A = this->inf_car.getAp();
  
  UInt n_obs = this->inf_car.getN_obs();
  UInt n_nodes = this->inf_car.getN_nodes();
  const SpMat * Psi = this->inf_car.getPsip();
  const SpMat * Psi_t = this->inf_car.getPsi_tp();
  UInt q = this->inf_car.getq(); 
  
  this->Lambda.resize(n_obs,n_obs);

  if(this->inf_car.getRegData()->getNumberOfRegions()>0){
    this->Lambda = (MatrixXr::Identity(n_obs,n_obs) - (*Psi)*((*E_inv).block(0,0, n_nodes, n_nodes)*(*Psi_t)*(A->asDiagonal()))); // I - Psi(Psi^T A Psi + P)^-1 Psi^T A
  }else{
    this->Lambda = (MatrixXr::Identity(n_obs,n_obs) - (*Psi)*((*E_inv).block(0,0, n_nodes, n_nodes)*(*Psi_t))); // I - Psi(Psi^T Psi + P)^-1 Psi^T
  }
  this->is_Lambda_computed = true;
  
  return; 
};
