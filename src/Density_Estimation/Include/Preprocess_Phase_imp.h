#ifndef __PREPROCESS_PHASE_IMP_H__
#define __PREPROCESS_PHASE_IMP_H__


template<UInt ORDER, UInt mydim, UInt ndim>
Preprocess<ORDER, mydim, ndim>::Preprocess(const DataProblem<ORDER, mydim, ndim>& dp,
                                           const FunctionalProblem<ORDER, mydim, ndim>& fp):
  dataProblem_(dp), funcProblem_(fp){

    densityInit_ = DensityInitialization_factory<ORDER, mydim, ndim>::createInitializationSolver(dp, fp);

    fInit_.resize(dp.getNlambda());
    fillFInit();

}


template<UInt ORDER, UInt mydim, UInt ndim>
void Preprocess<ORDER, mydim, ndim>::fillFInit()
{

    for(UInt l = 0; l < dataProblem_.getNlambda(); ++l) {
        fInit_[l] = densityInit_-> chooseInitialization(dataProblem_.getLambda(l));
    }

}


template<UInt ORDER, UInt mydim, UInt ndim>
void NoCrossValidation<ORDER, mydim, ndim>::performPreprocessTask()
{

    this->bestLambda_ = this->dataProblem_.getLambda(0);
    this->gInit_ = (*(this->fInit_[0])).array().log();

}


template<UInt ORDER, UInt mydim, UInt ndim>
CrossValidation<ORDER, mydim, ndim>::CrossValidation(const DataProblem<ORDER, mydim, ndim>& dp,
                                                     const FunctionalProblem<ORDER, mydim, ndim>& fp,
                                                     std::shared_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> ma):
  Preprocess<ORDER, mydim, ndim>(dp, fp), minAlgo_(ma), error_(dp){

    K_folds_.resize(dp.dataSize());
    CV_errors_.resize(dp.getNlambda(), 0);
    g_sols_.resize(dp.getNlambda());

}


template<UInt ORDER, UInt mydim, UInt ndim>
std::pair<VectorXr, Real> CrossValidation<ORDER, mydim, ndim>::performCV()
{

    UInt N = this->dataProblem_.dataSize();
    UInt K = this->dataProblem_.getNfolds();

    for (UInt i = 0; i< N; ++i) {
        UInt length = ((i % K) <= (N % K))? (i % K)*(N/K +1) : (N % K) + (i % K)*(N/K);
        K_folds_[length + i/K] = i;
    }

    // Cycle on the folds
    for (UInt i = 0; i < K; ++i){

        if(this->dataProblem_.Print()) {
            Rprintf("X_valid is the fold number %d\n", i);
        }

        std::vector<UInt> x_valid, x_train;

        if (i < N % K) { // "Big" folds
            std::set_union(K_folds_.cbegin(), K_folds_.cbegin()+ i*(N/K +1), K_folds_.cbegin()+ (i + 1)*(N/K +1), K_folds_.cend(), std::back_inserter(x_train));
            std::copy(K_folds_.cbegin()+ i*(N/K +1), K_folds_.cbegin()+ (i + 1)*(N/K +1), std::back_inserter(x_valid));
        }
        else { // "Small" folds
            std::set_union(K_folds_.cbegin(), K_folds_.cbegin()+ (N % K) + i*(N/K), K_folds_.cbegin()+ (N % K) + (i+1)*(N/K) , K_folds_.cend(), std::back_inserter(x_train));
            std::copy(K_folds_.cbegin()+ (N % K) + i*(N/K), K_folds_.cbegin()+ (N % K) + (i+1)*(N/K), std::back_inserter(x_valid));
        }

        SpMat Psi_train = this->dataProblem_.computePsi(x_train);
        SpMat Psi_valid = this->dataProblem_.computePsi(x_valid);

        performCV_core(i, Psi_train, Psi_valid); // It fills g_sols, CV_errors_

    }

    UInt init_best_lambda = std::distance(CV_errors_.cbegin(), std::min_element(CV_errors_.cbegin(), CV_errors_.cend()));

    return std::pair<VectorXr, Real> (g_sols_[init_best_lambda], this->dataProblem_.getLambda(init_best_lambda));

}


template<UInt ORDER, UInt mydim, UInt ndim>
void CrossValidation<ORDER, mydim, ndim>::performPreprocessTask()
{

    std::tie(this->gInit_, this->bestLambda_) = performCV();

}


template<UInt ORDER, UInt mydim, UInt ndim>
void
SimplifiedCrossValidation<ORDER, mydim, ndim>::performCV_core(UInt fold, const SpMat& Psi_train, const SpMat& Psi_valid)
{

    if(this->dataProblem_.Print()) {
        Rprintf("lambda: %f\n", this->dataProblem_.getLambda(fold));
    }

    this->g_sols_[fold] = this->minAlgo_->apply_core(Psi_train, this->dataProblem_.getLambda(fold), (*(this->fInit_[fold])).array().log());

    this->CV_errors_[fold] = this->error_(Psi_valid, this->g_sols_[fold]);

}


template<UInt ORDER, UInt mydim, UInt ndim>
RightCrossValidation<ORDER, mydim, ndim>::RightCrossValidation(const DataProblem<ORDER, mydim, ndim>& dp,
                                                               const FunctionalProblem<ORDER, mydim, ndim>& fp,
                                                               std::shared_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> ma):
  CrossValidation<ORDER, mydim, ndim>(dp, fp, ma){

    best_loss_.resize(this->dataProblem_.getNlambda(), std::numeric_limits<double>::max());

}


template<UInt ORDER, UInt mydim, UInt ndim>
void RightCrossValidation<ORDER, mydim, ndim>::performCV_core(UInt fold, const SpMat& Psi_train, const SpMat& Psi_valid)
{

    for (UInt l = 0; l < this->dataProblem_.getNlambda(); ++l) {

        std::unique_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> minimizationAlgo = this->minAlgo_->clone();

        if(this->dataProblem_.Print()) {
            Rprintf("lambda: %f\n", this->dataProblem_.getLambda(l));
        }

        VectorXr sols;

        sols = minimizationAlgo->apply_core(Psi_train, this->dataProblem_.getLambda(l), (*(this->fInit_[l])).array().log());

        this->CV_errors_[l] += this->error_(Psi_valid, sols);

        Real loss = std::get<0>(this->funcProblem_.computeFunctional_g(sols, this->dataProblem_.getLambda(l), Psi_train));

        if(loss < best_loss_[l]) {
            best_loss_[l] = loss;
            this->g_sols_[l] = sols;
        }
    }

}

// -----------------------------------------------
// --------------- Preprocess_time ---------------
// -----------------------------------------------


template<UInt ORDER, UInt mydim, UInt ndim>
Preprocess_time<ORDER, mydim, ndim>::Preprocess_time(const DataProblem_time<ORDER, mydim, ndim>& dp,
                                                     const FunctionalProblem_time<ORDER, mydim, ndim>& fp):
  dataProblem_(dp), funcProblem_(fp)
{

    densityInit_ = DensityInitialization_factory_time<ORDER, mydim, ndim>::createInitializationSolver(dp, fp);

    fInit_.resize(dp.getNlambda()*dp.getNlambda_time());
    fillFInit();

}


template<UInt ORDER, UInt mydim, UInt ndim>
void Preprocess_time<ORDER, mydim, ndim>::fillFInit()
{

    for(UInt i = 0; i < dataProblem_.getNlambda(); ++i) {
        for(UInt j = 0; j < dataProblem_.getNlambda_time(); ++j) {
            fInit_[i*dataProblem_.getNlambda_time()+j] = densityInit_->chooseInitialization(dataProblem_.getLambda(i), dataProblem_.getLambda_time(j));
        }
    }
}


template<UInt ORDER, UInt mydim, UInt ndim>
void NoCrossValidation_time<ORDER, mydim, ndim>::performPreprocessTask()
{

    this->bestLambda_S = this->dataProblem_.getLambda(0);
    this->bestLambda_T = this->dataProblem_.getLambda_time(0);

    if(this->dataProblem_.Print()) {
        Rprintf("Best lambda_S: %f,\nBest lambda_T %f\n", this->bestLambda_S, this->bestLambda_T);
    }

    this->gInit_ = (*(this->fInit_[0])).array().log();

}


template<UInt ORDER, UInt mydim, UInt ndim>
CrossValidation_time<ORDER, mydim, ndim>::CrossValidation_time(const DataProblem_time<ORDER, mydim, ndim>& dp,
                                                               const FunctionalProblem_time<ORDER, mydim, ndim>& fp,
                                                               std::shared_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> ma):
  Preprocess_time<ORDER, mydim, ndim>(dp, fp), minAlgo_(ma), error_(dp){

    K_folds_.resize(dp.dataSize());
    CV_errors_.resize(dp.getNlambda()*dp.getNlambda_time(), 0);
    g_sols_.resize(dp.getNlambda()*dp.getNlambda_time());

}


template<UInt ORDER, UInt mydim, UInt ndim>
std::tuple<VectorXr, Real, Real> CrossValidation_time<ORDER, mydim, ndim>::performCV()
{

    UInt N = this->dataProblem_.dataSize();
    UInt K = this->dataProblem_.getNfolds();

    for (UInt i = 0; i < N; ++i) {
        UInt length = ((i % K) <= (N % K)) ? (i % K)*(N/K +1) : (N % K) + (i % K)*(N/K);
        K_folds_[length + i/K] = i;
    }

    // Cycle on the folds
    for (UInt i = 0; i < K; i++){

        if(this->dataProblem_.Print()) {
            Rprintf("X_valid is the fold number %d\n", i);
        }

        std::vector<UInt> x_valid, x_train;

        if (i < N % K) { // "Big" folds
            std::set_union(K_folds_.cbegin(), K_folds_.cbegin() + i*(N/K +1), K_folds_.cbegin() + (i + 1)*(N/K +1), K_folds_.cend(), std::back_inserter(x_train));
            std::copy(K_folds_.cbegin()+ i*(N/K +1), K_folds_.cbegin()+ (i + 1)*(N/K +1), std::back_inserter(x_valid));
        }
        else { // "Small" folds
            std::set_union(K_folds_.cbegin(), K_folds_.cbegin() + (N % K) + i*(N/K), K_folds_.cbegin()+ (N % K) + (i+1)*(N/K) , K_folds_.cend(), std::back_inserter(x_train));
            std::copy(K_folds_.cbegin()+ (N % K) + i*(N/K), K_folds_.cbegin()+ (N % K) + (i+1)*(N/K), std::back_inserter(x_valid));
        }

        SpMat Upsilon_train = this->dataProblem_.computeUpsilon(x_train);
        SpMat Upsilon_valid = this->dataProblem_.computeUpsilon(x_valid);

        performCV_core(i, Upsilon_train, Upsilon_valid); // It fills g_sols, CV_errors_

    }

    UInt init_best_lambda = std::distance(CV_errors_.cbegin(), std::min_element(CV_errors_.cbegin(), CV_errors_.cend()));
    UInt init_best_lambda_S = static_cast<UInt>(init_best_lambda / this->dataProblem_.getNlambda_time());
    UInt init_best_lambda_T = init_best_lambda - init_best_lambda_S * this ->dataProblem_.getNlambda_time();

    if(this->dataProblem_.Print()) {
        Rprintf("Best lambda_S: %f\nBest lambda_T: %f\n", this->dataProblem_.getLambda(init_best_lambda_S), this->dataProblem_.getLambda_time(init_best_lambda_T));
    }

    return std::make_tuple(g_sols_[init_best_lambda], this->dataProblem_.getLambda(init_best_lambda_S), this->dataProblem_.getLambda_time(init_best_lambda_T));
}


template<UInt ORDER, UInt mydim, UInt ndim>
void CrossValidation_time<ORDER, mydim, ndim>::performPreprocessTask()
{

    std::tie(this->gInit_, this->bestLambda_S, this->bestLambda_T) = performCV();

}


template<UInt ORDER, UInt mydim, UInt ndim>
void SimplifiedCrossValidation_time<ORDER, mydim, ndim>::performCV_core(UInt fold, const SpMat& Upsilon_train,
                                                                        const SpMat& Upsilon_valid)
{

    UInt fold_S = static_cast<UInt>(fold / this->dataProblem_.getNlambda_time());
    UInt fold_T = fold - fold_S * this ->dataProblem_.getNlambda_time();

    if(this->dataProblem_.Print()){
        Rprintf("lambda_S: %f\nlambda_T: %f\n", this->dataProblem_.getLambda(fold_S), this->dataProblem_.getLambda_time(fold_T));
    }

    this->g_sols_[fold] = this->minAlgo_->apply_core(Upsilon_train, this->dataProblem_.getLambda(fold_S), this->dataProblem_.getLambda_time(fold_T), (*(this->fInit_[fold])).array().log());

    this->CV_errors_[fold] = this->error_(Upsilon_valid, this->g_sols_[fold]);

}


template<UInt ORDER, UInt mydim, UInt ndim>
RightCrossValidation_time<ORDER, mydim, ndim>::RightCrossValidation_time(const DataProblem_time<ORDER, mydim, ndim>& dp,
                                                                         const FunctionalProblem_time<ORDER, mydim, ndim>& fp,
                                                                         std::shared_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> ma):
  CrossValidation_time<ORDER, mydim, ndim>(dp, fp, ma){

    best_loss_.resize(this->dataProblem_.getNlambda()*this->dataProblem_.getNlambda_time(), std::numeric_limits<double>::max());

}


template<UInt ORDER, UInt mydim, UInt ndim>
void RightCrossValidation_time<ORDER, mydim, ndim>::performCV_core(UInt fold, const SpMat& Upsilon_train,
                                                                   const SpMat& Upsilon_valid)
{

    for (UInt l = 0; l < this->dataProblem_.getNlambda() * this->dataProblem_.getNlambda_time(); ++l) {

        UInt l_S = static_cast<UInt>(l / this->dataProblem_.getNlambda_time());
        UInt l_T = l - l_S * this ->dataProblem_.getNlambda_time();

        std::unique_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> minimizationAlgo = this->minAlgo_->clone();

        if(this->dataProblem_.Print()) {
            Rprintf("lambda_S: %f\nlambda_T: %f\n", this->dataProblem_.getLambda(l_S), this->dataProblem_.getLambda_time(l_T));
        }

        VectorXr sols;

        sols = minimizationAlgo->apply_core(Upsilon_train, this->dataProblem_.getLambda(l_S), this->dataProblem_.getLambda_time(l_T), (*(this->fInit_[l])).array().log());

        this->CV_errors_[l] += this->error_(Upsilon_valid, sols);

        Real loss = std::get<0>(this->funcProblem_.computeFunctional_g(sols, this->dataProblem_.getLambda(l_S), this->dataProblem_.getLambda_time(l_T), Upsilon_train));

        if(loss < best_loss_[l]){
            best_loss_[l] = loss;
            this->g_sols_[l] = sols;
        }
    }

}


#endif
