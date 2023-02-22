#ifndef __OPTIMIZATION_ALGORITHM_IMP_H__
#define __OPTIMIZATION_ALGORITHM_IMP_H__


template<UInt ORDER, UInt mydim, UInt ndim>
MinimizationAlgorithm<ORDER, mydim, ndim>::MinimizationAlgorithm(const DataProblem<ORDER, mydim, ndim>& dp,
                                                                 const FunctionalProblem<ORDER, mydim, ndim>& fp,
                                                                 const std::string& d):
  dataProblem_(dp), funcProblem_(fp)
{

    direction_ = DescentDirection_factory<ORDER,  mydim,  ndim>::createDirectionSolver(dp, fp, d);

}


template<UInt ORDER, UInt mydim, UInt ndim>
MinimizationAlgorithm<ORDER, mydim, ndim>::MinimizationAlgorithm(const MinimizationAlgorithm<ORDER, mydim, ndim>& rhs):
  dataProblem_(rhs.dataProblem_), funcProblem_(rhs.funcProblem_)
{

    direction_ = rhs.direction_->clone();

};


template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> FixedStep<ORDER, mydim, ndim>::clone() const
{

  return fdaPDE::make_unique<FixedStep<ORDER, mydim, ndim>>(*this);

}


template<UInt ORDER, UInt mydim, UInt ndim>
VectorXr FixedStep<ORDER,mydim,ndim>::apply_core(const SpMat& Psi, Real lambda, const VectorXr& g) const
{

    // Termination criteria variables
    const Real toll1 = this->dataProblem_.getTol1(), toll2 = this->dataProblem_.getTol2();
    Real norm_grad, dloss = toll1+1, dllik = toll1+1, dpen = toll1+1;

    // To save the current point
    VectorXr g_curr;

    // Variables
    VectorXr grad, d;
    Real loss, loss_old, llik, llik_old, pen, pen_old;
    UInt i;

    for(UInt e = 0; e < this->dataProblem_.getNstepProposals(); ++e)
    {
        // Start always with the initial point
        g_curr = g;

        std::tie(loss, grad, llik, pen) = this->funcProblem_.computeFunctional_g(g_curr, lambda, Psi);
        norm_grad = std::sqrt(grad.dot(grad));

        if(this->dataProblem_.Print()){
            Rprintf("loss %f, llik %f, pen %f, norm_Lp %f\n", loss, llik, pen, norm_grad);
        }

        for(i = 0; i < this->dataProblem_.getNsimulations() && (dloss > toll1 || dllik > toll1 || dpen > toll1) && norm_grad > toll2; ++i) {

            // Termination criteria variables
            loss_old = loss;
            llik_old = llik;
            pen_old = pen;

            // Compute a descent direction
            d = this->direction_->computeDirection(g_curr, grad);

            // Update the point
            g_curr = g_curr + this->dataProblem_.getStepProposals(e)*d;

            // Update termination criteria variables
            std::tie(loss, grad, llik, pen) = this->funcProblem_.computeFunctional_g(g_curr, lambda, Psi);
            dloss = std::abs((loss - loss_old)/loss_old);
            dllik = std::abs((llik - llik_old)/llik_old);
            dpen = std::abs((pen - pen_old)/pen_old);
            norm_grad = std::sqrt(grad.dot(grad));

            // Check if the step is ok
            if((loss_old - loss) < 0){
                if(this->dataProblem_.Print()){
                    Rprintf("The loss function increases: not good. Try decreasing the optimization parameter.\n");
                }
                break;
            }

            if(this->dataProblem_.Print()){
                Rprintf("Iter %d, loss %f, llik %f, pen %f, norm_Lp %f\n", i+1, loss, llik, pen, norm_grad);
            }
        }

        this->direction_->resetParameters();

        if ((loss_old - loss) < 0){}
        else if(dloss <= toll1 && dllik <= toll1 && dpen <= toll1){
            if(this->dataProblem_.Print()){
                Rprintf("The algorithm reaches the tolerance in terms of the functional. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
            }
            return g_curr;
        }
        else if(norm_grad <= toll2){
            if(this->dataProblem_.Print()){
                Rprintf("The algorithm reaches the tolerance in terms of the slope. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
            }
            return g_curr;
        }
        else if(i == this->dataProblem_.getNsimulations()){
            if(this->dataProblem_.Print()){
                Rprintf("The algorithm reaches the maximum number of iterations. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
            }
            return g_curr;
        }
    }

    // If you arrive here you don't have a good gradient parameter
    Rprintf("ERROR: The loss function increases: not good. Try decreasing the optimization parameter.\n");
    //std::abort();
    return VectorXr::Constant(g.size(),0);

}


template<UInt ORDER, UInt mydim, UInt ndim>
VectorXr AdaptiveStep<ORDER,mydim,ndim>::apply_core(const SpMat& Psi, Real lambda, const VectorXr& g) const
{

    // Termination criteria variables
    const Real toll1 = this->dataProblem_.getTol1(), toll2 = this->dataProblem_.getTol2();
    Real norm_grad, dloss = toll1+1, dllik = toll1+1, dpen = toll1+1;

    // To save the current point
    VectorXr g_curr = g;

    // Variables
    VectorXr grad, d;
    Real loss, loss_old, llik, llik_old, pen, pen_old, step;
    UInt i;

    std::tie(loss, grad, llik, pen) = this->funcProblem_.computeFunctional_g(g_curr, lambda, Psi);
    norm_grad = std::sqrt(grad.dot(grad));

    if(this->dataProblem_.Print()) {
        Rprintf("loss %f, llik %f, pen %f, norm_Lp %f\n", loss, llik, pen, norm_grad);
    }

    for(i = 0; i < this->dataProblem_.getNsimulations() && (dloss > toll1 || dllik > toll1 || dpen > toll1) && norm_grad > toll2; ++i) {

        // Termination criteria variables
        loss_old = loss;
        llik_old = llik;
        pen_old = pen;

        // Compute a descent direction
        d = this->direction_->computeDirection(g_curr, grad);

        // Compute a step
        step = computeStep(g_curr, loss, grad, d, lambda, Psi);

        // Update the point
        g_curr = g_curr + step*d;

        // Update termination criteria variables
        std::tie(loss, grad, llik, pen) = this->funcProblem_.computeFunctional_g(g_curr, lambda, Psi);
        dloss = std::abs((loss - loss_old)/loss_old);
        dllik = std::abs((llik - llik_old)/llik_old);
        dpen = std::abs((pen - pen_old)/pen_old);
        norm_grad = std::sqrt(grad.dot(grad));

        if(this->dataProblem_.Print()){
            Rprintf("Iter %d, loss %f, llik %f, pen %f, norm_Lp %f\n", i+1, loss, llik, pen, norm_grad);
        }

    }

    this->direction_->resetParameters();

    if(dloss <= toll1 && dllik <= toll1 && dpen <= toll1) {
        if(this->dataProblem_.Print()){
            Rprintf("The algorithm reaches the tolerance in terms of the functional. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
        }
        return g_curr;
    }
    else if(norm_grad <= toll2) {
        if(this->dataProblem_.Print()){
            Rprintf("The algorithm reaches the tolerance in terms of the slope. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
        }
        return g_curr;
    }
    else {
        if(this->dataProblem_.Print()){
            Rprintf("The algorithm reaches the maximum number of iterations. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
        }
        return g_curr;
    }
}


template<UInt ORDER, UInt mydim, UInt ndim>
Real BacktrackingMethod<ORDER,mydim,ndim>::computeStep(const VectorXr& g, Real loss, const VectorXr& grad,
                                                       const VectorXr& dir, Real lambda, const SpMat& Psi) const
{

    Real ro = 0.5, alpha = 1/ro, c = 0.5;

    Real loss_new, llik_new, pen_new, slope, grad_dir;
    VectorXr grad_new, new_point;

    grad_dir =  grad.dot(dir);

    do{
        // Update step
        alpha *= ro;

        slope = c*alpha*(grad_dir);

        // Update the point
        new_point = g + alpha*dir;

        // Functional in the new point
        std::tie(loss_new, grad_new, llik_new, pen_new) = this->funcProblem_.computeFunctional_g(new_point, lambda, Psi);

    } while(loss_new > (loss + slope));

    return alpha;
}


template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> BacktrackingMethod<ORDER, mydim, ndim>::clone() const
{

  return fdaPDE::make_unique<BacktrackingMethod<ORDER, mydim, ndim>>(*this);

}


template<UInt ORDER, UInt mydim, UInt ndim>
Real WolfeMethod<ORDER,mydim,ndim>::computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir,
                                                Real lambda, const SpMat& Psi) const
{

    Real alpha = 1, alphamax = 0, alphamin = 0, c1 = 1e-4, c2 = 0.9;

    Real loss_new, llik_new, pen_new, slope, grad_dir;
    VectorXr grad_new, new_point;

    grad_dir = grad.dot(dir);
    slope = c1*alpha*grad_dir;

    // Update the point
    new_point = g + alpha*dir;

    // Functional in the new point
    std::tie(loss_new, grad_new, llik_new, pen_new) = this->funcProblem_.computeFunctional_g(new_point, lambda, Psi);

    bool again = true;

    while(again) {

        again = false;

        while(loss_new > (loss + slope)){
            // Update step
            alphamax = alpha;
            alpha = 0.5*(alphamin + alphamax);

            // Try with the new point
            new_point = g + alpha*dir;
            std::tie(loss_new, grad_new, llik_new, pen_new) = this->funcProblem_.computeFunctional_g(new_point, lambda, Psi);
            slope = c1*alpha*grad_dir;
        }

        if(grad_new.dot(dir) < c2*grad_dir && std::abs(grad_dir) > 1e-2){

            again = true;

            // Update step
            alphamin = alpha;
            alpha = alphamax==0 ? 2*alphamin : 0.5*(alphamin+alphamax);

            // Try with the new point
            new_point = g + alpha*dir;
            std::tie(loss_new, grad_new, llik_new, pen_new) = this->funcProblem_.computeFunctional_g(new_point, lambda, Psi);
            slope =  alpha*c1*grad_dir;
        }
    }

    return alpha;
}


template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> WolfeMethod<ORDER, mydim, ndim>::clone() const
{

  return fdaPDE::make_unique<WolfeMethod<ORDER, mydim, ndim>>(*this);

}

// ----------------------------------------------------------
// --------------- MinimizationAlgorithm_time ---------------
// ----------------------------------------------------------

template<UInt ORDER, UInt mydim, UInt ndim>
MinimizationAlgorithm_time<ORDER, mydim, ndim>::MinimizationAlgorithm_time(const DataProblem_time<ORDER, mydim, ndim>& dp,
                                                                           const FunctionalProblem_time<ORDER, mydim, ndim>& fp,
                                                                           const std::string& d):
  dataProblem_(dp), funcProblem_(fp)
{

    direction_ = DescentDirection_factory_time<ORDER,  mydim,  ndim>::createDirectionSolver(dp, fp, d);

}

template<UInt ORDER, UInt mydim, UInt ndim>
MinimizationAlgorithm_time<ORDER, mydim, ndim>::MinimizationAlgorithm_time(const MinimizationAlgorithm_time<ORDER, mydim, ndim>& rhs):
  dataProblem_(rhs.dataProblem_), funcProblem_(rhs.funcProblem_)
{

    direction_ = rhs.direction_->clone();
};

template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> FixedStep_time<ORDER, mydim, ndim>::clone() const
{

    return fdaPDE::make_unique<FixedStep_time<ORDER, mydim, ndim>>(*this);

}

template<UInt ORDER, UInt mydim, UInt ndim>
VectorXr
FixedStep_time<ORDER,mydim,ndim>::apply_core(const SpMat& Upsilon, Real lambda_S, Real lambda_T,const VectorXr& g) const{

    // Termination criteria variables
    const Real toll1 = this->dataProblem_.getTol1(), toll2 = this->dataProblem_.getTol2();
    Real norm_grad, dloss = toll1+1, dllik = toll1+1, dpen_S = toll1+1, dpen_T = toll1+1;

    // To save the current point
    VectorXr g_curr;

    // Variables
    VectorXr grad, d;
    Real loss, loss_old, llik, llik_old, pen_S, pen_T, pen_old_S, pen_old_T;
    UInt i;

    for(UInt e = 0; e < this->dataProblem_.getNstepProposals(); ++e) {
        // Start always with the initial point
        g_curr = g;

        std::tie(loss, grad, llik, pen_S, pen_T) = this->funcProblem_.computeFunctional_g(g_curr, lambda_S, lambda_T, Upsilon);
        norm_grad = std::sqrt(grad.dot(grad));

        if(this->dataProblem_.Print()){
            Rprintf("loss %f, llik %f, pen_S %f, pen_T %f, norm_Lp %f\n", loss, llik, pen_S, pen_T, norm_grad);
        }

        for(i = 0; i < this->dataProblem_.getNsimulations() && (dloss > toll1 || dllik > toll1 || dpen_S > toll1 || dpen_T > toll1) && norm_grad > toll2; ++i) {

            // Termination criteria variables
            loss_old = loss;
            llik_old = llik;
            pen_old_S = pen_S;
            pen_old_T = pen_T;

            // Compute a descent direction
            d = this->direction_->computeDirection(g_curr, grad);

            // Update the point
            g_curr = g_curr + this->dataProblem_.getStepProposals(e)*d;

            // Update termination criteria variables
            std::tie(loss, grad, llik, pen_S, pen_T) = this->funcProblem_.computeFunctional_g(g_curr, lambda_S, lambda_T, Upsilon);
            dloss = std::abs((loss - loss_old)/loss_old);
            dllik = std::abs((llik - llik_old)/llik_old);
            dpen_S = std::abs((pen_S - pen_old_S)/pen_old_S);
            dpen_T = std::abs((pen_T - pen_old_T)/pen_old_T);
            norm_grad = std::sqrt(grad.dot(grad));

            // Check if the step is ok
            if((loss_old - loss) < 0){
                if(this->dataProblem_.Print()){
                    Rprintf("The loss function increases: not good. Try decreasing the optimization parameter.\n");
                }
                break;
            }

            if(this->dataProblem_.Print()) {
                Rprintf("Iter %d, loss %f, llik %f, pen_S %f, pen_T %f, norm_Lp %f\n", i+1, loss, llik, pen_S, pen_T, norm_grad);
            }
        }

        this->direction_->resetParameters();

        if ((loss_old - loss) < 0){}
        else if(dloss <= toll1 && dllik <= toll1 && dpen_S <= toll1 && dpen_T <= toll1){
            if(this->dataProblem_.Print()){
                Rprintf("The algorithm reaches the tolerance in terms of the functional. Norm of Lp: %f, dloss: %f, dllik: %f, dpen_S: %f, dpen_T: %f\n", norm_grad, dloss, dllik, dpen_S, dpen_T);
            }
            return g_curr;
        }
        else if(norm_grad <= toll2){
            if(this->dataProblem_.Print()){
                Rprintf("The algorithm reaches the tolerance in terms of the slope. Norm of Lp: %f, dloss: %f, dllik: %f, dpen_S: %f, dpen_T: %f\n", norm_grad, dloss, dllik, dpen_S, dpen_T);
            }
            return g_curr;
        }
        else if(i == this->dataProblem_.getNsimulations()){
            if(this->dataProblem_.Print()){
                Rprintf("The algorithm reaches the maximum number of iterations. Norm of Lp: %f, dloss: %f, dllik: %f, dpen_S: %f, dpen_T: %f\n", norm_grad, dloss, dllik, dpen_S, dpen_T);
            }
            return g_curr;
        }
    }

    // If you arrive here you don't have a good gradient parameter
    Rprintf("ERROR: The loss function increases: not good. Try decreasing the optimization parameter\n");
    //std::abort();
    return VectorXr::Constant(g.size(),0);

}

template<UInt ORDER, UInt mydim, UInt ndim>
VectorXr AdaptiveStep_time<ORDER,mydim,ndim>::apply_core(const SpMat& Upsilon, Real lambda_S, Real lambda_T,
                                                         const VectorXr& g) const
{

    // Termination criteria variables
    const Real toll1 = this->dataProblem_.getTol1(), toll2 = this->dataProblem_.getTol2();
    Real norm_grad, dloss = toll1+1, dllik = toll1+1, dpen_S = toll1+1, dpen_T = toll1+1;

    // To save the current point
    VectorXr g_curr = g;

    // Variables
    VectorXr grad, d;
    Real loss, loss_old, llik, llik_old, pen_S, pen_T, pen_old_S, pen_old_T, step;
    UInt i;

    std::tie(loss, grad, llik, pen_S, pen_T) = this->funcProblem_.computeFunctional_g(g_curr, lambda_S, lambda_T, Upsilon);
    norm_grad = std::sqrt(grad.dot(grad));

    if(this->dataProblem_.Print()) {
        Rprintf("loss %f, llik %f, pen_S %f, pen_T %f, norm_Lp %f\n", loss, llik, pen_S, pen_T, norm_grad);
    }

    for(i = 0; i < this->dataProblem_.getNsimulations() && (dloss > toll1 || dllik > toll1 || dpen_S > toll1 || dpen_T > toll1) && norm_grad > toll2; ++i) {

        // Termination criteria variables
        loss_old = loss;
        llik_old = llik;
        pen_old_S = pen_S;
        pen_old_T = pen_T;

        // Compute a descent direction
        d = this->direction_->computeDirection(g_curr, grad);

        // Compute a step
        step = computeStep(g_curr, loss, grad, d, lambda_S, lambda_T, Upsilon);

        // Update the point
        g_curr = g_curr + step*d;

        // Update termination criteria variables
        std::tie(loss, grad, llik, pen_S, pen_T) = this->funcProblem_.computeFunctional_g(g_curr, lambda_S, lambda_T, Upsilon);
        dloss = std::abs((loss - loss_old)/loss_old);
        dllik = std::abs((llik - llik_old)/llik_old);
        dpen_S = std::abs((pen_S - pen_old_S)/pen_old_S);
        dpen_T = std::abs((pen_T - pen_old_T)/pen_old_T);
        norm_grad = std::sqrt(grad.dot(grad));

        if(this->dataProblem_.Print()) {
            Rprintf("Iter %d, loss %f, llik %f, pen_S %f, pen_T %f, norm_Lp %f\n", i+1, loss, llik, pen_S, pen_T, norm_grad);
        }

    }

    this->direction_->resetParameters();

    if(dloss <= toll1 && dllik <= toll1 && dpen_S <= toll1 && dpen_T <= toll1) {
        if(this->dataProblem_.Print()) {
            Rprintf("The algorithm reaches the tolerance in terms of the functional. Norm of Lp: %f, dloss: %f, dllik: %f, dpen_S: %f, dpen_T: %f\n", norm_grad, dloss, dllik, dpen_S, dpen_T);
        }
        return g_curr;
    }
    else if(norm_grad <= toll2) {
        if(this->dataProblem_.Print()){
            Rprintf("The algorithm reaches the tolerance in terms of the slope. Norm of Lp: %f, dloss: %f, dllik: %f, dpen_S: %f, dpen_T: %f\n", norm_grad, dloss, dllik, dpen_S, dpen_T);
        }
        return g_curr;
    }
    else{
        if(this->dataProblem_.Print()) {
            Rprintf("The algorithm reaches the maximum number of iterations. Norm of Lp: %f, dloss: %f, dllik: %f, dpen_S: %f, dpen_T: %f\n", norm_grad, dloss, dllik, dpen_S, dpen_T);
        }
        return g_curr;
    }
}

template<UInt ORDER, UInt mydim, UInt ndim>
Real
BacktrackingMethod_time<ORDER,mydim,ndim>::computeStep(const VectorXr& g, Real loss, const VectorXr& grad,
                                                       const VectorXr& dir, Real lambda_S, Real lambda_T,
                                                       const SpMat& Upsilon) const
{

    Real ro = 0.5, alpha = 1/ro, c = 0.5;

    Real loss_new, llik_new, pen_new_S, pen_new_T, slope, grad_dir;
    VectorXr grad_new, new_point;

    grad_dir = grad.dot(dir);

    do{
        // Update step
        alpha *= ro;

        slope = c*alpha*(grad_dir);

        // Update the point
        new_point = g + alpha*dir;

        // Functional in the new point
        std::tie(loss_new, grad_new, llik_new, pen_new_S, pen_new_T) = this->funcProblem_.computeFunctional_g(new_point, lambda_S, lambda_T, Upsilon);

    } while(loss_new > (loss + slope));

    return alpha;
}


template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> BacktrackingMethod_time<ORDER, mydim, ndim>::clone() const
{

    return fdaPDE::make_unique<BacktrackingMethod_time<ORDER, mydim, ndim>>(*this);

}

template<UInt ORDER, UInt mydim, UInt ndim>
Real WolfeMethod_time<ORDER,mydim,ndim>::computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir,
                                                     Real lambda_S, Real lambda_T, const SpMat& Upsilon) const
{

    Real alpha = 1, alphamax = 0, alphamin = 0, c1 = 1e-4, c2 = 0.9;

    Real loss_new, llik_new, pen_new_S, pen_new_T, slope, grad_dir;
    VectorXr grad_new, new_point;

    grad_dir = grad.dot(dir);
    slope = c1*alpha*grad_dir;


    // Update the point
    new_point = g + alpha*dir;

    // Functional in the new point
    std::tie(loss_new, grad_new, llik_new, pen_new_S, pen_new_T) = this->funcProblem_.computeFunctional_g(new_point, lambda_S, lambda_T, Upsilon);

    bool again = true;

    while(again) {

        again = false;

        while(loss_new > (loss + slope)){
            // Update step
            alphamax = alpha;
            alpha = 0.5*(alphamin + alphamax);

            // Try with the new point
            new_point = g + alpha*dir;
            std::tie(loss_new, grad_new, llik_new, pen_new_S, pen_new_T) = this->funcProblem_.computeFunctional_g(new_point, lambda_S, lambda_T, Upsilon);
            slope = c1*alpha*grad_dir;
        }

        if(grad_new.dot(dir) < c2*grad_dir && std::abs(grad_dir) > 1e-2){

            again = true;

            // Update step
            alphamin = alpha;
            alpha = alphamax==0 ? 2*alphamin : 0.5*(alphamin+alphamax);

            // Try with the new point
            new_point = g + alpha*dir;
            std::tie(loss_new, grad_new, llik_new, pen_new_S, pen_new_T) = this->funcProblem_.computeFunctional_g(new_point, lambda_S, lambda_T, Upsilon);
            slope =  alpha*c1*grad_dir;
        }
    }

    return alpha;
}


template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> WolfeMethod_time<ORDER, mydim, ndim>::clone() const
{

    return fdaPDE::make_unique<WolfeMethod_time<ORDER, mydim, ndim>>(*this);

}


#endif
