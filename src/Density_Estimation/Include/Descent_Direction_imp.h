#ifndef __DESCENT_DIRECTION_IMP_H__
#define __DESCENT_DIRECTION_IMP_H__

template<UInt ORDER, UInt mydim, UInt ndim, class T>
std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>> DirectionGradient<ORDER, mydim, ndim, T>::clone() const {

    return fdaPDE::make_unique<DirectionGradient<ORDER, mydim, ndim, T>>(*this);

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
VectorXr DirectionGradient<ORDER, mydim, ndim, T>::computeDirection(const VectorXr& g, const VectorXr& grad){

    return (- grad);

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
DirectionConjugateGradient<ORDER, mydim, ndim, T>::DirectionConjugateGradient(const DirectionConjugateGradient<ORDER, mydim, ndim, T> &rhs):
        DirectionBase<ORDER, mydim, ndim, T>(rhs) {

    first_iteration_ = true;
    restart = 0;
    betaFormula_ = rhs.betaFormula_;

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>> DirectionConjugateGradient<ORDER, mydim, ndim, T>::clone() const {

    return fdaPDE::make_unique<DirectionConjugateGradient<ORDER, mydim, ndim, T>>(*this);

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
void
DirectionConjugateGradient<ORDER, mydim, ndim, T>::resetParameters() {

    first_iteration_ = true;
    restart = 0;

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
VectorXr DirectionConjugateGradient<ORDER, mydim, ndim, T>::computeDirection(const VectorXr& g, const VectorXr& grad){

    // Moving direction
    VectorXr direction_;

    if(first_iteration_){
        direction_ = -grad;
        first_iteration_ = false;
    }
    else {
        // Compute beta
        Real beta;
        if (restart % 10 == 0){
            beta = 0.;
        }
        else if(betaFormula_ == 0){
            // Fletcher-Reeves formula
            beta = (grad.squaredNorm()) / (gradOld_.squaredNorm());
        }
        else if(betaFormula_ == 1){
            // Polak-Ribiére-Polyak formula
            beta = (grad.dot(grad - gradOld_)) / (gradOld_.squaredNorm());
        }
        else if(betaFormula_ == 2){
            // Hestenes-Stiefel formula
            beta = (grad.dot(grad - gradOld_)) / (directionOld_.dot(grad - gradOld_));
        }
        else if(betaFormula_ == 3){
            // Dai-Yuan formula
            beta = (grad.squaredNorm()) / (directionOld_.dot(grad - gradOld_));
        }
        else if(betaFormula_ == 4){
            // Conjugate-Descent formula
            beta = - (grad.squaredNorm()) / (directionOld_.dot(gradOld_));
        }
        else if(betaFormula_ == 5){
            // Liu-Storey formula
            beta = - (grad.dot(grad - gradOld_)) / (directionOld_.dot(gradOld_));
        }

        // Update the direction
        direction_ = -grad + beta * directionOld_;

    }

    gradOld_ = grad;
    directionOld_ = direction_;
    ++restart;

    return (direction_);

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
DirectionBFGS<ORDER, mydim, ndim, T>::DirectionBFGS(const DirectionBFGS<ORDER, mydim, ndim, T>& rhs):
DirectionBase<ORDER, mydim, ndim, T>(rhs) {

    updateH_ = false;
    HInit_ = rhs.HInit_;
    HOld_ = rhs.HInit_;

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>> DirectionBFGS<ORDER, mydim, ndim, T>::clone() const {

    return fdaPDE::make_unique<DirectionBFGS<ORDER, mydim, ndim, T>>(*this);

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
VectorXr DirectionBFGS<ORDER, mydim, ndim, T>::computeDirection(const VectorXr& g, const VectorXr& grad){

    if(updateH_){
        const VectorXr delta = g - gOld_;
        const VectorXr gamma = grad - gradOld_;

        const Real dg = delta.dot(gamma);
        const VectorXr Hg = HOld_*gamma;

        HOld_ = HOld_ + (1 + (gamma.dot(Hg))/dg)*(delta*delta.transpose())/dg -
                (Hg*delta.transpose() + delta*Hg.transpose())/dg;
    }

    gOld_ = g;
    gradOld_ = grad;

    if(!updateH_) updateH_ = true;

    return (-HOld_*grad);

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
void DirectionBFGS<ORDER, mydim, ndim, T>::resetParameters(){

    updateH_ = false;
    HOld_ = HInit_;

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
DirectionLBFGS<ORDER, mydim, ndim, T>::DirectionLBFGS(const T& fp, const UInt m):
DirectionBase<ORDER, mydim, ndim, T>(fp), m_(m), first_iteration_(true) {

    s_.resize(m);
    y_.resize(m);
    ys_.resize(m);
    alpha_.resize(m);
    ncorr_ = 0;
    ptr_ = m_; // This makes sure that ptr_ % m_ == 0 in the first step

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
DirectionLBFGS<ORDER, mydim, ndim, T>::DirectionLBFGS(const DirectionLBFGS<ORDER, mydim, ndim, T>& rhs):
DirectionBase<ORDER, mydim, ndim, T>(rhs) {

    first_iteration_ = true;
    m_ = rhs.m_;
    ncorr_ = 0;
    ptr_ = rhs.m_;
    s_.resize(rhs.m_);
    y_.resize(rhs.m_);
    ys_.resize(rhs.m_);
    alpha_.resize(rhs.m_);

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>> DirectionLBFGS<ORDER, mydim, ndim, T>::clone() const {

    return fdaPDE::make_unique<DirectionLBFGS<ORDER, mydim, ndim, T>>(*this);

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
VectorXr DirectionLBFGS<ORDER, mydim, ndim, T>::computeDirection(const VectorXr& g, const VectorXr& grad){

    // Moving direction
    VectorXr direction_;

    if(first_iteration_){
        direction_ = grad;
        first_iteration_ = false;
    }
    else{
        // Update s and y
        // s_{k+1} = x_{k+1} - x_k
        // y_{k+1} = g_{k+1} - g_k
        add_correction(g-gOld_, grad-gradOld_);

        // Compute the direction
        // Loop 1
        direction_ = grad;
        UInt j = ptr_ % m_;
        for(UInt i = 0; i < ncorr_; ++i)
        {
            j = (j + m_ - 1) % m_;
            alpha_[j] = (s_[j].dot(direction_)) / ys_[j];
            direction_ -= alpha_[j] * y_[j];
        }

        // Apply initial H0
        direction_ *= gamma_;

        // Loop 2
        for(UInt i = 0; i < ncorr_; ++i)
        {
            const Real beta = (y_[j].dot(direction_)) / ys_[j];
            direction_ += (alpha_[j] - beta) * s_[j];
            j = (j + 1) % m_;
        }

    }

    gOld_ = g;
    gradOld_ = grad;

    return (-direction_);

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
void DirectionLBFGS<ORDER, mydim, ndim, T>::add_correction(const VectorXr& s, const VectorXr& y) {

    const UInt loc = ptr_ % m_;

    s_[loc] = s;
    y_[loc] = y;

    // ys = y's = 1/rho
    const Real ys = s_[loc].dot(y_[loc]);
    ys_[loc] = ys;

    gamma_ = ys / y.squaredNorm();

    if(ncorr_ < m_)
        ++ncorr_;

    ptr_ = loc + 1;

}


template<UInt ORDER, UInt mydim, UInt ndim, class T>
void DirectionLBFGS<ORDER, mydim, ndim, T>::resetParameters(){

    first_iteration_ = true;
    ncorr_ = 0;
    ptr_ = m_;

}

#endif
