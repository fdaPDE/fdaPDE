#ifndef __DESCENT_DIRECTION_IMP_H__
#define __DESCENT_DIRECTION_IMP_H__

template<UInt ORDER, UInt mydim, UInt ndim, typename T>
std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>>
DirectionGradient<ORDER, mydim, ndim, T>::clone() const {

  return make_unique<DirectionGradient<ORDER, mydim, ndim, T>>(*this);

}

template<UInt ORDER, UInt mydim, UInt ndim, typename T>
VectorXr
DirectionGradient<ORDER, mydim, ndim, T>::computeDirection(const VectorXr& g, const VectorXr& grad){

  return (- grad);
}

template<UInt ORDER, UInt mydim, UInt ndim, typename T>
DirectionBFGS<ORDER, mydim, ndim, T>::DirectionBFGS(const DirectionBFGS<ORDER, mydim, ndim, T>& rhs):
DirectionBase<ORDER, mydim, ndim, T>(rhs) {

  updateH_ = false;
  HInit_ = rhs.HInit_;
  HOld_ = rhs.HInit_;

}

template<UInt ORDER, UInt mydim, UInt ndim, typename T>
std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>>
DirectionBFGS<ORDER, mydim, ndim, T>::clone() const {

  return make_unique<DirectionBFGS<ORDER, mydim, ndim, T>>(*this);

}

template<UInt ORDER, UInt mydim, UInt ndim, typename T>
VectorXr
DirectionBFGS<ORDER, mydim, ndim, T>::computeDirection(const VectorXr& g, const VectorXr& grad){

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

template<UInt ORDER, UInt mydim, UInt ndim, typename T>
void
DirectionBFGS<ORDER, mydim, ndim, T>::resetParameters(){
  updateH_ = false;
  HOld_ = HInit_;
}

#endif
