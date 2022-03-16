#ifndef __DESCENT_DIRECTION_H__
#define __DESCENT_DIRECTION_H__

#include "../../Global_Utilities/Include/Make_Unique.h"

// This file contains the direction search technique useful for the optimization algorithm of the Density Estimation problem

/*! @brief An abstract class for computing the descent direction. The father is pure virtual; the right direction
is computed inside the children according to a proper method chosen.
*/
template<UInt ORDER, UInt mydim, UInt ndim, class T = FunctionalProblem<ORDER, mydim, ndim>>
class DirectionBase{
  protected:
    // to give generality if you want to add other children
    // T is a typename for FunctionalProblem (spatial problem, by default) or for FunctionalProblem_time (spatio-temporal problem)
    const T& funcProblem_;

  public:
    //! A constructor
    DirectionBase(const T& fp): funcProblem_(fp){};
    //! A destructor.
    virtual ~DirectionBase(){};
    //! A pure virtual clone method.
    virtual std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>> clone() const = 0;
    //! A pure virtual method to compute the descent direction.
    virtual VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) = 0;
    //! A pure virtual method to reset all the old parameters.
    virtual void resetParameters() = 0;
};


//! @brief A class for computing the gradient descent direction.
template<UInt ORDER, UInt mydim, UInt ndim, class T = FunctionalProblem<ORDER, mydim, ndim>>
class DirectionGradient : public DirectionBase<ORDER, mydim, ndim, T>{
  public:
    //! A delegating constructor.
    DirectionGradient(const T& fp):
    DirectionBase<ORDER, mydim, ndim, T>(fp){};
    //! Clone method overridden.
    std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>> clone() const override;
    //! A method to compute the gradient descent direction.
    VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) override;
    //! A method to reset all the old parameters. In the gradient method they aren't.
    void resetParameters() override {};

};


//! @brief A class for computing the conjugate gradient descent direction.
template<UInt ORDER, UInt mydim, UInt ndim, class T = FunctionalProblem<ORDER, mydim, ndim>>
class DirectionConjugateGradient : public DirectionBase<ORDER, mydim, ndim, T>{
private:
    // Old gradient and old direction
    VectorXr gradOld_, directionOld_;
    // Flag for the first iteration
    bool first_iteration_;
    // Flag to choose the update formula for beta
    UInt betaFormula_;
    // Counter for restart
    UInt restart;
public:
    //! A delegating constructor.
    DirectionConjugateGradient(const T& fp, const UInt beta):
            DirectionBase<ORDER, mydim, ndim, T>(fp), betaFormula_(beta), first_iteration_(true), restart(0) {};
    //! A copy constructor: it just creates a DirectionConjugateGradient object with the same features of rhs but it initializes the flag.
    DirectionConjugateGradient(const DirectionConjugateGradient<ORDER, mydim, ndim, T>& rhs);
    //! Clone method overridden.
    std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>> clone() const override;
    //! A method to compute the conjugate gradient descent direction.
    VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) override;
    //! A method to reset all the old parameters. In the conjugate gradient method they aren't.
    void resetParameters() override;

};


//! @brief A class for computing the BFGS descent direction.
template<UInt ORDER, UInt mydim, UInt ndim, class T = FunctionalProblem<ORDER, mydim, ndim>>
class DirectionBFGS : public DirectionBase<ORDER, mydim, ndim, T>{
  private:
    MatrixXr HInit_, HOld_;
    VectorXr gOld_, gradOld_;
    bool updateH_;

  public:
    //! A constructor.
    DirectionBFGS(const T& fp, UInt k):
    DirectionBase<ORDER, mydim, ndim, T>(fp), HInit_(MatrixXr::Identity(k, k)), HOld_(MatrixXr::Identity(k, k)), updateH_(false){};
    //! A copy constructor: it just creates a DirectionBFGS object with the same features of rhs but it initializes the matrices.
    DirectionBFGS(const DirectionBFGS<ORDER, mydim, ndim, T>& rhs);
    //! Clone method overridden.
    std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>> clone() const override;
    //! A method to compute the BFGS descent direction.
    VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) override;
    //! A method to reset all the old parameters.
    void resetParameters() override;

};


//! @brief A class for computing the L-BFGS descent direction.
template<UInt ORDER, UInt mydim, UInt ndim, class T = FunctionalProblem<ORDER, mydim, ndim>>
class DirectionLBFGS : public DirectionBase<ORDER, mydim, ndim, T>{
private:
    // Maximum number of correction vectors
    UInt m_;
    // Old g and old gradient
    VectorXr gOld_, gradOld_;
    // gamma_ * I is the initial approximation to the Hessian matrix
    Real gamma_;
    // History of s and y vectors
    std::vector<VectorXr> s_, y_;
    // History of the y's values
    std::vector<Real> ys_;
    // Temporary values used in computing H * v
    std::vector<Real> alpha_;
    // Number of correction vectors in the history, ncorr_ <= m_
    UInt ncorr_;
    // A Pointer to locate the most recent history, 1 <= ptr_ <= m_
    // Details: s and y vectors are stored in cyclic order.
    //          For example, if the current s-vector is stored in s_[m-1],
    //          then in the next iteration s_[0] will be overwritten.
    //          s_[ptr_-1] points to the most recent history,
    //          and s_[ptr_ % m_] points to the most distant one.
    UInt ptr_;
    // Flag for the first iteration
    bool first_iteration_;

public:
    //! A constructor.
    DirectionLBFGS(const T& fp, const UInt m);

    //! A copy constructor: it just creates a DirectionLBFGS object with the same features of rhs but it initializes the members.
    DirectionLBFGS(const DirectionLBFGS<ORDER, mydim, ndim, T>& rhs);
    //! Clone method overridden.
    std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>> clone() const override;
    //! A method to compute the L-BFGS descent direction.
    VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) override;
    //! A method to reset all the old parameters.
    void resetParameters() override;
    //! Add correction vectors
    void add_correction(const VectorXr& s, const VectorXr& y);

};

#include "Descent_Direction_imp.h"
#endif
