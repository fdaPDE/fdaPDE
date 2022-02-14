#ifndef __FUNCTIONAL_PROBLEM_H__
#define __FUNCTIONAL_PROBLEM_H__

#include "Data_Problem.h"

// This file implements the functional of the Density Estimation problem

//! @brief A class to store methods regarding the functional of the problem.
template<UInt ORDER, UInt mydim, UInt ndim>
class FunctionalProblem{
  private:
    using Integrator = typename DensityIntegratorHelper::Integrator<mydim>;
    static constexpr UInt EL_NNODES = how_many_nodes(ORDER,mydim);

    // A member to acess data problem methods
    const DataProblem<ORDER, mydim, ndim>& dataProblem_;

    //! A method to compute the integrals of the functional.
    std::pair<Real, VectorXr> computeIntegrals(const VectorXr& g) const;

  public:
    //! A constructor
    FunctionalProblem(const DataProblem<ORDER, mydim, ndim>& dp): dataProblem_(dp){};
    //! A method to compute the functional for the g-function. Output: loss, gradient, llik, penterm.
    std::tuple<Real, VectorXr, Real, Real> computeFunctional_g(const VectorXr& g, Real lambda, const SpMat& Psi) const;
    //! A method to compute the log-likelihood and the penalization term for the f-function.
    std::pair<Real,Real> computeLlikPen_f(const VectorXr& f) const;

};

//! @brief A class to store methods regarding the functional of the spatio-temporal problem.
template<UInt ORDER, UInt mydim, UInt ndim>
class FunctionalProblem_time {
private:
    using Integrator = typename DensityIntegratorHelper::Integrator<mydim>;
    //using Integrator_t = IntegratorGaussP5;
    using Integrator_t = IntegratorGaussP9;
    static constexpr UInt EL_NNODES = how_many_nodes(ORDER, mydim);

    //! A member to access DataProblem_time methods.
    const DataProblem_time<ORDER, mydim, ndim>& dataProblem_time_;

    //! A method to compute the integrals of the functional.
    std::pair<Real, VectorXr> computeIntegrals(const VectorXr& g) const;

public:
    //! A constructor.
    FunctionalProblem_time(const DataProblem_time<ORDER, mydim, ndim>& dp_t) : dataProblem_time_(dp_t){};
    //! A method to compute the functional for the g-function. Output: loss, gradient, llik, penterms.
    std::tuple<Real, VectorXr, Real, Real, Real> computeFunctional_g(const VectorXr& g, Real lambda, Real lambda_T,
                                                                     const SpMat& Upsilon) const;
    //! A method to compute the log-likelihood and the penalization terms for the f-function.
    std::tuple<Real, Real, Real> computeLlikPen_f(const VectorXr& f) const;

};

#include "Functional_Problem_imp.h"

#endif
