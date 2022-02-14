#ifndef __OPTIMIZATION_ALGORITHM_H__
#define __OPTIMIZATION_ALGORITHM_H__

#include <memory>
#include "../../Global_Utilities/Include/Make_Unique.h"
#include "Descent_Direction.h"
#include "Descent_Direction_Factory.h"

// This file contains info of the optimization algorithm of the Density Estimation problem

//! @brief An abstract base class to perform the minimization algorithm.
template<UInt ORDER, UInt mydim, UInt ndim>
class MinimizationAlgorithm{
  protected:
    // A member to access data problem methods
    const DataProblem<ORDER, mydim, ndim>& dataProblem_;
    // A member to access functional  methods
    const FunctionalProblem<ORDER, mydim, ndim>& funcProblem_;
    // A pointer to the object which computes the descent direction
    std::unique_ptr<DirectionBase<ORDER, mydim, ndim>> direction_;

  public:
    //! A constructor.
    MinimizationAlgorithm(const DataProblem<ORDER, mydim, ndim>& dp, const FunctionalProblem<ORDER, mydim, ndim>& fp,
                          const std::string& d);
    //! A destructor.
    virtual ~MinimizationAlgorithm(){};
    //! A copy constructor.
    MinimizationAlgorithm(const MinimizationAlgorithm<ORDER, mydim, ndim>& rhs);
    //! A pure virtual clone method.
    virtual std::unique_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> clone() const = 0;
    //! A pure virtual method to perform the minimization task.
    virtual VectorXr apply_core(const SpMat& Psi, Real lambda, const VectorXr& g) const = 0;

};


//! @brief A class to perform the minimization algorithm when the step parameter is fixed among all the iterations.
template<UInt ORDER, UInt mydim, UInt ndim>
class FixedStep : public MinimizationAlgorithm<ORDER, mydim, ndim>{
  public:
    //! A delegating constructor.
    FixedStep(const DataProblem<ORDER, mydim, ndim>& dp, const FunctionalProblem<ORDER, mydim, ndim>& fp,
              const std::string& d):
      MinimizationAlgorithm<ORDER, mydim, ndim>(dp, fp, d){};

    //! A copy constructor.
    FixedStep(const FixedStep<ORDER, mydim, ndim>& rhs):
      MinimizationAlgorithm<ORDER, mydim, ndim>(rhs){};
    //! Clone method overridden.
    std::unique_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> clone() const override;
    //! A method to perform the minimization algorithm when the step parameter is fixed among all the iterations.
    VectorXr apply_core(const SpMat& Psi, Real lambda, const VectorXr& g) const override;

};


//! @brief An abstract class to perform the minimization algorithm when the step is computed for each iteration.
template<UInt ORDER, UInt mydim, UInt ndim>
class AdaptiveStep : public MinimizationAlgorithm<ORDER, mydim, ndim>{
  protected:
    //! A copy constructor.
    AdaptiveStep(const AdaptiveStep<ORDER, mydim, ndim>& rhs):
      MinimizationAlgorithm<ORDER, mydim, ndim>(rhs){};
    //! A pure virtual method to compute the step.
    virtual Real computeStep (const VectorXr& g,  Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda,
                              const SpMat& Psi) const = 0;

  public:
    //! A delegating constructor.
    AdaptiveStep(const DataProblem<ORDER, mydim, ndim>& dp, const FunctionalProblem<ORDER, mydim, ndim>& fp,
                 const std::string& d):
      MinimizationAlgorithm<ORDER, mydim, ndim>(dp, fp, d){};
    //! A pure virtual clone method.
    virtual std::unique_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> clone() const = 0;
    //! A method to perform the minimization algorithm when the step is computed for each iteration.
    VectorXr apply_core(const SpMat& Psi, Real lambda, const VectorXr& g) const override;

};


//! @brief A class to handle the Backtracking Method.
template<UInt ORDER, UInt mydim, UInt ndim>
class BacktrackingMethod : public AdaptiveStep<ORDER, mydim, ndim>{
  private:
    //! A method to compute the step using the Backtracking Method.
    Real computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda,
                     const SpMat& Psi) const override;
  public:
    //! A delegating constructor.
    BacktrackingMethod(const DataProblem<ORDER, mydim, ndim>& dp, const FunctionalProblem<ORDER, mydim, ndim>& fp,
                       const std::string& d):
      AdaptiveStep<ORDER, mydim, ndim>(dp, fp, d){};

    //! A copy constructor.
    BacktrackingMethod(const BacktrackingMethod<ORDER, mydim, ndim>& rhs):
      AdaptiveStep<ORDER, mydim, ndim>(rhs){};
    //! Clone method overridden.
    std::unique_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> clone() const override;

};


//! @brief A class to handle the Wolfe Method.
template<UInt ORDER, UInt mydim, UInt ndim>
class WolfeMethod : public AdaptiveStep<ORDER, mydim, ndim>{
  private:
    //! A method to compute the step using the Wolfe Method.
    Real computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda,
                     const SpMat& Psi) const override;
  public:
    //! A delegating constructor.
    WolfeMethod(const DataProblem<ORDER, mydim, ndim>& dp, const FunctionalProblem<ORDER, mydim, ndim>& fp,
                const std::string& d):
      AdaptiveStep<ORDER, mydim, ndim>(dp, fp, d){};

    //! A copy constructor.
    WolfeMethod(const WolfeMethod<ORDER, mydim, ndim>& rhs):
      AdaptiveStep<ORDER, mydim, ndim>(rhs){};
    //! Clone method overridden.
    std::unique_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> clone() const override;

};


//! @brief An abstract base class to perform the minimization algorithm (spatio-temporal setting).
template<UInt ORDER, UInt mydim, UInt ndim>
class MinimizationAlgorithm_time{
protected:
    // A member to access data problem methods
    const DataProblem_time<ORDER, mydim, ndim>& dataProblem_;
    // A member to access functional  methods
    const FunctionalProblem_time<ORDER, mydim, ndim>& funcProblem_;
    // A pointer to the object which computes the descent direction
    std::unique_ptr<DirectionBase<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>> direction_;

public:
    //! A constructor.
    MinimizationAlgorithm_time(const DataProblem_time<ORDER, mydim, ndim>& dp,
                               const FunctionalProblem_time<ORDER, mydim, ndim>& fp, const std::string& d);
    //! A destructor.
    virtual ~MinimizationAlgorithm_time(){};
    //! A copy constructor.
    MinimizationAlgorithm_time(const MinimizationAlgorithm_time<ORDER, mydim, ndim>& rhs);
    //! A pure virtual clone method.
    virtual std::unique_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> clone() const = 0;
    //! A pure virtual method to perform the minimization task.
    virtual VectorXr apply_core(const SpMat& Upsilon, Real lambda_S, Real lambda_T, const VectorXr& g) const = 0;

};


//! @brief A class to perform the minimization algorithm when the step parameter is fixed among all the iterations
//! (spatio-temporal setting).
template<UInt ORDER, UInt mydim, UInt ndim>
class FixedStep_time : public MinimizationAlgorithm_time<ORDER, mydim, ndim>{
public:
    //! A delegating constructor.
    FixedStep_time(const DataProblem_time<ORDER, mydim, ndim>& dp, const FunctionalProblem_time<ORDER, mydim, ndim>& fp,
                   const std::string& d):
      MinimizationAlgorithm_time<ORDER, mydim, ndim>(dp, fp, d){};

    //! A copy constructor.
    FixedStep_time(const FixedStep_time<ORDER, mydim, ndim>& rhs):
      MinimizationAlgorithm_time<ORDER, mydim, ndim>(rhs){};
    //! Clone method overridden.
    std::unique_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> clone() const override;
    //! A method to perform the minimization algorithm when the step parameter is fixed among all the iterations.
    VectorXr apply_core(const SpMat& Upsilon, Real lambda_S, Real lambda_T, const VectorXr& g) const override;

};

//! @brief An abstract class to perform the minimization algorithm when the step is computed for each iteration
//! (spatio-temporal setting).
template<UInt ORDER, UInt mydim, UInt ndim>
class AdaptiveStep_time : public MinimizationAlgorithm_time<ORDER, mydim, ndim>{
protected:
    //! A copy constructor.
    AdaptiveStep_time(const AdaptiveStep_time<ORDER, mydim, ndim>& rhs):
      MinimizationAlgorithm_time<ORDER, mydim, ndim>(rhs){};
    //! A pure virtual method to compute the step.
    virtual Real computeStep (const VectorXr& g,  Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda_S,
                              Real lambda_T, const SpMat& Upsilon) const = 0;

public:
    //! A delegating constructor.
    AdaptiveStep_time(const DataProblem_time<ORDER, mydim, ndim>& dp, const FunctionalProblem_time<ORDER, mydim, ndim>& fp,
                      const std::string& d):
      MinimizationAlgorithm_time<ORDER, mydim, ndim>(dp, fp, d){};
    //! A pure virtual clone method.
    virtual std::unique_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> clone() const = 0;
    //! A method to perform the minimization algorithm when the step is computed for each iteration.
    VectorXr apply_core(const SpMat& Upsilon, Real lambda_S, Real lambda_T, const VectorXr& g) const override;

};


//! @brief A class to handle the Backtracking Method (spatio-temporal setting).
template<UInt ORDER, UInt mydim, UInt ndim>
class BacktrackingMethod_time : public AdaptiveStep_time<ORDER, mydim, ndim>{
private:
    //! A method to compute the step using the Backtracking Method.
    Real computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda_S, Real lambda_T,
                     const SpMat& Upsilon) const override;
public:
    //! A delegating constructor.
    BacktrackingMethod_time(const DataProblem_time<ORDER, mydim, ndim>& dp,
                            const FunctionalProblem_time<ORDER, mydim, ndim>& fp,
                            const std::string& d):
            AdaptiveStep_time<ORDER, mydim, ndim>(dp, fp, d){};

    //! A copy constructor.
    BacktrackingMethod_time(const BacktrackingMethod_time<ORDER, mydim, ndim>& rhs):
            AdaptiveStep_time<ORDER, mydim, ndim>(rhs){};
    //! Clone method overridden.
    std::unique_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> clone() const override;

};


//! @brief A class to handle the Wolfe Method (spatio-temporal setting).
template<UInt ORDER, UInt mydim, UInt ndim>
class WolfeMethod_time : public AdaptiveStep_time<ORDER, mydim, ndim>{
private:
    //! A method to compute the step using the Wolfe Method.
    Real computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda_S, Real lambda_T,
                     const SpMat& Upsilon) const override;
public:
    //! A delegating constructor.
    WolfeMethod_time(const DataProblem_time<ORDER, mydim, ndim>& dp,
                     const FunctionalProblem_time<ORDER, mydim, ndim>& fp,
                     const std::string& d):
            AdaptiveStep_time<ORDER, mydim, ndim>(dp, fp, d){};

    //! A copy constructor.
    WolfeMethod_time(const WolfeMethod_time<ORDER, mydim, ndim>& rhs):
            AdaptiveStep_time<ORDER, mydim, ndim>(rhs){};
    //! Clone method overridden.
    std::unique_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> clone() const override;

};


#include "Optimization_Algorithm_imp.h"

#endif
