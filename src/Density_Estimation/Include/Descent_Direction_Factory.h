#ifndef __DESCENT_DIRECTION_FACTORY_H__
#define __DESCENT_DIRECTION_FACTORY_H__

#include <memory>
#include "../../Global_Utilities/Include/Make_Unique.h"

//! @brief A Factory class: a class for the choice of the direction for the optimization algorithm.
template<UInt ORDER, UInt mydim, UInt ndim>
class DescentDirection_factory
{
	public:
	//! A method that builds a pointer to the right object for the direction choice, taking as parameters a string and others objects needed for constructor.
	static std::unique_ptr<DirectionBase<ORDER,  mydim,  ndim>>
    createDirectionSolver(const DataProblem<ORDER, mydim, ndim>& dp,
    const FunctionalProblem<ORDER, mydim, ndim>& fp, const std::string& d)
	{
		if (d=="Gradient")
			return fdaPDE::make_unique<DirectionGradient<ORDER,mydim,ndim>>(fp);
		else if (d=="BFGS")
			return fdaPDE::make_unique<DirectionBFGS<ORDER,mydim,ndim>>(fp, dp.getNumNodes());
		else{

			Rprintf("Unknown direction option - using gradient direction");

			return fdaPDE::make_unique<DirectionGradient<ORDER,mydim,ndim>>(fp);
		}
	}

};

//! @brief A Factory class: a class for the choice of the direction for the optimization algorithm (spatio-temporal setting).
template<UInt ORDER, UInt mydim, UInt ndim>
class DescentDirection_factory_time{
public:
    //! A method that builds a pointer to the right object for the direction choice, taking as parameters a string and others objects needed for constructor.
    static std::unique_ptr<DirectionBase<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>
    createDirectionSolver(const DataProblem_time<ORDER, mydim, ndim>& dp,
                          const FunctionalProblem_time<ORDER, mydim, ndim>& fp, const std::string& d)
    {
        if (d=="Gradient")
            return fdaPDE::make_unique<DirectionGradient<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp);
        else if (d=="ConjugateGradientFR")
            return fdaPDE::make_unique<DirectionConjugateGradient<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp, 0);
        else if (d=="ConjugateGradientPRP")
            return fdaPDE::make_unique<DirectionConjugateGradient<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp, 1);
        else if (d=="ConjugateGradientHS")
            return fdaPDE::make_unique<DirectionConjugateGradient<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp, 2);
        else if (d=="ConjugateGradientDY")
            return fdaPDE::make_unique<DirectionConjugateGradient<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp, 3);
        else if (d=="ConjugateGradientCD")
            return fdaPDE::make_unique<DirectionConjugateGradient<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp, 4);
        else if (d=="ConjugateGradientLS")
            return fdaPDE::make_unique<DirectionConjugateGradient<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp, 5);
        else if (d=="BFGS")
            return fdaPDE::make_unique<DirectionBFGS<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp, dp.getNumNodes()*dp.getSplineNumber());
        else if (d=="L-BFGS5")
            return fdaPDE::make_unique<DirectionLBFGS<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp, 5);
        else if (d=="L-BFGS10")
            return fdaPDE::make_unique<DirectionLBFGS<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp, 10);
        else{
            Rprintf("Unknown direction option - using gradient direction");
            return fdaPDE::make_unique<DirectionGradient<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp);
        }
    }
};

#endif
