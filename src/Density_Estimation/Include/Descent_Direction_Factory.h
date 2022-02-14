#ifndef __DESCENT_DIRECTION_FACTORY_H__
#define __DESCENT_DIRECTION_FACTORY_H__

#include <memory>
#include "../../Global_Utilities/Include/Make_Unique.h"

//! brief@ A Factory class: a class for the choice of the step method for the optimization algorithm.
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
			return make_unique<DirectionGradient<ORDER,mydim,ndim>>(fp);
		else if (d=="BFGS")
			return make_unique<DirectionBFGS<ORDER,mydim,ndim>>(fp, dp.getNumNodes());
		else{

			Rprintf("Unknown direction option - using gradient direction");

			return make_unique<DirectionGradient<ORDER,mydim,ndim>>(fp);
		}
	}

};

//! brief@ A Factory class: a class for the choice of the step method for the optimization algorithm.
template<UInt ORDER, UInt mydim, UInt ndim>
class DescentDirection_factory_time{
public:
    //! A method that builds a pointer to the right object for the direction choice, taking as parameters a string and others objects needed for constructor.
    static std::unique_ptr<DirectionBase<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>
    createDirectionSolver(const DataProblem_time<ORDER, mydim, ndim>& dp,
                          const FunctionalProblem_time<ORDER, mydim, ndim>& fp, const std::string& d)
    {
        if (d=="Gradient") {
            std::cout << "Gradient direction" << std::endl;
            return make_unique<DirectionGradient<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp);
        }
        else if (d=="BFGS") {
            std::cout << "BFGS direction" << std::endl;
            return make_unique<DirectionBFGS<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp, dp.getNumNodes()*dp.getSplineNumber());
        }
        else{
            Rprintf("Unknown direction option - using gradient direction");
            return make_unique<DirectionGradient<ORDER, mydim, ndim, FunctionalProblem_time<ORDER, mydim, ndim>>>(fp);
        }
    }
};

#endif
