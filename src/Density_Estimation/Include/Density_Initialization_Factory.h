#ifndef __DENSITY_INITIALIZATION_FACTORY_H__
#define __DENSITY_INITIALIZATION_FACTORY_H__

#include <memory>
#include "../../Global_Utilities/Include/Make_Unique.h"

//! @brief A Factory class: a class for the choice of the density initialization.
template<UInt ORDER, UInt mydim, UInt ndim>
class DensityInitialization_factory{
public:
    //! A method that builds a pointer to the right object for the initialization choice.
	static std::unique_ptr<DensityInitialization<ORDER, mydim, ndim>> createInitializationSolver(const DataProblem<ORDER, mydim, ndim>& dp,
                                                                                                 const FunctionalProblem<ORDER, mydim, ndim>& fp)
    {

        if(!dp.isFvecEmpty())
    		return fdaPDE::make_unique<UserInitialization<ORDER, mydim, ndim>>(dp);
    	else
    		return fdaPDE::make_unique<HeatProcess<ORDER, mydim, ndim>>(dp, fp);

    }

};


//! @brief A Factory class: a class for the choice of the density initialization (spatio-temporal setting).
template<UInt ORDER, UInt mydim, UInt ndim>
class DensityInitialization_factory_time{
public:
    //! A method that builds a pointer to the right object for the initialization choice (spatio-temporal setting).
    static std::unique_ptr<DensityInitialization_time<ORDER,  mydim,  ndim>>
    createInitializationSolver(const DataProblem_time<ORDER, mydim, ndim>& dp,
                               const FunctionalProblem_time<ORDER, mydim, ndim>& fp)
    {

        if(!dp.isFvecEmpty())
            return fdaPDE::make_unique<UserInitialization_time<ORDER, mydim, ndim>>(dp);
        else
            return fdaPDE::make_unique<HeatProcess_time<ORDER, mydim, ndim>>(dp, fp);

    }

};


#endif
