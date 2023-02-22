#ifndef __PREPROCESS_FACTORY_H__
#define __PREPROCESS_FACTORY_H__

#include <memory>
#include "../../Global_Utilities/Include/Make_Unique.h"

//! @brief A Factory class: a class for the choice of the cross-validation method.
template<UInt ORDER, UInt mydim, UInt ndim>
class Preprocess_factory {
public:
	//! A method that builds a pointer to the right object for the cross-validation method choice, taking as parameters
	//! a string and others objects needed for constructor.
	static std::unique_ptr<Preprocess<ORDER,  mydim,  ndim>>
    createPreprocessSolver(const DataProblem<ORDER, mydim, ndim>& dp, const FunctionalProblem<ORDER, mydim, ndim>& fp,
                           std::shared_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> ma, const std::string& p)
    {
        if(p=="RightCV")
		    return fdaPDE::make_unique<RightCrossValidation<ORDER, mydim, ndim>>(dp, fp, ma);
        else if(p=="SimplifiedCV")
			return fdaPDE::make_unique<SimplifiedCrossValidation<ORDER, mydim, ndim>>(dp, fp, ma);
		else if(p=="NoCrossValidation")
            return fdaPDE::make_unique<NoCrossValidation<ORDER, mydim, ndim>>(dp, fp);
        else {
            Rprintf("Unknown preprocess option - using right cross validation\n");
            return fdaPDE::make_unique<RightCrossValidation<ORDER, mydim, ndim>>(dp, fp, ma);
        }
    }

};

//! @brief A Factory class: a class for the choice of the cross-validation method (spatio-temporal setting).
template<UInt ORDER, UInt mydim, UInt ndim>
class Preprocess_factory_time {
public:
    //! A method that builds a pointer to the right object for the cross-validation method choice, taking as parameters
    //! a string and others objects needed for constructor.
    static std::unique_ptr<Preprocess_time<ORDER, mydim, ndim>>
    createPreprocessSolver(const DataProblem_time<ORDER, mydim, ndim>& dp, const FunctionalProblem_time<ORDER, mydim, ndim>& fp,
                           std::shared_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> ma, const std::string& p)
    {
        if(p=="RightCV")
            return fdaPDE::make_unique_time<RightCrossValidation_time<ORDER, mydim, ndim>>(dp, fp, ma);
        else if(p=="SimplifiedCV")
            return fdaPDE::make_unique_time<SimplifiedCrossValidation_time<ORDER, mydim, ndim>>(dp, fp, ma);
        else if(p=="NoCrossValidation")
            return fdaPDE::make_unique_time<NoCrossValidation_time<ORDER, mydim, ndim>>(dp, fp);
        else {
            Rprintf("Unknown preprocess option - using right cross validation\n");
            return fdaPDE::make_unique_time<RightCrossValidation_time<ORDER, mydim, ndim>>(dp, fp, ma);
        }

    }

};

#endif
