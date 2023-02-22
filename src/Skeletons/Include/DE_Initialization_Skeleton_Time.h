#ifndef __DE_INITIALIZATION_SKELETON_TIME_H__
#define __DE_INITIALIZATION_SKELETON_TIME_H__

#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"
#include "../../Global_Utilities/Include/Solver_Definitions.h"

//Density Estimation
#include "../../Density_Estimation/Include/Data_Problem.h"
#include "../../Density_Estimation/Include/Functional_Problem.h"
#include "../../Density_Estimation/Include/Optimization_Algorithm.h"
#include "../../Density_Estimation/Include/Optimization_Algorithm_Factory.h"
#include "../../Density_Estimation/Include/FE_Density_Estimation.h"


template<UInt ORDER, UInt mydim, UInt ndim>
SEXP DE_init_skeleton_time(SEXP Rdata, SEXP Rdata_time, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda,
                           SEXP Rlambda_time, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2, SEXP Rprint,
                           SEXP Rmesh, SEXP Rmesh_time, SEXP Rsearch, SEXP RisTimeDiscrete, SEXP RflagMass, SEXP RflagLumped,
                           const std::string& init, UInt init_fold)
{

    // Convert Rmesh_time into a vector
    UInt diml = Rf_length(Rmesh_time);
    std::vector<Real> mesh_time;
    mesh_time.reserve(diml);
    for(UInt i = 0; i < diml; ++i)
    {
        mesh_time.push_back(REAL(Rmesh_time)[i]);
    }

    // Construct data problem object
    DataProblem_time<ORDER, mydim, ndim> dataProblem(Rdata, Rdata_time, Rorder, Rfvec, RheatStep, RheatIter, Rlambda,
                                                     Rlambda_time, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint,
                                                     Rsearch, Rmesh, mesh_time, RisTimeDiscrete, RflagMass, RflagLumped);

    // Construct functional problem object
    FunctionalProblem_time<ORDER, mydim, ndim> functionalProblem(dataProblem);

    if(init == "Heat"){

        // Construct densityInit object
        std::unique_ptr<DensityInitialization_time<ORDER, mydim, ndim>> densityInit = fdaPDE::make_unique<HeatProcess_time<ORDER, mydim, ndim>>(dataProblem, functionalProblem);

        // Fill fInit
        std::vector<VectorXr> fInit(dataProblem.getNlambda() * dataProblem.getNlambda_time());
        for(UInt l = 0; l < dataProblem.getNlambda(); ++l){
            for(UInt k = 0; k < dataProblem.getNlambda_time(); ++k){
                fInit[k+l*dataProblem.getNlambda_time()] = *(densityInit->chooseInitialization(dataProblem.getLambda(l), dataProblem.getLambda_time(k)));
            }
        }

        // Copy result in R memory
        SEXP result = NILSXP;
        result = PROTECT(Rf_allocVector(VECSXP, 1));
        SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, ((fInit[0])).size(), fInit.size()));

        Real *rans = REAL(VECTOR_ELT(result, 0));
        for(UInt j = 0; j < fInit.size(); ++j)
        {
            for(UInt i = 0; i < (fInit[0]).size(); ++i)
                rans[i + (fInit[0]).size()*j] = (fInit[j])[i];
        }

        UNPROTECT(1);

        return(result);
    }

    else if(init=="CV"){

        // Construct densityInit object
        std::unique_ptr<Heat_CV_time<ORDER, mydim, ndim>> densityInit = fdaPDE::make_unique<Heat_CV_time<ORDER, mydim, ndim>>(dataProblem, functionalProblem, init_fold);

        // Fill fInit
        VectorXr fInit;
        fInit = *(densityInit->chooseInitialization(0, 0));

        // Copy result in R memory
        SEXP result = NILSXP;
        result = PROTECT(Rf_allocVector(VECSXP, 1));
        SET_VECTOR_ELT(result, 0, Rf_allocVector(REALSXP, fInit.size()));

        Real *rans = REAL(VECTOR_ELT(result, 0));
        for(UInt i = 0; i < fInit.size(); ++i)
        {
            rans[i] = fInit[i];
        }

        UNPROTECT(1);

        return(result);
    }
    else{

        #ifdef R_VERSION_
        Rprintf("Invalid initialization");
        #endif

        return NILSXP;
    }

}


#endif
