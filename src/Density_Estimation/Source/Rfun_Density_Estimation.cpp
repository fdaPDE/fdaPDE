#include "../../FdaPDE.h"
#include "../../Skeletons/Include/DE_Skeleton.h"
#include "../../Skeletons/Include/DE_Skeleton_Time.h"
#include "../../Skeletons/Include/DE_Initialization_Skeleton.h"
#include "../../Skeletons/Include/DE_Initialization_Skeleton_Time.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include "../../FE_Assemblers_Solvers/Include/Integration.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"
#include "../../Global_Utilities/Include/Solver_Definitions.h"

//Density Estimation
#include "../Include/Data_Problem.h"
#include "../Include/Functional_Problem.h"
#include "../Include/Optimization_Algorithm.h"
#include "../Include/Optimization_Algorithm_Factory.h"
#include "../Include/FE_Density_Estimation.h"


extern "C" {

        //! This function manages the various options for DE-PDE algorithm.
        /*!
        	This function is called from R code.
        	\param Rdata an R-matrix containing the data.
        	\param Rmesh an R-object containing the output mesh from Trilibrary
        	\param Rorder an R-integer containing the order of the approximating basis.
        	\param Rmydim an R-integer containing the dimension of the problem we are considering.
        	\param Rndim an R-integer containing the dimension of the space in which the locations are.
        	\param Rfvec an R-vector containing the initial solution coefficients given by the user.
        	\param RheatStep an R-double containing the step for the heat equation initialization.
        	\param RheatIter an R-integer containing the number of iterations to perform the heat equation initialization.
        	\param Rlambda an R-vector containing the penalization terms.
        	\param Rnfolds an R-integer specifying the number of folds for cross validation.
        	\param Rnsim an R-integer specifying the number of iterations to use in the optimization algorithm.
        	\param RstepProposals an R-vector containing the step parameters useful for the descent algorithm.
        	\param Rtol1 an R-double specifying the tolerance to use for the termination criterion based on the percentage differences.
        	\param Rtol2 an R-double specifying the tolerance to use for the termination criterion based on the norm of the gradient.
        	\param Rprint and R-integer specifying if print on console.
        	\param RstepMethod an R-string containing the method to use to choose the step in the optimization algorithm.
        	\param RdirectionMethod an R-string containing the descent direction to use in the optimization algorithm.
        	\param RpreprocessMethod an R-string containing the cross-validation method to use.
        	\param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).

        	\return R-list containing solutions.
        */

        SEXP Density_Estimation(SEXP Rdata, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP Rfvec, SEXP RheatStep,
                                SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1,
                                SEXP Rtol2, SEXP Rprint, SEXP RstepMethod, SEXP RdirectionMethod, SEXP RpreprocessMethod,
                                SEXP Rsearch)
        {

            UInt order= INTEGER(Rorder)[0];
            UInt mydim=INTEGER(Rmydim)[0];
        	UInt ndim=INTEGER(Rndim)[0];

        	std::string step_method=CHAR(STRING_ELT(RstepMethod, 0));
        	std::string direction_method=CHAR(STRING_ELT(RdirectionMethod, 0));
        	std::string preprocess_method=CHAR(STRING_ELT(RpreprocessMethod, 0));

            if(order== 1 && mydim==2 && ndim==2)
        		return(DE_skeleton<1, 2, 2>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
        	else if(order== 2 && mydim==2 && ndim==2)
        		return(DE_skeleton<2, 2, 2>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
        	else if(order== 1 && mydim==2 && ndim==3)
        		return(DE_skeleton<1, 2, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
        	else if(order== 2 && mydim==2 && ndim==3)
        		return(DE_skeleton<2, 2, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
        	else if(order == 1 && mydim==3 && ndim==3)
        		return(DE_skeleton<1, 3, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
        	else if(order == 2 && mydim==3 && ndim==3)
                return(DE_skeleton<2, 3, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));

        	return(NILSXP);

        }


        SEXP Density_Initialization(SEXP Rdata, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP Rfvec,
                                    SEXP RheatStep, SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim,
                                    SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch, SEXP Rinit,
                                    SEXP Rinit_fold)
        {

        	UInt order = INTEGER(Rorder)[0];
            UInt mydim = INTEGER(Rmydim)[0];
        	UInt ndim = INTEGER(Rndim)[0];

        	UInt init_fold = INTEGER(Rinit_fold)[0];

        	std::string init = CHAR(STRING_ELT(Rinit, 0));

            if(order == 1 && mydim == 2 && ndim == 2)
        		return(DE_init_skeleton<1, 2, 2>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
        	else if(order== 2 && mydim == 2 && ndim == 2)
        		return(DE_init_skeleton<2, 2, 2>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
        	else if(order == 1 && mydim == 2 && ndim == 3)
        		return(DE_init_skeleton<1, 2, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
        	else if(order == 2 && mydim == 2 && ndim == 3)
        		return(DE_init_skeleton<2, 2, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
        	else if(order == 1 && mydim == 3 && ndim == 3)
        		return(DE_init_skeleton<1, 3, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
            else if(order == 2 && mydim == 3 && ndim == 3)
                return(DE_init_skeleton<2, 3, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));

        	return(NILSXP);

        }

        //! This function manages the various options for STDE-PDE algorithm.
        /*!
            This function is called from R code.
            \param Rdata an R-matrix containing the spatial data (locations).
            \param Rdata_time an R-vector containing the temporal data (times).
            \param Rmesh an R-object containing the output spatial mesh from Trilibrary.
            \param Rmesh_time an R-vector containing the temporal mesh (with its nodes in increasing order).
            \param Rorder an R-integer containing the order of the approximating basis.
            \param Rmydim an R-integer containing the dimension of the problem we are considering.
            \param Rndim an R-integer containing the dimension of the space where locations are.
            \param Rfvec an R-vector containing the initial solution coefficients given by the user.
            \param RheatStep an R-double containing the step for the heat equation initialization.
            \param RheatIter an R-integer containing the number of iterations to perform the heat equation initialization.
            \param Rlambda an R-vector containing the penalization terms in space.
            \param Rlambda_time an R-vector containing the penalization terms in time.
            \param Rnfolds an R-integer specifying the number of folds for cross validation.
            \param Rnsim an R-integer specifying the number of iterations to use in the optimization algorithm.
            \param RstepProposals an R-vector containing the step parameters useful for the descent algorithm.
            \param Rtol1 an R-double specifying the tolerance to use for the termination criterion based on the percentage differences.
            \param Rtol2 an R-double specifying the tolerance to use for the termination criterion based on the norm of the gradient.
            \param Rprint and R-integer specifying whether to print on console.
            \param RstepMethod an R-string containing the method to use to choose the step in the optimization algorithm.
            \param RdirectionMethod an R-string containing the descent direction to use in the optimization algorithm.
            \param RpreprocessMethod an R-string containing the cross-validation method to use.
            \param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).
            \param RisTimeDiscrete an R-integer specifying the time data type: discrete (with many duplicates) or continuous.
            \param RflagMass an R-integer specifying whether to consider full mass matrices or identity mass matrices.
            \param RflagLumped an R-integer specifying whether to perform the mass lumping technique.

            \return R-list containing solutions.
        */
        SEXP Density_Estimation_time(SEXP Rdata, SEXP Rdata_time, SEXP Rmesh, SEXP Rmesh_time, SEXP Rorder, SEXP Rmydim,
                                     SEXP Rndim, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda,
                                     SEXP Rlambda_time, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1,
                                     SEXP Rtol2, SEXP Rprint, SEXP RstepMethod, SEXP RdirectionMethod,
                                     SEXP RpreprocessMethod, SEXP Rsearch, SEXP RisTimeDiscrete, SEXP RflagMass,
                                     SEXP RflagLumped)
        {

            UInt order = INTEGER(Rorder)[0];
            UInt mydim = INTEGER(Rmydim)[0];
            UInt ndim = INTEGER(Rndim)[0];

            std::string step_method = CHAR(STRING_ELT(RstepMethod, 0));
            std::string direction_method = CHAR(STRING_ELT(RdirectionMethod, 0));
            std::string preprocess_method = CHAR(STRING_ELT(RpreprocessMethod, 0));

            if(order == 1 && mydim == 2 && ndim == 2)
                return(DE_skeleton_time<1, 2, 2>(Rdata, Rdata_time, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rlambda_time, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rmesh_time, Rsearch, RisTimeDiscrete, RflagMass, RflagLumped, step_method, direction_method, preprocess_method));
            else if(order == 2 && mydim == 2 && ndim == 2)
                return(DE_skeleton_time<2, 2, 2>(Rdata, Rdata_time, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rlambda_time, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rmesh_time, Rsearch, RisTimeDiscrete, RflagMass, RflagLumped, step_method, direction_method, preprocess_method));
            else if(order == 1 && mydim == 2 && ndim == 3)
                return(DE_skeleton_time<1, 2, 3>(Rdata, Rdata_time, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rlambda_time, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rmesh_time, Rsearch, RisTimeDiscrete, RflagMass, RflagLumped, step_method, direction_method, preprocess_method));
            else if(order == 2 && mydim == 2 && ndim == 3)
                return(DE_skeleton_time<2, 2, 3>(Rdata, Rdata_time, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rlambda_time, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rmesh_time, Rsearch, RisTimeDiscrete, RflagMass, RflagLumped, step_method, direction_method, preprocess_method));
            else if(order == 1 && mydim == 3 && ndim == 3)
                return(DE_skeleton_time<1, 3, 3>(Rdata, Rdata_time, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rlambda_time, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rmesh_time, Rsearch, RisTimeDiscrete, RflagMass, RflagLumped, step_method, direction_method, preprocess_method));
            else if(order == 2 && mydim == 3 && ndim == 3)
                return(DE_skeleton_time<2, 3, 3>(Rdata, Rdata_time, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rlambda_time, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rmesh_time, Rsearch, RisTimeDiscrete, RflagMass, RflagLumped, step_method, direction_method, preprocess_method));

            return(NILSXP);

        }


        SEXP Density_Initialization_time(SEXP Rdata, SEXP Rdata_time, SEXP Rmesh, SEXP Rmesh_time, SEXP Rorder,
                                         SEXP Rmydim, SEXP Rndim, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda,
                                         SEXP Rlambda_time, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2,
                                         SEXP Rprint, SEXP Rsearch, SEXP RisTimeDiscrete, SEXP RflagMass, SEXP RflagLumped,
                                         SEXP Rinit, SEXP Rinit_fold)

        {

            UInt order = INTEGER(Rorder)[0];
            UInt mydim = INTEGER(Rmydim)[0];
            UInt ndim = INTEGER(Rndim)[0];

            UInt init_fold = INTEGER(Rinit_fold)[0];

            std::string init = CHAR(STRING_ELT(Rinit, 0));

            if(order == 1 && mydim == 2 && ndim == 2)
                return(DE_init_skeleton_time<1, 2, 2>(Rdata, Rdata_time, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rlambda_time, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rmesh_time, Rsearch, RisTimeDiscrete, RflagMass, RflagLumped, init, init_fold));
            else if(order == 2 && mydim == 2 && ndim == 2)
                return(DE_init_skeleton_time<2, 2, 2>(Rdata, Rdata_time, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rlambda_time, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rmesh_time, Rsearch, RisTimeDiscrete, RflagMass, RflagLumped, init, init_fold));
            else if(order == 1 && mydim == 2 && ndim == 3)
                return(DE_init_skeleton_time<1, 2, 3>(Rdata, Rdata_time, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rlambda_time, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rmesh_time, Rsearch, RisTimeDiscrete, RflagMass, RflagLumped, init, init_fold));
            else if(order == 2 && mydim == 2 && ndim == 3)
                return(DE_init_skeleton_time<2, 2, 3>(Rdata, Rdata_time, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rlambda_time, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rmesh_time, Rsearch, RisTimeDiscrete, RflagMass, RflagLumped, init, init_fold));
            else if(order == 1 && mydim == 3 && ndim == 3)
                return(DE_init_skeleton_time<1, 3, 3>(Rdata, Rdata_time, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rlambda_time, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rmesh_time, Rsearch, RisTimeDiscrete, RflagMass, RflagLumped, init, init_fold));
            else if(order == 2 && mydim == 3 && ndim == 3)
                return(DE_init_skeleton_time<2, 3, 3>(Rdata, Rdata_time, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rlambda_time, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rmesh_time, Rsearch, RisTimeDiscrete, RflagMass, RflagLumped, init, init_fold));

            return(NILSXP);

        }

}
