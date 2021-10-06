#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh.h"
#include "../Include/Auxiliary_Mesh_Skeletons.h"
#include "../../FE_Assemblers_Solvers/Include/Projection.h"

extern "C"{

    SEXP get_meshHandler(SEXP Rmesh_, SEXP order_, SEXP ndim_, SEXP mydim_, SEXP Rpoints_) {
        UInt order = INTEGER(order_)[0];
        UInt mydim = INTEGER(mydim_)[0];
        UInt ndim = INTEGER(ndim_)[0];

        if (order == 1 && mydim == 1 && ndim == 2)
            return (Auxiliary_Mesh_Skeleton<1, 1, 2>(Rmesh_, Rpoints_));
        else if (order == 2 && mydim == 1 && ndim == 2)
            return (Auxiliary_Mesh_Skeleton<2, 1, 2>(Rmesh_, Rpoints_));

        return (NILSXP);
    }

    SEXP reading_RObject(SEXP Rnodes, SEXP Redges){

        RNumericMatrix nodes(Rnodes);
        Element<3,1,2> myElement(1,std::array<Point<2> , 3>(
                {Point<2>({nodes(0,0),nodes(0,1)}),
                 Point<2>({nodes(1,0),nodes(1,1)}),
                 Point<2>({nodes(2,0),nodes(2,1)}) }));

        SEXP result;
        PROTECT(result = Rf_allocMatrix(REALSXP,1,1));

        RNumericMatrix tmp(result);
        tmp[0] = myElement.getMeasure();

        UNPROTECT(1);
        return result;
    }


    SEXP eval_FEM_fd_Auxiliary(SEXP Rmesh, SEXP Rlocations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rorder, SEXP Rfast, SEXP Rmydim, SEXP Rndim, SEXP Rsearch, SEXP RbaryLocations){

        UInt order = INTEGER(Rorder)[0];
        UInt mydim = INTEGER(Rmydim)[0];
        UInt ndim  = INTEGER(Rndim)[0];

        if(order==1 && mydim==1 && ndim==2)
            return Eval_FEM_fd_Skeleton<1,1,2>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==1 && ndim==2)
            return Eval_FEM_fd_Skeleton<2,1,2>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);

        return NILSXP;
    }

}