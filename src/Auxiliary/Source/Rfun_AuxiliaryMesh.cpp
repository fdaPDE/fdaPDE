#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh.h"
#include "../Include/Auxiliary_Mesh_Skeletons.h"
#include "../../FE_Assemblers_Solvers/Include/Projection.h"
#include "../../Mesh/Include/Domain.h"

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

    SEXP reading_RObject(SEXP Rmatrix, SEXP Rflag) {

        // flag = 0 -> leggo matrici di double
        // flag = 1 -> leggo matrici interi
        UInt flag = INTEGER(Rflag)[0];
        SEXP result;
        PROTECT(result = Rf_allocMatrix(INTSXP, 1, 2));
        RIntegerMatrix result_(result);
        if (flag == 0) {
            RNumericMatrix matrix_(Rmatrix);
            result_[0] = matrix_.nrows();
            result_[1] = matrix_.ncols();
        } else if (flag == 1) {
            RIntegerMatrix matrix_(Rmatrix);
            result_[0] = matrix_.nrows();
            result_[1] = matrix_.ncols();
        } else
            return NILSXP;

        UNPROTECT(1);
        return result;
    }

    SEXP reading_RIntMatrixMatrix(SEXP Rmatrix){

        RIntMatrixMatrix matrix_(Rmatrix);
        SEXP result;
        PROTECT(result=Rf_allocMatrix(INTSXP,matrix_.nrows()*matrix_.ncols(),2));

        RIntegerMatrix result_(result);

        for(UInt i=0; i<matrix_.nrows()*matrix_.ncols(); ++i) {
            result_(i, 0) = matrix_[i].nrows();
            result_(i,1)  = matrix_[i].ncols();
        }
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

    SEXP eval_FEM_fd_Auxiliary_new(SEXP Rmesh, SEXP Rlocations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rorder, SEXP Rfast, SEXP Rmydim, SEXP Rndim, SEXP Rsearch, SEXP RbaryLocations){

    UInt order = INTEGER(Rorder)[0];
    UInt mydim = INTEGER(Rmydim)[0];
    UInt ndim  = INTEGER(Rndim)[0];

    if(order==1 && mydim==1 && ndim==2)
        return Eval_FEM_fd_Skeleton_new<1,1,2>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==2 && mydim==1 && ndim==2)
        return Eval_FEM_fd_Skeleton_new<2,1,2>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==1 && mydim==2 && ndim==2)
        return Eval_FEM_fd_Skeleton_new<1,2,2>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==2 && mydim==2 && ndim==2)
        return Eval_FEM_fd_Skeleton_new<2,2,2>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==1 && mydim==2 && ndim==3)
        return Eval_FEM_fd_Skeleton_new<1,2,3>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==2 && mydim==2 && ndim==3)
        return Eval_FEM_fd_Skeleton_new<2,2,3>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==1 && mydim==3 && ndim==3)
        return Eval_FEM_fd_Skeleton_new<1,3,3>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==2 && mydim==3 && ndim==3)
        return Eval_FEM_fd_Skeleton_new<2,3,3>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);

    return NILSXP;
    }

}