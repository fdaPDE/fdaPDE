#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh.h"
#include "../Include/Auxiliary_Mesh_Skeletons.h"
#include "../../FE_Assemblers_Solvers/Include/Projection.h"
#include "../../Mesh/Include/Domain.h"

extern "C"{

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

    SEXP eval_FEM_time_Auxiliary_new(SEXP Rmesh, SEXP Rmesh_time, SEXP Rlocations, SEXP Rtime_locations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rorder, SEXP Rfast, SEXP Rflag_parabolic, SEXP Rmydim, SEXP Rndim, SEXP Rsearch, SEXP RbaryLocations)
    {
        UInt order = INTEGER(Rorder)[0];
        UInt mydim = INTEGER(Rmydim)[0];
        UInt ndim  = INTEGER(Rndim)[0];
        UInt flag_par = INTEGER(Rflag_parabolic)[0]; // UInt DEGREE =  flag_par ? 1 : 3

        if(order==1 && mydim==2 && ndim==2 && flag_par==1)
            return Eval_FEM_time_skeleton_new<1,2,2,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==1 && mydim==2 && ndim==2 && flag_par!=1)
            return Eval_FEM_time_skeleton_new<1,2,2,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==2 && ndim==2 && flag_par==1)
            return Eval_FEM_time_skeleton_new<2,2,2,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==2 && ndim==2 && flag_par!=1)
            return Eval_FEM_time_skeleton_new<2,2,2,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==1 && mydim==2 && ndim==3 && flag_par==1)
            return Eval_FEM_time_skeleton_new<1,2,3,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==1 && mydim==2 && ndim==3 && flag_par!=1)
            return Eval_FEM_time_skeleton_new<1,2,3,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==2 && ndim==3 && flag_par==1)
            return Eval_FEM_time_skeleton_new<2,2,3,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==2 && ndim==2 && flag_par!=1)
            return Eval_FEM_time_skeleton_new<1,2,3,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==1 && mydim==3 && ndim==3 && flag_par==1)
            return Eval_FEM_time_skeleton_new<1,3,3,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==1 && mydim==2 && ndim==2 && flag_par!=1)
            return Eval_FEM_time_skeleton_new<1,3,3,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==3 && ndim==3 && flag_par==1)
            return Eval_FEM_time_skeleton_new<2,3,3,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==2 && ndim==2 && flag_par!=1)
            return Eval_FEM_time_skeleton_new<2,3,3,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);

        return NILSXP;


    }

    SEXP isInside(SEXP Rmesh,SEXP Rpoints,SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP Rsearch, SEXP Rredundancy)
    {
        UInt order = INTEGER(Rorder)[0];
        UInt mydim = INTEGER(Rmydim)[0];
        UInt ndim = INTEGER(Rndim)[0];

        if( order==1 && mydim==2 && ndim==2)
            return isInside_skeleton<1,2,2>(Rmesh, Rpoints, Rsearch, Rredundancy);
        else if(order==2 && mydim==2 && ndim==2)
            return isInside_skeleton<2,2,2>(Rmesh, Rpoints, Rsearch, Rredundancy);
        else if(order==1 && mydim==2 && ndim==3)
            return isInside_skeleton<1,2,3>(Rmesh, Rpoints, Rsearch, Rredundancy);
        else if(order==2 && mydim==2 && ndim==3)
            return isInside_skeleton<2,2,3>(Rmesh, Rpoints, Rsearch, Rredundancy);
        else if(order==1 && mydim==3 && ndim==3)
            return isInside_skeleton<1,3,3>(Rmesh, Rpoints, Rsearch, Rredundancy);
        else if(order==2 && mydim==3 && ndim==3)
            return isInside_skeleton<2,3,3>(Rmesh, Rpoints, Rsearch, Rredundancy);
        else if(order==1 && mydim==1 && ndim==2)
            return isInside_skeleton<1,1,2>(Rmesh, Rpoints, Rsearch, Rredundancy);
        else if(order==2 && mydim==1 && ndim==2)
            return isInside_skeleton<2,1,2>(Rmesh, Rpoints, Rsearch, Rredundancy);

        return NILSXP;
    }

    SEXP CPP_integrate_f(SEXP Rmesh, SEXP Rsearch, SEXP Rcoeff){
        UInt search_ = INTEGER(Rsearch)[0];
        const RNumericMatrix coef(Rcoeff);
        MeshHandler<1,1,2> mesh_(Rmesh,search_);

        SEXP result=PROTECT(Rf_allocVector(REALSXP, 1));
        Real* integral = REAL(result);
        (*integral) = 0. ;

        MeshHandler<1,1,2>::meshElement current_element;
        Eigen::Matrix<Real,2,1> coefficients;

        for(UInt e=0; e<mesh_.num_elements(); ++e){
            current_element = mesh_.getElement(e);
            for (UInt i=0; i<2; ++i)
                coefficients[i]=coef[current_element[i].getId()];
            (*integral) += current_element.integrate(coefficients);
        }

        UNPROTECT(1);
        return result;
    }

}