#ifndef __AUXILIARY_MESH_SKELETONS_H
#define __AUXILIARY_MESH_SKELETONS_H
#include "../../Mesh/Include/Mesh.h"
#include "../../FdaPDE.h"
#include "../../FE_Assemblers_Solvers/Include/Projection.h"
#include "../../FE_Assemblers_Solvers/Include/Evaluator.h"
#include "Evaluator_New.h"
#include "../../Global_Utilities/Include/Timing.h"

template<UInt ORDER, UInt mydim, UInt ndim>
SEXP Auxiliary_Mesh_Skeleton(SEXP Rmesh, SEXP Rpoints) {

    RNumericMatrix PointsToProject(Rpoints);

    UInt num_points = PointsToProject.nrows();
    MeshHandler<ORDER, mydim, ndim> meshHandler(Rmesh);
    UInt num_elements = meshHandler.num_elements();

    SEXP result;
    PROTECT(result = Rf_allocVector(VECSXP, 1 + num_points * 2));

    SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, num_elements, 1));
    for (UInt i = 1; i < 2 * num_points; i += 2) {
        SET_VECTOR_ELT(result, i, Rf_allocMatrix(INTSXP, num_elements, 1)); //isPointInside column
        SET_VECTOR_ELT(result, i + 1, Rf_allocMatrix(REALSXP, num_elements, 2)); //coord of projected point
    }

    RNumericMatrix MeasuresColumn(VECTOR_ELT(result, 0));

    for (UInt i = 0; i < num_elements; ++i) {
        typename MeshHandler<ORDER, mydim, ndim>::meshElement curr = meshHandler.getElement(i);
        MeasuresColumn[i] = curr.getMeasure();
    }

    for (UInt j = 0; j < num_points; ++j) {
        RIntegerMatrix PointIsInsideColumn(VECTOR_ELT(result, 1 + 2 * j));
        RNumericMatrix CoordsColumn(VECTOR_ELT(result, 2 + 2 * j));
        Point<2> curr_point({PointsToProject(j, 0), PointsToProject(j, 1)});
        for(UInt i=0;i<num_elements;++i){
            typename MeshHandler<ORDER, mydim, ndim>::meshElement curr = meshHandler.getElement(i);
            PointIsInsideColumn[i] = curr.isPointInside(curr_point);
            Point<2> curr_proj = curr.computeProjection(curr_point);
            CoordsColumn(i, 0) = curr_proj[0];
            CoordsColumn(i, 1) = curr_proj[1];
        }
    }

    UNPROTECT(1);
    return(result);

}

template<UInt ORDER, UInt mydim, UInt ndim>
SEXP Eval_FEM_fd_Skeleton(SEXP Rmesh, SEXP Rlocations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rfast, SEXP Rsearch, SEXP RbaryLocations){

    //NB anche se non viene passato DEVI essere certo di creare un RbaryLocations "NULLO" RbaryLocations è una lista di liste
    // in modo tale da poter utilizzare gli RObjects!!! (al più sono matrici che hanno righe o colonne = 0
    // RbaryLocations[0] contiene Rlocations (esattamente uguali a Rlocations!!!!)
    // RbaryLocations[1] contiene elements_ids
    // RbaryLocations[2] contiene coordinate dei baricentri degli elementi contenuti in elements_ids

    RNumericMatrix barycenters( VECTOR_ELT(RbaryLocations,2));
    RIntegerMatrix id_element( VECTOR_ELT(RbaryLocations,1));
    RIntegerMatrix incidenceMatrix( RincidenceMatrix );
    RNumericMatrix locations(Rlocations);

    UInt n_X = locations.nrows();
    UInt nRegions = incidenceMatrix.nrows();
    RNumericMatrix coef(Rcoef);
    UInt search;
    bool fast;

    fast  = INTEGER(Rfast)[0];
    search  = INTEGER(Rsearch)[0];
    MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, search);
    Evaluator<ORDER, mydim, ndim> evaluator(mesh);

    SEXP result;

    if(n_X >0) {
        PROTECT(result=Rf_allocMatrix(REALSXP,n_X,1));
        RNumericMatrix result_(result);

        std::vector<bool> isinside(n_X);
        if (barycenters.nrows() == 0) { //doesn't have location information
            evaluator.eval(locations, coef, fast, result_, isinside);
        } else { //have location information
            evaluator.evalWithInfo(locations, coef, fast, result_, isinside, id_element, barycenters);
        }

        for (int i = 0; i < n_X; ++i) {
            if (!(isinside[i])) {
                result_[i] = NA_REAL;
            }
        }
    }
    else{
        PROTECT(result=Rf_allocMatrix(REALSXP, nRegions,1));
        RNumericMatrix result_(result);
        evaluator.integrate(incidenceMatrix,coef,result_);

    }

    UNPROTECT(1);
    return result;

}

template<UInt ORDER, UInt mydim, UInt ndim>
SEXP Eval_FEM_fd_Skeleton_new(SEXP Rmesh, SEXP Rlocations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rfast, SEXP Rsearch, SEXP RbaryLocations){

    //NB anche se non viene passato DEVI essere certo di creare un RbaryLocations "NULLO" RbaryLocations è una lista di liste
    // in modo tale da poter utilizzare gli RObjects!!! (al più sono matrici che hanno righe o colonne = 0
    // RbaryLocations[0] contiene Rlocations (esattamente uguali a Rlocations!!!!)
    // RbaryLocations[1] contiene elements_ids
    // RbaryLocations[2] contiene coordinate dei baricentri degli elementi contenuti in elements_ids

    RNumericMatrix barycenters( VECTOR_ELT(RbaryLocations,2));
    RIntegerMatrix id_element( VECTOR_ELT(RbaryLocations,1));
    RIntegerMatrix incidenceMatrix( RincidenceMatrix );
    RNumericMatrix locations(Rlocations);

    UInt n_X = locations.nrows();
    UInt nRegions = incidenceMatrix.nrows();
    RNumericMatrix coef(Rcoef);
    UInt search;
    bool fast;

    fast  = INTEGER(Rfast)[0];
    search  = INTEGER(Rsearch)[0];
    MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, search);
    //NB unico cambio rispetto allo skeleton di prima
    Evaluator_New<ORDER, mydim, ndim> evaluator(mesh);

    SEXP result;

    if(n_X >0) {
        PROTECT(result=Rf_allocMatrix(REALSXP,n_X,1));
        RNumericMatrix result_(result);

        std::vector<bool> isinside(n_X);
        if (barycenters.nrows() == 0) { //doesn't have location information
            evaluator.eval(locations, coef, fast, result_, isinside);
        } else { //have location information
            evaluator.evalWithInfo(locations, coef, fast, result_, isinside, id_element, barycenters);
        }

        for (int i = 0; i < n_X; ++i) {
            if (!(isinside[i])) {
                result_[i] = NA_REAL;
            }
        }
    }
    else{
        PROTECT(result=Rf_allocMatrix(REALSXP, nRegions,1));
        RNumericMatrix result_(result);
        evaluator.integrate(incidenceMatrix,coef,result_);

    }

    UNPROTECT(1);
    return result;

}

template<UInt ORDER, UInt mydim, UInt ndim>
SEXP TimingSearch_Skeleton(SEXP Rmesh, SEXP Rpoints){

    MeshHandler<ORDER,mydim,ndim> meshHandler(Rmesh,2); //ADtree
    RNumericMatrix points(Rpoints);
    UInt num_points = points.nrows();

    SEXP result;
    PROTECT(result=Rf_allocMatrix(REALSXP,num_points,2));
    RNumericMatrix result_(result);
    timer Timer;

    for(UInt i=0; i<num_points; ++i) {
        std::array <Real, ndim> coords;
        for (UInt n = 0; n < ndim; ++n)
            coords[n] = points(i, n);

        Point<ndim> curr_point(coords);

        //Naive
        Timer.start();
        typename MeshHandler<ORDER, mydim, ndim>::meshElement curr1 = meshHandler.findLocationNaive(curr_point);
        timespec delta_naive = Timer.stop();
        result_(i,0) = double(delta_naive.tv_nsec);

        //Tree
        Timer.start();
        typename MeshHandler<ORDER, mydim, ndim>::meshElement curr2 = meshHandler.findLocationTree(curr_point);
        timespec delta_tree = Timer.stop();
        result_(i,1) = double(delta_tree.tv_nsec);
    }

    UNPROTECT(1);
    return result;
}

#endif //__AUXILIARY_MESH_SKELETONS_H
