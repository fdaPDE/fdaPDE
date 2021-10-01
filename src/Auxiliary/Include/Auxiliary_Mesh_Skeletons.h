#ifndef __AUXILIARY_MESH_SKELETONS_H
#define __AUXILIARY_MESH_SKELETONS_H
#include "../../Mesh/Include/Mesh.h"
#include "../../FdaPDE.h"

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


#endif //__AUXILIARY_MESH_SKELETONS_H
