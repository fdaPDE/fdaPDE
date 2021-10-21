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
/*
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
*/

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

template<UInt ORDER, UInt mydim, UInt ndim, UInt DEGREE>
SEXP Eval_FEM_time_skeleton_new (SEXP Rmesh, SEXP Rmesh_time, SEXP Rlocations, SEXP Rtime_locations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rfast, SEXP Rsearch, SEXP RbaryLocations)
{/*
    UInt n = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
    UInt ns = INTEGER(Rf_getAttrib(VECTOR_ELT(Rmesh, 0), R_DimSymbol))[0];
    UInt nt = Rf_length(Rmesh_time);
    UInt nRegions = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[0];
*/
    RNumericMatrix coef(Rcoef); // ns x M ?
    RNumericMatrix locations(Rlocations);
    RIntegerMatrix incidenceMatrix(RincidenceMatrix);

    UInt n = locations.nrows();
    UInt ns = INTEGER(Rf_getAttrib(VECTOR_ELT(Rmesh, 0), R_DimSymbol))[0];
    UInt nt = Rf_length(Rmesh_time);
    UInt nRegions = incidenceMatrix.nrows();
    UInt nElements = incidenceMatrix.ncols();

    Real *mesh_time, *t; //questi diventano RNumericMatrix... Nope servono a Spline e sticazzi...
    mesh_time = REAL(Rmesh_time);
    t = REAL(Rtime_locations);

    UInt M = nt + DEGREE - 1;
    //SpMat phi(n,M);
    //UInt N = nRegions==0 ? n : nRegions;
    UInt N = nRegions==0 ? n : nRegions;
    SpMat phi(N,M);
    Spline<DEGREE,DEGREE-1>spline(mesh_time,nt);
    Real value;
    for (UInt i = 0; i < N; ++i)
    {
        for (UInt j = 0; j < M; ++j)
        {
            value = spline.BasisFunction(j, t[i]);
            if (value!=0)
            {
                phi.coeffRef(i,j) = value;
            }
        }
    }
    phi.makeCompressed();

    SEXP result;
    PROTECT(result=Rf_allocVector(REALSXP, N));  // 0 Attenzione

    SEXP Rcoef_0;
    PROTECT(Rcoef_0=Rf_allocMatrix(REALSXP, ns,1)); // 1
    RNumericMatrix coef_0(Rcoef_0);

    //!evaluates the solution on the given points location at the first
    //!node of the time mesh to initialize the array of results and retrieve the points out of mesh (NA)
    for(UInt j=0; j<ns; ++j)
    {
        coef_0[j] = coef[j];
    }
    //temp non va protetto?
    SEXP temp = Eval_FEM_fd_Skeleton_new<ORDER,mydim,ndim>(Rmesh,Rlocations, RincidenceMatrix,Rcoef_0, Rfast,Rsearch, RbaryLocations);
    UNPROTECT(1); //UNPROTECT Rcoef_0 ?!

    // nb) temp è una matrice a priori ma non dovrebbero esserci problemi
    for(UInt k=0; k < N; k++) {
        REAL(result)[k] = REAL(temp)[k];
        if (!ISNA(REAL(result)[k]))
            REAL(result)[k] = REAL(result)[k] * phi.coeff(k, 0);
    }

    //! loop over time b-splines basis and evaluate the solution only on the points that have
    //! the coefficient corresponding to that basis different from 0
    std::vector<UInt> indices;
    for(UInt i=1; i<M; ++i)
    {

        SEXP Rcoef_i;
        PROTECT(Rcoef_i=Rf_allocMatrix(REALSXP, ns,1)); // 1
        RNumericMatrix coef_i(Rcoef_i);
        for(UInt j=0; j<ns; ++j)
        {
            coef_i[j] = coef[i*ns+j];
        }
        // if n > 0
        UInt num_locations_i = 0; //non entro se n=0
        for(UInt k=0; k<n; k++) {
            if (phi.coeff(k, i) != 0 && !ISNA(REAL(result)[k])) {
                ++num_locations_i;
                indices.push_back(k);
            }
        }

        SEXP Rlocations_i;
        PROTECT(Rlocations_i = Rf_allocMatrix(REALSXP,num_locations_i,ndim)); //2
        RNumericMatrix locations_i(Rlocations_i);
        UInt count_locations_i = 0;
        for(UInt k=0; k<indices.size(); ++k){
            for(UInt l=0; l<ndim;++l)
             locations_i(count_locations_i,l) = locations(indices[k],l);
            ++count_locations_i;
        }

        SEXP RincidenceMatrix_i;
        PROTECT(RincidenceMatrix_i=Rf_allocMatrix(INTSXP,phi.col(i).nonZeros(), nElements)); //3
        RIntegerMatrix incidenceMatrix_i(RincidenceMatrix_i);
        UInt count_incidence_i = 0; // non entro se nRegions=0
        for (UInt k=0; k<nRegions; ++k)
        {
            if(phi.coeff(k,i)!=0 && !ISNA(REAL(result)[k]))
            {
                for (UInt j=0; j<nElements; j++)
                {
                    incidenceMatrix_i(count_incidence_i , j) = incidenceMatrix(k,j);
                }
                indices.push_back(k);
                ++count_incidence_i;
            }
        }
        temp = Eval_FEM_fd_Skeleton_new<ORDER,mydim,ndim>(Rmesh,Rlocations_i, RincidenceMatrix_i,Rcoef_i, Rfast,Rsearch, RbaryLocations);
        UNPROTECT(3); // UNPROTECT Rcoef_i, Rlocations_i RincidenceMatrix_i

        for(UInt k=0; k<indices.size(); ++k)
        {
            REAL(result)[indices[k]] = REAL(result)[indices[k]] + REAL(temp)[k]*phi.coeff(indices[k],i);
        }
        indices.clear();

    }

    UNPROTECT(1); //UNPROTECT result;
    return result;


}

#endif //__AUXILIARY_MESH_SKELETONS_H
