#ifndef __FE_SKELETON_H__
#define __FE_SKELETON_H__

#include "../../FdaPDE.h"
#include "Projection.h"
#include "../../Mesh/Include/Mesh.h"

template<UInt ORDER, UInt mydim, UInt ndim>
SEXP tree_mesh_skeleton(SEXP Rmesh) {
    MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, 2);

    //Copy result in R memory
    SEXP result = NILSXP;
    result = PROTECT(Rf_allocVector(VECSXP, 5));


    //SEND TREE INFORMATION TO R
    SET_VECTOR_ELT(result, 0, Rf_allocVector(INTSXP, 1)); //tree_header information
    int *rans = INTEGER(VECTOR_ELT(result, 0));
    rans[0] = mesh.getTree().gettreeheader().gettreelev();

    SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain origin
    Real *rans1 = REAL(VECTOR_ELT(result, 1));
    for(UInt i = 0; i < ndim*2; i++)
        rans1[i] = mesh.getTree().gettreeheader().domainorig(i);

    SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain scale
    Real *rans2 = REAL(VECTOR_ELT(result, 2));
    for(UInt i = 0; i < ndim*2; i++)
        rans2[i] = mesh.getTree().gettreeheader().domainscal(i);


    UInt num_tree_nodes = mesh.num_elements()+1; //Be careful! This is not equal to number of elements
    SET_VECTOR_ELT(result, 3, Rf_allocMatrix(INTSXP, num_tree_nodes, 3)); //treenode information
    int *rans3 = INTEGER(VECTOR_ELT(result, 3));
    for(UInt i = 0; i < num_tree_nodes; i++)
        rans3[i] = mesh.getTree().gettreenode(i).getid();

    for(UInt i = 0; i < num_tree_nodes; i++)
        rans3[i + num_tree_nodes*1] = mesh.getTree().gettreenode(i).getchild(0);

    for(UInt i = 0; i < num_tree_nodes; i++)
        rans3[i + num_tree_nodes*2] = mesh.getTree().gettreenode(i).getchild(1);

    SET_VECTOR_ELT(result, 4, Rf_allocMatrix(REALSXP, num_tree_nodes, ndim*2)); //treenode box coordinate
    Real *rans4 = REAL(VECTOR_ELT(result, 4));
    for(UInt j = 0; j < ndim*2; j++)
    {
        for(UInt i = 0; i < num_tree_nodes; i++)
            rans4[i + num_tree_nodes*j] = mesh.getTree().gettreenode(i).getbox().get()[j];
    }


    UNPROTECT(1);
    return(result);
}

template<UInt ORDER,UInt mydim,UInt ndim>
SEXP point_projection_skeleton(SEXP Rmesh, SEXP Rlocations){
    //RECIEVE PROJECTION INFORMATION FROM R
    RNumericMatrix locations(Rlocations);
    UInt n_X = locations.nrows();

    // Cast all computation parameters
    std::vector<Point<ndim> > deData_(n_X); // the points to be projected
    std::vector<Point<ndim> > prjData_(n_X); // the projected points

    std::array<Real,ndim> coords;
    for(UInt i=0; i<n_X; ++i) {
        for(UInt n=0; n<ndim; ++n)
            coords[n] = locations(i,n);
        deData_[i] = Point<ndim>(coords);

    }
    SEXP result;

    if (n_X>0) //pointwise data
    {
        PROTECT(result = Rf_allocMatrix(REALSXP, n_X, ndim));
        MeshHandler<ORDER,mydim,ndim> mesh(Rmesh);
        projection<ORDER,mydim,ndim> projector(mesh, deData_);
        prjData_ = projector.computeProjection();

        RNumericMatrix res(result);
        for(UInt i=0; i<n_X; ++i){
            for(UInt n=0; n<ndim; ++n)
                res(i,n) = prjData_[i][n];
        }

        UNPROTECT(1);
        return(result);
    }

    return(NILSXP);
}
/*
template<UInt ORDER, UInt mydim, UInt ndim>
SEXP Eval_FEM_fd_skeleton(SEXP Rmesh, SEXP Rlocations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rfast, SEXP Rsearch, SEXP RbaryLocations){

    //NB anche se non viene passato DEVI essere certo di creare un RbaryLocations "NULLO" RbaryLocations è una lista di liste
    // in modo tale da poter utilizzare gli RObjects!!! (al più sono matrici che hanno righe o colonne = 0
    // RbaryLocations[0] contiene Rlocations (esattamente uguali a Rlocations!!!!)
    // RbaryLocations[1] contiene elements_ids
    // RbaryLocations[2] contiene coordinate dei baricentri degli elementi contenuti in elements_id
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
#endif
