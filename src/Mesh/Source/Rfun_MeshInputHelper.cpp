#include "../../FdaPDE.h"

#include "../Include/Mesh_Input_Helper.h"

#include <array>

extern "C" {


SEXP CPP_SurfaceMeshHelper(SEXP Rtriangles, SEXP Rnodes){

  static constexpr std::array<UInt, 6> EDGES_ORDERING = {1,2,0,2,0,1};


  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 4));

  {
    simplex_container<2> edges_list(Rtriangles, Rnodes, EDGES_ORDERING);
    edges_list.assemble_subs(result, 0);
    edges_list.mark_boundary(result, 1);
    edges_list.compute_neighbors(result, 3);
  }

  mark_boundary_nodes(result, Rnodes, 2, 0, 1);

	UNPROTECT(1);

  return result;
}


SEXP CPP_SurfaceMeshOrder2(SEXP Rtriangles, SEXP Rnodes){

  static constexpr std::array<UInt, 6> EDGES_ORDERING = {1,2,0,2,0,1};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 6));

  {
    simplex_container<2> edges_list(Rtriangles, Rnodes, EDGES_ORDERING);
    edges_list.assemble_subs(result, 0);
    edges_list.mark_boundary(result, 1);
    edges_list.compute_neighbors(result, 3);
    edges_list.order2extend(result, 5);
  }

  mark_boundary_nodes(result, Rnodes, 2, 0, 1);
  compute_midpoints(result, Rnodes, 4, 0);

	UNPROTECT(1);

  return result;
}

SEXP CPP_TriangleMeshSplit(SEXP Rtriangles, SEXP Rnodes){

  static constexpr std::array<UInt, 6> EDGES_ORDERING = {1,2,0,2,0,1};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 2));

  {
    simplex_container<2> edges_list(Rtriangles, Rnodes, EDGES_ORDERING);
    split(result, Rtriangles, 0, edges_list);
    compute_midpoints(result, Rnodes, 1, edges_list);
  }

	UNPROTECT(1);

  return result;
}

SEXP CPP_TriangleMeshSplitOrder2(SEXP Rtriangles, SEXP Rnodes){

  static constexpr std::array<UInt, 6> EDGES_ORDERING = {1,2,0,2,0,1};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 1));

  {
    simplex_container<2> edges_list(Rtriangles, Rnodes, EDGES_ORDERING);
    split(result, Rtriangles, 0, edges_list);
  }

  UNPROTECT(1);

  return result;
}



SEXP CPP_VolumeMeshHelper(SEXP Rtetrahedrons, SEXP Rnodes){

  static constexpr std::array<UInt, 12> FACES_ORDERING = {1,2,3,0,2,3,0,1,3,0,1,2};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 4));

  {
    simplex_container<3> faces_list(Rtetrahedrons, Rnodes, FACES_ORDERING);
    faces_list.assemble_subs(result, 0);
    faces_list.mark_boundary(result, 1);
    faces_list.compute_neighbors(result, 3);
  }

  mark_boundary_nodes(result, Rnodes, 2, 0, 1);

	UNPROTECT(1);

  return result;

}

SEXP CPP_VolumeMeshOrder2(SEXP Rtetrahedrons, SEXP Rnodes){

  static constexpr std::array<UInt, 12> FACES_ORDERING = {1,2,3,0,2,3,0,1,3,0,1,2};
  static constexpr std::array<UInt, 12> EDGES_ORDERING = {0,1,0,2,0,3,1,2,2,3,1,3};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 6));

  {
    simplex_container<3> faces_list(Rtetrahedrons, Rnodes, FACES_ORDERING);
    faces_list.assemble_subs(result, 0);
    faces_list.mark_boundary(result, 1);
    faces_list.compute_neighbors(result, 3);
  }

  mark_boundary_nodes(result, Rnodes, 2, 0, 1);

  {
    simplex_container<2> edges_list(Rtetrahedrons, Rnodes, EDGES_ORDERING);
    edges_list.order2extend(result, 5);
    compute_midpoints(result, Rnodes, 4, edges_list);
  }

	UNPROTECT(1);

  return result;

}


SEXP CPP_TetraMeshSplit(SEXP Rtetrahedrons, SEXP Rnodes){

  static constexpr std::array<UInt, 12> EDGES_ORDERING = {0,1,0,2,0,3,1,2,2,3,1,3};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 2));

  {
    simplex_container<2> edges_list(Rtetrahedrons, Rnodes, EDGES_ORDERING);
    split3D(result, Rtetrahedrons, 0, edges_list);
    compute_midpoints(result, Rnodes, 1, edges_list);
  }

	UNPROTECT(1);

  return result;
}

SEXP CPP_TetraMeshSplitOrder2(SEXP Rtetrahedrons, SEXP Rnodes){

  static constexpr std::array<UInt, 12> EDGES_ORDERING = {0,1,0,2,0,3,1,2,2,3,1,3};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 1));

  {
    simplex_container<2> edges_list(Rtetrahedrons, Rnodes, EDGES_ORDERING);
    split3D(result, Rtetrahedrons, 0, edges_list);
  }

  UNPROTECT(1);

  return result;
}

SEXP CPP_EdgeMeshHelper(SEXP Redges, SEXP Rnodes){
    
    static constexpr std::array<UInt, 2> NODES_ORDERING = {1,0};
    
    
    SEXP result = NILSXP;
    result = PROTECT(Rf_allocVector(VECSXP, 4));
    
    {
      simplex_container<1> nodes_list(Redges, Rnodes, NODES_ORDERING);
      nodes_list.assemble_subs(result, 0);
      nodes_list.mark_boundary(result, 1);
      // compute_neighbors fills the index-th and the (index+1)-th positions in resultf
      nodes_list.compute_neighbors(result, 2);
    }
    
    UNPROTECT(1);
    
    return result;
}

SEXP CPP_EdgeMeshOrder2(SEXP Redges, SEXP Rnodes){

    static constexpr std::array<UInt, 2> NODES_ORDERING = {1,0};


    SEXP result = NILSXP;
    result = PROTECT(Rf_allocVector(VECSXP, 6));

    {
        simplex_container<1> nodes_list(Redges, Rnodes, NODES_ORDERING);
        nodes_list.assemble_subs(result, 0);
        nodes_list.mark_boundary(result, 1);
        //compute_neighbors fills index 2 and index 3
        nodes_list.compute_neighbors(result, 2);
        compute_midpoints(result, Rnodes,Redges,4);

        //midpoints global numbering
        SET_VECTOR_ELT(result,5, Rf_allocMatrix(INTSXP, nodes_list.get_num_elements(),1 ));
        RIntegerMatrix midpoints(VECTOR_ELT(result,5));
        UInt num_points = nodes_list.get_num_points();
        for(UInt i=0; i<nodes_list.get_num_elements(); ++i, ++num_points)
            midpoints[i]=num_points;
    }


    UNPROTECT(1);

    return result;
}

SEXP CPP_EdgeMeshSplit(SEXP Redges, SEXP Rnodes){
  
  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 2));
  
  {
    //split1D -> 0: splitted_edges
    //        -> 1: vectors of global midpoints numbers
    split1D(result, Rnodes, Redges, 0);
    // midpoints coordinates
    compute_midpoints(result, Rnodes, Redges, 1);
  }
  
  UNPROTECT(1);
  
  return result;
  }

}