#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh.h"
#include "../Include/Auxiliary_Mesh_Skeletons.h"

extern "C"{

    SEXP get_meshHandler(SEXP Rmesh_, SEXP order_, SEXP ndim_, SEXP mydim_, SEXP Rpoints_){
        UInt order = INTEGER(order_)[0];
        UInt mydim = INTEGER(mydim_)[0];
        UInt ndim = INTEGER(ndim_)[0];

        if(order==1 && mydim==1 && ndim==2)
            return (Auxiliary_Mesh_Skeleton<1,1,2>(Rmesh_,Rpoints_));
        else if(order==2 && mydim==1 && ndim==2)
            return (Auxiliary_Mesh_Skeleton<2,1,2>(Rmesh_,Rpoints_));

        return(NILSXP);

    }

}