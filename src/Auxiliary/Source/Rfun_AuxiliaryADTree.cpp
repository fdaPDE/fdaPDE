
#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../Mesh/Include/Domain.h"
#include "../../Mesh/Include/Bounding_Box.h"
#include "../../Mesh/Include/Tree_Header.h"
#include "../../Mesh/Include/Tree_Node.h"
#include "../../Mesh/Include/AD_Tree.h"
#include "../Include/Auxiliary_Mesh_Skeletons.h"


extern "C"{

    SEXP TimingSearch(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP Rpoints){
        UInt order = INTEGER(Rorder)[0];
        UInt mydim = INTEGER(Rmydim)[0];
        UInt ndim = INTEGER(Rndim)[0];

        if(order==1 && mydim==1 && ndim==2)
            return TimingSearch_Skeleton<1,1,2>(Rmesh,Rpoints);
        else if(order==2 && mydim==1 && ndim==2)
            return TimingSearch_Skeleton<2,1,2>(Rmesh,Rpoints);
        else if( order==1 && mydim==2 && ndim==2)
            return TimingSearch_Skeleton<1,2,2>(Rmesh,Rpoints);
        else if( order==2 && mydim==2 && ndim==2)
            return TimingSearch_Skeleton<2,2,2>(Rmesh,Rpoints);
        else if( order==1 && mydim==2 && ndim==3)
            return TimingSearch_Skeleton<1,2,3>(Rmesh,Rpoints);
        else if( order==2 && mydim==2 && ndim==3)
            return TimingSearch_Skeleton<2,2,3>(Rmesh,Rpoints);
        else if( order==1 && mydim==3 && ndim==3)
            return TimingSearch_Skeleton<1,3,3>(Rmesh,Rpoints);
        else if( order==2 && mydim==3 && ndim==3)
            return TimingSearch_Skeleton<2,3,3>(Rmesh,Rpoints);

        return NILSXP;
    }


}
