#ifndef __PROJECTION_H__
#define __PROJECTION_H__

#include "../../Mesh/Include/Mesh.h"

template<UInt ORDER, UInt mydim, UInt ndim>
class projection{
    /*static_assert((ORDER==1 || ORDER==2) &&
    (mydim==1 || mydim==2 || mydim==3) &&
    mydim <= ndim,
    "ERROR! TRYING TO INSTANTIATE PROJECTION WITH WRONG NUMBER OF PARAMETERS See projection.h");
    */
private:
    const MeshHandler<ORDER,mydim,ndim>& mesh_;
    const std::vector<Point<ndim> > & deData_; // the points to be projected
    const UInt num_points;

    template<bool isManifold=(ndim!=mydim)>
    typename std::enable_if<isManifold, std::vector<UInt> >::type
    computeNodePatch(UInt id_node) const;

public:
    projection(const MeshHandler<ORDER,mydim,ndim>& m, const std::vector<Point<ndim> > & d): mesh_(m), deData_(d), num_points(d.size()) {};

    template<bool isManifold=(ndim!=mydim)>
    typename std::enable_if<isManifold, std::vector<Point<ndim>> >::type
    computeProjection();

    template<bool isManifold=(ndim!=mydim)>
    typename std::enable_if<!isManifold, std::vector<Point<ndim>> >::type
    computeProjection();
};

#include "Projection_imp.h"

#endif
