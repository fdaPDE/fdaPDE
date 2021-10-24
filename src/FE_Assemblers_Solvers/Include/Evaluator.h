#ifndef __EVALUATOR_H__
#define __EVALUATOR_H__

#include <iostream>
#include <algorithm>

#include "../../FdaPDE.h"
#include "Finite_Element.h"
#include "Matrix_Assembler.h"
#include "../../Mesh/Include/Mesh.h"

template<UInt ORDER, UInt mydim, UInt ndim>
class Evaluator
{
public:
    //! A constructor. It initializes the constructor given a mesh object.
    Evaluator(const MeshHandler<ORDER,mydim,ndim>& mesh): mesh_(mesh){};

    //! A member that computes the evaluation of a Point in a mesh, given the bases' coefficients.
    /*!
    \param locations a RNumericMatrix (num_points x ndim) containing the coordinates of points to evaluate.
    \param coef a RNumericMatrix (num_basis x 1) containing the coefficients of the solution, the value in position i
    is associated to the basis \phi(i).
    \param redundancy a boolean that specifies if the algorithm is completely based on the walking
            algorithm (can miss locations in case of non convex structures)
    \param result a RNumericMatrix (num_points x 1) to an already protected memory space, where the evaluations
    will be stored
    */
    template<bool isManifold=(ndim!=mydim)>
    typename std::enable_if<isManifold,void>::type
    eval(const RNumericMatrix& locations,const RNumericMatrix& coef, bool redundancy,RNumericMatrix& result, std::vector<bool>& isinside);

    template<bool isManifold=(ndim!=mydim)>
    typename std::enable_if<!isManifold,void>::type
    eval(const RNumericMatrix& locations,const RNumericMatrix& coef, bool redundancy,RNumericMatrix& result, std::vector<bool>& isinside);

    void evalWithInfo(const RNumericMatrix& locations,const RNumericMatrix& coef, bool redundancy,RNumericMatrix& result, std::vector<bool>& isinside, const RIntegerMatrix& element_id,const RNumericMatrix& barycenters);

    //! A member that computes the integral over regions divided by the measure of the region in a mesh,
    //!  given the bases' coefficients.
    /*!
    \param incidenceMatrix a nRegions*nElements RIntegerMatrix telling which triangle compose each region.
    \param coef a pointer to the vector of coefficients of the solution, the value in position i
    is associated to the basis \phi(i).
    \param result a RNumericMatrix to an already protected memory space, where the evaluations
    will be stored.
    */
    void integrate(const RIntegerMatrix& incidenceMatrix, const RNumericMatrix& coef, RNumericMatrix& result);

private:
    const MeshHandler<ORDER,mydim,ndim> & mesh_;
};



#include "Evaluator_imp.h"

#endif
