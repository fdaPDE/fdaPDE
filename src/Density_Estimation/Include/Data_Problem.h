#ifndef __DATA_PROBLEM_H__
#define __DATA_PROBLEM_H__

#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <utility>
#include <numeric>
//#include <omp.h>
#include "../../FdaPDE.h"
#include "DE_Data.h"
#include "../../FE_Assemblers_Solvers/Include/Projection.h"
#include "../../FE_Assemblers_Solvers/Include/Integration.h"
#include "../../FE_Assemblers_Solvers/Include/Kronecker_Product.h"
#include "../../FE_Assemblers_Solvers/Include/Spline.h"
#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"

// This file contains data information for the Density Estimation problem

//! @brief A class to store common (spatial) data for the problem.
template<UInt ORDER, UInt mydim, UInt ndim>
class DataProblem{
protected:
    using Integrator = typename DensityIntegratorHelper::Integrator<mydim>;
    static constexpr UInt EL_NNODES = how_many_nodes(ORDER,mydim);
    DEData<ndim> deData_;
    MeshHandler<ORDER, mydim, ndim> mesh_;
    SpMat R0_, R1_, GlobalPsi_;
    MatrixXr P_;
    Eigen::Matrix<Real, Integrator::NNODES, EL_NNODES> PsiQuad_;

    //! A method to compute the finite element matrices.
    void fillFEMatrices();
    //! A method to compute the matrix which evaluates the basis function at the quadrature EL_NNODES.
    void fillPsiQuad();

public:
    //! A constructor: it delegates DEData and MeshHandler constructors.
    DataProblem(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds,
                SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch, SEXP Rmesh,
                bool isTime = 0);

    //! A method to compute the integral of a function (over the spatial domain).
    Real FEintegrate(const VectorXr& f) const {return (R0_*f).sum();}
    //! A method to compute the integral of the square of a function (over the spatial domain).
    Real FEintegrate_square(const VectorXr& f) const {return f.dot(R0_*f);}
    //! A method to compute the integral of the exponential of a function (over the spatial domain).
    Real FEintegrate_exponential(const VectorXr& g) const;
    //! A method to compute the matrix which evaluates the basis function at the data points.
    SpMat computePsi(const std::vector<UInt>& indices) const;

    // Getters
	//! A method to access the data. It calls the same method of DEData class.
    const std::vector<Point<ndim>>& data() const {return deData_.data();}
    //! A method returning a datum. It calls the same method of DEData class.
    const Point<ndim>& data(UInt i) const {return deData_.data(i);}
    //! A method returning the number of observations. It calls the same method of DEData class.
    UInt dataSize() const {return deData_.dataSize();}
	//! A method returning the the input order. It calls the same method of DEData class.
    UInt getOrder() const {return deData_.getOrder();}
	//! A method returning the initial coefficients for the density. It calls the same method of DEData class.
    VectorXr getFvec() const {return deData_.getFvec();}
	//! A method returning a bool which says if there is a user's initial density. It calls the same method of DEData class.
    bool isFvecEmpty() const {return deData_.isFvecEmpty();}
    //! A method returning the heat diffusion process alpha parameter. It calls the same method of DEData class.
    Real getHeatStep() const {return deData_.getHeatStep();}
    //! A method returning the number of iterations for the heat diffusion process. It calls the same method of DEData class.
    UInt getHeatIter() const {return deData_.getHeatIter();}
	//! A method returning the penalization parameters (in space). It calls the same method of DEData class.
    Real getLambda(UInt i) const {return deData_.getLambda(i);}
	//! A method returning the number of lambdas (in space). It calls the same method of DEData class.
    UInt getNlambda()  const {return deData_.getNlambda();}
	//! A method returning the number of folds for CV. It calls the same method of DEData class.
    UInt getNfolds()  const {return deData_.getNfolds();}
	//! A method returning the number of iterations to use in the optimization algorithm. It calls the same method of DEData class.
    UInt getNsimulations() const {return deData_.getNsimulations();}
	//! A method returning the number of parameters for fixed step methods. It calls the same method of DEData class.
    UInt getNstepProposals() const {return deData_.getNstepProposals();}
	//! A method returning a parameter for fixed step methods. It calls the same method of DEData class.
    Real getStepProposals(UInt i) const {return deData_.getStepProposals(i);}
    //! A method returning the tolerance for optimization algorithm first termination criterion. It calls the same method of DEData class.
    Real getTol1() const {return deData_.getTol1();}
    //! A method returning the tolerance for optimization algorithm second termination criterion. It calls the same method of DEData class.
    Real getTol2() const {return deData_.getTol2();}
    //! A method returning the boolean print member. It calls the same method of DEData class.
    bool Print() const {return deData_.Print();}
    //! A method returning the integer that specifies the search algorithm type.
    UInt getSearch() const {return deData_.getSearch();}

    // Getters for mesh
    //! A method returning the mesh.
    const MeshHandler<ORDER, mydim, ndim>& getMesh() const {return mesh_;}
    // Getters for specific mesh features
    //! A method returning the number of mesh EL_NNODES. It calls the same method of MeshHandler class.
    UInt getNumNodes() const {return mesh_.num_nodes();}
    //! A method returning the number of mesh elements. It calls the same method of MeshHandler class.
    UInt getNumElements() const {return mesh_.num_elements();}
    //! A method returning a node. It calls the same method of MeshHandler class.
    Point<ndim> getPoint(Id id) const {return mesh_.getPoint(id);}
    //! A method returning an element. It calls the same method of MeshHandler class.
    Element<EL_NNODES,mydim,ndim> getElement(Id id) const {return mesh_.getElement(id);}
    //! A method returning the element in which the point in input is located. It calls the same method of MeshHandler class.
    Element<EL_NNODES,mydim,ndim> findLocation(const Point<ndim>& point) const {return mesh_.findLocation(point);}

    // Getters for matrices
    //! A method returning the penalty matrix.
    MatrixXr getP() const {return P_;}
    //! A method returning the mass matrix.
    SpMat getMass() const {return R0_;}
    //! A method returning the stiffness matrix.
    SpMat getStiffness() const {return R1_;}
    //! A method returning the PsiQuad_ matrix.
    const Eigen::Matrix<Real, Integrator::NNODES, EL_NNODES>& getPsiQuad() const {return PsiQuad_;}
    //! A method returning the GlobalPsi_ matrix.
    const SpMat& getGlobalPsi() const {return GlobalPsi_;}
};


//! @brief A class to store common (time) data for the problem.
template<UInt ORDER, UInt mydim, UInt ndim>
class DataProblem_time : public DataProblem<ORDER, mydim, ndim>{
private:
    using Integrator = typename DensityIntegratorHelper::Integrator<mydim>;

    // Data member related to time discretization.
    static const UInt SPLINE_DEGREE = 3;
    static const UInt ORDER_DERIVATIVE = 2;
    using Integrator_t = IntegratorGaussP9;
    Spline<SPLINE_DEGREE, ORDER_DERIVATIVE> spline_;

    // Time data (stored in a DEData_time object).
    DEData_time deData_time_;
    // Time mesh.
    std::vector<Real> mesh_time_;
    // Matrix containing evaluations of the spline basis functions at different time instants; time mass matrix.
    SpMat GlobalPhi_, K0_;
    // Time and space penalty matrices.
    SpMat Pt_, Ps_;
    // Matrix obtained by computing the Kronecker product of GlobalPhi_ and GlobalPsi_ matrices.
    SpMat Upsilon_;
    // Data structure used during the initialization procedure via discretized heat diffusion process.
    std::vector<std::vector<UInt>> data_Heat_;
    // Flags related to penalty matrices.
    bool flagMass_, flagLumped_;

    //! A method to compute the matrix which evaluates the temporal basis functions at the time instants.
    void fillGlobalPhi();
    //! A method to compute the time mass matrix.
    void fillTimeMass();
    //! A method to compute the time penalty matrix.
    void fillTimeSecondDerivative();
    //! A method to build space and time mass matrices according to the mass lumping procedure.
    SpMat makeLumped(const SpMat& mass) const;
    //! A method to compute the space penalty matrix.
    void fillPenaltySpace();
    //! A method to compute the time penalty matrix.
    void fillPenaltyTime();
    //! A method to generate the data structure used during the initialization procedure via discretized heat diffusion process.
    void setDataHeat();

public:
    //! A constructor: it delegates DataProblem, DEData_time and Spline constructors.
    DataProblem_time(SEXP Rdata, SEXP Rdata_time, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda,
                     SEXP Rlambda_time, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2, SEXP Rprint,
                     SEXP Rsearch, SEXP Rmesh, const std::vector<Real>& mesh_time, SEXP RisTimeDiscrete,
                     SEXP RflagMass, SEXP RflagLumped, bool isTime = 1);

    //! A method to compute the integral of a function (over the temporal domain).
    Real FEintegrate_time(const VectorXr& f) const {return (kroneckerProduct(getTimeMass(), this->getMass())*f).sum();}

    //! A method to compute the matrix which evaluates the temporal basis functions at a given time_node (for the
    //! discretization of the integral of the exponential).
    MatrixXr fillPhiQuad(UInt time_node) const;

    //! A method to compute the Upsilon_ matrix as the Kronecker product between GlobalPhi_ and GlobalPsi_.
    //! In discrete time settings, this method builds Upsilon_indices_ for the efficient extraction of Upsilon_ rows
    //! (for CV preprocessing).
    SpMat computeUpsilon(const SpMat& phi, const SpMat& psi) const;
    //! A method to compute the Upsilon_ matrix by considering only locations and times in the positions stored in
    //! indices. This method is needed for CV preprocessing (only points that are in the considered fold are used).
    SpMat computeUpsilon(const std::vector<UInt>& indices) const;

    //! A method returning the matrix needed for the penalizing term in space.
    const SpMat computePen_s() const {return Ps_;}
    //! A method returning the matrix needed for the penalizing term in time.
    const SpMat computePen_t() const {return Pt_;}

    // Getters
    //! A method to access the data. It calls the same method of DEData_time class.
    const std::vector<Real>& data_time() const {return deData_time_.data();}
    //! A method returning a datum.
    Real data_time(UInt i) const {return data_time()[i];}
    //! A method returning the spline degree.
    UInt getSplineDegree() const {return SPLINE_DEGREE;}
    //! A method returning the total number of B-splines basis functions.
    UInt getSplineNumber () const {return spline_.num_knots()-SPLINE_DEGREE-1;}
    //! A method returning the penalization parameters (in time). It calls the same method of DEData_time class.
    Real getLambda_time(UInt i) const {return deData_time_.getLambda_time(i);}
    //! A method returning the number of lambdas (in time). It calls the same method of DEData_time class.
    UInt getNlambda_time()  const {return deData_time_.getNlambda_time();}
    //! A method returning the vector of the point indices that are active for the j-th B-spline.
    const std::vector<UInt>& getDataIndex_Heat(UInt j) const {return data_Heat_[j];}

    // Getters for matrices
    //! A method returning the time penalty matrix.
    SpMat getPt() const {return Pt_;}
    //! A method returning the matrix containing the evaluations of the spline basis functions.
    const SpMat& getGlobalPhi() const {return GlobalPhi_;}
    //! A method returning the Upsilon_ matrix.
    const SpMat& getUpsilon() const {return Upsilon_;}
    //! A method returning the time mass matrix.
    SpMat getTimeMass() const {return K0_;}

    // Getters for time mesh
    //! A method returning the time mesh.
    const std::vector<Real>& getMesh_time() const {return mesh_time_;}
    //! A method returning the number of time mesh nodes.
    UInt getNumNodes_time() const {return mesh_time_.size();}
};

#include "Data_Problem_imp.h"

#endif
