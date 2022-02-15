#ifndef __DATA_PROBLEM_IMP_H__
#define __DATA_PROBLEM_IMP_H__

template<UInt ORDER, UInt mydim, UInt ndim>
DataProblem<ORDER, mydim, ndim>::DataProblem(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter,
                                             SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1,
                                             SEXP Rtol2, SEXP Rprint, SEXP Rsearch, SEXP Rmesh, bool isTime):
  deData_(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rsearch),
  mesh_(Rmesh, INTEGER(Rsearch)[0]){

    std::vector<Point<ndim>>& data = deData_.data();

    // PROJECTION
    if((mydim == 2 && ndim == 3) || (mydim == 1 && ndim == 2)){
      Rprintf("##### DATA PROJECTION #####\n");
      projection<ORDER, mydim, ndim> projection(mesh_, data);
      data = projection.computeProjection();
    }
    
    // REMOVE POINTS NOT IN THE DOMAIN
    if(!isTime) {
        for(auto it = data.begin(); it != data.end(); ){
            Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.findLocation(data[it - data.begin()]);
            if(tri_activated.getId() == Identifier::NVAL)
            {
                it = data.erase(it);
                Rprintf("WARNING: an observation is not in the domain. It is removed and the algorithm proceeds.\n");
            }
            else {
                ++it;
            }
        }
    }

    // FILL SPACE MATRICES
    fillFEMatrices();
    fillPsiQuad();

    if(!isTime) {
        std::vector <UInt> v(deData_.dataSize());
        std::iota(v.begin(), v.end(), 0);
        GlobalPsi_ = computePsi(v);
    }
}


template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem<ORDER, mydim, ndim>::fillFEMatrices(){

  //fill R0 and R1
  FiniteElement<ORDER, mydim, ndim> fe;
  typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
  typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
  Assembler::operKernel(mass, mesh_, fe, R0_);
  Assembler::operKernel(stiff, mesh_, fe, R1_);

  //fill P
  Eigen::SparseLU<SpMat> solver;
	solver.compute(R0_);
	auto X2 = solver.solve(R1_);
	P_ = R1_.transpose()* X2;
}


template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem<ORDER, mydim, ndim>::fillPsiQuad(){
	for(UInt i=0; i<Integrator::NNODES; ++i)
	   PsiQuad_.row(i)=reference_eval_point<EL_NNODES, mydim>(Integrator::NODES[i]);
}


template<UInt ORDER, UInt mydim, UInt ndim>
Real DataProblem<ORDER, mydim, ndim>::FEintegrate_exponential(const VectorXr& g) const{

  using EigenMap2WEIGHTS = Eigen::Map<const Eigen::Matrix<Real, Integrator::NNODES, 1> >;

  Real total_sum = 0.;

  for(UInt triangle=0; triangle<mesh_.num_elements(); ++triangle){

    Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.getElement(triangle);

// (3) -------------------------------------------------
    Eigen::Matrix<Real,EL_NNODES,1> sub_g;
    for (UInt i=0; i<EL_NNODES; i++){
      sub_g[i]=g[tri_activated[i].getId()];
    }

// (4) -------------------------------------------------
    Eigen::Matrix<Real,Integrator::NNODES,1> expg = (PsiQuad_*sub_g).array().exp();

    total_sum+=expg.dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0]))*tri_activated.getMeasure();

  }

  return total_sum;
}


template<UInt ORDER, UInt mydim, UInt ndim>
SpMat
DataProblem<ORDER, mydim, ndim>::computePsi(const std::vector<UInt>& indices) const{

    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    UInt nnodes = mesh_.num_nodes();
    UInt nlocations = indices.size();
	SpMat psi(nlocations, nnodes);

	std::vector<coeff> triplets;
	triplets.reserve(EL_NNODES*nlocations);

	for(auto it = indices.cbegin(); it != indices.cend(); ++it)
	{
        Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.findLocation(deData_.data(*it));

		if(tri_activated.getId() == Identifier::NVAL)
		{
			Rprintf("WARNING: the following observation is not in the domain\n");
            //operator<<(std::cout, deData_.data(*it));
		}
        else
        {
			for(UInt node = 0; node < EL_NNODES ; ++node)
			{
				Real evaluator = tri_activated.evaluate_point(deData_.data(*it), Eigen::Matrix<Real,EL_NNODES,1>::Unit(node));
				triplets.emplace_back(it-indices.cbegin(), tri_activated[node].getId(), evaluator);
			}
		}
	}

	psi.setFromTriplets(triplets.begin(),triplets.end());

	psi.prune(tolerance);
	psi.makeCompressed();

	return psi;
}

// ------------------------------------------------
// --------------- DataProblem_time ---------------
// ------------------------------------------------

template<UInt ORDER, UInt mydim, UInt ndim>
DataProblem_time<ORDER, mydim, ndim>::DataProblem_time(SEXP Rdata, SEXP Rdata_time, SEXP Rorder, SEXP Rfvec, SEXP RheatStep,
                                                       SEXP RheatIter, SEXP Rlambda, SEXP Rlambda_time, SEXP Rnfolds, SEXP Rnsim,
                                                       SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch,
                                                       SEXP Rmesh, const std::vector<Real>& mesh_time, SEXP RisTimeDiscrete,
                                                       SEXP RflagMass, SEXP RflagLumped, bool isTime):
  DataProblem<ORDER, mydim, ndim>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals,
                                    Rtol1, Rtol2, Rprint, Rsearch, Rmesh, isTime),
  deData_time_(Rdata_time, Rlambda_time), mesh_time_(mesh_time), spline_(mesh_time) {

    flagMass_ = INTEGER(RflagMass)[0];
    flagLumped_ = INTEGER(RflagLumped)[0];

    std::vector<Point<ndim>>& data_ = this->deData_.data();
    std::vector<Real>& data_time_ = deData_time_.data();
    const Real t_min = mesh_time_.front();
    const Real t_max = mesh_time_.back();

    // REMOVE POINTS NOT IN THE DOMAIN
    for (auto it = data_.begin(); it != data_.end();) {
        Element<this->EL_NNODES, mydim, ndim> tri_activated = this->mesh_.findLocation(data_[it - data_.begin()]);
        if (tri_activated.getId() == Identifier::NVAL || (data_time_[it - data_.begin()] < t_min || data_time_[it - data_.begin()] > t_max)) {
            data_time_.erase(data_time_.begin() + (it - data_.begin()));
            it = data_.erase(it);
            Rprintf("WARNING: an observation is not in the domain. It is removed and the algorithm proceeds.\n");
        } else {
            ++it;
        }
    }

    Rprintf("WARNING: %d observations used in the algorithm.\n", data_.size());

    // FILL SPACE MATRIX
    std::vector <UInt> v(this->deData_.dataSize());
    std::iota(v.begin(), v.end(), 0);
    this->GlobalPsi_ = this->computePsi(v);

    // DISCRETE TIME DATA
    if(static_cast<bool> (INTEGER(RisTimeDiscrete)[0])) {
        deData_time_.setTimes2Locations();
        //deData_time_.printTimes2Locations(std::cout);

        // DATA STRUCTURE FOR EFFICIENT UPSILON COMPUTATION (USEFUL DURING CV PREPROCESSING)
        Upsilon_indices_.resize(deData_time_.getNTimes());
    }

    // DATA STRUCTURE FOR HEAT INITIALIZATION
    if(this->isFvecEmpty())
        setDataHeat();

    // FILL TIME MATRICES
    fillGlobalPhi();
    fillTimeMass();
    fillTimeSecondDerivative();

    // COMPUTE SPACE AND TIME PENALTY MATRICES
    fillPenaltySpace();
    fillPenaltyTime();

    // ASSEMBLE SPACE-TIME MATRICES
    Upsilon_ = computeUpsilon(GlobalPhi_, this->GlobalPsi_);
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillGlobalPhi()
{
    const UInt M = getSplineNumber();
    const UInt m = deData_time_.getNTimes();

    GlobalPhi_.resize(m, M);
    Real value;

    for(UInt i = 0; i < m; ++i)
    {
        for(UInt j = 0; j < M; ++j)
        {
            value = spline_.BasisFunction(j, this->deData_time_.time(i));
            if(value != 0)
            {
                GlobalPhi_.coeffRef(i,j) = value;
            }
        }
    }

    GlobalPhi_.makeCompressed();
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillTimeMass(void)
{
    Assembler::operKernel(spline_, K0_);
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillTimeSecondDerivative(void)
{
    Assembler::operKernel(spline_, Pt_);
}

template<UInt ORDER, UInt mydim, UInt ndim>
SpMat DataProblem_time<ORDER, mydim, ndim>::makeLumped(const SpMat& mass) const
{
    VectorXr diag = mass * VectorXr::Ones(mass.cols());
    SpMat lumped_mass(diag.asDiagonal());

    return mass;
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillPenaltySpace(void)
{
    // Update R1_
    SpMat R1_temp;
    UInt N_r = this->R1_.rows(), N_c = this->R1_.cols();
    UInt M_r = K0_.rows(), M_c = K0_.cols();

    SpMat K0temp(K0_);
    if(!flagMass_)
        K0temp.setIdentity();

    R1_temp = kroneckerProduct(K0temp, this->getStiffness());
    R1_temp.makeCompressed();

    // Update R0
    SpMat R0_temp;
    R0_temp = kroneckerProduct(K0temp, this->getMass());
    R0_temp.makeCompressed();

    if(flagLumped_)
        R0_temp = makeLumped(R0_temp);

    // Compute Space Penalty
    Ps_.resize(N_r * M_r, N_c * M_c);
    Eigen::SparseLU<SpMat> factorized_R0(R0_temp);
    Ps_ = (R1_temp).transpose()*factorized_R0.solve(R1_temp);     // Ps_ == R1_^t * R0_^{-1} * R1_
    Ps_.makeCompressed();
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillPenaltyTime() {
    SpMat mass_temp(this->getMass());
    if(!flagMass_)
        mass_temp.setIdentity();
    Pt_ = kroneckerProduct(getPt(),mass_temp);
    Pt_.makeCompressed();
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::setDataHeat()
{
    const UInt M = getSplineNumber();
    data_Heat_.resize(M);

    //! ### POSSIBLE PARALLELIZATION via openMP ###
    for (int i = 0; i < deData_time_.getNTimes(); ++i) {
        for (int j = 0; j < M; ++j) {
            if(spline_.BasisFunction(j, deData_time_.time(i)) != 0) //alternative: std::abs(spline_.BasisFunction(j, data_time_[i])) >= tol)
                data_Heat_[j].push_back(i);
        }
    }
}

template<UInt ORDER, UInt mydim, UInt ndim>
MatrixXr DataProblem_time<ORDER, mydim, ndim>::fillPhiQuad(UInt time_node) const
{
    MatrixXr phi;
    phi.resize(Integrator_t::NNODES, SPLINE_DEGREE+1);
    Real t_a = mesh_time_[time_node], t_b = mesh_time_[time_node+1];
    std::array<Real, Integrator_t::NNODES> ref_nodes;
    for(UInt k = 0; k < Integrator_t::NNODES; ++k)
        ref_nodes[k] = ((t_b-t_a) * Integrator_t::NODES[k] + t_a + t_b) / 2;
    for(UInt j = 0; j < phi.cols(); j++){
        for(UInt i = 0; i < phi.rows(); i++)
            phi(i,j) = spline_.BasisFunction(time_node+j, ref_nodes[i]);
    }
    return phi;
}


template<UInt ORDER, UInt mydim, UInt ndim>
SpMat DataProblem_time<ORDER, mydim, ndim>::computeUpsilon(const SpMat& phi, const SpMat& psi)
{
    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    const UInt phi_r = phi.rows();
    const UInt phi_c = phi.cols();
    const UInt psi_c = psi.cols();

    if (deData_time_.getNTimes() != deData_time_.dataSize() & this->Print()) {
        //Rprintf("WARNING: %d temporal duplicates.\n", deData_time_.dataSize() - deData_time_.getNTimes());
        Rprintf("%d distinct time instants.\n", deData_time_.getNTimes());
    }

    std::vector<coeff> Upsilon_tripletList;
    Upsilon_tripletList.reserve(deData_time_.dataSize() * phi_c * psi_c);

    if(deData_time_.getNTimes() != deData_time_.dataSize()) // time duplicates: phi_r < psi_r
    {
        UInt global_row_counter = 0;
        for(UInt i = 0; i < phi_r; ++i) {
            const std::vector<UInt>& v = deData_time_.getTimes2Locations(i);
            for(UInt j : v) {
                Upsilon_indices_[j] = global_row_counter;
                SpMat localKProd_(1, phi_c * psi_c);
                localKProd_ = kroneckerProduct(phi.row(i), psi.row(j));
                for (UInt idx = 0; idx < localKProd_.outerSize(); ++idx)
                    Upsilon_tripletList.emplace_back(global_row_counter,idx,localKProd_.coeff(0, idx));
                ++global_row_counter;
            }
        }
    }
    else // NO time duplicates: phi_r = psi_r = #observations
    {
        for(UInt i = 0; i < phi_r; ++i) {
            SpMat localKProd_(1, phi_c * psi_c);
            localKProd_ = kroneckerProduct(phi.row(i), psi.row(i));
            for (UInt idx = 0; idx < localKProd_.outerSize(); ++idx)
                Upsilon_tripletList.emplace_back(i,idx,localKProd_.coeff(0, idx));
        }
    }

    SpMat upsilon(deData_time_.dataSize(), phi_c * psi_c);
    upsilon.setFromTriplets(Upsilon_tripletList.begin(), Upsilon_tripletList.end());

    upsilon.prune(tolerance);
    upsilon.makeCompressed();

    return upsilon;
}

template<UInt ORDER, UInt mydim, UInt ndim>
SpMat DataProblem_time<ORDER, mydim, ndim>::computeUpsilon(const std::vector<UInt>& indices) const
{
    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    Eigen::SparseMatrix<Real, Eigen::RowMajor> upsilon(indices.size(), Upsilon_.cols());

    if(deData_time_.getNTimes() != deData_time_.dataSize()) // time duplicates
    {
        for(UInt i = 0; i < indices.size(); ++i) {
            upsilon.row(i) = Upsilon_.row(Upsilon_indices_[indices[i]]);
        }
    }
    else // NO time duplicates
    {
        for(UInt i = 0; i < indices.size(); ++i) {
            upsilon.row(i) = Upsilon_.row(indices[i]);
        }
    }

    upsilon.prune(tolerance);
    upsilon.makeCompressed();

    return upsilon;

}

/*
// ALTERNATIVE VERSION (WITHOUT Upsilon_indices_)
template<UInt ORDER, UInt mydim, UInt ndim>
SpMat DataProblem_time<ORDER, mydim, ndim>::computeUpsilon(const std::vector<UInt>& indices) const
{
    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    SpMat psi = this->computePsi(indices);
    SpMat phi = GlobalPhi_;

    const UInt phi_r = phi.rows();
    const UInt phi_c = phi.cols();
    const UInt psi_c = psi.cols();

    std::vector<coeff> Upsilon_tripletList;
    Upsilon_tripletList.reserve(indices.size() * phi_c * psi_c);

    if(deData_time_.getNTimes() != deData_time_.dataSize()) // time duplicates
    {
        UInt global_row_counter = 0;
        for(UInt i = 0; i < phi_r; ++i) {
            for(UInt j = 0; j < deData_time_.getTimes2Locations(i).size(); ++j) {
                auto it = std::find(indices.begin(), indices.end(), deData_time_.getTimes2Locations(i)[j]);
                //auto it = std::lower_bound(indices.begin(), indices.end(), deData_time_.getTimes2Locations(i)[j]);
                if(it != indices.end()) {
                    //if(*it == deData_time_.getTimes2Locations(i)[j]) {
                    SpMat localKProd_(1, phi_c * psi_c);
                    localKProd_ = kroneckerProduct(phi.row(i), psi.row(indices[it-indices.begin()]));
                    for (UInt idx = 0; idx < localKProd_.outerSize(); ++idx)
                        Upsilon_tripletList.emplace_back(global_row_counter,idx,localKProd_.coeff(0, idx));
                }
                ++global_row_counter;
            }
        }
    }
    else // NO time duplicates
    {
        for(UInt i = 0; i < indices.size(); ++i) {
            SpMat localKProd_(1, phi_c * psi_c);
            localKProd_ = kroneckerProduct(phi.row(indices[i]), psi.row(i));
            for (UInt idx = 0; idx < localKProd_.outerSize(); ++idx)
                Upsilon_tripletList.emplace_back(i,idx,localKProd_.coeff(0, idx));
        }
    }

    SpMat upsilon(indices.size(), phi_c * psi_c);
    upsilon.setFromTriplets(Upsilon_tripletList.begin(), Upsilon_tripletList.end());

    upsilon.prune(tolerance);
    upsilon.makeCompressed();

    return upsilon;

}
*/

#endif
