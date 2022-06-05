#ifndef __MIXED_FE_FPCA_IMP_H__
#define __MIXED_FE_FPCA_IMP_H__

#include <iostream>
#include<iterator>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>
#include "../../Global_Utilities/Include/Timing.h"
#include <fstream>

#include "R_ext/Print.h"

template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase::computeDelta(const MeshHandler<ORDER, mydim, ndim> & mesh)
{
	UInt nRegions = fpcaData_.getNumberOfRegions();
	Delta_.resize(nRegions,1);
	for (int i=0; i<nRegions; i++)
	{
		Delta_(i)=0;
		for (int j=0; j<fpcaData_.getIncidenceMatrix().cols(); j++)
		{
			if (fpcaData_.getIncidenceMatrix()(i,j) == 1)
			{
				Delta_(i)+=mesh.elementMeasure(j);
			}
		}
	}
}

template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase::computeBasisEvaluations(const MeshHandler<ORDER, mydim, ndim> & mesh)
{
	//std::cout<<"Data Matrix Computation by Basis Evaluation.."<<std::endl;
	UInt nnodes = mesh.num_nodes();
	UInt nlocations = fpcaData_.getNumberofObservations();
	static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 100 * eps;

	Psi_.resize(nlocations, nnodes);
	if (fpcaData_.isLocationsByNodes()) //pointwise data
	{
		std::vector<coeff> tripletAll;
		auto k = fpcaData_.getObservationsIndices();

		tripletAll.reserve(k.size());
		for (int i = 0; i< k.size(); ++i)
		{
			tripletAll.push_back(coeff(i,k[i],1.0));
		}
		Psi_.setFromTriplets(tripletAll.begin(),tripletAll.end());
		Psi_.makeCompressed();
	}
	else if (fpcaData_.isLocationsByBarycenter() && (fpcaData_.getNumberOfRegions()==0))
	{
    	static constexpr UInt EL_NNODES = how_many_nodes(ORDER,mydim);

		for(UInt i=0; i<nlocations;i++)
		{

			Element<EL_NNODES, mydim, ndim> tri_activated = mesh.getElement(fpcaData_.getElementId(i));

			if(tri_activated.getId() == Identifier::NVAL)
			{
				Rprintf("WARNING: Observation %d is not in the domain, remove point and re-perform smoothing\n", i+1);
			}
			else
			{
				for(UInt node = 0; node < EL_NNODES ; ++node)
				{
					Real evaluator = fpcaData_.getBarycenter(i,node);
					Psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		} //end of for loop

		Psi_.prune(tolerance);
		Psi_.makeCompressed();
	}
	else if ((!fpcaData_.isLocationsByBarycenter()) && fpcaData_.getNumberOfRegions()==0)
	{
    	static constexpr UInt EL_NNODES = how_many_nodes(ORDER,mydim);

		this->barycenters_.resize(nlocations, EL_NNODES);
		this->element_ids_.resize(nlocations);
		for(UInt i=0; i<nlocations;i++)
		{

			Element<EL_NNODES, mydim, ndim> tri_activated = mesh.findLocation(fpcaData_.template getLocations<ndim>(i));

			if(tri_activated.getId() == Identifier::NVAL)
			{
				Rprintf("WARNING: Observation %d is not in the domain, remove point and re-perform smoothing\n", i+1);
			}
			else
			{
				element_ids_(i)=tri_activated.getId();
				for(UInt node = 0; node < EL_NNODES ; ++node)
				{
					Real evaluator = tri_activated.evaluate_point(fpcaData_.template getLocations<ndim>(i), Eigen::Matrix<Real,EL_NNODES,1>::Unit(node));
					barycenters_(i,node)=evaluator;
					Psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		} //end of for loop

		Psi_.prune(tolerance);
		Psi_.makeCompressed();
	}
	else //areal data
	{
    	static constexpr UInt EL_NNODES = how_many_nodes(ORDER,mydim);

		Real *tab; //Psi_i
		tab = (Real*) malloc(sizeof(Real)*nnodes);
		for(UInt i=0; i<nlocations;i++) //nlocations = number of regions
		{
			for (UInt k=0; k<nnodes; k++) {tab[k]=0;}
			for (UInt j=0; j<mesh.num_elements(); j++)
			{
				if (fpcaData_.getIncidenceMatrix()(i,j) == 1) //element j is in region i
				{
					Element<EL_NNODES, mydim, ndim> tri = mesh.getElement(j); //can also be a tetrahedron
					for (UInt k=0; k<EL_NNODES; k++)
					{
						tab[tri[k].getId()] += tri.integrate(Eigen::Matrix<Real,EL_NNODES,1>::Unit(k)); // integral over tri of psi_k	
					}				
				}
			}
			for (int k=0; k<nnodes; k++)
			{
				if (tab[k] != 0)
				{
					Psi_.insert(i,k) = tab[k]/Delta_(i);
				}
			}
		}
		free(tab);
		Psi_.makeCompressed();
	}
}

template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase::SetAndFixParameters(const MeshHandler<ORDER, mydim, ndim> & mesh)
{
	FiniteElement<ORDER, mydim, ndim> fe;
	this->nnodes_ = mesh.num_nodes();

	this-> template computeDelta<ORDER,mydim,ndim>(mesh);
	this-> template computeBasisEvaluations<ORDER,mydim,ndim>(mesh); //compute Psi
	computeDataMatrix(DMat_); //NW block

	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	Assembler::operKernel(stiff, mesh, fe, AMat_);
	Assembler::operKernel(mass, mesh, fe, MMat_);


	/*const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,Eigen::DontAlignCols,", ","\n");

	std::string mass_name("Mass3D.csv");
	std::string stiff_name("Stiff3D.csv");

	std::ofstream file_mass(mass_name.c_str());
	std::ofstream file_stiff(stiff_name.c_str());

	file_mass << MatrixXr(MMat_).format(CSVFormat);
	file_stiff << MatrixXr(AMat_).format(CSVFormat);*/

	scores_mat_.resize(fpcaData_.getNPC());
	loadings_mat_.resize(fpcaData_.getNPC());
	lambda_PC_.resize(fpcaData_.getNPC());


	datamatrixResiduals_ = fpcaData_.getDatamatrix();
	solution_.resize(fpcaData_.getLambda().size());
}



#endif
