#ifndef __DE_DATA_IMP_H__
#define __DE_DATA_IMP_H__

template<UInt ndim>
DEData<ndim>::DEData(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals,
  SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch)
{
  setData(Rdata);

  order_ = INTEGER(Rorder)[0];

  setFvec(Rfvec);

  heatStep_ = REAL(RheatStep)[0];

  heatIter_ = INTEGER(RheatIter)[0];

  setLambda(Rlambda);

  Nfolds_ = INTEGER(Rnfolds)[0];

  nsim_ = INTEGER(Rnsim)[0];

  setStepProposals(RstepProposals);

  tol1_ = REAL(Rtol1)[0];

  tol2_ = REAL(Rtol2)[0];

  print_ = INTEGER(Rprint)[0];

  search_ = INTEGER(Rsearch)[0];

}

template<UInt ndim>
DEData<ndim>::DEData(const std::vector<Point<ndim> >& data, const UInt& order, const VectorXr& fvec, Real heatStep, UInt heatIter, const std::vector<Real>& lambda,
               const UInt& nfolds, const UInt& nsim, const std::vector<Real>& stepProposals, Real tol1, Real tol2,
               bool print, UInt search):
                data_(data), order_(order), fvec_(fvec), heatStep_(heatStep), heatIter_(heatIter), lambda_(lambda), Nfolds_(nfolds),
                nsim_(nsim), stepProposals_(stepProposals), tol1_(tol1), tol2_(tol2), print_(print), search_(search)
{
}

template<UInt ndim>
void DEData<ndim>::setData(SEXP Rdata)
{
  const RNumericMatrix data(Rdata);
  UInt n_=data.nrows();
  if(n_ > 0){
    data_.reserve(n_);
    for(int i = 0; i < n_; ++i)
      data_.emplace_back(i, data);
  }
}

template<UInt ndim>
void DEData<ndim>::setFvec(SEXP Rfvec)
{
  UInt dimc = Rf_length(Rfvec);
  fvec_.resize(dimc);
  for(UInt i = 0; i < dimc; ++i)
  {
      fvec_[i] = REAL(Rfvec)[i];
  }
}

template<UInt ndim>
void DEData<ndim>::setLambda(SEXP Rlambda)
{
  UInt diml = Rf_length(Rlambda);
  lambda_.reserve(diml);
  for(UInt i = 0; i < diml; ++i)
  {
      lambda_.push_back(REAL(Rlambda)[i]);
  }
}

template<UInt ndim>
void DEData<ndim>::setStepProposals(SEXP RstepProposals)
{
  UInt dimPG = Rf_length(RstepProposals);
  stepProposals_.reserve(dimPG);
  for(UInt i = 0; i < dimPG; ++i)
  {
      stepProposals_.push_back(REAL(RstepProposals)[i]);
  }
}

template<UInt ndim>
void DEData<ndim>::printData(std::ostream & out) const
{
  for(int i = 0; i < data_.size(); ++i)
	{
		out<<data_[i]<<std::endl;
	}
}

// -------------------------------------------
// --------------- DEData_time ---------------
// -------------------------------------------

DEData_time::DEData_time(SEXP Rdata_time, SEXP Rlambda_time)
{
    setDataTime(Rdata_time);

    setLambdaTime(Rlambda_time);

}

DEData_time::DEData_time(const std::vector<Real>& data_time, const std::vector<Real>& lambda_time):
        data_time_(data_time), lambda_time_(lambda_time)
{
}

void DEData_time::setDataTime(SEXP Rdata_time)
{
    UInt diml = Rf_length(Rdata_time);
    data_time_.reserve(diml);
    for(UInt i = 0; i < diml; ++i)
    {
        data_time_.push_back(REAL(Rdata_time)[i]);
    }
}

void DEData_time::setLambdaTime(SEXP Rlambda_time)
{
    UInt diml = Rf_length(Rlambda_time);
    lambda_time_.reserve(diml);
    for(UInt i = 0; i < diml; ++i)
    {
        lambda_time_.push_back(REAL(Rlambda_time)[i]);
    }
}

void DEData_time::setTimes2Locations()
{
    // Extraction of non-duplicated times
    std::set<Real> set_times(data_time_.cbegin(), data_time_.cend());
    times_.resize(set_times.size());
    (void)std::copy(set_times.begin(), set_times.end(), times_.begin());
    set_times.clear();

    // Creation of the data structure containing for each time index i, the location indices that are observed at time i
    // (ordered as data appear in the clean dataset).
    Times2Locations_.resize(times_.size());
    for(UInt t = 0; t < times_.size(); ++t)
    {
        for(UInt d = 0; d < data_time_.size(); ++d) {
            if(data_time_[d] == times_[t])
                Times2Locations_[t].push_back(d);
        }
    }
}

void DEData_time::printData(std::ostream& out) const
{
    for(int i = 0; i < data_time_.size(); ++i)
    {
        out << data_time_[i] << std::endl;
    }
}

void DEData_time::printTimes2Locations(std::ostream& out) const
{
    for(std::vector<std::vector<UInt>>::const_iterator it = Times2Locations_.cbegin(); it != Times2Locations_.cend(); ++it) {
        out << "time index: " << it-Times2Locations_.cbegin() << '\t' << "[location index]:";
        for(UInt i : *it)
            out << " " << i;
        out << std::endl;
    }
}

#endif
