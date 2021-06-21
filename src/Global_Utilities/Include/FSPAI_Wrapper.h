#ifndef __FSPAI_WRAPPER_H__
#define __FSPAI_WRAPPER_H__

#include <iostream>
#include "../../fdaPDE.h"
#include <unsupported/Eigen/SparseExtra>
#include <cstdio>

//! Wrapper for FSPAI methods used to approximate the inverse of a sparse matrix
/*
  \param A the sparse matrix to be inverted 
  \param A_inv the inverse of A that will be computed by FSPAI
  \return bool
*/

bool FSPAI_Wrapper(const SpMat & A, SpMat & A_inv);

#endif
