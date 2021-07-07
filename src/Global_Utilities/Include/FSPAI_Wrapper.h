#ifndef __FSPAI_WRAPPER_H__
#define __FSPAI_WRAPPER_H__

#include <iostream>
#include "../../FdaPDE.h"
#include <unsupported/Eigen/SparseExtra>
#include <cstdio>


struct FSPAI_data{

    // From imput // TO BE RETRIEVED
    std::string tol_Inverse     = "0.05";                        // Controls the quality of approximatin, default 0.05 //TO BE FIXED FROM EXTERNAL
    std::string max_Step_Col    = "10";                          // Max number of improvement steps per columns
    std::string max_New_Nz      = "10";                          // Max number of new nonzero candidates per step
    std::string out_File;                                        // Temporary file on which will be written the inverse matrix
    std::string sol             = "0";                           // Type of solver used for the system Lx=rhs, 0 means no solver is used
    
  };


//! Wrapper for FSPAI methods used to approximate the inverse of a sparse matrix
/*
  \param A the sparse matrix to be inverted 
  \param A_inv the inverse of A that will be computed by FSPAI
  \return bool
*/

bool FSPAI_Wrapper(const SpMat & A, SpMat & A_inv);



#endif
