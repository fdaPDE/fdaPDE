#ifndef __KRONECKER_PRODUCT_H__
#define __KRONECKER_PRODUCT_H__

#include "../../FdaPDE.h"

SpMat kroneckerProduct(const SpMat&, const SpMat&);

MatrixXr kroneckerProduct_Matrix (const MatrixXr&, const MatrixXr&);

#endif
