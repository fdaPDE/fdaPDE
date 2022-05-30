/*
    =======================================================================
    =======================================================================
    ==                                                                   ==
    ==  FSPAI:  Factorized SPAI algorithm to compute a Factorized SParse ==
    ==          Approximate Inverse matrix for symmetric positive        ==
    ==          definite systems.                                        ==
    ==                                                                   ==
    ==  Copyright (C)  2011 by                                           ==
    ==                 Matous Sedlacek <sedlacek@in.tum.de>              ==
    ==                 Chair of Scientific Computing -- Informatics V    ==
    ==                 Technische Universität München                    ==
    ==                                                                   ==
    ==  This file is part of FSPAI.                                      ==
    ==                                                                   ==
    ==  FSPAI is free software: you can redistribute it and/or           ==
    ==  modify it under the terms of the GNU Lesser General Public       ==
    ==  License as published by the Free Software Foundation, either     ==
    ==  version 3 of the License, or (at your option) any later version. ==
    ==                                                                   ==
    ==  FSPAI is distributed in the hope that it will be useful,         ==
    ==  but WITHOUT ANY WARRANTY; without even the implied warranty of   ==
    ==  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    ==
    ==  GNU Lesser General Public License for more details.              ==
    ==                                                                   ==
    ==  You should have received a copy of the GNU Lesser General        ==
    ==  Public License along with FSPAI.                                 ==
    ==  If not, see <http://www.gnu.org/licenses/>.                      ==
    ==                                                                   ==
    =======================================================================
    =======================================================================
*/
///////////////////////////////////////////////////////////////
///     \brief The external lapack functions
///
///     The lapack library provides huge contingent
///     of functions for computing mathematical operations.
///     The are optimized for speed and memory. Using these
///     library functions makes it possible to efficiently
///     run the mathematical operations needed for best
///     performance and to concentrate on the code around.
///     See the lapack documentation for details.
///////////////////////////////////////////////////////////////
extern "C"
{
    ///////////////////////////////////////////////////////////////////////
    ///     \brief Solving real SPD system of linear equations.
    ///
    ///     DPOSV computes the solution to a real system of linear equations
    ///     A * X = B, where A is an N-by-N symmetric positive definite matrix
    ///     and X and B are N-by-NRHS matrices.
    ///
    ///     \see http://www.netlib.org/lapack/double/dposv.f
    ///
    ///     \param UPLO 'U':  A is upper triangular;'L':  A is lower triangular.
    ///     \param N  The number of columns of the matrix A.  N >= 0.
    ///     \param NHRS The number of columns of the matrix B.  NRHS >= 0.
    ///     \param A  The M-by-N matrix A.
    ///     \param LDA The leading dimension of the array A.  LDA >= max(1,N).
    ///     \param B The right hand side matrix B.
    ///     \param LDB The leading dimension of the array B.  LDB >= max(1,N).
    ///     \param INFO  Function return value
    ////////////////////////////////////////////////////////////////////////
    void    dposv_(const char* UPLO,
                   int*        N,
                   int*        NHRS,
                   double      A[],
                   int*        LDA,
                   double      B[],
                   int*        LDB,
                   int*        INFO);

    ///////////////////////////////////////////////////////////////////////
    ///     \brief Solving dot product of double precision arrays.
    ///
    ///     \see http://www.netlib.org/blas/ddot.f
    ///
    ///     \param N  The dimension of the vectors.  N >= 0.
    ///     \param DX  The first vector.
    ///     \param INCX Increment between elements of x.
    ///     \param DY The second vector.
    ///     \param INCY Increment between elements of y.
    ///     \return The dot product value
    ////////////////////////////////////////////////////////////////////////
    double  ddot_(  int*    N,
                    double  DX[],
                    int*    INCX,
                    double  DY[],
                    int*    INCY);

    ///////////////////////////////////////////////////////////////////////
    ///     \brief Copying double precision vector.
    ///
    ///     \see http://www.netlib.org/blas/dcopy.f
    ///
    ///     \param N  The dimension of the vectors.  N >= 0.
    ///     \param DX  The first vector.
    ///     \param INCX Increment between elements of x.
    ///     \param DY The second vector.
    ///     \param INCY Increment between elements of y.
    ////////////////////////////////////////////////////////////////////////
    void    dcopy_( int*    N,
                    double  DX[],
                    int*    INCX,
                    double  DY[],
                    int*    INCY);

    ///////////////////////////////////////////////////////////////////////
    ///     \brief Scales a double precision vector by a constant.
    ///
    ///     \see http://www.netlib.org/blas/dscal.f
    ///
    ///     \param N  The dimension of the vector.  N >= 0.
    ///     \param DA  Scalar value to perform scaling with.
    ///     \param DX The vector to be scaled.
    ///     \param INCX Increment between elements of x.
    ////////////////////////////////////////////////////////////////////////
    void    dscal_( int*    N,
                    double* DA,
                    double  DX[],
                    int*    INCX);

    ///////////////////////////////////////////////////////////////////////
    ///     \brief Constant times a vector plus a vector (double precision)
    ///
    ///     \see http://www.netlib.org/blas/daxpy.f
    ///
    ///     \param N  The dimension of the vectors.  N >= 0.
    ///     \param DA The constant value to be multiplied with.
    ///     \param DX  The first vector.
    ///     \param INCX Increment between elements of x.
    ///     \param DY The second vector.
    ///     \param INCY Increment between elements of y.
    ////////////////////////////////////////////////////////////////////////
    void    daxpy_( int*    N,
                    double* DA,
                    double  DX[],
                    int*    INCX,
                    double  DY[],
                    int*    INCY);

    //============================================================================
    //============================================================================
    //=============== Template specifications for COMPLEX matrices ===============
    //============================================================================
    //============================================================================

    ///////////////////////////////////////////////////////////////////////
    ///     \brief Solving complex HPD system of linear equations
    ///
    ///     ZPOSV computes the solution to a real system of linear equations
    ///     A * X = B, where A is an N-by-N Hermitian positive definite matrix
    ///     and X and B are N-by-NRHS matrices.
    ///
    ///     \see http://www.netlib.org/lapack/double/zposv.f
    ///
    ///     \param UPLO 'U':  A is upper triangular;'L':  A is lower triangular.
    ///     \param N  The number of columns of the matrix A.  N >= 0.
    ///     \param NHRS The number of columns of the matrix B.  NRHS >= 0.
    ///     \param A  The M-by-N matrix A.
    ///     \param LDA The leading dimension of the array A.  LDA >= max(1,N).
    ///     \param B The right hand side matrix B.
    ///     \param LDB The leading dimension of the array B.  LDB >= max(1,N).
    ///     \param INFO  Function return value
    ////////////////////////////////////////////////////////////////////////
    void    zposv_(const char* UPLO,
                   int*        N,
                   int*        NHRS,
                   COMPLEX     A[],
                   int*        LDA,
                   COMPLEX     B[],
                   int*        LDB,
                   int*        INFO);

    ///////////////////////////////////////////////////////////////////////
    ///     \brief Solving dot product of COMPLEX precision arrays.
    ///
    ///     \see http://www.netlib.org/blas/zdotu.f
    ///
    ///     \param N  The dimension of the vectors.  N >= 0.
    ///     \param DX  The first vector.
    ///     \param INCX Increment between elements of x.
    ///     \param DY The second vector.
    ///     \param INCY Increment between elements of y.
    ///     \return The dot product value
    ////////////////////////////////////////////////////////////////////////
    COMPLEX zdotu_( int*    N,
                    COMPLEX DX[],
                    int*    INCX,
                    COMPLEX DY[],
                    int*    INCY);

    ///////////////////////////////////////////////////////////////////////
    ///     \brief Copying COMPLEX precision vector.
    ///
    ///     \see http://www.netlib.org/blas/zcopy.f
    ///
    ///     \param N  The dimension of the vectors.  N >= 0.
    ///     \param DX  The first vector.
    ///     \param INCX Increment between elements of x.
    ///     \param DY The second vector.
    ///     \param INCY Increment between elements of y.
    ////////////////////////////////////////////////////////////////////////
    void    zcopy_( int*    N,
                    COMPLEX DX[],
                    int*    INCX,
                    COMPLEX DY[],
                    int*    INCY);

    ///////////////////////////////////////////////////////////////////////
    ///     \brief Scales a COMPLEX vector by a constant.
    ///
    ///     \see http://www.netlib.org/blas/zscal.f
    ///
    ///     \param N  The dimension of the vector.  N >= 0.
    ///     \param DA  Scalar value to perform scaling with.
    ///     \param DX The vector to be scaled.
    ///     \param INCX Increment between elements of x.
    ////////////////////////////////////////////////////////////////////////
    void    zscal_( int*     N,
                    COMPLEX* DA,
                    COMPLEX  DX[],
                    int*    INCX);

    ///////////////////////////////////////////////////////////////////////
    ///     \brief Constant times a vector plus a vector (COMPLEX)
    ///
    ///     \see http://www.netlib.org/blas/zaxpy.f
    ///
    ///     \param N  The dimension of the vectors.  N >= 0.
    ///     \param DA The constant value to be multiplied with.
    ///     \param DX  The first vector.
    ///     \param INCX Increment between elements of x.
    ///     \param DY The second vector.
    ///     \param INCY Increment between elements of y.
    ////////////////////////////////////////////////////////////////////////
    void    zaxpy_( int*     N,
                    COMPLEX* DA,
                    COMPLEX  DX[],
                    int*     INCX,
                    COMPLEX  DY[],
                    int*     INCY);
}
