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

#include "../Include/Fspai_Sub.h"
#include "../Include/Lapack.h"

//============================================================================
//============================================================================
//================ Template specifications for double matrices ===============
//============================================================================
//============================================================================

template<>  void
Fspai_Sub<double>::Print_Red_Mtx
(   const double*   red_mtx,
    const int       n,
    const int       m)
{
    /*std::cout << "\n\tred_mtx:  \n\t\t";
    for (int i = 0; i < m*n; i++)
        std::cout << red_mtx[i] << " ";
    std::cout << std::endl;

    for (int i = 0; i < m; i++)
    {
        std::cout << "\n\t\t";
        for (int j = 0; j < n; j++)
            std::cout << red_mtx[i + j * m] << " ";
    }
    std::cout << "\n" << std::endl;*/
    return;
}



template<>  void
Fspai_Sub<double>::Print_Diags
(   const int dim )
{
    /*
    std::cout << "\tdiag_elements: ";
    for(int i= 0; i < dim; i++)
        std::cout << diag_elements[i].diag_el << " ";
    std::cout << std::endl;
    */
    return;

}



template<>  void
Fspai_Sub<double>::Init_Vec
(   double*&  vec,
    const int dim )
{
    for(int i= 0; i < dim; i++)
        vec[i] = 0.0;
}



template<>  void
Fspai_Sub<double>::Init_Diags
(   const int dim )
{
    for(int i = 0; i < dim; i++)
    {
       diag_elements[i].set = false;
       diag_elements[i].diag_el = 0.0;
    }
}



template<>  void
Fspai_Sub<double>::Init_Diagonal
(   double& A_kk )
{
    A_kk = 0.0;
}



template<> void
Fspai_Sub<double>::Reset_Sumvec
(   const int r_idx )
{
    sumvec[r_idx] = 0.0;
}



template<> void
Fspai_Sub<double>::Solve_HPD_System
(   double*     red_mtx,
    double*     red_col,
    int&        n,
    const char* uplo,
    int&        nrhs,
    int&        info)
{
    int lda = std::max(1, n);
    int ldb = std::max(1, n);

    // computing the solution to a real system of linear equations
    // where A is an N-by-N symmetric positive definite matrix
    dposv_( uplo,
            &n,
            &nrhs,
            red_mtx,
            &lda,
            red_col,
            &ldb,
            &info);
}



template<> double
Fspai_Sub<double>::Dot_Product
(   const double*   vec1,
    const double*   vec2,
    const int       dim)
{
    double out = 0.0;
    for(int i = 0; i < dim; i++)
        out += (vec1[i]*vec2[i]);
    return out;
}



template<> double
Fspai_Sub<double>::Field_Sqrt
(   const double number)
{
    return sqrt(number);
}

//============================================================================
//============================================================================
//================ Template specifications for COMPLEX matrices ==============
//============================================================================
//============================================================================


template<>  void
Fspai_Sub<COMPLEX>::Print_Red_Mtx
(   const COMPLEX*  red_mtx,
    const int       n,
    const int       m)
{
    /*std::cout << "\n\tred_mtx:  \n\t\t";
    for (int i = 0; i < m*n; i++)
        std::cout << red_mtx[i].real
                  << " " << red_mtx[i].imag << " | ";
    std::cout << std::endl;

    for (int i = 0; i < m; i++)
    {
        std::cout << "\n\t\t";
        for (int j = 0; j < n; j++)
            std::cout << red_mtx[i + j * m].real
                      << " " << red_mtx[i+ j * m].imag << " | ";
    }
    std::cout << "\n" << std::endl;
    */
    return;
}



template<>  void
Fspai_Sub<COMPLEX>::Print_Diags
(   const int dim )
{
    /*
    std::cout << "diag_elements: ";
    for(int i= 0; i < dim; i++)
        std::cout << diag_elements[i].diag_el.real << " "
                  << diag_elements[i].diag_el.imag << " ";
    std::cout << std::endl;
    */
    return;
}



template<>  void
Fspai_Sub<COMPLEX>::Init_Vec
(   COMPLEX*& vec,
    const int dim )
{
    for(int i= 0; i < dim; i++)
    {
        vec[i].imag = 0.0;
        vec[i].real = 0.0;
    }
}



template<>  void
Fspai_Sub<COMPLEX>::Init_Diags
(   const int dim )
{
    for(int i = 0; i < dim; i++)
    {
       diag_elements[i].set = false;
       diag_elements[i].diag_el.real = 0.0;
       diag_elements[i].diag_el.imag = 0.0;
    }
}



template<>  void
Fspai_Sub<COMPLEX>::Init_Diagonal
(   COMPLEX& A_kk )
{
    A_kk.real = 0.0;
    A_kk.imag = 0.0;
}



template<> void
Fspai_Sub<COMPLEX>::Reset_Sumvec
(   const int r_idx )
{
    sumvec[r_idx].real = 0.0;
    sumvec[r_idx].imag = 0.0;
}



template<> void
Fspai_Sub<COMPLEX>::Solve_HPD_System
(   COMPLEX*    red_mtx,
    COMPLEX*    red_col,
    int&        n,
    const char* uplo,
    int&        nrhs,
    int&        info)
{
    int lda = std::max(1, n);
    int ldb = std::max(1, n);

    // computing the solution to a complex system of linear equations
    // where A is an N-by-N Hermitian positive definite matrix
    zposv_( uplo,
            &n,
            &nrhs,
            red_mtx,
            &lda,
            red_col,
            &ldb,
            &info);
}



template<> COMPLEX*
Fspai_Sub<COMPLEX>::Conjugate_Complex_Vec
(   const COMPLEX*  vec1,
    const int       dim)
{
    COMPLEX* ccvec = new COMPLEX[dim];
    for(int i = 0; i < dim; i++)
    {
        ccvec[i].real = vec1[i].real;
        ccvec[i].imag = (-1.0)*vec1[i].imag;
    }
    return ccvec;
}



template<> COMPLEX
Fspai_Sub<COMPLEX>::Dot_Product
(   const COMPLEX*    vec1,
    const COMPLEX*    vec2,
    const int  dim)
{
    COMPLEX out;
    out.real = 0.0;
    out.imag = 0.0;

    // Note that in the complex case we
    // have to use the conjugate complex vector for
    // the product A(Jktilde,k)^H*y_k
    COMPLEX* ccvec = Conjugate_Complex_Vec(vec1,dim);

    for(int i = 0; i < dim; i++)
        out = out + (ccvec[i]*vec2[i]);
    delete [] ccvec;
    return out;
}



template<> COMPLEX
Fspai_Sub<COMPLEX>::Field_Sqrt
(   const COMPLEX number)
{
    COMPLEX out;
    double r   = sqrt(number.real*number.real+number.imag*number.imag);
    double phi = acos(number.real/r);
    out.real   = sqrt(r)*cos(phi);
    out.imag   = sqrt(r)*sin(phi);
    return out;
}
