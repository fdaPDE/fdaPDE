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

#ifndef ENV_MPI

#include "../Include/PCG.h"

//============================================================================
//============================================================================
//================ Template specifications for double matrices ===============
//============================================================================
//============================================================================

template<> void
PCG<double>::Print_Vec
(   const double* vec)
{
    std::cout << std::endl << "Vector: " << std::endl;
    for(int i = 0; i < mtx->my_nbr_cols; i++)
        std::cout << vec[i] << std::endl;
    std::cout << std::endl;
}



template<> void
PCG<double>::Init_Vec
(   double*&  vec,
    const int dim )
{
    memset(vec,0,dim*sizeof(double) );
}



template<> void
PCG<double>::RHS_Ones
(   )
{
    for(int i = 0; i < mtx->my_nbr_cols; i++)
        b[i] = 1.0;
}



template<> void
PCG<double>::RHS_Rand
(   )
{
    double r;
    // seed random generator once
    srand((unsigned)time(0));
    for(int i = 0; i < mtx->my_nbr_cols; i++)
    {
        r = static_cast<double> (rand())/static_cast<double>(RAND_MAX);
        b[i] = r;
    }
}



template<> void
PCG<double>::Write_Header
(   FILE* f  ) const
{
    fprintf(f, "matrix array real general\n");
}



template<> void
PCG<double>::Write_Line
(   const double&  val,
    FILE*          f ) const
{
    // Write data
    fprintf(f, "%.13e\n", val);
}



template<> double
PCG<double>::Inner_Prod
(   int     dim,
    double* x,
    double* y )
{
    double  result;
    int     one = 1;
    result = ddot_(&dim, x, &one, y, &one);
    return result;
}



template<> void
PCG<double>::Copy_Vector
(   int     dim,
    double* x,
    double* y)
{
    int one = 1;
    dcopy_(&dim, x, &one, y, &one);
}



template<> void
PCG<double>::Scale_Vector
(   int     dim,
    double  alpha,
    double  *x)
{
    int one = 1;
    dscal_(&dim, &alpha, x, &one);
}



template<> void
PCG<double>::Axpy
(   int     dim,
    double  alpha,
    double  *x,
    double  *y)
{
    int one = 1;
    daxpy_(&dim, &alpha, x, &one, y, &one);
}



template<> void
PCG<double>::MV_Prod
(   double*  loc_vec,
    double*& res_vec)
{
    int     row_len,
            *row_idcs_buf  = NULL;
    double  sum,
            *row_buf       = NULL;

    for(int row = 0; row < mtx->my_nbr_cols; row++)
    {
        comm_handler->Get_Mtx_Col(
                      mtx->my_start_idx+row,
                      row_len,
                      row_idcs_buf,
                      row_buf);
        sum = 0.0;
        for(int c = 0; c < row_len; c++)
            sum += (row_buf[c] * loc_vec[row_idcs_buf[c]]);
        res_vec[row] = sum;
    }
}



template<> void
PCG<double>::FSPAI_Apply
(   double*& res_vec,
    double*  loc_vec)
{
    int     col_len,
            rowidx,
            *col_idcs_buf = NULL;

    double  *col_buf      = NULL;

    // As this code is serial and matrices till approx.
    // size 10^6 may be solved in meaningull time memset
    // ok here.
    memset(res_vec,0,mtx->my_nbr_cols*sizeof(double));

    // Applying preconditioner L_M*loc_vec
    for(int col = 0; col < mtx->my_nbr_cols; col++)
    {
        // Get preconditioner column -> only local columns are possible.
        col_len      = precond->c_lines->len_cols[col];
        col_idcs_buf = precond->c_lines->col_idcs[col];
        col_buf      = precond->c_lines->mtx_vals[col];

        for(int c = 0; c < col_len; c++)
        {
            rowidx = col_idcs_buf[c];
            res_vec[rowidx] += col_buf[c] * loc_vec[col];
        }
    }
}



template<> void
PCG<double>::FSPAI_ApplyTrans
(   double*&      res_vec,
    double* loc_vec)
{
    int     col_len,
            *col_idcs_buf = NULL;
    double  sum,
            *col_buf       = NULL;

    // Applying preconditioner L_M^T*r
    for(int col = 0; col < mtx->my_nbr_cols; col++)
    {
        // Get preconditioner column. L^T only local columns are necessary.
        col_len      = precond->c_lines->len_cols[col];
        col_idcs_buf = precond->c_lines->col_idcs[col];
        col_buf      = precond->c_lines->mtx_vals[col];

        sum = 0.0;
        for(int c = 0; c < col_len; c++)
            sum += col_buf[c] * loc_vec[col_idcs_buf[c]];
        res_vec[col] = sum;
    }
}

//============================================================================
//============================================================================
//=============== Template specifications for COMPLEX matrices ===============
//============================================================================
//============================================================================

template<> void
PCG<COMPLEX>::Print_Vec(const COMPLEX* vec)
{
    std::cout << std::endl << "Vector: " << std::endl;
    for(int i = 0; i < mtx->my_nbr_cols; i++)
        std::cout << vec[i].real << " " << vec[i].imag << std::endl;;
    std::cout << std::endl;
}



template<> void
PCG<COMPLEX>::Init_Vec
(   COMPLEX*& vec,
    const int dim )
{
    for(int i= 0; i < dim; i++)
    {
        vec[i].imag = 0.0;
        vec[i].real = 0.0;
    }
}



template<> void
PCG<COMPLEX>::RHS_Ones()
{
    COMPLEX tmp;
    tmp.real = 1.0;
    tmp.imag = 0.0;
    for(int i = 0; i < mtx->my_nbr_cols; i++)
        b[i] = tmp;
}



template<> void
PCG<COMPLEX>::RHS_Rand
(   )
{
    double r1, r2;
    // seed random generator once
    srand((unsigned)time(0));
    for(int i = 0; i < mtx->my_nbr_cols; i++)
    {
        r1 = static_cast<double> (rand())/static_cast<double>(RAND_MAX);
        r2 = static_cast<double> (rand())/static_cast<double>(RAND_MAX);
        b[i].real = r1;
        b[i].imag = r2;
    }
}



template<> void
PCG<COMPLEX>::Write_Header
(   FILE* f  ) const
{
    fprintf(f, "matrix array complex general\n");
}



template<> void
PCG<COMPLEX>::Write_Line
(   const COMPLEX&  val,
    FILE*           f ) const
{
    // Write data
    fprintf(f, "%.13e %.13e\n", val.real, val.imag);
}



template<> COMPLEX*
PCG<COMPLEX>::Conjugate_Complex_Vec
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



template<> double
PCG<COMPLEX>::Inner_Prod
(   int dim,
    COMPLEX* x,
    COMPLEX* y )
{
    double  result;
    int     one = 1;
    COMPLEX* ccvec = Conjugate_Complex_Vec(x, dim);
    COMPLEX tmp = zdotu_(&dim, ccvec, &one, y, &one);
    result = tmp.real;
    delete [] ccvec;
    return result;
}



template<> void
PCG<COMPLEX>::Copy_Vector
(   int      dim,
    COMPLEX* x,
    COMPLEX* y)
{
    int one = 1;
    zcopy_(&dim, x, &one, y, &one);
}



template<> void
PCG<COMPLEX>::Scale_Vector
(   int      dim,
    double   alpha,
    COMPLEX  *x)
{
    int one = 1;
    COMPLEX calpha;
    calpha.real = alpha;
    calpha.imag = 0.0;
    zscal_(&dim, &calpha, x, &one);
}



template<> void
PCG<COMPLEX>::Axpy
(   int      dim,
    double   alpha,
    COMPLEX  *x,
    COMPLEX  *y)
{
    int one = 1;
    COMPLEX calpha;
    calpha.real = alpha;
    calpha.imag = 0.0;
    zaxpy_(&dim, &calpha, x, &one, y, &one);
}



template<> void
PCG<COMPLEX>::MV_Prod
(   COMPLEX*  loc_vec,
    COMPLEX*& res_vec)
{
    int     row_len,
            *row_idcs_buf = NULL;
    COMPLEX sum,
            tmp,
            *row_buf       = NULL;

    // compute own part of MV-Prod.
    for(int row = 0; row < mtx->my_nbr_cols; row++)
    {
        comm_handler->Get_Mtx_Col(
                      mtx->my_start_idx+row,
                      row_len,
                      row_idcs_buf,
                      row_buf);
        sum.real = 0.0;
        sum.imag = 0.0;
        for(int c = 0; c < row_len; c++)
        {
            tmp = row_buf[c];
            // as A is conjugate complex back-conjugation is necessary
            tmp.imag = -1.0*tmp.imag;
            sum = sum + (tmp * loc_vec[row_idcs_buf[c]]);
        }
        res_vec[row] = sum;
    }
}



template<> void
PCG<COMPLEX>::FSPAI_Apply
(   COMPLEX*& res_vec,
     COMPLEX* loc_vec)
{
    int     col_len,
            rowidx,
            *col_idcs_buf = NULL;

    COMPLEX  *col_buf      = NULL;

    // As this code is serial and matrices till approx.
    // size 10^6 may be solved in meaningull time memset
    // ok here.
    Init_Vec(res_vec, mtx->my_nbr_cols);

    // Applying preconditioner L_M*loc_vec
    for(int col = 0; col < mtx->my_nbr_cols; col++)
    {
        // Get preconditioner column -> only local columns are possible.
        col_len      = precond->c_lines->len_cols[col];
        col_idcs_buf = precond->c_lines->col_idcs[col];
        col_buf      = precond->c_lines->mtx_vals[col];

        for(int c = 0; c < col_len; c++)
        {
            rowidx = col_idcs_buf[c];
            res_vec[rowidx] = res_vec[rowidx] + (col_buf[c] * loc_vec[col]);
        }
    }
}



template<> void
PCG<COMPLEX>::FSPAI_ApplyTrans
(   COMPLEX*& res_vec,
    COMPLEX*  loc_vec)
{
    int     col_len,
            *col_idcs_buf = NULL;
    COMPLEX sum,
            tmp,
            *col_buf      = NULL;

    // Applying preconditioner L_M^H*r
    for(int col = 0; col < mtx->my_nbr_cols; col++)
    {
        // Get preconditioner column. L^T only local columns are necessary.
        col_len      = precond->c_lines->len_cols[col];
        col_idcs_buf = precond->c_lines->col_idcs[col];
        col_buf      = precond->c_lines->mtx_vals[col];

        sum.real = 0.0;
        sum.imag = 0.0;
        for(int c = 0; c < col_len; c++)
        {
            tmp = col_buf[c];
            tmp.imag = -1.0*tmp.imag;   // due to L^H
            sum = sum + (tmp * loc_vec[col_idcs_buf[c]]);
        }
        res_vec[col] = sum;
    }
}

#endif
