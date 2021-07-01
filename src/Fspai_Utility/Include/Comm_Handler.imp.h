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


template <class T_Field> int
Comm_Handler<T_Field>::LB_Get_Precond_Col()
{
    int     pe,
            index,
            col;

    // Do I still have local columns to process?
    if (P->next_col < P->my_nbr_cols)
    {
        pe = P->my_id;
        index = P->next_col;
        col = index + P->start_indices[pe];
        P->next_col++;

        return col; // success
    }
    return -1;
}



template <class T_Field> void
Comm_Handler<T_Field>::Get_Mtx_Col
(   int         col,
    int&        col_len,
    int*&       col_idcs_buf,
    T_Field*&   col_buf)
{
    col_len      = mtx->c_lines->len_cols[col];
    col_idcs_buf = mtx->c_lines->col_idcs[col];
    col_buf      = mtx->c_lines->mtx_vals[col];
}



template <class T_Field> void
Comm_Handler<T_Field>::Get_P_Col
(   const int   col,
    Index_Set*& J   )
{
    int pe,idx;
    Index_Set       *J_pattern = NULL;

    pe          = P->pe[col];
    idx         = col - P->start_indices[pe];
    J_pattern   = P->j_sets[idx];

    // Copy the vlues into the big working Index_Set J
    memcpy(J->idcs,
           J_pattern->idcs,
           J_pattern->len * sizeof(int));

    J->len = J_pattern->len;
}



template <class T_Field> void
Comm_Handler<T_Field>::Insert_Precond_Solution
(   const T_Field*      Lk,
    const int           col,
    const Index_Set*    Jk_tilde,
    const Index_Set*    Jk)
{
    int         pe,
                idx,
                nnz,
                *col_idcs_buf = NULL;
    T_Field     *col_buf = NULL;

    pe                              = mtx->pe[col];
    idx                             = col - mtx->start_indices[pe];
    nnz                             = Jk->len;
    precond->c_lines->len_cols[idx] = nnz;
    col_buf                         = new T_Field[nnz];
    col_idcs_buf                    = new int[nnz];
    precond->c_lines->mtx_vals[idx] = col_buf;
    precond->c_lines->col_idcs[idx] = col_idcs_buf;

    for ( int i = 0; i < nnz; i++)
    {
        col_idcs_buf[i]             = Jk->idcs[i];
        col_buf[i]                  = Lk[i];
    }
}

