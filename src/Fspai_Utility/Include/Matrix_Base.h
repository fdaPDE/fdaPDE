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

#ifndef MATRIX_BASE_H_
#define MATRIX_BASE_H_


//////////////////////////////////////
///     \class Matrix_Base
///     \brief Base class of Matrix
//////////////////////////////////////
class Matrix_Base
{
    public:

        /// The pe's Id within MPI environment
        int     my_id;

        /// The number of pe's within MPI environment
        int     num_procs;

        /// The number of columns to solve by this pe
        int     my_nbr_cols;

        /// The start index of the first column to solve
        /// within whole input matrix
        int     my_start_idx;

        /// Number of columns of whole system matrix
        int     n;

        /// Number of rows of whole system matrix
        int     m;

        /// my_nbr_cols of each pe
        int     *all_nbr_cols;

        /// my_start_index of each pe
        int     *start_indices;

        /// The len of all columns on every pe -
        /// will be filled with Allgatherv
        int     *len_all_cols;

        /// The maximum number of nnz of all columns
        /// on this pe
        int     max_nnz;

        /// nnz on this pe
        int     my_nnz;

        /// processor assignment for every row
        /// and column
        int     *pe;

        /// Remote buffer for transferring the column
        /// indices between the pe's
        int     *remote_col_idcs_buf;

        /// The next column to process
        int     next_col;
};

#endif
