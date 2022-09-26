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

#ifndef COMPRESSED_LINES_H_
#define COMPRESSED_LINES_H_

//C++ includings
#include <iostream>


/////////////////////////////////////////////////////////
///     \class Compressed_Lines
///     \brief Implementing the compressed
///            column storage (CCS)
///
///     Each matrix has its own compressed lines
///     structure which holds all matrix data.
///     Data is stored in compressed column storage (CCS)
///     with additional information. This is useful for
///     the communication methods.
/////////////////////////////////////////////////////////
template <class T_Field>
class Compressed_Lines
{
    public:

        /// Empty Constructor
        Compressed_Lines<T_Field>() { };

        /// Constructor
        Compressed_Lines<T_Field>(const int nbr_cols);

        /// Destructor
        ~Compressed_Lines<T_Field>();

        //Member variables

        /// 2D array holds column-arrays named col_buf
        T_Field** mtx_vals;

        /// This is the buffer array which holds the values of the specific
        /// column of the matrix, it is accessed via the pointers of mtx_vals.
        T_Field* col_buf;

        /// 2D array holds column specific pointers to col_idcs_buf
        int **col_idcs;

        /// Buffer array which holds the row indices of the
        /// specific nnz elementes of each column of the matrix.
        /// It is accessed via the pointers of col_idcs.
        int *col_idcs_buf;

        /// Number of nonzeros in each column
        int *len_cols;

    private:

        /// memory size of buffers
        int nbr_cols;
};

#include "Compressed_Lines.imp.h"

#endif
