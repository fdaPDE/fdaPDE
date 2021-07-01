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


template <class T_Field>
Compressed_Lines<T_Field>::Compressed_Lines
( const int nbr_cols_ )
{
    mtx_vals        = NULL;
    col_buf         = NULL;
    col_idcs        = NULL;
    col_idcs_buf    = NULL;
    len_cols        = NULL;
    nbr_cols        = nbr_cols_;

    mtx_vals = new T_Field*[nbr_cols];
    for (int i = 0; i < nbr_cols; i++)
        mtx_vals[i] = NULL;

    col_idcs = new int*[nbr_cols];
    for (int i = 0; i < nbr_cols; i++)
        col_idcs[i] = NULL;

    len_cols    = new int[nbr_cols];
    memset(len_cols, 0, nbr_cols * sizeof(int));
}


template <class T_Field>
Compressed_Lines<T_Field>::~Compressed_Lines( )
{
    if (col_buf)
        delete [] col_buf;
    else
        for (int i = 0; i < nbr_cols; i++)
            if (mtx_vals[i])
                delete [] mtx_vals[i];

    if (col_idcs_buf)
        delete [] col_idcs_buf;
    else
        for (int i = 0; i < nbr_cols; i++)
            if (col_idcs[i])
                delete [] col_idcs[i];

    if (mtx_vals)   delete [] mtx_vals;
    if (col_idcs)   delete [] col_idcs;
    if (len_cols)   delete [] len_cols;
}
