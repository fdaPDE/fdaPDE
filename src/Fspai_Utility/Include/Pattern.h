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


#ifndef PATTERN_H_
#define PATTERN_H_

//file includings
#include "Index_Set.h"
#include "Distribution.h"
#include "Timer.h"
#include "Pe_Exception.h"

//C/C++ includings
#include <stdio.h>
#include <stdlib.h>

//////////////////////////////////////////
///     \brief Row Column structure
///            representing a row column
///            format for pattern
//////////////////////////////////////////
struct RC
{
    int i;
    int j;
};


////////////////////////////////////////////
///     \brief This comparator is necessary
///            for sorting the RC arrays in
///            ascending order
////////////////////////////////////////////
struct Idx_Comparator
{
    template <class T>
    bool operator()(const T& a, const T& b)
    {
        if ( (a.j <  b.j) ||
            ((a.j == b.j) && (a.i < b.i)) ||
            ((a.j == b.j) && (a.i == b.i)) )
            return  true;

        return false;
    }
};


//////////////////////////////////////////
///     \class Pattern
///     \brief The Pattern data structure
///
///     Every pe has its own local chunk
///     of pattern data he is processing
///     on. If a pe needs data which his
///     matrix does not contain, he has
///     to request it from a remote pe.
///     The pattern consists mainly of a
///     2D integer array containing the
///     row column indices within the
///     pattern matrix
//////////////////////////////////////////
class Pattern
{
    public:
        /// Empty Constructor
        Pattern();

        /// Destructor
        ~Pattern();

        //Member variables

        /// 2D array containing index
        /// sets of each pattern column
        Index_Set   **j_sets;

        /// The number of columns/index sets
        /// this pattern contains
        int         len;

        /// The start index of the chunk pattern
        /// of this pe within the whole pattern
        int         my_start_idx;

        /// The maximum number of indices per pattern
        /// column
        int         max_nnz;

         /// The number of columns to solve by this pe
        int         my_nbr_cols;

        /// The next column to process
        int         next_col;

        /// my_start_index of each pe
        int         *start_indices;

        /// processor assignment for every column
        int         *pe;

        /// my_nbr_cols of each pe
        int         *all_nbr_cols;

        /// The pe's Id within MPI environment
        int         my_id;

        /// The number of pe's within MPI environment
        int         num_procs;

        //Methods

        //////////////////////////////////////////////////////
        ///     \brief  Printing all pattern data
        //////////////////////////////////////////////////////
        void    Print_Pattern_Data();

        //////////////////////////////////////////////////////
        ///     \brief  Building pattern out of pattern file
        ///
        ///     Building pattern data structure out of pattern
        ///     file.
        ///
        ///     \param env_handler Interface to environment
        ///                        handler specific methods.
        ///     \param file Path to file to build pattern from.
        ///     \param dim  Number of columns of input matrix.
        //////////////////////////////////////////////////////
        void    Arbitrary_Pattern
                (    ENV_Handler& env_handler,
                     char         *file,
                     const int    dim);

        //////////////////////////////////////////////////////
        ///     \brief  Building pattern if user wants
        ///             unit matrix as input pattern
        ///
        ///     \param env_handler Interface to environment
        ///                        handler specific methods.
        ///     \param mtx_dim Number of columns of input system.
        ///     \param mtx_my_nbr_cols Number of columns pattern
        ///                            columns on this pe.
        ///     \param mtx_start_idx Start index of this pe's
        ///                      local work chunk within whole
        ///                      pattern.
        //////////////////////////////////////////////////////
        void    Diagonal_Pattern
                (    ENV_Handler& env_handler,
                     const int    mtx_dim,
                     const int    mtx_my_nbr_cols,
                     const int    mtx_my_start_idx );
    private:

        //////////////////////////////////////////////////////
        ///     \brief  Parsing file and filling arrays with
        ///             row/column elements specific to this
        ///             pe.
        ///
        ///     \param f Path to the pattern file matrix
        //////////////////////////////////////////////////////
        void    Read_Data( FILE *f );

        //////////////////////////////////////////////////////
        ///     \brief  Mapping read out pattern data to
        ///             Pattern data structure
        ///
        ///     \param env_handler Interface to environment
        ///                        handler specific methods.
        //////////////////////////////////////////////////////
        void    Data_To_Pattern( ENV_Handler& env_handler );

        //////////////////////////////////////////////////////
        ///     \brief  Counting nnz in each column
        ///
        ///     \param in array containing the RC elements
        ///               read out from the pattern file
        ///     \param size size of the in array
        ///     \param out array containing the number of
        ///                nnzs of each column within the
        ///                pattern
        //////////////////////////////////////////////////////
        void    Count_NNZ_Cols
                (   RC const    *in,
                    size_t      size,
                    int         *out);

        //////////////////////////////////////////////////////
        ///     \brief  Counting nnz in each row
        ///
        ///     \param in array containing the RC elements
        ///               read out from the pattern file
        ///     \param size size of the in array
        ///     \param out array containing the number of
        ///                nnzs of each row within the
        ///                pattern
        //////////////////////////////////////////////////////
        void    Count_NNZ_Rows
                (   RC const    *in,
                    size_t      size,
                    int         *out);

        //////////////////////////////////////////////////////
        ///     \brief  Parsing file and filling the cols
        ///             structure with matrix elements.
        ///
        ///     \param f Path to the pattern file matrix
        //////////////////////////////////////////////////////
        void    Data_To_Memory( FILE *f );

        //////////////////////////////////////////////////////
        ///     \brief  Counting number of nonzeros of pattern
        ///             for memory allocation size.
        ///
        ///     \param f Path to the pattern file matrix
        //////////////////////////////////////////////////////
        void    Count_NNZ( FILE* f );

        // Member
        RC      *cols;

        int     file_cols,
                file_rows,
                file_nnz,
                nnz_cols,
                split_idx,
                split_pe,
                start_idx;

        bool    symmetric;
};

#endif
