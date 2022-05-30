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

#ifndef READ_MATRIX_
#define READ_MATRIX_

//file includings
#include "Pattern.h"
#include "Matrix.h"
#include "MMio.h"
#include "Distribution.h"
#include "ENV_Handler.h"
#include "Comm_Handler.h"
#include "Pe_Exception.h"

//C/C++ includings
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <vector>


////////////////////////////////////////////
///     \brief Simple row column value struct
////////////////////////////////////////////
struct RCV
{
    int i;          //row
    int j;          //column
    double val;     //value
};


////////////////////////////////////////////
///     \brief Row Column Complex struct
////////////////////////////////////////////
struct RCC
{
    int i;
    int j;
    COMPLEX c;
};


////////////////////////////////////////////
///     \brief This comparator is necessary
///            for sorting the RCC and RCV
///            column arrays in acsending
///            order
////////////////////////////////////////////
struct Col_Comparator
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


////////////////////////////////////////////
///     \brief This comparator is necessary
///            for sorting the RCC and RCV
///            row arrays in acsending
///            order
////////////////////////////////////////////
struct Row_Comparator
{
    template <class T>
    bool operator()(const T& a, const T& b)
    {
        if ( (a.i <  b.i) ||
            ((a.i == b.i) && (a.j < b.j)) ||
            ((a.i == b.i) && (a.j == b.j)) )
            return  true;

        return false;
    }
};


//////////////////////////////////////////
///     \class Matrix_Reader
///
///     \brief Provides the mapping from
///            a matrix file in matrix market
///            format to the matrix
///            data structure.
//////////////////////////////////////////
template <class T_Field, class T_Format>
class Matrix_Reader
{
    public:

        /// Constructor
        Matrix_Reader<T_Field, T_Format>();

        /// Destructor
        ~Matrix_Reader<T_Field, T_Format>();

        //===============================================================
        //=================== template based methods ====================
        //===============================================================

        ///////////////////////////////////////////////////////////
        ///     \brief  Parses the matrix file and fills the
        ///             matrix structure.
        ///
        ///     Step list:
        ///     * check header
        ///     * determine chunk size with distribution in MPI environment
        ///     * parse matrix file and read data
        ///     * fill data into matrix data structure
        ///
        ///     \param mtx The matrix data structure to be filled
        ///     \param matrix_file path to matrix file
        ///     \param mmio Object with specific typecode
        ///     \param env_handler Interface to environment
        ///                        specific methods.
        ///////////////////////////////////////////////////////////
        void    To_Memory
                (   Matrix<T_Field> *&  mtx,
                    char*               matrix_file,
                    const MMio          mmio,
                    ENV_Handler&        env_handler);

        ///////////////////////////////////////////////////////////
        ///     \brief  Counts the number of nonzeros of the system
        ///             and the probing matrix this pe will have.
        ///
        ///     Before allocating memory the space is computed by
        ///     this method. In case of symmetric matrices given in
        ///     matrix market the necessary space will be double size.
        ///
        ///     \param f The file to read the nnz from
        ///     \param my_nnz_col The current number of neccessary
        ///                       columns.
        ///////////////////////////////////////////////////////////
        void    Count_NNZ
                (   FILE    *f,
                    int&    my_nnz_col);

    private:

        //================================================================
        //===================== real matrix methods ======================
        //================================================================

        ///////////////////////////////////////////////////////////
        ///     \brief  Filling one double element to cols array.
        ///
        ///     \param row row index
        ///     \param col column index
        ///     \param val value of matrix element
        ///     \param inc position in memory array
        ///////////////////////////////////////////////////////////
        void    Fill_Field_Date
                (   const int    col,
                    const int    row,
                    const double val,
                    int&         inc );

        //================================================================
        //==================== complex matrix methods ====================
        //================================================================

        ///////////////////////////////////////////////////////////
        ///     \brief  Filling one complex element to cols array.
        ///
        ///     \param row row index
        ///     \param col column index
        ///     \param real real part of complex number
        ///     \param imag imaginary part of complex number
        ///     \param cc_val hermitian value if system is hermitian
        ///     \param inc position in memory array
        ///////////////////////////////////////////////////////////
        void    Fill_Field_Date
                (   const int    row,
                    const int    col,
                    const double real,
                    const double imag,
                    const double cc_val,
                    int&         inc );

        //===============================================================
        //=================== template based methods ====================
        //===============================================================

        ///////////////////////////////////////////////////////////
        ///     \brief  Filling lines structure with columns
        ///             specific matrix values.
        ///
        ///     \param len_col Length of current column
        ///     \param lines Compressed Lines structure to be filled
        ///     \param col The current column
        ///     \param pos Position in data array
        ///////////////////////////////////////////////////////////
        void    Fill_Lines
                (   const int                   len_col,
                    Compressed_Lines<T_Field>*& lines,
                    const int                   col,
                    int&                        pos );

        ///////////////////////////////////////////////////////////
        ///     \brief  Parses the matrix file and fills data to
        ///             memory.
        ///
        ///      In the external Matrix Market format indices are
        ///      1-based.
        ///      All internal representations are 0-based.
        ///      This procedure fills the cols array.
        ///      The values are stored in i, j and val. It is a
        ///      structure. -> Array of elements - each a structure
        ///      containing the nnz-elements. Each nnz-element is
        ///      i, j, val.
        ///
        ///     \param f file pointer in matrix file
        ///////////////////////////////////////////////////////////
        void     Data_To_Memory
                 (  FILE    *f);

        ///////////////////////////////////////////////////////////
        ///     \brief  Parses the matrix file and invokes method to
        ///             fill temporary array with matrix data.
        ///
        ///      In the external Matrix Market format indices are
        ///      1-based.
        ///      The file is passed twice. First to get the final
        ///      size the arrays will have, and second to fill
        ///      in the values.
        ///
        ///     \param f file pointer in matrix file
        ///////////////////////////////////////////////////////////
        void     Read_Data
                 (   FILE *f );

        ///////////////////////////////////////////////////////////
        ///     \brief  Filling matrix datastructure with data
        ///
        ///     For staying memory efficient, no row indices arrays
        ///     are filled, if the input matrix is symmetric.
        ///
        ///     \param mtx The matrix structure to be filled
        ///     \param env_handler Interface to environment
        ///                        specific methods.
        ///////////////////////////////////////////////////////////
       void     Data_To_Matrix
                (   Matrix<T_Field>*&  mtx,
                    const ENV_Handler& env_handler);

        ///////////////////////////////////////////////////////////
        ///     \brief  Counting nnz in each column
        ///
        ///     \param in array containing the RCV elements
        ///               read out from the matrix file
        ///     \param size size of the in array
        ///     \param out array containing the number of
        ///                nnzs of each column within the
        ///                pattern
        ///////////////////////////////////////////////////////////
        void    Count_NNZ_Cols
                (   T_Format const* in,
                    size_t          size,
                    int*            out);

        ///////////////////////////////////////////////////////////
        ///     \brief  Counting nnz in each row
        ///
        ///     \param in array containing the RCC elements
        ///               read out from the matrix file
        ///     \param size size of the in array
        ///     \param out array containing the number of
        ///                nnzs of each row within the
        ///                pattern
        ///////////////////////////////////////////////////////////
        void     Count_NNZ_Rows
                 (  T_Format const* in,
                    size_t          size,
                    int*            out);

        // Member variables

        int         file_cols,
                    file_rows,
                    file_nnz,
                    start_idx,
                    nbr_rows,
                    nbr_cols,
                    nnz,
                    nnz_cols,
                    split_idx,
                    split_pe,
                    my_nbr_cols,
                    my_nbr_rows;

        bool        hermitian;

        /// RCV/RCC array of matrix data
        T_Format*   cols;

        /// holding column flags whether matrix is set on diagonal
        bool        *diagonal_elements;
};

#include "Matrix_Reader.imp.h"

#endif
