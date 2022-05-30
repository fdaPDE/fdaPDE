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

#ifndef MMIO_H_
#define MMIO_H_

//C++ includings
#include <iostream>

// file includings
#include "Macros.h"

/// MatrixMarket typecode
typedef char MM_typecode[4];

////////////////////////////////////////////
///     \brief Matrix field data input
////////////////////////////////////////////
enum
{
    real,
    complex
};


///////////////////////////////////////////
///     \class MMio
///     \brief Matrix Market input/output
///     is responsible that the input
///     files, containing matrix data, are
///     all in correct matrix market
///     format.
///////////////////////////////////////////
class MMio
{
    public:

        /////////////////////////////////////////////////////////
        ///     \brief Skipping matrix file header
        ///
        ///     \param Path to the matrix file
        /////////////////////////////////////////////////////////
        void    Skip_Header
                ( FILE *f );

        /////////////////////////////////////////////////////////
        ///     \brief  Open matrix file and read the matrix
        ///             header input
        ///
        ///     \param f Path to the matrix file
        ///     \param type Type of the input matrix
        /////////////////////////////////////////////////////////
        void    Parse_Header
                ( FILE    *f,
                  int     &type );

        /////////////////////////////////////////////////////////
        ///     \brief  Determines the matrix field
        ///             ( data input )
        ///
        ///     \param matrix_file path to matrix file
        ///     \return Type of the matrix
        /////////////////////////////////////////////////////////
        int     Get_Matrix_Type
                (   char *matrix_file );

        /////////////////////////////////////////////////////////
        ///     \brief Reading the banner of a
        ///            matrix market file.
        ///
        ///     If anything is not correct -> abort
        ///     and return with error. If alright,
        ///     initialize the typecode by means of the
        ///     banner.
        ///
        ///     \param f File pointer of file, reading
        ///              the banner from
        ///     \param matcode The matcode to be set
        ///     \return return/error code for this
        ///             method
        /////////////////////////////////////////////////////////
        int     MM_Read_Banner
                (   FILE *f,
                    MM_typecode *matcode );

        /////////////////////////////////////////////////////////
        ///     \brief Reading the dimension and nnz
        ///            line of a matrix market file.
        ///
        ///     \param f File pointer of file reading
        ///              the line from
        ///     \param M The m-dimension of the matrix
        ///              to be set
        ///     \param N The n-dimension of the matrix
        ///              to be set
        ///     \param nz The number of nnz's of the
        ///               matrix to be set
        /////////////////////////////////////////////////////////
        void    MM_Read_Mtx_Crd_Size
                (   FILE *&f,
                    int &M,
                    int &N,
                    int &nz )                               const;

        /////////////////////////////////////////////////////////
        ///     \brief Reading the dimension and nnz
        ///            line of a matrix market pattern
        ///            file.
        ///
        ///     \param f File pointer of pattern file
        ///              reading the line from
        ///     \param M The m-dimension of the pattern
        ///              to be set
        ///     \param N The n-dimension of the pattern
        ///              to be set
        ///     \param nz The number of nnz's of the
        ///               pattern to be set
        /////////////////////////////////////////////////////////
        void    MM_Read_Pattern_Crd_Size
                (   FILE *&f,
                    int &M,
                    int &N,
                    int &nz )                               const;

        /////////////////////////////////////////////////////////
        ///     \brief  Checks whether pattern file can be opened
        ///             and whether it is a valid pattern file.
        ///
        ///     \param Path to the matrix file
        /////////////////////////////////////////////////////////
        void    Parse_Pattern_Header
                (   FILE *f );

        /////////////////////////////////////////////////////////
        ///     \brief  Clears matrix market typecode.
        /////////////////////////////////////////////////////////
        void    Clear_Typecode
                ( );

        /////////////////////////////////////////////////////////
        ///     \brief  Prints the matrix market typecode.
        /////////////////////////////////////////////////////////
        void    Print_Typecode
                ( )                                         const;

        /////////////////////////////////////////////////////////
        ///     \brief  Whether the typecode/object has symmetric
        ///             input or not.
        /////////////////////////////////////////////////////////
        bool    Is_Symmetric
                ( )                                         const;

        /// MatrixMarket typecode
        MM_typecode matcode;
};

#endif
