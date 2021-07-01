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

#ifndef MACROS_H_
#define MACROS_H_


// Macros for colorful shell output

#define COLOR_NORMAL     "\033[0m"
#define COLOR_RED        "\033[0;31m"
#define COLOR_GREEN      "\033[0;32m"
#define COLOR_YELLOW     "\033[0;33m"
#define COLOR_BLUE       "\033[0;34m"
#define COLOR_PURPLE     "\033[0;35m"
#define COLOR_CYAN       "\033[0;36m"
#define COLOR_GRAY       "\033[0;37m"

#define COLOR_B_RED      "\033[1;31m"
#define COLOR_B_GREEN    "\033[1;32m"
#define COLOR_B_YELLOW   "\033[1;33m"
#define COLOR_B_BLUE     "\033[1;34m"
#define COLOR_B_PURPLE   "\033[1;35m"
#define COLOR_B_CYAN     "\033[1;36m"
#define COLOR_B_WHITE    "\033[1;37m"

#define BGCOLOR_RED        "\033[0;41m"
#define BGCOLOR_GREEN      "\033[0;42m"
#define BGCOLOR_YELLOW     "\033[0;43m"
#define BGCOLOR_BLUE       "\033[0;44m"
#define BGCOLOR_PURPLE     "\033[0;45m"
#define BGCOLOR_CYAN       "\033[0;46m"
#define BGCOLOR_GRAY       "\033[0;47m"

#define BGCOLOR_B_RED      "\033[1;41m"
#define BGCOLOR_B_GREEN    "\033[1;42m"
#define BGCOLOR_B_YELLOW   "\033[1;43m"
#define BGCOLOR_B_BLUE     "\033[1;44m"
#define BGCOLOR_B_PURPLE   "\033[1;45m"
#define BGCOLOR_B_CYAN     "\033[1;46m"
#define BGCOLOR_B_WHITE    "\033[1;47m"

// Return value macros

#define EXIT_FAILURE 1


// Matrix Market macros

#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64


// MM_typecode query fucntions

#define MM_Is_Matrix(typecode)      ((typecode)[0]=='M')

#define MM_Is_Sparse(typecode)      ((typecode)[1]=='C')
#define MM_Is_Coordinate(typecode)  ((typecode)[1]=='C')
#define MM_Is_Dense(typecode)       ((typecode)[1]=='A')
#define MM_Is_Array(typecode)       ((typecode)[1]=='A')

#define MM_Is_Complex(typecode)     ((typecode)[2]=='C')
#define MM_Is_Real(typecode)        ((typecode)[2]=='R')
#define MM_Is_Pattern(typecode)     ((typecode)[2]=='P')
#define MM_Is_Integer(typecode)     ((typecode)[2]=='I')

#define MM_Is_Symmetric(typecode)   ((typecode)[3]=='S')
#define MM_Is_General(typecode)     ((typecode)[3]=='G')
#define MM_Is_Skew(typecode)        ((typecode)[3]=='K')
#define MM_Is_Hermitian(typecode)   ((typecode)[3]=='H')


// MM_typecode modify fucntions

#define MM_Set_Matrix(typecode)     ((*typecode)[0]='M')
#define MM_Set_Coordinate(typecode) ((*typecode)[1]='C')
#define MM_Set_Array(typecode)      ((*typecode)[1]='A')
#define MM_Set_Dense(typecode)      MM_Set_Array(typecode)
#define MM_Set_Sparse(typecode)     MM_Set_Coordinate(typecode)

#define MM_Set_Complex(typecode)    ((*typecode)[2]='C')
#define MM_Set_Real(typecode)       ((*typecode)[2]='R')
#define MM_Set_Pattern(typecode)    ((*typecode)[2]='P')
#define MM_Set_Integer(typecode)    ((*typecode)[2]='I')

#define MM_Set_Symmetric(typecode)  ((*typecode)[3]='S')
#define MM_Set_General(typecode)    ((*typecode)[3]='G')
#define MM_Set_Skew(typecode)       ((*typecode)[3]='K')
#define MM_Set_Hermitian(typecode)  ((*typecode)[3]='H')

#define MM_Clear_Typecode(typecode) ((*typecode)[0]=(*typecode)[1]= \
                                    (*typecode)[2]=' ',(*typecode)[3]=' ')

#define MM_Initialize_Typecode(typecode) MM_Clear_Typecode(typecode)


// Matrix Market error codes

#define MM_COULD_NOT_READ_FILE  11
#define MM_PREMATURE_EOF        12
#define MM_NOT_MTX              13
#define MM_NO_HEADER            14
#define MM_UNSUPPORTED_TYPE     15
#define MM_LINE_TOO_LONG        16
#define MM_COULD_NOT_WRITE_FILE 17


/******************** Matrix Market internal definitions ****************

   MM_matrix_typecode: 4-character sequence

                    ojbect      sparse/     data        storage
                                dense       type        scheme

   string position:  [0]        [1]         [2]         [3]

   Matrix typecode:  M(atrix)  C(oord)      R(eal)      G(eneral)
                                A(array)    C(omplex)   H(ermitian)
                                            P(attern)   S(ymmetric)
                                            I(nteger)   K(kew)

 ***********************************************************************/

#define MM_MTX_STR      "MATRIX"
#define MM_ARRAY_STR    "ARRAY"
#define MM_DENSE_STR    "ARRAY"
#define MM_COORDINATE_STR "COORDINATE"
#define MM_SPARSE_STR   "COORDINATE"
#define MM_COMPLEX_STR  "COMPLEX"
#define MM_REAL_STR     "REAL"
#define MM_INT_STR      "INTEGER"
#define MM_GENERAL_STR  "GENERAL"
#define MM_SYMM_STR     "SYMMETRIC"
#define MM_HERM_STR     "HERMITIAN"
#define MM_SKEW_STR     "SKEW-SYMMETRIC"
#define MM_PATTERN_STR  "PATTERN"



// Communication IDs for parallel Environment

// Messages with these tags are handled by Comm_Handler class

/// communication id for sending
/// request to remote pe
#define GET_RC_TAG          1

/// communication id for requesting
/// a column from a remote pe
#define REQUEST_MCOL_TAG    2

/// communication id for putting
/// remotly computed preconditioner
/// column into own cache
#define PUT_MCOL_TAG        3

/// communication id that this pe
/// has finished
#define IM_DONE_TAG         4

/// communication id that every
/// pe has finished
#define DONE_SIGNAL_TAG     5

/// communication id if error
/// occured within MPI communication
#define EXIT_TAG            6

/// Communication id for reqeusting
/// a column of the start pattern
#define GET_P_TAG           7

// Messages with these tags aren't

/// communication id for requesting
/// column indices
#define SEND_COLS_TAG       8

/// communication id for requesting
/// column values
#define SEND_VALS_TAG       9

/// communication id for
#define SEND_MCOL_TAG       10

/// communication id for sending
/// column indices of remotly computed
/// preconditioner column
#define PUT_MCOL_INDS_TAG    11

/// communication id for sending
/// column values of remotly computed
/// preconditioner column
#define PUT_MCOL_VALS_TAG    12

/// communication id for sending
/// the start pattern indices
#define SEND_PCOL_TAG        13

/// communication id for sending
/// the start pattern length
#define SEND_PLEN_TAG        14

#endif /* MACROS_H_ */
