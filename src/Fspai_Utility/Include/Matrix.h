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


#ifndef MATRIX_H_
#define MATRIX_H_

//C++ includings
#include <stdexcept>

//file includings
#include "Pattern.h"
#include "Compressed_Lines.h"
#include "Matrix_Base.h"
#include "Timer.h"
#include "ENV_Handler.h"
#include "Pe_Exception.h"


//////////////////////////////////////////
///     \brief Simple structure representing
///            a complex number
//////////////////////////////////////////
struct COMPLEX
{
    /// real part of complex number
    double real;
    /// imaginary part of complex number
    double imag;
};


/// Operator overloading for "*" between
/// complex numbers
inline const COMPLEX
operator*(const COMPLEX& a, const COMPLEX& b)
{
    COMPLEX out;
    out.real = a.real * b.real - a.imag * b.imag;
    out.imag = a.real * b.imag + a.imag * b.real;
    return out;
}


/// Operator overloading for "*" between
/// complex number and scalar
inline const COMPLEX
operator*(const COMPLEX& a, const double& b)
{
    COMPLEX tmp;
    tmp.real = b;
    tmp.imag = 0.0;
    return a*tmp;
}


/// Operator overloading for "+" between
/// complex numbers
inline const COMPLEX
operator+(const COMPLEX& a, const COMPLEX& b)
{
    COMPLEX out;
    out.real = a.real + b.real;
    out.imag = a.imag + b.imag;
    return out;
}


/// Operator overloading for "-" between
/// complex numbers
inline const COMPLEX
operator-(const COMPLEX& a, const COMPLEX& b)
{
    COMPLEX out;
    out.real = a.real - b.real;
    out.imag = a.imag - b.imag;
    return out;
}


/// Operator overloading for "/" between
/// complex numbers
inline const COMPLEX
operator/(const COMPLEX& a, const COMPLEX& b)
{
    COMPLEX out;
    double div = b.real*b.real+b.imag*b.imag;
    out.real = (a.real*b.real+a.imag*b.imag) / div;
    out.imag = (a.imag*b.real-a.real*b.imag) / div;
    return out;
}


/// Operator overloading for "/" between
/// complex numbers
inline const COMPLEX
operator/(const COMPLEX& a, const double& b)
{
    COMPLEX tmp;
    tmp.real = b;
    tmp.imag = 0.0;
    return a / tmp;
}


/// Operator overloading for "/" between
/// complex numbers
inline const COMPLEX
operator/(const double& a, const COMPLEX& b)
{
    COMPLEX tmp;
    tmp.real = a;
    tmp.imag = 0.0;
    return tmp / b;
}


//////////////////////////////////////////
///     \class Matrix
///     \brief The Matrix data structure
///
///     Every pe has its own local chunk
///     of matrix data he is processing
///     on. If a pe needs data which his
///     matrix does not contain, he has
///     to request it from the remote pe.
///     The matrix consists mainly of the
///     compressed lines which store the
///     matrix relevant data.
//////////////////////////////////////////
template <class T_Field>
class Matrix : public Matrix_Base
{
    public:

        /// Empty Constructor
        Matrix<T_Field>() { };

        /////////////////////////////////////////////////////////
        ///     \brief  Constructor
        ///
        ///     \param env_handler Environment handler
        ///     \param nbr_cols Number of columns of matrix
        ///     \param nbr_rows Number of rows of matrix
        /////////////////////////////////////////////////////////
        Matrix<T_Field>
        (  const ENV_Handler& env_handler,
            const int         nbr_cols,
            const int         nbr_rows );

        /// Destructor
        ~Matrix<T_Field>();

        //Member variables

        /// Compressed column storage (CCS)
        /// This datastructure contains the
        /// matrix relevant data on this pe
        Compressed_Lines<T_Field> *c_lines;

        /// Remote buffer for transferring
        /// Values between the pe's
        T_Field                   *remote_col_buf;

        // Methods

        /////////////////////////////////////////////////////////
        ///     \brief  Printing all matrix data
        ///
        ///     This method prints all data in shortest
        ///     screen output.
        /////////////////////////////////////////////////////////
        void        Print_Matrix_Data()                    const;

        /////////////////////////////////////////////////////////
        ///     \brief  Printing the matrix values in
        ///             human readable format
        ///
        ///     This method prints only the matrix values
        ///     in human readable format.
        /////////////////////////////////////////////////////////
        void        Print_Matrix_Human_Readable
                    ( )                                     const;

        /////////////////////////////////////////////////////////
        ///     \brief  Writing one single line to file during
        ///             writing matrix to file.
        ///
        ///     \param c Column index of value.
        ///     \param r Row index of value.
        ///     \param f File to write data to.
        /////////////////////////////////////////////////////////
        void        Write_Line
                    (   int c,
                        int r,
                        FILE* f)                            const;

        /////////////////////////////////////////////////////////
        ///     \brief  Writing header of preconditioner matrix.
        ///
        ///     \param f File to write data to.
        /////////////////////////////////////////////////////////
        void        Write_Header
                    ( FILE* f )                             const;

        //===============================================================
        //============== Template methods - see Matrix.imp ==============
        //===============================================================

        /////////////////////////////////////////////////////////
        ///     \brief  Initializing the preconditioner matrix
        ///             with default values
        ///
        ///     \param  mtx The system matrix as access object to
        ///                 some members.
        ///     \param  env_handler Interface to environment
        ///                        specific methods.
        /////////////////////////////////////////////////////////
        void        Init_Preconditioner
                    (   Matrix<T_Field>* mtx,
                        ENV_Handler&     env_handler);

        /////////////////////////////////////////////////////////
        ///     \brief  Converts given matrix to a pattern.
        ///
        ///     \param P Pattern object to be filled.
        ///     \param env_handler Interface to environment
        ///            specific methods.
        /////////////////////////////////////////////////////////
        void        To_Pattern
                    (   Pattern*           P,
                        const ENV_Handler& env_handler );

        /////////////////////////////////////////////////////////
        ///     \brief  Printing the real matrix values
        ///             in matrix market format into a
        ///             file.
        ///
        ///     This method prints only the matrix values
        ///     into a file in matrix market format.
        ///     This format is a simple row-column-value
        ///     format (e.g.:) 2 1 9.003
        ///
        ///     \param env_handler Interface to environment
        ///                        specific methods.
        ///     \param file The file where the matrix
        ///                 should be printed to.
        /////////////////////////////////////////////////////////
        void        To_File
                    ( ENV_Handler& env_handler,
                      const char*  file )                  const;

    private:

        /////////////////////////////////////////////////////////
        ///     \brief  Counting the nz elements of this matrix.
        ///
        ///     \param env_handler Interface to environment
        ///                        specific methods.
        ///     \return Number of nz within this matrix.
        /////////////////////////////////////////////////////////
        int         Count_NNZ
                    ( ENV_Handler& env_handler )            const;
};

#include "Matrix.imp.h"

#endif
