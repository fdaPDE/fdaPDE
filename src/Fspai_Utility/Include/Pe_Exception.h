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

#ifndef PE_EXCEPTION_H_
#define PE_EXCEPTION_H_

// C/C++ includings
#include <stdexcept>

///////////////////////////////////////////
///     \class Pe_Exception
///     \brief Own Exception class
///
///     Providing own exception class for
///     catching processor specific errors
///     in parallel environment. E.g. if
///     a pe with ID > 0 throws an exception
///     the error message cannot be streamed
///     to shell within the environment-global
///     runtime_error catch block, because this
///     is used for error messages which occur
///     on every pe. Therefore a pe specific
///     catch block for messages which are pe
///     specific is added.
///////////////////////////////////////////
class Pe_Exception : public std::runtime_error
{
    public:
        Pe_Exception(const std::string& err) : std::runtime_error(err) { };
};

#endif /* PE_EXCEPTION_H_ */
