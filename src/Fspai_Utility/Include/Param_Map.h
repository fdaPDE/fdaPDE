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

#ifndef PARAM_MAP_H_
#define PARAM_MAP_H_

// C++ includings
#include <map>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


////////////////////////////////////////
///     \class Param_Map
///     \brief Dictionary holding all
///            shell parameters
////////////////////////////////////////
class Param_Map : public std::map<std::string, std::string>
{
    public:

        /////////////////////////////////////////////////
        ///     \brief  Gets a input string and converts
        ///             it to the right type argument.
        ///
        ///     \param par The parameter to be converted
        ///     \param arg The converted parameter
        ///     \return Whether the conversion was
        ///             successful
        /////////////////////////////////////////////////
        template<class T_Type>
        int     Get_Arg_Typed
                (   std::string par,
                    T_Type&     arg ) ;

        /////////////////////////////////////////////////
        ///     \brief  Converts input shell parameter to
        ///             correct output type
        ///
        ///     \param in Input parameter
        ///     \param out Converted input parameter
        ///     \return Whether the conversion was successful
        /////////////////////////////////////////////////
        template<class T_Type>
        int     Converter
                (   std::string in,
                    T_Type&     out ) ;
};

#include "Param_Map.imp.h"

#endif

