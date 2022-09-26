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

#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

// C/C++ includings
#include <math.h>

//file includings
#include "ENV_Handler.h"


///////////////////////////////////////////
///     \class Distribution
///     \brief This class distributes
///            the work chunks to all pe's.
///////////////////////////////////////////
class Distribution
{
    public:
        ////////////////////////////////////////////////
        ///     \brief Distributing the work chunks to
        ///            all pe's.
        ///
        ///     Every cluster node gets the number
        ///     of columns it will have to compute the
        ///     preconditioner solution for. If x is the
        ///     dimension of the input matrix and n the
        ///     the number of pe's than the first m = x-n
        ///     pe's will have to compute one column more
        ///     than the remaining n - m cluster nodes.
        ///     This way there is a maximum work chunk
        ///     difference of one between some nodes.
        ///
        ///     \param env_handler Environment handler interface
        ///     \param cols number of columns of matrix
        ///     \param my_nbr_cols number of columns this
        ///                pe will have to solve
        ///     \param split_pe index where the work chunk
        ///                     difference occur the first
        ///                     time within matrix
        ///     \param split_idx index where the work chunk
        ///                      is splitted.
        ///     \param start_idx index the work chunk
        ///                      starts within matrix
        ////////////////////////////////////////////////
        void    Basic_Distribution
                (   ENV_Handler&            env_handler,
                    int                     cols,
                    int&                    my_nbr_cols,
                    int&                    split_pe,
                    int&                    split_idx,
                    int&                    start_idx);
};

#endif
