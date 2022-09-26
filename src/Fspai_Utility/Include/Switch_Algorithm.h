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


#ifndef SWITCH_ALGORITHM_H_
#define SWITCH_ALGORITHM_H_

// file includings
#include "Macros.h"
#include "Fspai_Base.h"
#include "Fspai_Unrestrained.h"
#include "Comm_Handler.h"

// C/C++ includings
#include <iostream>


////////////////////////////////////////////
///     \brief Which algorithm to switch to
///            by given algorithm level
////////////////////////////////////////////
enum
{
    unrestrained
};


//////////////////////////////////////////////
///     \class Switch_Algorithm
///     \brief This class invokes the
///            correct FSPAI algorithm by
///            means of the user specified
///            algorithm level.
///
///     This class is responsible to invoke
///     the correct derived Fspai algorihm
///     class which derives from Fspai_Base.
///
//////////////////////////////////////////////
template <class T_Field>
class Switch_Algorithm
{
    public:
        /////////////////////////////////////////////////////////
        ///     \brief  Get the Fspai algorithm by means of
        ///             the algorithm level.
        ///
        ///     Due to the algorithm level the user invoked,
        ///     the specific Fspai algorithm will be instantiated.
        ///
        ///     \param alg_level Algorithm level
        ///     \param mtx_ Local matrix chunk on this pe
        ///     \param precond_ Local preconditioner chunk on this
        ///                     pe
        ///     \param P_ Local pattern chunk on this pe
        ///     \param env_handler Interface to environment
        ///                        specific methods.
        ///     \param hash_param Size of hash table to be used.
        ///     \param epsilon_param Epsilon as tolerance for
        ///            residual
        ///     \param updates_param Number of update steps
        ///     \param max_idcs_param Max. number of indices to be
        ///            augmented per update step
        ///     \param use_mean_param Whether to use the mean value
        ///            bound when augmenting indices per step
        ///     \return From Fspai_Base derived class which
        ///             implements the virtual methods.
        /////////////////////////////////////////////////////////
        Fspai_Base<T_Field>*    Get_Algorithm
                                (   const int         alg_level,
                                    Matrix<T_Field>*  mtx,
                                    Matrix<T_Field>*& precond,
                                    Pattern*          P,
                                    ENV_Handler&      env_handler,
                                    const int         hash_param,
                                    const double      epsilon_param,
                                    const int         updates_param,
                                    const int         max_idcs_param,
                                    const bool        use_mean_param);
};

#include "Switch_Algorithm.imp.h"

#endif
