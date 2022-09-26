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

#ifndef FSPAI_UNRESTRAINED_H
#define FSPAI_UNRESTRAINED_H

// file includings
#include "Fspai_Base.h"

//C/C++ includings
#include <iostream>
#include <ctime>


/////////////////////////////////////////////////////////////
///     \class Fspai_Unrestrained
///     \brief  Implementing the standard  FSPAI algorithm.
///
///     This is the FSPAI algorithm implementation as
///     presented by T. Huckle.
/////////////////////////////////////////////////////////////
template <class T_Field>
class Fspai_Unrestrained : public Fspai_Base<T_Field>
{
    public:

        // Inheriting template members
        using Fspai_Base<T_Field>::mtx;
        using Fspai_Base<T_Field>::precond;
        using Fspai_Base<T_Field>::P;
        using Fspai_Base<T_Field>::residual_norm;
        using Fspai_Base<T_Field>::Jk;
        using Fspai_Base<T_Field>::fspai_sub;
        using Fspai_Base<T_Field>::comm_handler;
        using Fspai_Base<T_Field>::hash_param;
        using Fspai_Base<T_Field>::updates_param;
        using Fspai_Base<T_Field>::max_idcs_param;
        using Fspai_Base<T_Field>::use_mean_param;
        using Fspai_Base<T_Field>::epsilon_param;

        /////////////////////////////////////////////////////////
        ///     \brief  Constructor
        ///
        ///     \param mtx_ Local matrix chunk on this pe
        ///     \param precond_ Local preconditioner chunk on this
        ///                     pe
        ///     \param P_ Local pattern chunk on this pe
        ///     \param env_handler Environment handler
        ///     \param hash_param Size of hash table
        ///     \param epsilon_param Epsilon as tolerance for
        ///            residual
        ///     \param updates_param Number of update steps
        ///     \param max_idcs_param Max. number of indices to be
        ///            augmented per update step
        ///     \param use_mean_param Whether to use the mean value
        ///            bound when augmenting indices per ste
        /////////////////////////////////////////////////////////
        Fspai_Unrestrained
        (   Matrix<T_Field>*  mtx_,
            Matrix<T_Field>*& precond_,
            Pattern*          P_,
            ENV_Handler&      env_handler,
            const int         hash_param,
            const double      epsilon_param,
            const int         updates_param,
            const int         max_idcs_param,
            const bool        use_mean_param ) :
                Fspai_Base<T_Field>(
                   mtx_, precond_, P_, env_handler,
                   hash_param, epsilon_param, updates_param,
                   max_idcs_param, use_mean_param  ) { }

    protected:

        /////////////////////////////////////////////////////////
        ///     \brief  Computing FSPAI for one preconditioner
        ///             column.
        ///     \param col Column to comupute preconditioner for
        /////////////////////////////////////////////////////////
        /*virtual*/
        void    Fspai_Column ( const int col );


        /////////////////////////////////////////////////////////
        ///     \brief  Computing Approximate Inverse for L_k
        ///
        ///     \param Lk_tilde Buffer where the solution will be
        ///                     written to
        ///     \param red_mtx According to pattern set Jk_tilde
        ///                    reduced submatrix A(Jk_tilde,Jk_tilde)
        ///     \param red_col According to pattern set Jk_tilde
        ///                    reduced column A(Jk_tilde,k)
        ///     \param dim Dimension of subsystem
        ///     \param col Currenct processed column k
        ///     \param A_kk Diagonal element of A = A_kk
        /////////////////////////////////////////////////////////
        void    Approximate_Inverse
                ( T_Field*&     Lk,
                  T_Field*      red_mtx,
                  T_Field*      red_col,
                  const int     dim,
                  const int     col,
                  const T_Field A_kk);
};

#include "Fspai_Unrestrained.imp.h"

#endif
