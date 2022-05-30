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

#ifndef ENV_MPI

#   ifndef COMM_HANDLER_H_
#   define COMM_HANDLER_H_

// C++ includings
#include <iostream>
#include <stdlib.h>
#include <string.h>

// file includings
#include "ENV_Handler.h"
#include "Matrix.h"
#include "Pattern.h"


///////////////////////////////////////////////////////////
///     \class Comm_Handler
///     \brief Handling communication due to specified
///            environment. Responsible for the MPI
///            communication between the remote cluster
///            nodes. Every FSPAI Algorithm has a
///            communication handler which includes
///            an environment handler. For the sequential
///            version no communication becomes necessary.
/////////////////////////////////////////////////////////
template<class T_Field>
class Comm_Handler
{
    public:

        /// Empty constructor
        Comm_Handler<T_Field>() { }

        /////////////////////////////////////////////////////////
        ///     \brief Constructor
        ///
        ///     \param env_handler Environment handler.
        ///     \param hash_size size of the hash table to be used.
        ///     \param mtx_ Local matrix data
        ///     \param precond_ Local preconditioner chunk
        ///     \param P_ Local pattern chunk
        /////////////////////////////////////////////////////////
        Comm_Handler<T_Field>
            (   ENV_Handler&      env_handler,
                const int         hash_size,
                Matrix<T_Field>*  mtx_,
                Matrix<T_Field>*& precond_,
                Pattern*          P_) :
                    mtx( mtx_ ), precond( precond_ ), P(P_) {}

        /// Destructor
        ~Comm_Handler<T_Field>() {}

        /////////////////////////////////////////////////////
        ///     \brief Provided for parallel implementation
        ///            which uses this method to drive on
        ///            communication between pes.
        /////////////////////////////////////////////////////
        void    Communicate( ) { }

        ////////////////////////////////////////////////////
        ///     \brief  Load/Work balancing mechanism.
        ///
        ///     Responsible for balancing the work
        ///     between all pes. Get the new column
        ///     index for which a preconditioner solution
        ///     should be computed for. In the sequential
        ///     version no load/work balancing necessary.
        ///
        ///     \return Column index to be work balanced.
        /////////////////////////////////////////////////////
        int     LB_Get_Precond_Col();

        ///////////////////////////////////////////////////////
        ///     \brief Requesting column data of the input
        ///            matrix.
        ///
        ///     If the requested column data is local, get
        ///     it from local cache. If not, look into the
        ///     hash table if this column was requested before.
        ///     If not, send the request to the specific pe and
        ///     receive the data. Finally store this request
        ///     into hash table for later requests. Thus the
        ///     MPI traffic can be reduced. In the sequential
        ///     version only the local cache version will be
        ///     used.
        ///
        ///     \param col The comlumn which is requested
        ///     \param col_len The length of the requested column
        ///     \param col_idcs_buf Column indices of requested
        ///                         column
        ///     \param col_buf  Values of requested column
        ///////////////////////////////////////////////////////
        void    Get_Mtx_Col
                (   int         col,
                    int&        col_len,
                    int*&       col_idcs_buf,
                    T_Field*&   col_buf);

        /////////////////////////////////////////////////////
        ///     \brief Finalizes the communication for this pe.
        ///
        ///     After processing the local work chunk, the pe
        ///     signalizes that it has finished by invoking
        ///     this method. The corresponding sequential
        ///     method is empty.
        /////////////////////////////////////////////////////
        void    Finalize_Communication( )  { }

        /////////////////////////////////////////////////////
        ///     \brief Some pe needs a start pattern column
        ///            from a remote pe.
        ///
        ///     This is part of the load balancing mechanism.
        ///     If some pe has finished his work, he
        ///     previously requested a preconditioner column
        ///     and needs now the specific start pattern
        ///     column.
        ///
        ///     \param col The start pattern column to request
        ///     \param J The index set which is to be set
        /////////////////////////////////////////////////////
        void    Get_P_Col
                (   const int col,
                    Index_Set*& J );

        /////////////////////////////////////////////////////
        ///     \brief Inserting computed preconditioner
        ///            column into local cache
        ///
        ///     \param Lk_tilde Solution of preconditioner column
        ///            without L_kk
        ///     \param col The computed column to insert
        ///     \param Jk_tilde Index set of pattern without k
        ///     \param Jk Index set of pattern
        /////////////////////////////////////////////////////
        void    Insert_Precond_Solution
                (   const T_Field*   Lk,
                    const int        col,
                    const Index_Set* Jk_tilde,
                    const Index_Set* Jk);

    private:

        /// Matrix
        Matrix<T_Field>*  mtx;

        /// Solution/Preconditioner
        Matrix<T_Field>*& precond;

        /// Pattern
        Pattern*          P;
};

#   include "Comm_Handler.imp.h"

#   endif /* COMM_HANDLER_H_ */

#endif

