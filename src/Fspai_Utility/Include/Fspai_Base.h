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

#ifndef FSPAI_BASE_H_
#define FSPAI_BASE_H_


// file includings
#include "Fspai_Sub.h"
#include "Comm_Handler.h"
#include "ENV_Handler.h"


///////////////////////////////////////////////////
///     \class Fspai_Base
///     \brief This is the base class of
///            all FSPAI algorithms.
///
///     It provides the base loop whose
///     implementation depends on the derived
///     algorithms.
///     Furthermore it is checked, whether some pes
///     communicate or wait for remote data. Idle
///     pes are not waiting, but invoking the
///     communication method within Com_Handler.
///////////////////////////////////////////////////
template <class T_Field>
class Fspai_Base
{
    public:
        /////////////////////////////////////////////////////////
        ///     \brief  Constructor
        ///
        ///     \param mtx_ The system matrix to be preconditioned
        ///     \param precond_ The local preconditioner chunk
        ///     \param P_ The local start pattern chunk on this pe
    	///		\param env_handler Environment handler
    	///		\param hash_param Size of the hash table to be
    	///						  used.
        ///     \param epsilon_param Epsilon as tolerance for
        ///            residual
        ///     \param updates_param Number of update steps
        ///     \param max_idcs_param Max. number of indices to be
        ///            augmented per update step
        ///     \param use_mean_param Whether to use the mean value
        ///            bound when augmenting indices per step
        /////////////////////////////////////////////////////////
        Fspai_Base
        (   Matrix<T_Field>* 	mtx_,
            Matrix<T_Field>*& 	precond_,
            Pattern*         	P_,
            ENV_Handler& 		env_handler,
            const int           hash_param,
            const double        epsilon_param,
            const int           updates_param,
            const int           max_idcs_param,
            const bool          use_mean_param);

        /// Destructor
        virtual ~Fspai_Base<T_Field>() {
            if(fspai_sub) 		delete fspai_sub;
            if(comm_handler) 	delete comm_handler;
            };

        /////////////////////////////////////////////////////////
        ///     \brief  The base loop for all Fspai algorithms.
        ///
        ///     \param env_handler Interface to environment specific
        ///                        methods.
        /////////////////////////////////////////////////////////
        void            Fspai_Algorithm
                        ( ENV_Handler& env_handler );

        /// Interface to communication Handler
        Comm_Handler<T_Field>*    comm_handler;

    protected:
        /////////////////////////////////////////////////////////
        ///     \brief  Computing FSPAI for one preconditioner
        ///             column
        /////////////////////////////////////////////////////////
        virtual void    Fspai_Column( const int col )       = 0;

        /// Member pointer to matrix chunk
        const Matrix<T_Field>*    mtx;

        /// Member pointer to preconditioner chunk
        Matrix<T_Field>*&         precond;

        Pattern*                  P;
        /// Large first bound for "approximating" infinity
        double                    residual_norm;

        /// Index set J containing the local pattern indices
        Index_Set                 *Jk;

        /// Interface to FSPAI subroutines
        Fspai_Sub<T_Field>*       fspai_sub;

        /// How big the hash table is
        int                       hash_param;

        /// Epsilon tolerance (stopping criterion)
        double                    epsilon_param;

        /// Number of pattern updates
        int                       updates_param;

        /// Number of augmenting indices per step
        int                       max_idcs_param;

        /// Whether to use mean value as upper bound for augmenting
        /// indices in tau calculation
        bool                      use_mean_param;
};

#include "Fspai_Base.imp.h"

#endif
