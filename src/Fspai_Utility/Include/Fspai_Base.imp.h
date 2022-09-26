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


template <class T_Field>
Fspai_Base<T_Field>::Fspai_Base
(   Matrix<T_Field>*  mtx_,
    Matrix<T_Field>*& precond_,
    Pattern*          P_ ,
    ENV_Handler&      env_handler,
    const int         hash_param_,
    const double      epsilon_param_,
    const int         updates_param_,
    const int         max_idcs_param_,
    const bool        use_mean_param_)
        :   mtx( mtx_ ),
            precond( precond_ ),
            P(P_),
            hash_param(hash_param_),
            epsilon_param(epsilon_param_),
            updates_param(updates_param_),
            max_idcs_param(max_idcs_param_),
            use_mean_param(use_mean_param_)
{
    residual_norm = 1e+50;
    Jk            = NULL;
    comm_handler  = new Comm_Handler<T_Field>(env_handler, hash_param,
                                             mtx_, precond_,P_);
    fspai_sub     = new Fspai_Sub<T_Field>(comm_handler,mtx->n);
}



template <class T_Field>  void
Fspai_Base<T_Field>::Fspai_Algorithm
(   ENV_Handler& env_handler )
{
    int		col = -1;
    Timer   timer = Timer();

    // Start time measurement
    env_handler.Barrier();
    timer.Start( env_handler );

    //Iterating over each column
    while( true )
    {
        col = comm_handler->LB_Get_Precond_Col();
        if (col < 0) break;
        Fspai_Column( col );
    }

    comm_handler->Finalize_Communication();
    env_handler.Barrier();

    // Stop time measurement
    timer.Stop( env_handler );
    timer.Report( env_handler );
}
