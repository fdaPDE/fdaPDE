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


template <class T_Field> void
Distribution<T_Field>::Basic_Distribution
(   ENV_Handler&            env_handler,
    int                     cols,
    int&                    my_nbr_cols,
    int&                    split_pe,
    int&                    split_idx,
    int&                    start_idx,
    const Comm_Handler<T_Field>& comm_handler)
{
    int     small_chunk,
            big_chunk,
            rem,
            sum_nbr_cols,
            num_procs,
            my_id;

    env_handler.Get_Environment_Params( num_procs, my_id );
    env_handler.Barrier();

    small_chunk = static_cast< int >
        ( floor( static_cast< double > ( cols ) / num_procs ) );
    big_chunk   = static_cast< int >
        ( ceil( static_cast< double > ( cols ) / num_procs ) );
    rem         = static_cast< int >
        (fmod(static_cast< double > ( cols ),
              static_cast< double > ( num_procs ) ) );

    if ( my_id < rem )
        my_nbr_cols = big_chunk;
    else
        my_nbr_cols = small_chunk;

    // Computes sum of work chunks
    comm_handler.Dist_Scan( my_nbr_cols, sum_nbr_cols );

    split_pe = rem;
    split_idx = big_chunk * rem;
    start_idx = sum_nbr_cols - my_nbr_cols;

    // Adjust my_nbr_cols in last processor
    // if n is not a multiple of chuck size
    if ( my_id == ( num_procs-1 ) )
        if ( sum_nbr_cols != cols )
            my_nbr_cols -=  ( sum_nbr_cols - cols );
}
