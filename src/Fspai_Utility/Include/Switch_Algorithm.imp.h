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


template <class T_Field>  Fspai_Base<T_Field>*
Switch_Algorithm<T_Field>::Get_Algorithm
(  const int         alg_level,
   Matrix<T_Field>*  mtx,
   Matrix<T_Field>*& precond,
   Pattern*          P,
   ENV_Handler&      env_handler,
   const int         hash_param,
   const double      epsilon_param,
   const int         updates_param,
   const int         max_idcs_param,
   const bool        use_mean_param)
{
    if ( mtx->my_id == 0 )   std::cout << "\t    Optimizations:\n";
    switch( alg_level )
    {
        case unrestrained: // Unrestrained FSPAI without any dictionary
            if ( mtx->my_id == 0 )
            {
                std::cout << "\t    CACHE:\t";
                std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;
                std::cout << "\t    HASH:\t";
                std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;

                // todo: Abstrakte klasse für Lapack und CSparse gemeinsam!!!!
                std::cout << "\t    CHOLESKY DATA MODE:\t";
                std::cout << COLOR_RED << "Lapack\n" << COLOR_NORMAL << std::endl;
            }
            return new Fspai_Unrestrained<T_Field>(
                    mtx, precond, P, env_handler, hash_param, epsilon_param,
                    updates_param, max_idcs_param, use_mean_param);

        // no default case
    }
    //dirty - suppress warning message
    return NULL;
}
