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

////////////////////////////////////////////////////////////
// This is the main page of the doxygen documentation
////////////////////////////////////////////////////////////
/// \mainpage FSPAI (Factorized Sparse Approximate Inverses)
/// \section section1 ABOUT
/// FSPAI is the implementation of the known FSPAI algorithm
/// introduced by T. Huckle. \n It generates a preconditioner for large sparse
/// and ill-conditioned symmetric positive \n definite systems
/// of linear equations. FSPAI is inherently parallel and the  computed \n
/// preconditioner approximates the inverse of the Cholesky
/// factor of the system matrix. \n It is the factorized version of
/// the SPAI algorithm. The \e FSPAI was implemented at the \n research unit \n\n
/// \e Informatik \e V -- \e Scientific \e Computing \e in \e Computer \e Science \n
/// \e Technische \e Universität \e München. \n \n \n
/// \b DESIGNED \b BY \n \n
/// Matous Sedlacek <sedlacek@in.tum.de> \n \n \n
/// \b RELEASED \b 2011 \n \n
/// FSPAI is published under the LGPL in year 2011. \n \n \n
/// \b LICENSE \n \n
/// FSPAI: Factorized Sparse Approximate Inverses \n
/// Copyright © 2011 Matous Sedlacek \n
/// Scientific Computing in Computer Science -- Informatics V \n
/// Technische Universität München \n \n
/// FSPAI is free software: you can redistribute it and/or modify it under the \n
/// terms of the GNU Lesser General Public License as published by the Free Software \n
/// Foundation, either version 3 of the License, or (at your option) any later version. \n
/// \n
/// FSPAI is distributed in the hope that it will be useful, but WITHOUT ANY \n
/// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR \n
/// A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details. \n
/// \n
/// You should have received a copy of the GNU Lesser General Public License along with \n
/// FSPAI. If not, see http://www.gnu.org/licenses/. \n
/// \n
/// If you obtain any results with FSPAI we would appreciate that you refer to FSPAI. \n
/////////////////////////////////////////////////////////////

#ifndef COMMAND_PARSER_H_
#define COMMAND_PARSER_H_

//file includings
#include "Param_Map.h"
#include "ENV_Handler.h"

// C/C++ includings
#include <string>


////////////////////////////////////////
///     \class Command_Parser
///     \brief This class does the
///            mapping from shell input
///            to the program variables.
////////////////////////////////////////
class Command_Parser
{
    public:
        /// Constructor
        Command_Parser
        (   ENV_Handler& env_handler );

        /////////////////////////////////////////////////
        ///     \brief  Reads parameters from shell to
        ///             member variables
        ///
        ///     \param argc Number of shell parameters
        ///     \param argv Shell parameters
        ///     \param env_handler Environment handler
        /////////////////////////////////////////////////
        void    Read_Parameters
                (   int          argc,
                    char         *argv[],
                    ENV_Handler& env_handler );

        // Member variables

                /// Epsilon tolerance (stopping criterion)
        double  epsilon_param,

                /// Solver tolerance
                solver_tol,

                // ParaSails threshold
                para_threshold,

                // ParaSails filter parameter
                para_filter;

                /// Whether to write FSPAI to file or not
        int     write_param,

                /// How big the hash table will be for remote data
                hash_param,

                /// Pattern parameter as switch for pattern generation
                pattern_param,

                /// Which Fspai algorithm will be invoked
                alg_level,

                /// Number of pattern updates
                updates_param,

                /// Solver parameter: 0-> no solver, 1->FSPAI solver,
                /// 2->unprec. solver & FSPAI solver
                solver_param,

                /// number of maximum iterations in solver
                solver_maxit,

                /// Number of augmenting indices per step
                max_idcs_param,

                /// Rhs to use in solver
                rhs_param,

                // ParaSails parameter for levels
                para_max_levels;

                /// Path to pattern file
        char    *pattern_file,

                /// Path to matrix file
                *matrix_file,

                /// Path of output file
                output_file[],

                /// Path of solution output file
                solver_file[];

                /// Whether to use mean value as upper bound for augmenting
                /// indices in tau calculation
        bool    use_mean_param;

    private:

        // Methods

        /////////////////////////////////////////////////
        ///     \brief  Filling parammap with shell
        ///             parameters.
        ///
        ///     \param argc Number of shell parameters
        ///     \param argv Shell parameters
        /////////////////////////////////////////////////
        void    Fill_Param_Map( int argc, char *argv[] );

        /////////////////////////////////////////////////
        ///     \brief  Maps parammap to member variables.
        ///
        ///     The entries of the already filled
        ///     parammap are type converted and mapped to
        ///     the specific member variables.
        /////////////////////////////////////////////////
        void    Set_Params( );

        /////////////////////////////////////////////////
        ///     \brief  Checks whether user made valid
        ///             input.
        /////////////////////////////////////////////////
        void    Check_Params(ENV_Handler& env_handler );

        // Member variables

        // Environment specific help string
        std::string     options;

        /// Containing all environment known parameters
        std::string     known_params;

        /// Holds all shell parameters
        Param_Map       parameters;
};

#endif
