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
#include "../Include/FSPAI_Solver_Wrapper.h"

//C++/MPI includings
#include <iostream>
#include <stdexcept>

//file includings
#include "../Include/ENV_Handler.h"
#include "../Include/Macros.h"
#include "../Include/Pe_Exception.h"
#include "../Include/Command_Parser.h"
#include "../Include/Timer.h"
#include "../Include/Type_Base_Handler.h"
#include "../Include/Double_Handler.h"
#include "../Include/Complex_Handler.h"
#include "../Include/MMio.h"
#include "../Include/Pattern.h"



int FSPAI_Solver_Wrapper(int argc, char *argv[])
{
    int                 num_procs,
                        my_id;
    Timer               timer;
    Type_Base_Handler*  type_handler = NULL;
    MMio                mmio;
    ENV_Handler         env_handler = ENV_Handler();
    Pattern             *P = NULL;

    // Initializing environment-based parameters
    env_handler.Initialize_Environment( argc, argv );

    // Get environment parameters
    env_handler.Get_Environment_Params( num_procs, my_id );

    // No exception class because only message string will be different
    try
    {
        if( my_id == 0 )
        {
            std::cout <<
            "\n\t==============================================================\n"
            "\t= \t\t\t\t\t\t\t     = \n"
            "\t= FSPAI v1.0  Copyright (C) 2010, 2011  Matous Sedlacek\t     =\n"
            "\t= \tChair of Scientific Computing -- Informatics V\t     =\n"
            "\t= \tTechnische Universität München\t\t\t     =\n"
            "\t= FSPAI v1.0 comes with ABSOLUTELY NO WARRANTY. This is free =\n"
            "\t= software, and you are welcome to redistribute it under     =\n"
            "\t= certain conditions; see the documentation or LICENSE file. ="
            "\n\t= \t\t\t\t\t\t\t     ="
            "\n\t==============================================================\n"
            "\t=====================   STARTING FSPAI   =====================\n\n"
            "\t-> on " << num_procs << " processor(s) " << std::endl;
            std::cout.flush();
        }

        // Getting the user input parameters from shell
        if( my_id == 0 )
        {
            std::cout << "\n\t* Reading input parameters...\t\t" << std::endl;
            std::cout.flush();
        }
        Command_Parser comm_parser = Command_Parser( env_handler );
        comm_parser.Read_Parameters( argc, argv, env_handler );

        // Start global time measurement
        timer.Start( env_handler );

        // Runtime polymorphism for matrix (real or complex)
        switch( mmio.Get_Matrix_Type( comm_parser.matrix_file ) )
        {
            case real:
            {
                type_handler = new Double_Handler(env_handler,
                                                  comm_parser.hash_param);
                break;
            }
            case complex:
            {
                // Checking for solver parameter as solution of
                // complex matrices in parallel environment not possible
                // via hypre-package
                env_handler.Check_Solver_Type(comm_parser.solver_param);
                type_handler = new Complex_Handler(env_handler,
                                                   comm_parser.hash_param);
                break;
            }
        }
        // Reading data input and generating the matrix object
        if( my_id == 0 )
        {
            type_handler->Print_Mtx_Type();
            std::cout << "\n\t* Reading matrix data...\t\t ";
            std::cout.flush();
        }
        type_handler->Mtx_To_Memory( mmio, comm_parser.matrix_file );

        // Reading pattern file and generating pattern
        if( my_id == 0 )
        {
            std::cout << "\n\t* Generating pattern data...\t\t ";
            std::cout.flush();
        }
        P = new Pattern();
        switch( comm_parser.pattern_param )
        {
            case 0: //Own pattern file
                P->Arbitrary_Pattern( env_handler, comm_parser.pattern_file,
                                      type_handler->Get_Mtx_Dimension( ) );
                break;
            case 1: //Diagonal pattern
                P->Diagonal_Pattern( env_handler,type_handler->Get_Mtx_Dimension( ),
                                     type_handler->Get_Mtx_MyCols( ),
                                     type_handler->Get_Mtx_StartIdx( ));
                break;
            case 2: //Lower triangular pattern of the system matrix
                type_handler->To_Pattern( P );
                break;
            //default case will not occur
        }

        // Checking algorithm level and getting the
        // requested SPAI algorithm
        if( my_id == 0 )
            std::cout << "\n\t* Checking algorithm level... "
                      << std::endl;
        type_handler->Set_Fspai_Algorithm( comm_parser.alg_level, P,
										   comm_parser.hash_param,
										   comm_parser.epsilon_param,
										   comm_parser.updates_param,
										   comm_parser.max_idcs_param,
										   comm_parser.use_mean_param);

        // Compute FSPAI with requested FSPAI algorithm
        if ( my_id == 0 )
        {
            std::cout << "\t* Computing FSPAI..." << std::endl;
            std::cout << "\t    Pattern Updates:\t"
                      << comm_parser.updates_param << std::endl;
            std::cout << "\t    Idcs. per Update:\t"
                      << comm_parser.max_idcs_param << std::endl;
            std::cout << "\t    Epsilon Tol.:\t"
                      << comm_parser.epsilon_param << std::endl;
            std::cout << "\t    Mean-Value:\t\t";
            if(comm_parser.use_mean_param) std::cout << "yes\t\t ";
            else                           std::cout << "no\t\t ";

            std::cout.flush();
        }
        type_handler->Invoke_Fspai( );

        // Invoke solver if requested
        if( comm_parser.solver_param )
        {
            type_handler->Init_Solver(comm_parser.solver_tol,
                                      comm_parser.solver_maxit,
                                      comm_parser.para_max_levels,
                                      comm_parser.para_threshold,
                                      comm_parser.para_filter);
            type_handler->Invoke_PCG(comm_parser.rhs_param,
                                     comm_parser.solver_param,
                                     my_id);
            type_handler->Solution_To_File(comm_parser.solver_file);
        }

        // Write preconditioner to file if requested
        if ( comm_parser.write_param )
        {
            if( my_id == 0 )
            {
                std::cout << "\n\t* Writing FSPAI to file " +
                    std::string( comm_parser.output_file ) + "...   ";
                std::cout.flush();
            }
            type_handler->Fspai_To_File( comm_parser.output_file );
        }

        // Stop time measurement
        if (my_id == 0)
            std::cout << "\t\t\t\t_____________________________________\n"
                           "\n\t\t\t\tTotal time: \t ";

        timer.Stop( env_handler );
        timer.Report( env_handler );
    }
    catch(Pe_Exception& ex)
    {
        delete type_handler;
        delete P;
        std::cerr << ex.what() << "\n\t============================"
                "==================================\n" << std::endl;
        env_handler.Abort_Environment();
        return EXIT_FAILURE;
    }
    catch ( std::runtime_error& ex )
    {
        delete type_handler;
        delete P;
        if ( my_id == 0 )
        {
            std::cerr << ex.what() << "\n\t============================"
                    "==================================\n" << std::endl;
            env_handler.Abort_Environment();
        }
        return EXIT_FAILURE;
    }
    catch( std::bad_alloc &ex )
    {
        delete type_handler;
        delete P;
        if ( my_id == 0 )
        {
            std::cerr << "\n\tPROGRAM ABORTED: " << ex.what() <<
                "\n\t\t  ERROR: Not enough memory.\n\n"
                "\t==================================="
                "===========================\n" << std::endl;
            env_handler.Abort_Environment();
        }
        return EXIT_FAILURE;
    }

    delete P;
    delete type_handler;
    if( my_id == 0 )
    {
        std::cout << "\n\n\t==================   SUCCESSFULLY FINISHED   "
        "=================\n\t=========================================="
        "====================\n" << std::endl;
    }

    // Closing environment
    env_handler.Barrier();
    env_handler.Finalize();
	return EXIT_SUCCESS;
}
