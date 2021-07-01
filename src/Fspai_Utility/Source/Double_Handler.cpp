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

// file includings
#include "../Include/Double_Handler.h"


Double_Handler::Double_Handler
(	ENV_Handler& env_handler_,
	const int 	 hash_size)
{
    mtx 					= NULL;
    precond 				= NULL;
    env_handler 			= env_handler_;
    Fspai_base_algorithm 	= NULL;
    solver                  = NULL;
}



Double_Handler::~Double_Handler( )
{
    if(mtx)                     delete mtx;
    if(precond)                 delete precond;
    if(Fspai_base_algorithm)    delete Fspai_base_algorithm;
    if(solver)                  delete solver;
}



void
Double_Handler::Mtx_To_Memory
(   const MMio	mmio,
    char*       matrix_file )
{
    Matrix_Reader<double, RCV>().To_Memory
        ( mtx, matrix_file,
          mmio, env_handler );
}



void
Double_Handler::Print_Mtx_Type
( ) const
{
    std::cout << "\t  Matrix is REAL!" << std::endl;
}



void
Double_Handler::Print_Matrix_Human_Readable
( ) const
{
    mtx->Print_Matrix_Human_Readable( );
}



void
Double_Handler::Print_Matrix_Data
( ) const
{
    mtx->Print_Matrix_Data();
}



void
Double_Handler::Print_Precond_Human_Readable
( ) const
{
    precond->Print_Matrix_Human_Readable( );
}



void
Double_Handler::Print_Precond_Data
( ) const
{
    precond->Print_Matrix_Data();
}



int
Double_Handler::Get_Mtx_Dimension
( ) const
{
    return mtx->n;
}



int
Double_Handler::Get_Mtx_MyCols
( ) const
{
    return mtx->my_nbr_cols;
}



int
Double_Handler::Get_Mtx_StartIdx
( ) const
{
    return mtx->my_start_idx;
}



void Double_Handler::To_Pattern
(   Pattern* P) const
{
    mtx->To_Pattern(P, env_handler);
}



void Double_Handler::Set_Fspai_Algorithm
(   const int       alg_level,
    Pattern*        P,
    const int 	    hash_param,
    const double    epsilon_param,
    const int       updates_param,
    const int       max_idcs_param,
    const bool      use_mean_param)
{
    Fspai_base_algorithm  =
        Switch_Algorithm<double>().Get_Algorithm(
                alg_level, mtx, precond, P, env_handler,
                hash_param, epsilon_param, updates_param,
                max_idcs_param, use_mean_param);
}



void Double_Handler::Invoke_Fspai
(  )
{
    precond = new Matrix<double>( env_handler, mtx->n, mtx->n );
    precond->Init_Preconditioner( mtx, env_handler );
    // Invoke FSPAI Algorithm
    Fspai_base_algorithm->Fspai_Algorithm( env_handler );
}



void Double_Handler::Init_Solver
(   const double tol,
    const int    maxit,
    const int    para_max_levels,
    const double para_threshold,
    const double para_filter)
{
    solver = new PCG<double>(env_handler,
            Fspai_base_algorithm->comm_handler, mtx, tol, maxit,
            para_max_levels, para_threshold, para_filter);
}



void Double_Handler::Invoke_PCG
(   const int  rhs_choice,
    const int  solver_param,
    const int  my_id)
{
    solver->Create_RHS(rhs_choice, mtx->my_nbr_cols);
    switch(solver_param)
    {
        case 1:
            if( my_id == 0 )
            {
                std::cout << "\n\t* Solving system with PCG "
                "using FSPAI... ";
                std::cout.flush();
            }
            solver->Set_Precond(precond);
            solver->Init_Solution(mtx->my_nbr_cols);
            solver->Init_ParTmp(mtx->my_nbr_cols);
            solver->Solve_System();
            break;
        case 2:
            if( my_id == 0 )
            {
                std::cout << "\n\t* Solving system with "
                "unprecond. CG...   ";
                std::cout.flush();
            }
            solver->Set_Precond(NULL);
            solver->Init_Solution(mtx->my_nbr_cols);
            solver->Init_ParTmp(mtx->my_nbr_cols);
            solver->Solve_System();
            if( my_id == 0 )
            {
                std::cout << "\n\t* Solving system with PCG using "
                "FSPAI... ";
                std::cout.flush();
            }
            solver->Set_Precond(precond);
            solver->Init_Solution(mtx->my_nbr_cols);
            solver->Init_ParTmp(mtx->my_nbr_cols);
            solver->Solve_System();
            break;
        case 3:
            if( my_id == 0 )
            {
                std::cout << "\n\t* Solving system with PCG using "
                "FSPAI... ";
                std::cout.flush();
            }
            solver->Set_Precond(precond);
            solver->Init_Solution(mtx->my_nbr_cols);
            solver->Init_ParTmp(mtx->my_nbr_cols);
            solver->Solve_System();
            // ParaSails
            solver->Set_Precond(NULL);
            solver->Init_Solution(mtx->my_nbr_cols);
            solver->Init_ParTmp(mtx->my_nbr_cols);
            solver->Solve_Precond_ParaSails();
            break;
    }
}



void Double_Handler::Solution_To_File
(   const char* output_file )
{
    solver->Solution_To_File( output_file );
}



void Double_Handler::Fspai_To_File
(   const char* output_file )
{
    precond->To_File( env_handler, output_file );
}
