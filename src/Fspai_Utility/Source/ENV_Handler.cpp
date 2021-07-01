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

// file includings
#include "../Include/ENV_Handler.h"


void
ENV_Handler::Get_Environment_Params
(   int &num_procs_,
    int &my_id_ ) const
{
    num_procs_  = num_procs;
    my_id_      = my_id;
}



void
ENV_Handler::Initialize_Environment
(   int     argc,
    char*   argv[] )
{
    num_procs = 1;
    my_id = 0;
}



void
ENV_Handler::Set_Help
(   std::string& known_params,
    std::string& options ) const
{
    known_params =
        "-ep -ns -mn -um -wp -out -sol -sol_tol -sol_maxit -sol_out -rhs";
    options      =
        "\n"
        "\t  -ep:  epsilon parameter for FSPAI\n"
        "\t\tdefault is ep = 0.4\n"
        "\n"
        "\t  -mn:  maximum number of new nonzero candidate per step\n"
        "\t\tdefault is mn = 5\n"
        "\n"
        "\t  -ns:  maximum number of improvement steps per row in FSPAI\n"
        "\t\tdefault is ns = 5\n"
        "\n"
        "\t  -wp:  write preconditioner in matrix market format to file precond.mtx\n"
        "\t\tdefault is vb = 1\n"
        "\n"
        "\t  -um:  Use mean value as bound for augmenting indices.\n"
        "\t\tdefault is um = 1.\n"
        "\n"
        "\t  -sol: Parameter whether PCG solver is to be invoked or not.\n"
        "\t\t0: don't use PCG solver\n"
        "\t\t1: use PCG solver with computed FSPAI\n"
        "\t\t2: use unpreconditioned CG solver and PCG solver\n"
        "\t\t\twith computed FSPAI\n"
        "\t\tdefault is sol = 1.\n"
        "\n"
        "\t  -sol_tol: Solver tolerance.\n"
        "\t\tdefault is sol_tol = 1e-06.\n"
        "\n"
        "\t  -sol_maxit: Maximum number of iterations for solver.\n"
        "\t\tdefault is sol_maxit = 1000.\n"
        "\n"
        "\t  -sol_out: Output file name of solution.\n"
        "\t\tStore solution of system here.\n"
        "\t\tdefault is sol_out = ./solution.mtx\n"
        "\n"
        "\t  -out: Output file name.\n"
        "\t\tStore preconditioner here.\n"
        "\t\tdefault is out = ./precond.mtx\n"
        "\n"
        "\t  -rhs: Right-hand side parameter for solver\n"
        "\t\t0: use random vector as rhs\n"
        "\t\t1: use ones vector as rhs\n"
        "\t\tdefault is rhs = 1\n";
}



void
ENV_Handler::Get_Time
(   double&  curr_time,
    timespec tp ) const
{
    clock_gettime(CLOCK_REALTIME, &tp);
    curr_time = 1e9 * tp.tv_sec + tp.tv_nsec;
}



void
ENV_Handler::Time_Diff
(   double& max,
    double  sum_time ) const
{
    max = sum_time * 1e-9;
}



void*
ENV_Handler::Get_Communicator()
{
    return NULL;
}



void
ENV_Handler::Dist_Scan
(   int& my_nbr_cols,
    int& sum_nbr_cols ) const
{
    sum_nbr_cols = my_nbr_cols;
}



void
ENV_Handler::Dist_Local_Chunks
(   int   my_nbr_cols,
    int*& all_nbr_cols ) const
{
    all_nbr_cols[0] = my_nbr_cols;
}



void
ENV_Handler::Dist_Max
(   int  max,
    int& max_nnz ) const
{
    max_nnz = max;
}



void
ENV_Handler::Dist_Col_Length
(   int*        len_cols,
    const int   my_nbr_cols,
    int*&       len_all_cols,
    int*        all_nbr_cols,
    int*        start_indices ) const
{
    memcpy( len_all_cols, len_cols,
            my_nbr_cols * sizeof(int) );
}



void
ENV_Handler::Sum_NNZ
(   int     my_nnz,
    int&    nnz ) const
{
    nnz = my_nnz;
}



void
ENV_Handler::Check_Solver_Param
(   const int solver_param ) const
{
    std::stringstream   out_str;
    if(solver_param > 2)
    {
        out_str << solver_param;
        throw std::runtime_error(
            "\n\t\tIllegal solver parameter -sol "
            + out_str.str() + "\n"
            "\t\tIn a sequential environment only the computation\n"
            "\t\tof the FSPAI preconditioner is possible.\n"
            "\t\tThere is no other preconditioner implementation\n"
            "\t\tavailable for now.\n"
            "\t\tPlease use the parameter -sol 1 or -sol 2");
    }
}


#endif
