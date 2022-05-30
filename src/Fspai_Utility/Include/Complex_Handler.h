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

#ifndef COMPLEX_HANDLER_H_
#define COMPLEX_HANDLER_H_

// file includings
#include "Type_Base_Handler.h"
#include "ENV_Handler.h"
#include "Matrix.h"
#include "Matrix_Reader.h"
#include "Fspai_Base.h"
#include "Switch_Algorithm.h"
#include "PCG.h"

// C++/MPI includings
#include <iostream>


//////////////////////////////////////////
///     \class Complex_Handler
///     \brief Complex matrix field interface
///
///     Complex matrix interface.
//////////////////////////////////////////
class Complex_Handler : public Type_Base_Handler
{
    public:
        /////////////////////////////////////////////////////////
        ///     \brief Constructor
        ///
        ///     \param env_handler Environment handler
        ///     \param hash_size size of the hash table to be used.
        /////////////////////////////////////////////////////////
        Complex_Handler
        (	ENV_Handler& env_handler,
            const int 	 hash_size);

        /// Destructor
        ~Complex_Handler();

        /////////////////////////////////////////////////////////
        ///     \brief Parses matrix file, distributes all local
        ///            work chunks among all pes and generates
        ///            local Matrix data.
        ///
        ///     \param mmio MatrixMarket I/O object for having
        ///                 access to members.
        ///     \param matrix_file Matrix file path
        /////////////////////////////////////////////////////////
        void    Mtx_To_Memory
                (   const MMio	mmio,
                    char*       matrix_file );

        /////////////////////////////////////////////////////////
        ///     \brief  Printing matrix field to shell.
        /////////////////////////////////////////////////////////
        void    Print_Mtx_Type
                ( )                                         const;

        /////////////////////////////////////////////////////////
        ///     \brief  Prints matrix in human readable matrix
        ///             format.
        /////////////////////////////////////////////////////////
        void    Print_Matrix_Human_Readable
                ( )                                         const;

        /////////////////////////////////////////////////////////
        ///     \brief  Prints all matrix data.
        /////////////////////////////////////////////////////////
        void    Print_Matrix_Data
                ( )                                         const;

        /////////////////////////////////////////////////////////
        ///     \brief  Prints preconditioner in human readable
        ///             matrix format.
        /////////////////////////////////////////////////////////
        void    Print_Precond_Human_Readable
                ( )                                         const;

        /////////////////////////////////////////////////////////
        ///     \brief  Prints all preconditioner data.
        /////////////////////////////////////////////////////////
        void    Print_Precond_Data
                ( )                                         const;

        /////////////////////////////////////////////////////////
        ///     \brief  Getting number of columns/rows of matrix.
        ///
        ///     \return Number of columns of the system.
        /////////////////////////////////////////////////////////
        int     Get_Mtx_Dimension
                ( )                                         const;

        /////////////////////////////////////////////////////////
        ///     \brief  Getting number of local columns on this
        ///             pe.
        ///
        ///     \return Number of columns on this pe.
        /////////////////////////////////////////////////////////
        int     Get_Mtx_MyCols
                ( )                                         const;

        /////////////////////////////////////////////////////////
        ///     \brief  Getting index of local work chunk on this
        ///             pe.
        ///
        ///     \return The index of the column where the chunk
        ///             starts
        /////////////////////////////////////////////////////////
        int     Get_Mtx_StartIdx
                ( )                                         const;

        /////////////////////////////////////////////////////////
        ///     \brief  Converts a given matrix to a pattern.
        ///
        ///     \param P Pattern object to be filled
        /////////////////////////////////////////////////////////
        void    To_Pattern
                (Pattern *P)								const;

        /////////////////////////////////////////////////////////
        ///     \brief  Sets the FSPAI algorithm.
        ///
        ///     \param alg_level User given input which decides
        ///                      which FSPAI algorithm to use.
        ///     \param P Local pattern chunk on this pe
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
        void    Set_Fspai_Algorithm
                (   const int       alg_level,
                    Pattern*        P,
                    const int       hash_param,
                    const double    epsilon_param,
                    const int       updates_param,
                    const int       max_idcs_param,
                    const bool      use_mean_param);

        /////////////////////////////////////////////////////////
        ///     \brief  Initializes preconditioner and invokes
        ///             FSPAI algorithm.
        /////////////////////////////////////////////////////////
        void    Invoke_Fspai(  );

        /////////////////////////////////////////////////////////
        ///     \brief  Initializes solver.
        ///
        ///     \param tol solver tolerance
        ///     \param maxit maximum number of iterations
        ///     \param para_max_levels ParaSails level parameter
        ///     \param para_threshold ParaSails threshold parameter
        ///     \param para_filter ParaSails filter parameter
        /////////////////////////////////////////////////////////
        void    Init_Solver
                (   const double tol,
                    const int    maxit,
                    const int    para_max_levels,
                    const double para_threshold,
                    const double para_filter);

        /////////////////////////////////////////////////////////
        ///     \brief  Invokes (preconditioned) conjugate gradient
        ///             method.
        ///
        ///     \param rhs_choice Which rhs to choose
        ///     \param solver_param User parameter for solver
        ///     \param my_id Id of this PE.
        /////////////////////////////////////////////////////////
        void    Invoke_PCG
                (   const int  rhs_choice,
                    const int  solver_param,
                    const int  my_id);

        /////////////////////////////////////////////////////////
        ///     \brief  Writes system solution to file.
        ///
        ///     \param  output_file Path of solution output file.
        /////////////////////////////////////////////////////////
        void    Solution_To_File
                (   const char* output_file );

        /////////////////////////////////////////////////////////
        ///     \brief  Writes FSPAI matrix to specified file:
        ///
        ///     \param output_file Path of output file
        /////////////////////////////////////////////////////////
        void    Fspai_To_File
                (const char* output_file );

        // Member variables

        /// COMPLEX matrix
        Matrix<COMPLEX>*        mtx;

        /// COMPLEX preconditioner matrix
        Matrix<COMPLEX>*        precond;

        /// Interface to environment handler
        ENV_Handler             env_handler;

    protected:

        /// FSPAI Algorithm pointer
        Fspai_Base<COMPLEX>*    Fspai_base_algorithm;

        /// PCG solver pointer
        PCG<COMPLEX>*           solver;
};

#endif /* COMPLEX_HANDLER_H_ */
