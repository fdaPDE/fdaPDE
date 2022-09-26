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

#ifndef TYPE_BASE_HANDLER_H_
#define TYPE_BASE_HANDLER_H_

// file includings
#include "ENV_Handler.h"
#include "MMio.h"
#include "Pattern.h"


//////////////////////////////////////////
///     \class Type_Base_Handler
///     \brief Matrix Field interface
///
///     Type independant base interface
///     for access for field specific
///     Matrix methods.
//////////////////////////////////////////
class Type_Base_Handler
{
    public:
        /// Constructor
                        Type_Base_Handler() { }

        /// Destructor
        virtual         ~Type_Base_Handler() { }

        /////////////////////////////////////////////////////////
        ///     \brief Parses matrix file, distributes all local
        ///            work chunks among all pes and generates
        ///            local Matrix data.
        ///
        ///     \param mmio MatrixMarket I/O object for having
        ///                 access to members.
        ///     \param matrix_file Matrix file path
        /////////////////////////////////////////////////////////
        virtual void    Mtx_To_Memory
                        (   const MMio mmio,
                            char* matrix_file ) = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Printing matrix field to shell.
        /////////////////////////////////////////////////////////
        virtual void    Print_Mtx_Type
                        ( )                             const = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Prints matrix in human readable matrix
        ///             format.
        /////////////////////////////////////////////////////////
        virtual void    Print_Matrix_Human_Readable
                        ( )                             const = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Prints all Matrix data.
        /////////////////////////////////////////////////////////
        virtual void    Print_Matrix_Data
                        ( )                             const = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Prints preconditioner in human readable
        ///             matrix format.
        /////////////////////////////////////////////////////////
        virtual void    Print_Precond_Human_Readable
                ( )                                     const = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Prints all preconditioner data.
        /////////////////////////////////////////////////////////
        virtual void    Print_Precond_Data
                ( )                                     const = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Getting number of columns/rows of matrix.
        ///
        ///     \return Number of columns of the system.
        /////////////////////////////////////////////////////////
        virtual int     Get_Mtx_Dimension
                        ( )                             const = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Getting number of local columns on this
        ///             pe.
        ///
        ///     \return Number of columns on this pe.
        /////////////////////////////////////////////////////////
        virtual int     Get_Mtx_MyCols
                        ( )                             const = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Getting index of local work chunk on this
        ///             pe.
        ///
        ///     \return The index of the column where the chunk
        ///             starts
        /////////////////////////////////////////////////////////
        virtual int     Get_Mtx_StartIdx
                        ( )                             const = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Converts a given matrix to a pattern.
        ///
        ///     \param P Pattern object to be filled
        /////////////////////////////////////////////////////////
        virtual void    To_Pattern
                        (   Pattern* P  ) 				const = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Sets the FSPAI algorithm.
        ///
        ///     \param alg_level User given input which decides
        ///                      which FSPAI algorithm to use.
        ///     \param P Local start pattern chunk on this pe.
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
        virtual  void   Set_Fspai_Algorithm
                        (   const int       alg_level,
                            Pattern*        P ,
                            const int       hash_param,
                            const double    epsilon_param,
                            const int       updates_param,
                            const int       max_idcs_param,
                            const bool      use_mean_param ) = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Initializes preconditioner and invokes
        ///             FSPAI algorithm.
        /////////////////////////////////////////////////////////
        virtual void    Invoke_Fspai(  )                     = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Initializes solver.
        ///
        ///     \param tol solver tolerance
        ///     \param maxit maximum number of iterations
        ///     \param para_max_levels ParaSails level parameter
        ///     \param para_threshold ParaSails threshold parameter
        ///     \param para_filter ParaSails filter parameter
        /////////////////////////////////////////////////////////
        virtual void    Init_Solver
                        (   const double tol,
                            const int    maxit,
                            const int    para_max_levels,
                            const double para_threshold,
                            const double para_filter)         = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Invokes (preconditioned) conjugate gradient
        ///             method.
        ///
        ///     \param rhs_choice Which rhs to choose
        ///     \param solver_param User parameter for solver
        ///     \param my_id Id of this PE.
        /////////////////////////////////////////////////////////
        virtual void    Invoke_PCG
                        (   const int  rhs_choice,
                            const int  solver_param,
                            const int my_id)                = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Writes system solution to file.
        ///
        ///     \param  output_file Path of solution output file.
        /////////////////////////////////////////////////////////
        virtual void    Solution_To_File
                        (   const char* output_file )        = 0;

        /////////////////////////////////////////////////////////
        ///     \brief  Writes FSPAI matrix to specified file:
        ///
        ///     \param output_file Path of output file
        /////////////////////////////////////////////////////////
        virtual void    Fspai_To_File
                        (   const char* output_file )        = 0;
};

#endif /* TYPE_BASE_HANDLER_ */
