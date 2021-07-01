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

#ifndef PCG_H_
#define PCG_H_

// file includings
#include "Matrix.h"
#include "Lapack.h"
#include "Comm_Handler.h"


////////////////////////////////////////////////////////
///     \class PCG
///     \brief This class represents the sequential
///            preconditioned conjugate gradient
///            algorithm. Note that the implementation
///            is quite similar to the CG
///            implemenation in the ParaSails package.
////////////////////////////////////////////////////////
template <class T_Field>
class PCG
{
    public:
        /////////////////////////////////////////////////////////
        ///     \brief  Constructor
        ///
        ///     \param env_handler_ Environment handler
        ///     \param comm_handler_ Pointer to communication
        ///                          handler
        ///     \param mtx_ Pointer to system matrix
        ///     \param tol_ Tolerance for solver
        ///     \param maxit_ Maximum number of iterations for
        ///                   solver
        ///     \param para_max_levels_ Number of levels in
        ///                             ParaSails preconditioner
        ///     \param para_threshold_ Threshold in ParaSails
        ///                            preconditioner
        ///     \param para_filter_ Filter parameter for ParaSails
        ///                         preconditioner
        /////////////////////////////////////////////////////////
        PCG<T_Field>( const ENV_Handler&     env_handler_,
                      Comm_Handler<T_Field>* comm_handler_,
                      Matrix<T_Field>*       mtx_,
                      const double           tol_,
                      const unsigned int     maxit_,
                      const int              para_max_levels_,
                      const double           para_threshold_,
                      const double           para_filter_);

        /// Destructor
        ~PCG<T_Field>( );

        /////////////////////////////////////////////////////////
        ///     \brief  Initializes solution vector with all 0.
        ///
        ///     \param dim Dimension of solution vector x.
        /////////////////////////////////////////////////////////
        void        Init_Solution( const int dim );

        /////////////////////////////////////////////////////////
        ///     \brief  Initializes temporary vector with all 0.
        ///
        ///     Dummy function in sequential version.
        ///
        ///     \param dim Dimension of temporary vector par_tmp.
        /////////////////////////////////////////////////////////
        void        Init_ParTmp( const int dim ){ }

        /////////////////////////////////////////////////////////
        ///     \brief  Creates right-hand side and builds it
        ///             according to chosen rhs_choice.
        ///
        ///     \param rhs_choice Choice of rhs.
        ///     \param dim Dimension of rhs vector
        ///                to be build with.
        /////////////////////////////////////////////////////////
        void        Create_RHS
                    (   const int rhs_choice,
                        const int dim );

        /////////////////////////////////////////////////////////
        ///     \brief  Builds rhs which is vector containing
        ///             random double values between 0.0 and 1.0.
        /////////////////////////////////////////////////////////
        void        RHS_Rand( );

        /////////////////////////////////////////////////////////
        ///     \brief  Builds rhs which is vector containing
        ///             all ones.
        /////////////////////////////////////////////////////////
        void        RHS_Ones( );

        /////////////////////////////////////////////////////////
        ///     \brief  Solves linear system with (P)CG.
        ///
        ///     Note that this implementation relies on the PCG
        ///     algorithm from ParaSails.
        /////////////////////////////////////////////////////////
        void        Solve_System( );

        /////////////////////////////////////////////////////////
        ///     \brief  Prints vector to std output
        ///
        ///     \param vec Vector to be printed.
        /////////////////////////////////////////////////////////
        void        Print_Vec(  const T_Field* vec );

        /////////////////////////////////////////////////////////
        ///     \brief  Setting preconditioner member.
        ///
        ///     \param precond_ Pointer to preconditioner to set
        ///                     member with.
        /////////////////////////////////////////////////////////
        void        Set_Precond(const Matrix<T_Field>* precond_);

        /////////////////////////////////////////////////////////
        ///     \brief  Writing system solution to file
        ///
        ///     \param file File name of file to write solution
        ///            to.
        /////////////////////////////////////////////////////////
        void        Solution_To_File(   const char* file );

        /////////////////////////////////////////////////////////
        ///     \brief  Dummy function for sequential environment
        /////////////////////////////////////////////////////////
        void        Solve_Precond_ParaSails(  ) {}

        // Members

        /// Solution vector
        T_Field*                  x;

        /// Right-hand side
        T_Field*                  b;

    private:

        /////////////////////////////////////////////////////////
        ///     \brief  Writing header to file.
        ///
        ///     \param file File pointer
        /////////////////////////////////////////////////////////
        void        Write_Header
                    (   FILE * file )                      const;

        /////////////////////////////////////////////////////////
        ///     \brief  Writing Line to file.
        ///
        ///     \param val Value to be written to file
        ///     \param f File pointer
        /////////////////////////////////////////////////////////
        void        Write_Line
                    (   const T_Field&  val,
                        FILE*           f )                const;

        /////////////////////////////////////////////////////////
        ///     \brief  Performing inner product / dot product.
        ///
        ///     \param dim Dimension of the vectors
        ///     \param x The x vector
        ///     \param y The y vector
        /////////////////////////////////////////////////////////
        double      Inner_Prod
                    (   int      dim,
                        T_Field* x,
                        T_Field* y );

        /////////////////////////////////////////////////////////
        ///     \brief  Copying a vector x->y.
        ///
        ///     \param dim Dimension of the vectors
        ///     \param x The x vector
        ///     \param x The y vector
        /////////////////////////////////////////////////////////
        void        Copy_Vector
                    (   int      dim,
                        T_Field* x,
                        T_Field* y);

        /////////////////////////////////////////////////////////
        ///     \brief  Scaling a vector x with alpha*x.
        ///
        ///     \param dim Dimension of the vectors
        ///     \param alpha constant alpha
        ///     \param x The x vector
        /////////////////////////////////////////////////////////
        void        Scale_Vector
                    (   int      dim,
                        double   alpha,
                        T_Field  *x);

        /////////////////////////////////////////////////////////
        ///     \brief  Performing alpha*x+y
        ///
        ///
        ///     \param dim Dimension of the vectors
        ///     \param alpha constant alpha
        ///     \param x The x vector
        ///     \param y The y vector
        /////////////////////////////////////////////////////////
        void        Axpy
                    (   int      dim,
                        double   alpha,
                        T_Field* x,
                        T_Field* y);

        /////////////////////////////////////////////////////////
        ///     \brief  Performing a Matrix-Vector product.
        ///
        ///     In a sequential environment loc_vec equals the
        ///     global vector.
        ///
        ///     \param res_vec Result vector
        ///     \param loc_vec Local vector
        /////////////////////////////////////////////////////////
        void        MV_Prod
                    (   T_Field*  loc_vec,
                        T_Field*& res_vec);

        /////////////////////////////////////////////////////////
        ///     \brief  Applying FSPAI preconditioner L to
        ///             L^T*loc_vec.
        ///
        ///     In a sequential environment loc_vec equals the
        ///     global vector.
        ///
        ///     \param res_vec Result vector
        ///     \param loc_vec Local vector
        /////////////////////////////////////////////////////////
        void        FSPAI_Apply(T_Field*& res_vec,
                                T_Field*  loc_vec);

        /////////////////////////////////////////////////////////
        ///     \brief  Applying FSPAI preconditioner L to
        ///             L*loc_vec.
        ///
        ///     In a sequential environment loc_vec equals the
        ///     global vector.
        ///
        ///     \param res_vec Result vector
        ///     \param loc_vec Local vector
        /////////////////////////////////////////////////////////
        void        FSPAI_ApplyTrans
                    (   T_Field*& res_vec,
                        T_Field*  loc_vec);

        /////////////////////////////////////////////////////////
        ///     \brief  Conjugating COMPLEX vector
        ///
        ///     \param vec Vector to be conjugated
        ///     \param dim Dimension of underlying model problem,
        ///                i.e., mtx->n
        /////////////////////////////////////////////////////////
        T_Field*    Conjugate_Complex_Vec
                    (   const T_Field*  vec1,
                        const int       dim);

        /////////////////////////////////////////////////////////
        ///     \brief  Initialize sumvec as it is of non POD
        ///             type.
        ///
        ///     \param vec Vector to be initialized
        ///     \param dim Dimension of underlying model problem,
        ///                i.e., mtx->n
        /////////////////////////////////////////////////////////
        void        Init_Vec
                    (   T_Field*& vec,
                        const int dim );

        // Private member

        /// Environment handler
        ENV_Handler               env_handler;

        /// Pointer to communication handler
        Comm_Handler<T_Field>*    comm_handler;

        /// Member pointer to matrix chunk
        const Matrix<T_Field>*    mtx;

        /// Member pointer to preconditioner chunk
        const Matrix<T_Field>*    precond;

        /// Solver tolerance
        const double              tol;

        /// Maximum number of iterations for solver
        const unsigned int        maxit;

        // Provided due to parallel implementation - not used in
        // sequential environment
        /// Maximum number of levels in ParaSails preconditioner
        int                       para_max_levels;

        /// Threshold parameter in ParaSails preconditioner
        double                    para_threshold;

        /// Filter parameter in ParaSails preconditioner
        double                    para_filter;
};

#include "PCG.imp.h"

#endif /* PCG_H_ */

#endif
