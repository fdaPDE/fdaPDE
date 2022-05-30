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

#ifndef FSPAI_SUB_H_
#define FSPAI_SUB_H_

// file includings
#include "Matrix.h"
#include "Comm_Handler.h"


//////////////////////////////////////////
///     \brief Representing a pair of
///            index and double
//////////////////////////////////////////
struct TAU_IDX
{
    int     idx;
    double  tau;
};


//////////////////////////////////////////
///     \brief How to sort the TAU_IDX
///            within an array - this is
///            ascending order
//////////////////////////////////////////
struct TAU_Comparator
{
    bool operator()(const TAU_IDX& a, const TAU_IDX& b)
    {
        return a.tau > b.tau;
    }
};


///////////////////////////////////////////////////
///     \class Fspai_Sub
///     \brief Implementing sub routines used
///            in all FSPAI algorithms
///////////////////////////////////////////////////
template <class T_Field>
class Fspai_Sub
{
    public:

        /// Empty constructor
        Fspai_Sub<T_Field>();

        /// Constructor with environment handler
        Fspai_Sub<T_Field>(Comm_Handler<T_Field>* comm_handler_, const int dim);

        /// Destructor
        ~Fspai_Sub<T_Field>();

        /////////////////////////////////////////////////////////
        ///     \brief  Reducing matrix to submatrix due to
        ///             computed index set J_k.
        ///
        ///     \param Jk Index set J_k from FSPAI algorithm
        ///     \return The reduced submatrix
        /////////////////////////////////////////////////////////
        T_Field*    Reduce_Mtx
                    (   const Index_Set* Jk);

        /////////////////////////////////////////////////////////
        ///     \brief  Extract column from matrix and reduce
        ///             it according to index set Jk_tilde.
        ///
        ///     Note: At the same time diagonal element is
        ///           extracted which is used later.
        ///
        ///     \param Jk_tilde Reducing index set
        ///     \param col column to be extracted from matrix
        ///     \param A_kk Diagonal element to be extracted
        ///     \return The extracted column from the matrix
        /////////////////////////////////////////////////////////
        T_Field*    Extract_Reduced_Column
                    (   const Index_Set* Jk_tilde,
                        const int        col,
                        T_Field&         A_kk);

        /////////////////////////////////////////////////////////
        ///     \brief  Computing index set Jk_tilde
        ///
        ///     \param Jk Index set to be deleted possible
        ///               diagonal element from.
        ///     \param col Column which is currently
        ///                computed.
        ///     \return The index set Jk_tilde.
        /////////////////////////////////////////////////////////
        Index_Set*  Compute_Jktilde
                    (   const Index_Set* Jk,
                        const int        col);

        /////////////////////////////////////////////////////////
        ///     \brief  Printing reduced submatrix
        ///
        ///     \param red_mtx Reduced matrix to be printed.
        ///     \param n Number of columns of submatrix.
        ///     \param m Number of rows of submatrix.
        /////////////////////////////////////////////////////////
        void        Print_Red_Mtx
                    (   const T_Field* red_mtx,
                        const int      n,
                        const int      m);

        /////////////////////////////////////////////////////////
        ///     \brief  Solve a symmetric/hermitian positive
        ///             definite linear system
        ///
        ///     \param red_mtx According to pattern set Jk_tilde
        ///                    reduced submatrix
        ///                    A(Jk_tilde,Jk_tilde)
        ///     \param red_col According to pattern set Jk_tilde
        ///                    reduced column A(Jk_tilde,k)
        ///     \param n The dimension of the linear system
        ///     \param uplo 'U':  A is upper triangular;
        ///                 'L':  A is lower triangular.
        ///     \param nrhs The number of columns of rhs
        ///     \param info Function return value
        /////////////////////////////////////////////////////////
        void        Solve_HPD_System
                    (   T_Field*    red_mtx,
                        T_Field*    red_col,
                        const int&        n,
                        const char* uplo,
                        int&        nrhs,
                        int&        info);

        /////////////////////////////////////////////////////////
        ///     \brief  Computes a dot prodcut between
        ///             two vectors
        ///
        ///     \param vec1 Vector 1
        ///     \param vec2 Vector 2
        ///     \param dim Dimension of the vectors.
        ///     \return Dot product of the two vectors
        /////////////////////////////////////////////////////////
        T_Field     Dot_Product
                    (   const T_Field*  vec1,
                        const T_Field*  vec2,
                        const int       dim);

        /////////////////////////////////////////////////////////
        ///     \brief  Computes the square root of a real
        ///             and complex number
        ///
        ///     \param number The number for which the square
        ///                   root is to be compupted.
        ///     \return Square root of number
        /////////////////////////////////////////////////////////
        T_Field     Field_Sqrt
                    (   const T_Field number);

        /////////////////////////////////////////////////////////
        ///     \brief  Get the diagonal element from A
        ///
        ///     \param col Column k from which the diagonal
        ///                element of A is to be returned
        ///     \return Diagonal element A_kk from A.
        /////////////////////////////////////////////////////////
        T_Field     Extract_A_kk
                    (   const int   col);

        /////////////////////////////////////////////////////////
        ///     \brief  Computes the conjugated complex vector.
        ///
        ///     \param vec1 The complex vector to be conjugated
        ///     \param dim The dimension of the vector
        ///     \return The conjugated complex vector
        /////////////////////////////////////////////////////////
        T_Field*    Conjugate_Complex_Vec
                    (   const T_Field*  vec1,
                        const int       dim);

        /////////////////////////////////////////////////////////
        ///     \brief  Updates current sparsity pattern with
        ///             new promising values.
        ///
        ///     \param Lk Current solution of kth column of L
        ///     \param Jk Current pattern of kth column of L
        ///     \param col Current column k of L
        ///     \param dim Dimension of matrix mtx
        ///     \param epsilon_param epsilon tolerance for each
        ///                          column.
        ///     \param max_idcs_param Maximum number of indices to
        ///                           be augmented per step
        ///                     computation.
        ///     \param use_mean_param Whether mean_value has to
        ///                           be used or not.
        ///     \param update flag whether pattern was updated or
        ///                   not.
        ///     \return New augmented index set
        /////////////////////////////////////////////////////////
        Index_Set*  Pattern_Update
                    (   const T_Field*  Lk,
                        Index_Set*      Jk,
                        const int       col,
                        const int       dim,
                        const double    epsilon_param,
                        const int       max_idcs_param,
                        const bool      use_mean_param,
                        bool&           update);

        /////////////////////////////////////////////////////////
        ///     \brief  Computes tau values for pattern update.
        ///
        ///     \param col Current column k of L
        ///     \param dim Dimension of matrix mtx
        ///     \param Lk Current solution of kth column of L
        ///     \param Jk Current pattern of kth column of L
        ///     \param mean_val Mean value to be set during tau
        ///                     computation.
        ///     \param use_mean_param Whether mean_value has to
        ///                           be used or not.
        ///     \param nbr_taus Number of tau values which were
        ///                     computed (which are > 0.0)
        ///     \return All tau values for current update
        /////////////////////////////////////////////////////////
        TAU_IDX*    Compute_Taus
                    (   const int           col,
                        const int           dim,
                        const T_Field*      Lk,
                        const Index_Set*    Jk,
                        double&             mean_val,
                        const bool          use_mean_param,
                        int&                nbr_taus);

        /////////////////////////////////////////////////////////
        ///     \brief  Resetting the bitvec vector to initial state.
        ///
        ///     For the computation of the shadow of the index set Jk,
        ///     the whole bitvector must be in initial state containing
        ///     only 0-elements. The algorithm Get_Shadow operates
        ///     only on specific positions of bitvec so a
        ///     complete memset is not necessary and to expensive for
        ///     large matrices beyond a size of 10^5.
        ///     The reset vector stores all used position in bitvec and
        ///     using the reset length the bitvec can be resetted easily.
        ///
        ///     \param reset_len The number of used position in
        ///                      current iteration.
        /////////////////////////////////////////////////////////
        void        Reset_bitvec( int reset_len );

        /////////////////////////////////////////////////////////
        ///     \brief  Initialize vector as it may be of non POD
        ///             type.
        ///
        ///     \param vec Vector to be initialized
        ///     \param dim Dimension of underlying model problem,
        ///                i.e., mtx->n
        /////////////////////////////////////////////////////////
        void        Init_Vec
                    (   T_Field*& vec,
                        const int dim );

    private:

        /////////////////////////////////////////////////////////
        ///     \brief  Computes the shadow on mtx for given
        ///             index set and sums up all values of
        ///             A(j,Jk)*Lk(Jk).
        ///
        ///     Computing subshadow from column k+1 because only
        ///     indices j > k has to be considered in shadow.
        ///     Furthermore, during the shadow computation all
        ///     values of A(j,Jk)*Lk(Jk) for j > k are stored
        ///     within array sumvec.
        ///
        ///     \param col The current column to begin shadow
        ///                computation from
        ///     \param Jk The index set to be computed the shadow
        ///               from
        ///     \param Lk Current approximate solution of kth col.
        /////////////////////////////////////////////////////////
        void        Compute_Sum_Shadow
                    (   const int        col,
                        const Index_Set* Jk,
                        const T_Field*   Lk);

        /////////////////////////////////////////////////////////
        ///     \brief  Testing whether an integer is subset in
        ///             at this bitvector position
        ///
        ///     \param bv The "bitset" integer vector
        ///     \param bit The integer to be tested
        ///     \return Whether the integer is included or not
        /////////////////////////////////////////////////////////
        int         Bit_Test
                    (   unsigned int bv,
                        int          bit );

        /////////////////////////////////////////////////////////
        ///     \brief  Setting bits of an integer into this
        ///             bitvector position.
        ///
        ///     \param bv The "bitset" integer vector
        ///     \param bit The integer to be set
        /////////////////////////////////////////////////////////
        void        Set_Bit
                    (   unsigned int *bv,
                        int          bit );

        /////////////////////////////////////////////////////////
        ///     \brief  Reset sumvec array at specified position.
        ///
        ///     \param r_idx Position in sumvec which has to be
        ///                  resetted.
        /////////////////////////////////////////////////////////
        void        Reset_Sumvec( const int r_idx );

        /////////////////////////////////////////////////////////
        ///     \brief  Computes tau value for given index j
        ///             during pattern updates.
        ///             tau_j = (|A(j,Jk)*Lk(Jk)|)^2/A_jj
        ///
        ///     This template-specific method is inline as
        ///     invoked heavily often.
        ///
        ///     \param sum The sum of A(j,Jk)*Lk(Jk) for j-th
        ///                index to compute tau of
        ///     \param diag The diagonal element of j-th column/row
        ///     \return The tau value for j-th index.
        /////////////////////////////////////////////////////////
        double      Compute_Tau
                    (   double& sum,
                        double& diag )
                    {   return ((sum * sum) / diag);  }

        /////////////////////////////////////////////////////////
        ///     \brief  Computes tau value for given index j
        ///             during pattern updates.
        ///             tau_j = (|A(j,Jk)*Lk(Jk)|)^2/A_jj
        ///
        ///     This template-specific method is inline as
        ///     invoked heavily often.
        ///
        ///     \param sum The sum of A(j,Jk)*Lk(Jk) for j-th
        ///                index to compute tau of
        ///     \param diag The diagonal element of j-th column/row
        ///     \return The tau value for j-th index.
        /////////////////////////////////////////////////////////
        double      Compute_Tau
                    (   COMPLEX& sum,
                        COMPLEX& diag )
                    {   return ((sum.real*sum.real +
                                 sum.imag*sum.imag) / diag.real);  }

        /////////////////////////////////////////////////////////
        ///     \brief  Print array of locally cached
        ///             diagonal elements.
        ///
        ///     \param dim Dimension of underlying model problem,
        ///                i.e., mtx->n
        /////////////////////////////////////////////////////////
        void        Print_Diags( const int dim );

        /////////////////////////////////////////////////////////
        ///     \brief  Initialize diagonal elements as it is of
        ///             non POD type.
        ///
        ///     \param dim Dimension of underlying model problem,
        ///                i.e., mtx->n
        /////////////////////////////////////////////////////////
        void        Init_Diags( const int dim );

        /////////////////////////////////////////////////////////
        ///     \brief  Initializes diagonal element. This method
        ///             is only provided to suppress compiler
        ///             warnings.
        ///
        ///     \param A_kk Diagonal element to be initialized
        /////////////////////////////////////////////////////////
        void        Init_Diagonal( T_Field& A_kk );

        /// Interface to communication handler
        Comm_Handler<T_Field>*  comm_handler;

        /// Bitvector holding the shadow values during shadow computation
        unsigned int*           bitvec;

        /// Reset vector used to reset positions from shadow computation
        unsigned int*           reset_vec;

        /// Holding current length of reset_vec
        int                     reset_len;

        /// Index set holding indices of shadow
        Index_Set*              shadow;

        /// Array holding sums of A(j,Jk)*Lk(Jk) for all j
        T_Field*                sumvec;

        /// \brief Local cache storing diagonal elements of A. It is filled by
        /// when remote columns are requested.
        struct  DIAG_ELEMENT {
            bool    set;
            T_Field diag_el;
        }*                      diag_elements;
};

#include "Fspai_Sub.imp.h"

#endif /* FSPAI_SUB_H_ */
