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


template <class T_Field>  void
Fspai_Unrestrained<T_Field>::Fspai_Column
( const int col )
{
    T_Field         *red_mtx    = NULL,
                    *red_col    = NULL,
                    *Lk         = NULL,
                     A_kk;
    Index_Set       *Jk_tilde   = NULL,
                    *Jk_updated = NULL;
    bool             update     = true;

    // J will be filled with the
    // new J values from pattern
    // and will be updated each augmenting
    // step.
    Jk = new Index_Set(mtx->n);

    // Extract pattern from pattern matrix
    // Create the Index Set J
    comm_handler->Get_P_Col(col, Jk);

    // Pattern Updates if updates_param > 0
    for (int step = -1; step < updates_param; step++)
    {
        // Compute Jk_tilde
        Jk_tilde = fspai_sub->Compute_Jktilde(Jk, col);

        // Drive on communication between pes
        comm_handler->Communicate();

        Lk = new T_Field[Jk->len];

        // Skip all computations if Jktilde is empty. In this case we
        // only have to compute the diagonal entry L_kk
        if(Jk_tilde->len == 0)
        {
            A_kk = fspai_sub->Extract_A_kk(col);
            Lk[0] = 1.0 / fspai_sub->Field_Sqrt(A_kk);
        }
        else
        {

            // Reduce Matrix to small size due to column Pattern Jk
            // A(Jk_tilde,Jk_tilde)
            red_mtx = fspai_sub->Reduce_Mtx(Jk_tilde);

            // Extract reduced column from Matrix A(Jk_tilde,k) AND
            // Extract diagonal element from A, i.e., A_kk
            red_col = fspai_sub->Extract_Reduced_Column(Jk_tilde, col, A_kk);

            // Computing the approximate inverse for column k, i.e.,
            // L_k(Jk_tilde) and L_kk
            Approximate_Inverse( Lk, red_mtx, red_col, Jk_tilde->len, col, A_kk);
        }
        // If no update steps requested finish computation
        if ( step == (updates_param - 1) )  break;

        // Augment pattern column
        Jk_updated = fspai_sub->Pattern_Update(Lk, Jk, col, mtx->n, epsilon_param,
                     max_idcs_param, use_mean_param, update);

        if(!update) break; // taus equal zero or epsilon tolerance reached
        else        Jk = Jk_updated;

        if (red_mtx)  delete [] red_mtx;
        if (red_col)  delete [] red_col;
        if (Lk)       delete [] Lk;
        if (Jk_tilde) delete Jk_tilde;
    }

    //Inserting solution column L_k into preconditioner L_M
    comm_handler->Insert_Precond_Solution(Lk, col, Jk_tilde, Jk);

    if (red_mtx)    delete [] red_mtx;
    if (red_col)    delete [] red_col;
    if (Lk)         delete [] Lk;
    if (Jk_tilde)   delete Jk_tilde;
    if (Jk)         delete Jk;
}



template <class T_Field>  void
Fspai_Unrestrained<T_Field>::Approximate_Inverse
(   T_Field*&     Lk,
    T_Field*      red_mtx,
    T_Field*      red_col,
    int           dim,
    const int     col,
    const T_Field A_kk)
{
    int           info = 0,
                  nrhs = 1;
    T_Field       tmp_Lkk,
                  L_kk,
                  Lk_tilde[dim];
    const char    *uplo = "L";

    // copy A(Jk_tilde,k) to Lk_solution because
    // Solve_HPD_System will overwrite red_col
    memcpy(Lk_tilde, red_col, dim * sizeof(T_Field));

    fspai_sub->Solve_HPD_System(red_mtx, red_col, dim, uplo, nrhs, info);
    // red_col contains solution of y_k = A(Jk_tilde,Jk_tilde)^{-1}A(Jk_tilde,k)

    // Computing L_kk = 1/(sqrt(A_kk-A(Jk_tilde,k)^T*y_k))
    L_kk = 1.0 / (fspai_sub->Field_Sqrt(A_kk -
                  fspai_sub->Dot_Product(Lk_tilde, red_col, dim)));

    Lk[0] = L_kk;

    // Computing L_k(Jk_tilde) = -L_kk*y_k;
    tmp_Lkk = L_kk * (-1.0);
    for(int i = 0; i < dim; i++)
        Lk[i+1] = red_col[i]*tmp_Lkk;
}

