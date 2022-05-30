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


template <class T_Field>
Fspai_Sub<T_Field>::Fspai_Sub
(  )
{
    bitvec          = NULL;
    reset_vec       = NULL;
    shadow          = NULL;
    diag_elements   = NULL;
    sumvec          = NULL;
}



template <class T_Field>
Fspai_Sub<T_Field>::Fspai_Sub
(   Comm_Handler<T_Field>* comm_handler_,
    const int              dim)
{
    comm_handler = comm_handler_;

    // reset vector must have mtx->n positions
    // because more positions can be stored at once.
    int bitvec_size =  dim;
    reset_vec   = new unsigned int[bitvec_size];

    // Only these positions are needed within
    // the bitvec. 5: log n bits of int
    bitvec_size = (bitvec_size >> 5) + 1;
    bitvec      = new unsigned int[bitvec_size];
    memset(bitvec, 0, bitvec_size * sizeof(unsigned int));

    shadow = new Index_Set(dim);
    shadow->len = 0;
    reset_len = 0;

    sumvec = new T_Field[dim];
    Init_Vec(sumvec, dim);

    diag_elements = new DIAG_ELEMENT[dim];
    Init_Diags(dim);
}



template <class T_Field>
Fspai_Sub<T_Field>::~Fspai_Sub()
{
    if (bitvec)        delete [] bitvec;
    if (reset_vec)     delete [] reset_vec;
    if (shadow)        delete shadow;
    if (diag_elements) delete [] diag_elements;
    if (sumvec)        delete [] sumvec;
}



template <class T_Field>
Index_Set* Fspai_Sub<T_Field>::Compute_Jktilde
(   const Index_Set* Jk,
    const int        col)
{
    Index_Set* Jk_tilde = new Index_Set(Jk->len);
    int        j,
               jtidx = 0;

    for(int jidx = 0; jidx < Jk->len; jidx++)
    {
        j = Jk->idcs[jidx];
        if (j != col)
            Jk_tilde->idcs[jtidx++] = j;
    }
    Jk_tilde->len = jtidx;
    return Jk_tilde;
}



template <class T_Field>
T_Field* Fspai_Sub<T_Field>::Reduce_Mtx
(   const Index_Set* Jk_tilde )
{
    int         jk_len = Jk_tilde->len,
                *col_idcs_buf = NULL,
                col_len = 0,
                jk_idx,
                nnz_idx,
                r_idx;
    T_Field     *red_mtx = new T_Field[jk_len * jk_len],
                *col_buf = NULL;

    memset(red_mtx, 0, jk_len * jk_len * sizeof(T_Field));

    for(int j = 0; j < jk_len; j++)
    {
        nnz_idx = 0;
        comm_handler->Get_Mtx_Col(
                Jk_tilde->idcs[j],
                col_len,
                col_idcs_buf,
                col_buf);

        for(int c = 0; c < col_len; c++)
        {
            if (nnz_idx >= jk_len) break;
            r_idx = col_idcs_buf[c];
            jk_idx = Jk_tilde->idcs[nnz_idx];
            if ( jk_idx < r_idx )
            {
                nnz_idx++; c--;
                continue;
            }
            if ( r_idx == jk_idx )
            {
                red_mtx[nnz_idx+j*jk_len] = col_buf[c];
                nnz_idx++;
            }
        }
    }
    return red_mtx;
}



template <class T_Field> T_Field*
Fspai_Sub<T_Field>::Extract_Reduced_Column
(   const Index_Set* Jk_tilde,
    const int        col,
    T_Field&         A_kk)
{
    int         jk_len = Jk_tilde->len,
                *col_idcs_buf = NULL,
                col_len = 0,
                nnz_idx = 0,
                col_idx,
                jk_idx;
    T_Field     *red_col = new T_Field[jk_len],
                *col_buf = NULL,
                diag;

    memset(red_col, 0, jk_len * sizeof(T_Field));

    comm_handler->Get_Mtx_Col(
                  col,
                  col_len,
                  col_idcs_buf,
                  col_buf);

    for (int jk = 0; jk < jk_len; jk++)
    {
        jk_idx = Jk_tilde->idcs[jk];
        for(int c = 0; c < col_len; c++)
        {
            col_idx = col_idcs_buf[c];
            if(col_idx == col)
            {
                diag = col_buf[c];
                diag_elements[col].set = true;
                diag_elements[col].diag_el = diag;
                A_kk = diag;
            }
            if(jk_idx < col_idcs_buf[c])
                break;
            if(jk_idx == col_idcs_buf[c])
                red_col[nnz_idx] = col_buf[c];
        }
        nnz_idx++;
    }
    return red_col;
}



template <class T_Field> T_Field
Fspai_Sub<T_Field>::Extract_A_kk
(   const int   col)
{
    int         *col_idcs_buf = NULL,
                col_len = 0;
    T_Field     *col_buf = NULL,
                diag,
                A_kk;

    // Suppressing compiler warnings and
    // initializing diagonal element.
    // It will be set for sure in the
    // following for-loop.
    Init_Diagonal(A_kk);

    comm_handler->Get_Mtx_Col(
                  col,
                  col_len,
                  col_idcs_buf,
                  col_buf);

    for(int c = 0; c < col_len; c++)
        if(col_idcs_buf[c] == col)
        {
            diag = col_buf[c];
            diag_elements[col].set = true;
            diag_elements[col].diag_el = diag;
            A_kk = diag;
        }
    return A_kk;
}



template <class T_Field>  void
Fspai_Sub<T_Field>::Compute_Sum_Shadow
(   const int           col,
    const Index_Set*    Jk,
    const T_Field* Lk)
{
    int           pos,
                  bit,
                  r_idx,
                  jk_idx,
                  col_len       = 0,
                  shadow_len    = 0,
                  reset_len     = 0,
                  *col_idcs_buf = NULL;
    T_Field       *col_buf = NULL;

    for (int j = 0; j < Jk->len; j++)
    {
        jk_idx = Jk->idcs[j];
        comm_handler->Get_Mtx_Col(
                      jk_idx,
                      col_len,
                      col_idcs_buf,
                      col_buf);

        for (int idx = 0; idx < col_len; idx++)
        {
            //Get the row index value
            r_idx = col_idcs_buf[idx];
            if(r_idx == jk_idx)
            {
                diag_elements[r_idx].set = true;
                diag_elements[r_idx].diag_el = col_buf[idx];
            }
            if(r_idx < col) continue;
            pos = r_idx >> 5;       //log n bits of int
            bit = r_idx % 32;       //n bits of int
            if (! Bit_Test( bitvec[pos], bit ))
            {
                Set_Bit( &bitvec[pos], bit );
                shadow->idcs[shadow_len++] = r_idx;
                reset_vec[reset_len++] = pos;
            }
            sumvec[r_idx] = sumvec[r_idx] + (col_buf[idx]*Lk[j]);
        }
    }
    // Sort I
    std::sort(shadow->idcs, shadow->idcs + shadow_len);
    shadow->len = shadow_len;

    // Avoiding memset of whole bitvec because
    // too expensive for large matrices. Just
    // reset the used positions.
    Reset_bitvec(reset_len);
}



template <class T_Field>  int
Fspai_Sub<T_Field>::Bit_Test
(   unsigned int bv,
    int          bit)
{
    return ( bv & (1 << bit) );
}



template <class T_Field>  void
Fspai_Sub<T_Field>::Set_Bit
(   unsigned int *bv,
    int           bit)
{
    *bv |= (1 << bit);
}



template<class T_Field> void
Fspai_Sub<T_Field>::Reset_bitvec
(   int            reset_len)
{
    for (int i = 0; i < reset_len; i++)
        bitvec[reset_vec[i]] = 0;
}



template <class T_Field> Index_Set*
Fspai_Sub<T_Field>::Pattern_Update
(   const T_Field*      Lk,
    Index_Set*          Jk,
    const int           col,
    const int           dim,
    const double        epsilon_param,
    const int           max_idcs_param,
    const bool          use_mean_param,
    bool&               update)
{
    TAU_IDX*    taus                = NULL;
    double      mean_val            = 0.0;
    int         nbr_taus            = 0,
                nbr_updated_idcs    = 0;
    Index_Set   *Jk_new             = NULL,
                *Jk_cand            = NULL;

    // Compute all tau values for j candidates
    taus = Compute_Taus(col, dim, Lk, Jk, mean_val,
                        use_mean_param, nbr_taus);

    if ( (nbr_taus == 0) || (taus[0].tau < epsilon_param) )
    {
        update = false;
        if (taus) delete [] taus;
        return NULL;
    }

    Jk_cand = new Index_Set(nbr_taus);

    // Augment with new indices
    // nbr_taus != 0 -> mean_val > 0;
    for (int el = 0;
            (el < max_idcs_param) &&
            (el < nbr_taus) &&
            (taus[el].tau >= mean_val);
         el++)
    {
        Jk_cand->idcs[el] = taus[el].idx;
        nbr_updated_idcs++;
    }
    Jk_cand->len = nbr_updated_idcs;

    if ( nbr_updated_idcs == 0)
    {
        update = false;
        if (taus) delete [] taus;
        if (Jk_cand) delete Jk_cand;
        return NULL;
    }

    //Sorting the new augmented Index_Set J
    std::sort(Jk_cand->idcs, Jk_cand->idcs + Jk_cand->len);

    // Union Jk, Jk_augmented
    Jk_new = Jk->Set_Union(Jk, Jk_cand);

    if (Jk)         delete Jk;
    if (taus)       delete [] taus;
    if (Jk_cand)    delete Jk_cand;
    return Jk_new;
}



template <class T_Field> TAU_IDX*
Fspai_Sub<T_Field>::Compute_Taus
(   const int         col,
    const int         dim,
    const T_Field*    Lk,
    const Index_Set*  Jk,
    double&           mean_val,
    const bool        use_mean_param,
    int&              nbr_taus)
{
    int               j,
                      jhat_len;
    T_Field           diag,
                      sum;
    double            tau;

    Compute_Sum_Shadow(col, Jk, Lk);
    jhat_len = shadow->len;

    // maximum number of taus == lenght(JHat);
    TAU_IDX* taus = new TAU_IDX[jhat_len];

    // For all j > k and j element from Jhat directly computing
    // A(j,Jk)*Lk(Jk).
    for(int l = 1; l < jhat_len; l++)
    {
        j = shadow->idcs[l];

        // j is larger than k as exluded in Get_Shadows
        // Check whether j is in Jk. If so go to next j
        if(Jk->Has_Idx(j))
        {
            Reset_Sumvec(j);
            continue;
        }

        // Extract diagonal element
        if(!diag_elements[j].set)
            diag = Extract_A_kk(j);
        else
            diag = diag_elements[j].diag_el;

        // Extract sum of A(j,Jk)*Lk(Jk).
        sum = sumvec[j];

        // Computing the dot product A(j,Jk)*Lk(Jk)
        tau = Compute_Tau(sum, diag);

        taus[nbr_taus].tau = tau;
        taus[nbr_taus++].idx = j;
        mean_val += tau;

        // Reset j-th position used sumvec
        Reset_Sumvec(j);
    }
    // Reset 0-th position in sumvec because
    // was skipped in previous for-loop
    Reset_Sumvec(shadow->idcs[0]);
    //Computing mean value of rho indices
    if (use_mean_param)
    {
        if(nbr_taus > 0)
            mean_val /= nbr_taus;
    } else  mean_val = 0.0;

    //Descending sort of taus
    std::sort(taus, taus + nbr_taus, TAU_Comparator());
    return taus;
}
