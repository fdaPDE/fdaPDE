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
Matrix<T_Field>::Matrix
(   const ENV_Handler& env_handler,
    const int          nbr_cols_,
    const int          nbr_rows_ )
{
    env_handler.Get_Environment_Params( num_procs, my_id );
    my_nbr_cols         = 0;
    my_start_idx        = 0;
    next_col            = 0;
    max_nnz             = 0;
    my_nnz              = 0;
    n                   = nbr_cols_;
    m                   = nbr_rows_;
    all_nbr_cols        = new int[num_procs];
    start_indices       = new int[num_procs];
    len_all_cols        = new int[nbr_cols_];
    pe                  = new int[nbr_cols_];
    remote_col_idcs_buf = NULL;
    remote_col_buf      = NULL;
}



template <class T_Field>
Matrix<T_Field>::~Matrix()
{
    if (c_lines)
        delete c_lines;
    if (all_nbr_cols)           delete [] all_nbr_cols;
    if (start_indices)          delete [] start_indices;
    if (pe)                     delete [] pe;
    if (len_all_cols)           delete [] len_all_cols;
    if (remote_col_buf)         delete [] remote_col_buf;
    if (remote_col_idcs_buf)    delete [] remote_col_idcs_buf;
}



template <class T_Field>  int
Matrix<T_Field>::Count_NNZ
(   ENV_Handler& env_handler ) const
{
    int nnz = 0,
        my_nnz = 0;

    for ( int i = 0; i < my_nbr_cols; i++ )
        my_nnz += c_lines->len_cols[i];

    env_handler.Sum_NNZ( my_nnz, nnz );

    return nnz;
}



template <class T_Field>  void
Matrix<T_Field>::Init_Preconditioner
(    Matrix<T_Field>*        mtx,
     ENV_Handler& env_handler )
{
    env_handler.Dist_Local_Chunks( mtx->my_nbr_cols, all_nbr_cols );

    // Filling start indices
    start_indices[0] = 0;
    for ( int pe_loc = 1; pe_loc < num_procs; pe_loc++ )
        start_indices[pe_loc] =  start_indices[pe_loc-1] +
                                 all_nbr_cols[pe_loc-1];

    my_nbr_cols = all_nbr_cols[my_id];
    my_start_idx = start_indices[my_id];

    c_lines = new Compressed_Lines<T_Field>( my_nbr_cols );
    memcpy( pe, mtx->pe, n * sizeof( int ) );
}



template <class T_Field>  void
Matrix<T_Field>::To_Pattern
(   Pattern*                P,
    const ENV_Handler& env_handler )
{
    Index_Set   *i_set   = NULL;
    int          clen,
                 start_index,
                 newlen,
                 *cbuf = NULL;

    env_handler.Get_Environment_Params( P->num_procs, P->my_id );

    P->all_nbr_cols  = new int[P->num_procs];
    P->start_indices = new int[P->num_procs];
    P->pe            = new int[n];
    P->j_sets        = new Index_Set*[my_nbr_cols];
    for (int i = 0; i < my_nbr_cols; i++)
        P->j_sets[i] = NULL;
    P->len           = n;

    env_handler.Dist_Local_Chunks( my_nbr_cols, P->all_nbr_cols );

    // Filling start indices
    P->start_indices[0] = 0;
    for (int pe_loc = 1; pe_loc < P->num_procs; pe_loc++)
        P->start_indices[pe_loc] = P->start_indices[pe_loc-1] +
                               P->all_nbr_cols[pe_loc-1];

    // filling pe array
    memset(P->pe, 0, n * sizeof(int));
    for (int pe_loc = 0; pe_loc < P->num_procs; pe_loc++)
    {
        start_index = P->start_indices[pe_loc];
        for (int i = 0; i < P->all_nbr_cols[pe_loc]; i++)
            P->pe[start_index + i] = pe_loc;
    }

    P->my_nbr_cols = P->all_nbr_cols[P->my_id];
    P->my_start_idx = P->start_indices[P->my_id];

    // Filling pattern data structure with indices
    // Copying only lower triangular part of the
    // system to pattern.
    for (int col = 0; col < my_nbr_cols; col++)
    {
        clen = c_lines->len_cols[col];
        cbuf = c_lines->col_idcs[col];

        newlen = 0;
        for (int r = 0; r < clen; r++)
        {
            if(cbuf[clen-r-1] < (my_start_idx + col))
                break;
            newlen++;
        }

        // Checking whether triangular pattern
        // has no zero column
        if ( newlen == 0 )
            throw Pe_Exception(
                "\n\n\tERROR:  Lower triangular part of system\n"
                "\t\thas zero column. Pattern will not be generated\n"
                "\t\tfrom this input. Please verify that lower\n"
                "\t\ttriangular part of system is not singular.\n");

        i_set = new Index_Set(newlen);
        cbuf += ( clen-newlen );
        memcpy(i_set->idcs,
               cbuf,
               newlen * sizeof(int));
        P->j_sets[col] = i_set;
    }

    // Get the maximum number of nnz per
    // column/row of all pes
    env_handler.Dist_Max( max_nnz, P->max_nnz );
}



template<class T_Field> void
Matrix<T_Field>::To_File
( ENV_Handler& env_handler,
 const char*   file ) const
{
    int         nnz, ierr;
    FILE        *f;
    const char  *mm_string = "%%MatrixMarket";
    char        fullname[1024],
                cat_cmd[1024],
                rm_cmd[1024];
    Timer       timer = Timer();

    // Start time measurement
    timer.Start( env_handler );

    if (num_procs > 1)
    {
        sprintf(fullname,   "%s_tmp%5.5d",      file, my_id);
        sprintf(cat_cmd,    "cat %s_tmp* > %s", file, file);
        sprintf(rm_cmd,     "rm -f %s_tmp*",    file);
    }
    else
        sprintf(fullname, "%s", file);

    if ( !( f = fopen( fullname,"w" ) ) )
        throw std::runtime_error(
            "\n\tERROR:  Failed writing preconditioner to file "
            + std::string(fullname) + "\n"
            "\n\t\tCheck your access rights!\n");

    nnz = Count_NNZ( env_handler );

    // write Matrix-Market header
    if (my_id == 0)
    {
        fprintf(f, "%s ", mm_string);
        Write_Header(f);
        fprintf(f, "%d %d %d\n", m, n, nnz);
        fflush(f);
    }

    for (int j = 0; j < my_nbr_cols; j++)
        for (int i = 0; i < c_lines->len_cols[j]; i++)
            Write_Line(j, i, f);

    fflush(f);
    fclose(f);

    env_handler.Barrier();

    if (num_procs > 1)
        if (my_id == 0)
            ierr = system(cat_cmd);

    env_handler.Barrier();

    if (num_procs > 1)
        if (my_id == 0)
            ierr = system(rm_cmd);

    // Stop and report time measurement
    timer.Stop( env_handler );
    timer.Report( env_handler );
}
