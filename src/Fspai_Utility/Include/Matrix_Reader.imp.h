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


template<class T_Field, class T_Format>
Matrix_Reader<T_Field, T_Format>::Matrix_Reader
( )
{
    file_cols   = 0;
    file_rows   = 0;
    file_nnz    = 0;
    start_idx   =  0;
    nbr_rows    = 0;
    nbr_cols    = 0;
    nnz         = 0;
    nnz_cols    = 0;
    split_idx   = -1;
    split_pe    = -1;
    my_nbr_cols = 0;
    my_nbr_rows = 0;
    hermitian   = false;
    cols        = NULL;
    diagonal_elements = NULL;
}



template<class T_Field, class T_Format>
Matrix_Reader<T_Field, T_Format>::~Matrix_Reader
( )
{
    if (cols)              delete [] cols;
    if (diagonal_elements) delete [] diagonal_elements;
}



template<class T_Field, class T_Format> void
Matrix_Reader<T_Field, T_Format>::To_Memory
(   Matrix<T_Field>*&   mtx,
    char*               matrix_file,
    const MMio          mmio,
    ENV_Handler&        env_handler)
{
    FILE              *f    = NULL;
    Timer             timer = Timer();

    // Start time measurement
    timer.Start( env_handler );

    if ( !( f=fopen( matrix_file,"r" ) ) )
        throw std::runtime_error(
            "\n\tERROR:  Could not open matrix file "
            + std::string(matrix_file) +
            " for read access!\n"
            "\n\t\tUse -h(elp) for details.\n");

    mmio.MM_Read_Mtx_Crd_Size( f, file_rows, file_cols, file_nnz );

    // matrix must be square
    if ( file_cols != file_rows )
        throw std::runtime_error(
            "\n\tERROR:  Could not read matrix "
            "size and nnz's properly!\n"
            "\n\t\tMatrix must be square\n");

    nbr_cols    = file_cols;
    nbr_rows    = file_rows;
    nnz         = file_nnz;
    my_nbr_cols = nbr_cols;
    my_nbr_rows = nbr_rows;

    if ( MM_Is_Hermitian( mmio.matcode ) )
        hermitian = true;

    // Determine the distribution of rows and
    // columns across processors.
    Distribution().Basic_Distribution( env_handler, nbr_cols, my_nbr_cols,
                             split_pe, split_idx, start_idx );

    Read_Data( f );
    fclose(f);

    Data_To_Matrix( mtx, env_handler );

    // Stop time measurement
    timer.Stop( env_handler );
    timer.Report( env_handler );
}



template<class T_Field, class T_Format> void
Matrix_Reader<T_Field, T_Format>::Read_Data
(   FILE *f )
{
    int     my_nnz_col  = 0;

    // First parse through files to precompute
    // memory space to be allocated
    Count_NNZ( f, my_nnz_col );
    nnz_cols = my_nnz_col;

    // Second pass through file
    rewind( f );
    //Skipping header
    MMio().Skip_Header( f );

    cols = new T_Format[my_nnz_col];

    Data_To_Memory( f );
}



template<class T_Field, class T_Format> void
Matrix_Reader<T_Field, T_Format>::Count_NNZ
(   FILE    *f,
    int&    my_nnz_col )
{
    int     row,
            col;
    char    line[128];

    for ( int i = 0; i < file_nnz; i++ )
    {
        if ( !fgets(line,128,f) ) break;

        // Change ',' to ' '
        for ( int c = 0; line[c]; c++ )
            if ( line[c] == ',' )
                line[c] = ' ';

        sscanf( line, "%d %d\n", &row, &col );

        row--;
        col--;

        if (( col >= start_idx ) && ( col < ( start_idx + my_nbr_cols )))
            my_nnz_col++;
        // if matrix is symmetric column
        // data structure must carry full matrix
        if ( row != col )
            if (( row >= start_idx ) && ( row < ( start_idx + my_nbr_cols )))
                my_nnz_col++;
    }
}



template<class T_Field, class T_Format> void
Matrix_Reader<T_Field, T_Format>::Data_To_Matrix
(   Matrix<T_Field> *& mtx,
    const ENV_Handler& env_handler )
{
    int     *nnz_per_col = NULL,
            start_idx,
            col_buf_size = 0,
            A_buf_size   = 0,
            len_col,
            gap          = 0,
            max          = 0,
            pos          = 0;

    mtx = new Matrix<T_Field>( env_handler, nbr_cols, nbr_rows );;
    Compressed_Lines<T_Field> *lines = NULL;

    env_handler.Dist_Local_Chunks( my_nbr_cols, mtx->all_nbr_cols );

    // Filling start indices
    mtx->start_indices[0] = 0;
    for ( int pe = 1; pe < mtx->num_procs; pe++ )
        mtx->start_indices[pe] =  mtx->start_indices[pe-1] +
                                  mtx->all_nbr_cols[pe-1];

    // filling pe array
    memset( mtx->pe, 0, file_cols * sizeof(int) );
    for ( int pe = 0; pe < mtx->num_procs; pe++ )
    {
        start_idx = mtx->start_indices[pe];
        for ( int i = 0; i < mtx->all_nbr_cols[pe]; i++ )
            mtx->pe[start_idx + i] = pe;
    }

    mtx->my_nbr_cols = mtx->all_nbr_cols[mtx->my_id];
    mtx->my_start_idx = mtx->start_indices[mtx->my_id];

    nnz_per_col = new int[file_cols];
    memset(nnz_per_col, 0, file_cols * sizeof( int ));

    //Sorting the cols array for next calculation
    //of nnz per column
    std::sort( cols, cols + nnz_cols, Col_Comparator() );
    // Count the number of nonzeros in each col
    Count_NNZ_Cols( cols,( size_t ) nnz_cols, nnz_per_col );

    mtx->my_nnz = nnz_cols;
    mtx->c_lines = new Compressed_Lines<T_Field>( mtx->my_nbr_cols );
    lines = mtx->c_lines;

    // Filling array with column length information
    // of all pes. Checking singularity of matrix
    for ( int i = 0; i < mtx->my_nbr_cols; i++ )
    {
        len_col = nnz_per_col[i];
        if ( len_col == 0 )
        {
            delete [] nnz_per_col;
            throw Pe_Exception(
                "\n\n\tERROR:  Matrix is singular!\n"
                "\t\tPlease use only matrices which\n"
                "\t\tare invertible.\n");
        }
        lines->len_cols[i] = len_col;
        col_buf_size += len_col;
        A_buf_size += len_col;
        if (diagonal_elements[i] == false)
        {
            delete [] nnz_per_col;
            throw Pe_Exception(
                "\n\n\tERROR:  Matrix has zero diagonal entry!\n"
                "\t\tPlease use only matrices which\n"
                "\t\thave values on the diagonal.\n");
        }
    }

    mtx->c_lines->col_idcs_buf = new int[col_buf_size];
    mtx->c_lines->col_buf      = new T_Field[A_buf_size];

    // Set pointers & Fill compressed lines
    for ( int col = 0; col < mtx->my_nbr_cols; col++ )
    {
        len_col = nnz_per_col[col];
        lines->mtx_vals[col] = &( lines->col_buf[gap] );
        lines->col_idcs[col] = &( lines->col_idcs_buf[gap] );
        gap += len_col;
        Fill_Lines( len_col, lines, col, pos );
    }

    // Get length of all cols
    env_handler.Dist_Col_Length( lines->len_cols,
                                 mtx->my_nbr_cols,
                                 mtx->len_all_cols,
                                 mtx->all_nbr_cols,
                                 mtx->start_indices);

    // Get the maximum number of nnz per
    // column/row of all pes
    for ( int i = 0; i < mtx->n; i++ )
        if ( nnz_per_col[i] > max )
            max = nnz_per_col[i];

    env_handler.Dist_Max( max, mtx->max_nnz );

    //Initializing the remote transfer buffers
    mtx->remote_col_buf = new T_Field[mtx->max_nnz];
    memset(mtx->remote_col_buf, 0, mtx->max_nnz * sizeof( T_Field ));
    mtx->remote_col_idcs_buf = new int[mtx->max_nnz];
    memset(mtx->remote_col_idcs_buf, 0, mtx->max_nnz * sizeof( int ));

    delete [] nnz_per_col;
}



template<class T_Field, class T_Format> void
Matrix_Reader<T_Field, T_Format>::Count_NNZ_Cols
(   T_Format const  *in,
    size_t          size,
    int             *out )
{
    // In case of zero column dont
    // increment the column lenght
    if (size == 0) return;

    //Columns with no nnz will be present with
    //no entry in out
    int diff = 0;
    ( *out )++;
    for (size_t i = 1; i < size; i++)
    {
        diff = in[i].j - in[i-1].j;
        if ( diff >= 1 )
            out += diff;
        ( *out )++;
    }
}



template<class T_Field, class T_Format> void
Matrix_Reader<T_Field, T_Format>::Count_NNZ_Rows
(   T_Format const  *in,
    size_t          size,
    int             *out )
{
    // In case of zero column dont
    // increment the column lenght
    if (size == 0) return;

    //Columns with no nnz will be present with
    //no entry in out
    int diff = 0;
    ( *out )++;
    for (size_t i = 1; i < size; i++)
    {
        diff = in[i].i - in[i-1].i;
        if ( diff >= 1 )
            out += diff;
        ( *out )++;
    }
}
