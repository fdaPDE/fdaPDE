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
#include "../Include/Pattern.h"
#include "../Include/MMio.h"

//C++ includings
#include <iostream>
#include <stdexcept>


Pattern::Pattern()
{
    num_procs       = 0;
    my_id           = -1;
    file_cols       = 0;
    file_rows       = 0;
    file_nnz        = 0;
    nnz_cols        = 0;
    split_idx       = -1;
    split_pe        = -1;
    start_idx       = -1;
    j_sets          = NULL;
    start_indices   = NULL;
    pe              = NULL;
    all_nbr_cols    = NULL;
    len             = 0;
    max_nnz         = 0;
    my_nbr_cols     = 0;
    next_col        = 0;
    j_sets          = NULL;
    cols            = NULL;
    symmetric       = false;
}



Pattern::~Pattern()
{
    if (j_sets)
    {
        for (int i = 0; i < my_nbr_cols; i++)
            if (j_sets[i]) delete j_sets[i];
        delete [] j_sets;
    }
    if (start_indices)  delete [] start_indices;
    if (pe)             delete [] pe;
    if (all_nbr_cols)   delete [] all_nbr_cols;
    if (cols)           delete [] cols;
}



void
Pattern::Arbitrary_Pattern
(   ENV_Handler&    env_handler,
    char            *pattern_file,
    const int       dim)
{
    FILE    *f = NULL;
    MMio    mmio;

    if ( !( f = fopen(pattern_file,"r") ) )
        throw std::runtime_error(
            "\n\n\tERROR:  Could not open pattern file " +
                std::string(pattern_file) +
            " for read access!\n"
            "\n\t\tUse -h(elp) for details.\n");

    mmio.Parse_Pattern_Header( f );
    symmetric  = mmio.Is_Symmetric();
    mmio.MM_Read_Pattern_Crd_Size( f, file_rows, file_cols, file_nnz );

    if ( file_rows != file_cols )
        throw std::runtime_error(
            "\n\tERROR:  Could not read pattern size "
            "and nnz's properly!\n"
            "\n\t\tPattern must be square.\n");

    if( file_rows != dim )
        throw std::runtime_error(
            "\n\tERROR:  Mismatch in matrix and pattern dimension!"
            "\n\t\tMatrix and pattern must have the same "
            "number of rows!\n");

    my_nbr_cols = file_cols;

    // Determine the distribution of rows and
    // columns across pocessors.
    Distribution().Basic_Distribution(
        env_handler, file_cols, my_nbr_cols,
        split_pe, split_idx, start_idx );

    Read_Data( f );
    fclose( f );
    Data_To_Pattern( env_handler );
}



void
Pattern::Read_Data
( FILE *f)
{
    int dummy1,dummy2,dummy3;

    Count_NNZ(f);

    // Second pass through file
    rewind(f);

    MMio().Parse_Pattern_Header( f );
    MMio().MM_Read_Pattern_Crd_Size( f, dummy1, dummy2, dummy3 );
    cols = new RC[nnz_cols];

    Data_To_Memory( f );
}



void Pattern::Count_NNZ
( FILE* f )
{
    int             row,
                    col;
    char            line[128],
                    *ierr;

    for ( int i = 0; i < file_nnz; i++ )
    {
        ierr = fgets( line, 128,f );

        // Change ',' to ' '
        for ( int ii = 0; line[ii]; ii++ )
            if (line[ii] == ',')
                line[ii] = ' ';

        sscanf( line, "%d %d \n", &row, &col );

        row--;
        col--;

        if (( col >= start_idx ) && ( col < ( start_idx + my_nbr_cols ) ))
            nnz_cols++;
    }
}



void Pattern::Data_To_Memory
( FILE *f )
{
    int     row,
            col,
            lencol = 0;
    char    line[128];

    for ( int i = 0; i < file_nnz; i++ )
    {
        if ( !fgets( line,128,f ) ) break;

        // Change ',' to ' '
        for ( int ii = 0; line[ii]; ii++ )
            if ( line[ii] == ',' )
                line[ii] = ' ';

        sscanf( line, "%d %d\n", &row, &col );

        row--;
        col--;

        if (( col >= start_idx ) && ( col < start_idx + my_nbr_cols ))
        {
            cols[lencol].i = row;
            cols[lencol].j = col;
            lencol++;
        }
    }
}



void
Pattern::Data_To_Pattern
(   ENV_Handler& env_handler )
{
    int         *nnz_per_col = NULL,
                idx = 0,
                leni,
                max = 0,
                start_index;

    Index_Set   *i_set = NULL;

    env_handler.Get_Environment_Params( num_procs,my_id );

    all_nbr_cols  = new int[num_procs];
    start_indices = new int[num_procs];
    pe            = new int[file_cols];
    j_sets        = new Index_Set*[my_nbr_cols];
    for ( int i = 0; i < my_nbr_cols; i++ )
        j_sets[i] = NULL;
    len           = file_cols;

    env_handler.Dist_Local_Chunks( my_nbr_cols, all_nbr_cols );

    // Filling start indices
    start_indices[0] = 0;
    for ( int pe_loc = 1; pe_loc < num_procs; pe_loc++ )
        start_indices[pe_loc] = start_indices[pe_loc-1] +
                                all_nbr_cols[pe_loc-1];

    // filling pe array
    memset( pe, 0, file_cols * sizeof( int ) );
    for ( int pe_loc = 0; pe_loc < num_procs; pe_loc++ )
    {
        start_index = start_indices[pe_loc];
        for ( int i = 0; i < all_nbr_cols[pe_loc]; i++ )
            pe[start_index + i] = pe_loc;
    }

    my_nbr_cols  = all_nbr_cols[my_id];
    my_start_idx = start_indices[my_id];
    nnz_per_col  = new int[my_nbr_cols];
    memset( nnz_per_col, 0, my_nbr_cols * sizeof( int ) );

    //Sorting the cols array for next calculation
    //of nnz per column
    std::sort(  cols,
                cols + nnz_cols,
                Idx_Comparator() );

    Count_NNZ_Cols( cols,
                    ( size_t ) nnz_cols,
                    nnz_per_col );

   // Filling pattern data structure
    for ( int col = 0; col < my_nbr_cols; col++ )
    {
        leni = nnz_per_col[col];
        // Checking whether triangular pattern
        // has no zero column
        if ( leni == 0 )
        {
            delete [] nnz_per_col;
            throw Pe_Exception(
                "\n\n\tERROR:  Pattern has zero column!\n"
                "\t\tPlease use only patterns which\n"
                "\t\thave at least one entry in each column.\n");
        }
        i_set = new Index_Set(leni);
        for ( int nnz = 0; nnz < leni; nnz++ )
            i_set->idcs[nnz] = cols[idx++].i;
        j_sets[col] = i_set;
    }

    // Get the maximum number of nnz per
    // column/row of all pes
    for ( int i = 0; i < my_nbr_cols; i++ )
        if ( nnz_per_col[i] > max )
            max = nnz_per_col[i];

    env_handler.Dist_Max( max, max_nnz );
    delete [] nnz_per_col;
}



void
Pattern::Count_NNZ_Cols
(   RC const    *in,
    size_t      size,
    int         *out )
{
    // In case of zero column dont
    // increment the column lenght
    if (size == 0) return;

    // Columns with no nnz will be
    // present with no entry in out
    int diff = 0;

    ( *out )++;
    for ( size_t i = 1; i < size; i++ )
    {
        diff = in[i].j - in[i-1].j;
        if ( diff >= 1 )
            out += diff;
        ( *out )++;
    }
}



void
Pattern::Count_NNZ_Rows
(   RC const *in,
    size_t   size,
    int      *out)
{
    // In case of zero column dont
    // increment the column lenght
    if (size == 0) return;

    // Columns with no nnz will be
    // present with no entry in out
    int diff = 0;

    ( *out )++;
    for ( size_t i = 1; i < size; i++ )
    {
        diff = in[i].i - in[i-1].i;
        if ( diff >= 1 )
            out += diff;
        ( *out )++;
    }
}



void
Pattern::Print_Pattern_Data
( )
{   /*
    std::cout << "\n\tPattern Data:\t\n"
              << "\n\tmy_id:\t\t" << my_id
              << "\n\tnum_procs:\t" << num_procs
              << "\n\tmy_nbr_cols:\t" << my_nbr_cols
              << "\n\tmy_start_idx:\t" << my_start_idx
              << std::endl;

    std::cout << "\tall_nbr_cols:\t";
    for ( int i = 0; i < num_procs; i++ )
            std::cout << all_nbr_cols[i] << " ";
    std::cout << std::endl;

    std::cout << "\tstart_indices:\t";
    for ( int i = 0; i < num_procs; i++ )
            std::cout << start_indices[i] << " ";
    std::cout << std::endl;

    std::cout << "\tpe:\t\t";
    for ( int i = 0; i < len; i++ )
            std::cout << pe[i] << " ";
    std::cout << std::endl;

    for ( int i = 0; i < my_nbr_cols; i++ )
    {
        int l = j_sets[i]->len;
        std::cout << "\n\tlen[" << i << "]: "
                  << l << std::endl;
        std::cout << "\t";
        std::cout << "j_set:  ";
        for ( int j = 0; j < l; j++ )
            std::cout << j_sets[i]->idcs[j] << " ";
        std::cout << std::endl;
        std::cout << "\t";
    }
   */
   return;
}



void
Pattern::Diagonal_Pattern
(   ENV_Handler& env_handler,
    const int    mtx_dim,
    const int    mtx_my_nbr_cols,
    const int    mtx_my_start_idx)
{
    int         start_index;
    Index_Set   *i_set;

    env_handler.Get_Environment_Params( num_procs,my_id );

    all_nbr_cols = new int[num_procs];
    start_indices= new int[num_procs];
    pe           = new int[mtx_dim];
    j_sets       = new Index_Set*[mtx_my_nbr_cols];
    for ( int i = 0; i < mtx_my_nbr_cols; i++ )
        j_sets[i] = NULL;
    len          = mtx_dim;

    env_handler.Dist_Local_Chunks( mtx_my_nbr_cols, all_nbr_cols );

    // Filling start indices
    start_indices[0] = 0;
    for ( int pe_loc = 1; pe_loc < num_procs; pe_loc++ )
        start_indices[pe_loc] = start_indices[pe_loc-1] +
                                all_nbr_cols[pe_loc-1];

    // filling pe array
    memset( pe, 0, mtx_dim * sizeof( int ) );
    for ( int pe_loc = 0; pe_loc < num_procs; pe_loc++ )
    {
        start_index = start_indices[pe_loc];
        for ( int i = 0; i < all_nbr_cols[pe_loc]; i++ )
            pe[start_index + i] = pe_loc;
    }

    my_nbr_cols = all_nbr_cols[my_id];
    my_start_idx = start_indices[my_id];

    // Diagonal pattern contains only index
    // sets of one element
    for( int col = 0; col < mtx_my_nbr_cols; col++ )
    {
        i_set = new Index_Set( 1 );
        i_set->idcs[0] = mtx_my_start_idx + col;
        j_sets[col] = i_set;
    }
    max_nnz = 1;
}
