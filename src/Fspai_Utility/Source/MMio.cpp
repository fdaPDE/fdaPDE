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

//file includings
#include "../Include/MMio.h"

//C++ includings
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdexcept>


int
MMio::Get_Matrix_Type
(   char *matrix_file )
{
    int     type;
    FILE    *f;

    if ( !( f = fopen( matrix_file,"r" ) ) )
        throw std::runtime_error(
            "\n\n\tERROR:  Could not open matrix file "
            + std::string( matrix_file ) +
            " for read access!\n"
            "\n\t\tUse -h(elp) for details.\n");

    Parse_Header( f, type );
    fclose( f );
    return type;
}



void
MMio::Parse_Header
(   FILE* f,
    int&  type)
{
    if ( MM_Read_Banner( f, &matcode ) != 0)
    {
        fclose( f );
        throw std::runtime_error(
            "\n\n\tERROR:  Could not read mm header properly!\n"
            "\n\t\tCould not process Matrix Market banner.\n");
    }

    if (! ( MM_Is_Real( matcode )  || MM_Is_Complex( matcode ) ) ||
           !MM_Is_Matrix( matcode ) ||
        !(  MM_Is_Symmetric( matcode ) || MM_Is_Hermitian( matcode ) ) )
    {
        fclose( f );
        throw std::runtime_error(
            "\n\n\tERROR:  Could not read mm header properly!"
            "\n\t\tMM matrix must be real or complex,"
            "\n\t\tcoordinate, and symmetric or hermitian.\n");
    }

    if ( MM_Is_Real( matcode ) )
         type = real;
    else type = complex;
}



void
MMio::Parse_Pattern_Header
(   FILE *f )
{
    if ( MM_Read_Banner( f, &matcode ) != 0 )
    {
        fclose( f );
        throw std::runtime_error(
            "\n\tERROR:  Could not read pattern "
                "header properly!\n"
            "\n\t\tCould not process Matrix Market "
                "pattern banner.\n");
    }

    if (!( MM_Is_Pattern( matcode )) ||
          !MM_Is_Matrix(  matcode ) ||
          !MM_Is_Sparse(  matcode ) ||
          !MM_Is_Symmetric( matcode ))
    {
        fclose(f);
        throw std::runtime_error(
            "\n\n\tERROR:  Could not read pattern "
                "header properly!\n"
            "\n\t\tMM matrix must be pattern, "
            "coordinate, and symmetric\n");
    }
}


bool MMio::Is_Symmetric
( ) const
{
    if ( MM_Is_Symmetric( matcode ) ) return true;
    else return false;
}



void
MMio::Skip_Header
( FILE *f )
{
    int     m,
            n,
            nnz,
            dummy1;

    Parse_Header( f, dummy1 );
    MM_Read_Mtx_Crd_Size( f, m, n, nnz );
}



void MMio::Clear_Typecode
( )
{
    MM_Clear_Typecode( &matcode );
}



void MMio::Print_Typecode
( ) const
{
    std::cout <<
    matcode[0] << " " <<
    matcode[1] << " " <<
    matcode[2] << " " <<
    matcode[3] << " " <<
    std::endl;
}



int
MMio::MM_Read_Banner
(   FILE*        f,
    MM_typecode* matcode )
{
    char line[MM_MAX_LINE_LENGTH];
    char banner[MM_MAX_TOKEN_LENGTH];
    char mtx[MM_MAX_TOKEN_LENGTH];
    char crd[MM_MAX_TOKEN_LENGTH];
    char data_type[MM_MAX_TOKEN_LENGTH];
    char storage_scheme[MM_MAX_TOKEN_LENGTH];
    char *p;

    // First init typecode with empty strings ,
    // typecode can carry 4 elements
    Clear_Typecode();

    if ( fgets( line, MM_MAX_LINE_LENGTH, f ) == NULL )
        return MM_PREMATURE_EOF;

    if ( sscanf( line, "%s %s %s %s %s",
                 banner,
                 mtx,
                 crd,
                 data_type,
                 storage_scheme) != 5 )
        return MM_PREMATURE_EOF;

    // convert to upper case
    for ( p=mtx; *p!='\0';              *p=toupper( *p ), p++ );
    for ( p=crd; *p!='\0';              *p=toupper( *p ), p++ );
    for ( p=data_type; *p!='\0';        *p=toupper( *p ), p++ );
    for ( p=storage_scheme; *p!='\0';   *p=toupper( *p ), p++ );

    // check for banner
    if ( strncmp( banner,
                  MatrixMarketBanner,
                  strlen( MatrixMarketBanner )) != 0 )
        return MM_NO_HEADER;


    // first field should be "mtx"
    if ( strcmp( mtx, MM_MTX_STR ) != 0 )
        return  MM_NOT_MTX;

    MM_Set_Matrix( matcode );

    // second field describes whether this is a
    // sparse matrix (in coordinate storage) or a dense array
    if ( strcmp( crd, MM_SPARSE_STR ) == 0 )
        MM_Set_Sparse( matcode );
    else
        if ( strcmp( crd, MM_DENSE_STR ) == 0 )
            MM_Set_Dense( matcode );

    // third field
    if ( strcmp( data_type, MM_REAL_STR ) == 0 )
        MM_Set_Real( matcode );
    else
        if ( strcmp( data_type, MM_COMPLEX_STR ) == 0 )
            MM_Set_Complex( matcode );
    else
        if ( strcmp( data_type, MM_PATTERN_STR ) == 0 )
            MM_Set_Pattern( matcode );
    else
        if ( strcmp( data_type, MM_INT_STR ) == 0 )
            MM_Set_Integer( matcode );


    // fourth field
    if ( strcmp( storage_scheme, MM_GENERAL_STR ) == 0 )
        MM_Set_General( matcode );
    else
        if ( strcmp( storage_scheme, MM_SYMM_STR ) == 0 )
            MM_Set_Symmetric( matcode );
    else
        if ( strcmp( storage_scheme, MM_HERM_STR ) == 0 )
            MM_Set_Hermitian( matcode );
    else
        if ( strcmp( storage_scheme, MM_SKEW_STR ) == 0 )
            MM_Set_Skew( matcode );

    //everything is correct
    return EXIT_SUCCESS;
}



void
MMio::MM_Read_Mtx_Crd_Size
(   FILE *&f,
    int  &M,
    int  &N,
    int  &nz ) const
{
    char line[MM_MAX_LINE_LENGTH];

    // set return null parameter values,
    // in case we exit with errors
    M = N = nz = 0;

    // now continue scanning until you reach the end-of-comments
    do
    {
        if ( fgets( line,MM_MAX_LINE_LENGTH,f ) == NULL )
        {
            fclose( f );
            throw std::runtime_error(
                "\n\tERROR:  Could not read matrix size "
                    "and nnz's properly!\n"
                "\n\t\tVerify that second line in *.mtx "
                    "file is:  m n nnz\n");
        }
    } while ( line[0] == '%' );


    if ( sscanf( line, "%d %d %d", &M, &N, &nz ) != 3 )
    {
        fclose( f );
        throw std::runtime_error(
            "\n\tERROR:  Could not read matrix size and "
                "nnz's properly!\n"
            "\n\t\tVerify that second line in *.mtx file "
                "is:  m n nnz\n");
    }
}



void
MMio::MM_Read_Pattern_Crd_Size
(   FILE *&f,
    int  &M,
    int  &N,
    int  &nz ) const
{

    char line[MM_MAX_LINE_LENGTH];

    // set return null parameter values,
    // in case we exit with errors
    M = N = nz = 0;

    // now continue scanning until you reach the end-of-comments
    do
    {
        if ( fgets( line,MM_MAX_LINE_LENGTH,f ) == NULL )
        {
            fclose( f );
            throw std::runtime_error(
                "\n\tERROR:  Could not read pattern size "
                    "and nnz's properly!\n"
                "\n\t\tVerify that second line in pattern "
                    "file is:  m n nnz\n");
        }

    } while ( line[0] == '%' );

    if ( sscanf( line, "%d %d %d", &M, &N, &nz ) != 3 )
    {
        fclose( f );
        throw std::runtime_error(
            "\n\tERROR:  Could not read pattern size and "
                "nnz's properly!\n"
            "\n\t\tVerify that second line in pattern "
                "file is:  m n nnz\n");
    }
}

