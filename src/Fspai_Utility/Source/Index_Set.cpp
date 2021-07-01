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
#include "../Include/Index_Set.h"

//C++ includings
#include <stdio.h>
#include <stdlib.h>
#include <iostream>


Index_Set::Index_Set
(   int len_a )
{
    idcs = NULL;
    len = len_a;
    idcs = new int[len];
}



Index_Set::~Index_Set
( )
{
    if ( idcs )
        delete [] idcs;
}



void
Index_Set::Print_Index_Set
(   char* str ) const
{
    std::cout <<"\t" << str <<  " ";
    for ( int i = 0; i < len; i++ )
            std::cout << idcs[i] << " ";
    std::cout << std::endl;
}



Index_Set*
Index_Set::Copy_Index_Set
(   Index_Set* in )
{
    Index_Set* out  = new Index_Set( in->len );
    memcpy( out->idcs, in->idcs, in->len * sizeof( int ) );
    return out;
}



Index_Set*
Index_Set::Set_Difference
(   Index_Set* is_a,
    Index_Set* is_b )
{
    int     a,
            b,
            pos = 0,
            idx_a = 0,
            idx_b = 0,
            diff,
            len_a = is_a->len,
            len_b = is_b->len;

    //is_a->len is dimension of matrix
    Index_Set*  is_diff = new Index_Set( len_a );

    // First sorting index set a (index set b is
    // already sorted, then parallel reading.
    // This way we get nearly O(nlogn) complexity
    std::sort( is_a->idcs, is_a->idcs + len_a );

    // "Parallel" reading of both index sets
    while ( idx_a < len_a && idx_b < len_b )
    {
        a = is_a->idcs[idx_a];
        b = is_b->idcs[idx_b];
        if ( a < b )
        {
            is_diff->idcs[pos++] = a;
            idx_a++;
        }
        else if ( a > b )
        {
            idx_b++;
        }
        else    // a == b
        {
            idx_a++;
            idx_b++;
        }
    }

    // copy possible remaining elements of
    // is_a into is_diff
    diff = len_a - idx_a;
    if ( diff > 0 )
    {
        memcpy( is_diff->idcs + pos,
                is_a->idcs + idx_a,
                diff * sizeof( int ) );
        pos += diff;
    }
    is_diff->len = pos;

    return is_diff;
}



Index_Set*
Index_Set::Set_Intersection
(   Index_Set*  is_a,
    Index_Set*  is_b )
{
    int     a,
            b,
            pos = 0,
            idx_a = 0,
            idx_b = 0,
            len_a = is_a->len,
            len_b = is_b->len;

    Index_Set*  is_int = NULL;


    // Size of new J_tilde cannot be larger
    // than the old size of J_tilde
    is_int = new Index_Set( len_a );

    while ( idx_a < len_a && idx_b < len_b )
    {
        a = is_a->idcs[idx_a];
        b = is_b->idcs[idx_b];
        if ( a < b )
            idx_a++;
        else if ( a > b )
            idx_b++;
        else  // a == b
        {
            is_int->idcs[pos++] = a;
            idx_a++;
            idx_b++;
        }
    }
    is_int->len = pos;

    // is_a is the old J_tilde heap index set.
    // delete it because no realloc is performed
    // Dirty - no input parameters should be
    // deleted
    delete is_a;
    return is_int;
}



Index_Set*
Index_Set::Set_Union
(   Index_Set*     is_a,
    Index_Set*     is_b )
{
    int     a,
            b,
            pos = 0,
            idx_a = 0,
            idx_b = 0,
            len_a = is_a->len,
            len_b = is_b->len,
            diff = 0;

    Index_Set *is_union = new Index_Set( len_a + len_b );

    while ( idx_a < len_a && idx_b < len_b )
    {
        a = is_a->idcs[idx_a];
        b = is_b->idcs[idx_b];
        if ( a < b )
        {
            is_union->idcs[pos++] = a;
            idx_a++;
        }
        else if ( a > b )
        {
            is_union->idcs[pos++] = b;
            idx_b++;
        }
        else    // a == b
        {
            is_union->idcs[pos++] = a;
            idx_a++;
            idx_b++;
        }
    }

    // one of the index set has finished,
    // remaining elements of the other index set
    // will now be copied.
    diff = len_a - idx_a;
    if ( diff )
    {
        memcpy( is_union->idcs + pos,
                is_a->idcs + idx_a,
                diff * sizeof(int));
        pos += diff;
    }
    else
    {
        diff = len_b - idx_b;
        memcpy( is_union->idcs + pos,
                is_b->idcs + idx_b,
                diff * sizeof( int ) );
        pos += diff;
    }
    is_union->len = pos;

    return is_union;
}



int
Index_Set::Get_El_Idx
(   const int  el )
{
    for ( int idx = 0; idx < len; idx++ )
        if ( el == idcs[idx] ) return idx;
    return -1;
}



bool
Index_Set::Has_Idx
(   const int  idx) const
{
    for (int i = 0; i < len; i++)
        if(idcs[i] == idx) return true;
    return false;
}
