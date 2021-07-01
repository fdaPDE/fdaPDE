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


template<class T_Field>
Hash_Table<T_Field>::Hash_Table
(   int size_m )
{
    // hash size should be
    // prime to 5 because of
    // linear rehashing
    switch(size_m)
    {
        case 0: size_m = 0;
        case 1: size_m = 101;
        case 2: size_m = 503;
        case 3: size_m = 2503;
        case 4: size_m = 12503;
        case 5: size_m = 62501;
        case 6: size_m = 104743;
    }

    if (! (size_m % jump))
    {
        std::cout << "\n\t  Warning in initializing hash table:"
                  << std::endl;
        std::cout << "\t\t Size of hash table and linear "
                "rehashing factor are not relative primes!\n"
                  << std::endl;
    }

    size            = size_m;
    col_idcs_table  = new int*[size];
    vals_table      = new T_Field*[size];
    for (int i = 0; i < size; i++)
    {
        col_idcs_table[i]   = NULL;
        vals_table[i]       = NULL;
    }
}



template<class T_Field>
Hash_Table<T_Field>::~Hash_Table
( )
{
    for (int i = 0; i < size; i++)
    {
        if (col_idcs_table[i])
            delete [] col_idcs_table[i];
        if (vals_table[i])
            delete [] vals_table[i];
    }
    if(col_idcs_table) delete [] col_idcs_table;
    if(vals_table) delete [] vals_table;
}



template <class T_Field>  bool
Hash_Table<T_Field>::Insert
(  int      idx,
   int      *col_idcs_buf,
   T_Field  *col_buf,
   int      col_len)
{
    int     loc,
            *htcol_idcs_buf;

    T_Field *htvals_buf;

    // If user does not want to use hash table
    // no look up has to be performed.
    // This may only occur if user set -hs 0
    if (! this) return false;

    loc = Find_Location(idx);

    // The column/row indices arrays have 2 elements more
    // because into first position the index of the
    // column/row is inserted, and on second position
    // there is the length of the column/row

    htcol_idcs_buf      = new int[col_len + 2];
    htcol_idcs_buf[0]   = idx;
    htcol_idcs_buf[1]   = col_len;
    memcpy(&htcol_idcs_buf[2], col_idcs_buf, col_len * sizeof(int));
    col_idcs_table[loc] = htcol_idcs_buf;
    htvals_buf          = new T_Field[col_len];
    memcpy(htvals_buf, col_buf, col_len * sizeof(T_Field));
    vals_table[loc]     = htvals_buf;

    return true;
}



template <class T_Field>  bool
Hash_Table<T_Field>::Look_Up
( int       idx,
  int       *&col_idcs_buf,
  T_Field   *&col_buf,
  int&      col_len,
  const int loc)
{
    // Hash table pointer is NUll, return.
    // This can only occur if user does not
    // want to use hash table an set -hs 0
    if (! this) return false;

    // found them
    col_idcs_buf    = col_idcs_table[loc];
    col_buf         = vals_table[loc];
    col_len         = col_idcs_buf[1];
    col_idcs_buf++;
    col_idcs_buf++;

    return true;
}



template <class T_Field>  int
Hash_Table<T_Field>::Find_Location
(   int idx )
{
    unsigned long   loc;
    int             *buf,
                    nbr_tries;

    loc         = idx % size;
    nbr_tries   = 1;
    buf         = col_idcs_table[loc];

    // Looking at loc position, if
    // the buffer is not set this
    // position is free.
    // If set, do linear rehash
    // max. 5 times.
    do
    {
        if (buf == NULL) return loc;
        else    // linear rehash
        {
            loc += jump;
            loc %= size;
            nbr_tries++;
            buf = col_idcs_table[loc];
        }
    }
    while (nbr_tries < 5);

    // The 5 jumps are done.
    // Delete what's there and return loc.
    if (buf != NULL)
    {
        delete [] buf;
        delete [] vals_table[loc];
    }
    return loc;
}



template <class T_Field>  int
Hash_Table<T_Field>::Get_Location
(   int idx )
{
    unsigned long   loc;
    int             nbr_tries,
                    *buf;

    // Looking at loc position, if
    // the buffer is not set this
    // position is free.
    // If set, do linear rehash
    // as long as hash table size
    // is not exceeded.
    if (! this) return -1;
    loc         = idx % size;
    nbr_tries   = 1;
    do
    {
        buf = col_idcs_table[loc];
        if (buf == NULL)    return -1;
        if (buf[0] == idx)  return loc;
        else    // linear rehash
        {
            loc += jump;
            loc %= size;
            nbr_tries++;
        }
    }
    while (nbr_tries < size);
    return -1;
}

