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

#ifndef GUARD_HASH_TABLE_H
#define GUARD_HASH_TABLE_H

/// prime number for linear rehashing
const int jump = 5;


///////////////////////////////////////////
///     \class Hash_Table
///     \brief This class represents the
///            hash table which stores all
///            remotely requested matrix columns.
///            This way the MPI traffic can
///            be reduced.
///////////////////////////////////////////
template <class T_Field>
class Hash_Table
{
    public:

        /// Empty Constructor
        Hash_Table() { };

        /// Constructor
        Hash_Table(int num);

        /// Destructor
        ~Hash_Table();

        // Members

        /// size of the hash table
        int         size;

        /// 2D array containing the column indices
        int         **col_idcs_table;

        /// 2D array containing the column values
        T_Field     **vals_table;


        // Methods
        //============================================================
        //========== Template methods - see Hash_Table.imp ===========
        //============================================================

        ///////////////////////////////////////////////
        ///     \brief Inserting column data into hash
        ///            table.
        ///
        ///     \param idx Index of column to be inserted
        ///     \param col_idcs_buf Column indices to be
        ///                         inserted
        ///     \param col_buf Column values to be
        ///                    inserted
        ///     \param col_len Lentgth of the column
        ///                    which has to be inserted
        ///     \return Data is not within hash table
        ///             if user set -hs 0, if hs > 0
        ///             then it is
        ///////////////////////////////////////////////
        bool    Insert
                (   int     idx,
                    int     *col_idcs_buf,
                    T_Field *col_buf,
                    int     col_len);

        ///////////////////////////////////////////////
        ///     \brief Look up if element is in hash table
        ///            and return data if so.
        ///
        ///     \param idx Index of column to be extraced
        ///     \param col_idcs_buf Column indices to be
        ///                         extraced
        ///     \param col_buf Column values to be
        ///                    extraced
        ///     \param col_len Lentgth of the column
        ///                    which has to be extraced
        ///     \param loc Location of column within hashtable.
        ///     \return Data is not within hash table
        ///             if user set -hs 0, if hs > 0
        ///             then it is
        ///////////////////////////////////////////////
        bool    Look_Up
                (   int         idx,
                    int         *&col_idcs_buf,
                    T_Field     *&col_buf,
                    int&        col_len,
                    const int   loc);

        ///////////////////////////////////////////////
        ///     \brief  Searching for position within
        ///             hash table.
        ///
        ///     This will attempt at as many jumps as
        ///     hash table size is not exceeded of
        ///     linear rehashing before returning a
        ///     location. If no location could be found
        ///     -1 will be returned - element is not
        ///     within hash table.
        ///
        ///     \param idx Index of column to be
        ///                found
        ///     \return The place where the column data
        ///             is stored within hash table
        //////////////////////////////////////////////
        int     Get_Location(int idx);

    private:

        ///////////////////////////////////////////////
        ///     \brief  Finds the location in which to
        ///             make an insertion.
        ///
        ///     This will attempt at most 5 jumps of
        ///     linear rehashing before returning a
        ///     location.  In that case, whatever was
        ///     already in the location will be
        ///     discarded.
        ///
        ///     \param idx Index of column to be
        ///                inserted
        ///     \return Location where to insert the
        ///             column
        //////////////////////////////////////////////
        int     Find_Location(int idx);
};

#include "Hash_Table.imp.h"

#endif
