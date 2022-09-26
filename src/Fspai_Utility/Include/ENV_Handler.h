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

#ifndef ENV_MPI

#	ifndef ENV_HANDLER_H_
#	define ENV_HANDLER_H_

// C/C++ includings
#include <iostream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <stdexcept>

/////////////////////////////////////////////////////////
///     \class ENV_Handler
///     \brief Interface methods to specified
///            environment. The environment handler is
///            built via compile-time polymorphism. If the
///            compile flag "ENV_MPI" is not used this
///            sequential handler will be used. Each method
///            has its corresponding parallel counterpart
///            in the parallel ENV_Handler class.
/////////////////////////////////////////////////////////
class ENV_Handler
{
	public:

	    /// Constructor
		ENV_Handler() {}

		/// Destructor
		~ENV_Handler() {}

        ////////////////////////////////////////////////////////
        ///     \brief  Sets the member environment parameters.
        ///
        ///     \param num_procs_ Number of processors used in
		///                       this environment.
        ///     \param my_id_ Environment identifier for this pe.
        ////////////////////////////////////////////////////////
        void    Get_Environment_Params
                (   int& num_procs_,
                    int& my_id_ )                       const ;

        ////////////////////////////////////////////////////////
        ///     \brief  Initializes the environment and sets the
        ///             environment specific parameters.
        ///
        ///     \param argc Number of shell parameters.
        ///     \param argv Shell parameters.
        ////////////////////////////////////////////////////////
        void    Initialize_Environment
                (   int   argc,
                    char* argv[] );

        ////////////////////////////////////////////////////////
        ///     \brief  Finalizes/Closes the open environment.
        ////////////////////////////////////////////////////////
        void    Finalize( )                             const { }

        ////////////////////////////////////////////////////////
        ///     \brief  Aborts the open environment.
        ///
        ///     This is used in case that an error occurs.
        ////////////////////////////////////////////////////////
        void    Abort_Environment( )                    const { }

        ////////////////////////////////////////////////////////
        ///     \brief  Environment barrier.
        ///
        ///     For a sequential environment no barrier is used.
        ////////////////////////////////////////////////////////
        void    Barrier()   {}

        ////////////////////////////////////////////////////////
        ///     \brief  Sets the environment specific help string.
        ///
        ///     Due to the compiled environment a different
        ///     help string (-help) has to be used as different
        ///     user parameters become available.
        ///
        ///     \param known_params The available parameters for
        ///                         this environment.
        ///     \param options The environment specific help
        ///                    string.
        ////////////////////////////////////////////////////////
        void    Set_Help
                ( std::string& known_params,
                  std::string& options )                const;

        ////////////////////////////////////////////////////////
        ///     \brief  Environment specific time measurement.
        ///
        ///     \param max Maximum measured time over all pes.
        ///     \param sum_time Sum of time of all pes in this
        ///                     environment.
        ////////////////////////////////////////////////////////
        void    Time_Diff
                (   double& max,
                    double  sum_time )                  const;

        ////////////////////////////////////////////////////////
        ///     \brief  Environment specific time measurement for
        ///             debug output.
        ///
        ///     This method is used for the parallel environment
        ///     to store different subtimes in the dgb-arrays.
        ///     See class Timer. It is mainly used for debugging
        ///     communication problems.
        ///
        ///     \param max Maximum measured time over all pes.
        ///     \param sum Sum of time of all pes in this
        ///                environment.
        ////////////////////////////////////////////////////////
        void    Dbg_Timers
                (   double     dbg_time,
                    double*&   sum )                    const { }

        ////////////////////////////////////////////////////////
        ///     \brief  Starts the current time on this pe in
        ///             the specific environment.
        ///
        ///     \param curr_time The current time on this pe.
        ///     \param tp Time structure.
        ////////////////////////////////////////////////////////
        void    Get_Time
                (  double&  curr_time,
                   timespec tp )                        const;

        ////////////////////////////////////////////////////////
        ///     \brief  Gets the environment specific communicator.
        ////////////////////////////////////////////////////////
        void*   Get_Communicator();

        ////////////////////////////////////////////////////////
        ///     \brief  Environment specific sum of number of
        ///             columns of every pe.
        ///
        ///     \param my_nbr_cols Number of columns to be
        ///                        preconditioned by this pe.
        ///     \param sum_nbr_cols Sum of my_nbr_cols of all pes.
        ////////////////////////////////////////////////////////
        void    Dist_Scan
                (   int& my_nbr_cols,
                    int& sum_nbr_cols )                 const;

        ////////////////////////////////////////////////////////
        ///     \brief Scatters my_nbr_cols and gathers all
        ///            local work chunks in all_nbr_cols in the
        ///            specific environment.
        ///
        ///     \param my_nbr_cols Number of columns to be
        ///                        preconditioned by this pe.
        ///     \param all_nbr_cols Array which holds
        ///                         my_nbr_cols of all pes.
        ////////////////////////////////////////////////////////
        void    Dist_Local_Chunks
                (   int   my_nbr_cols,
                    int*& all_nbr_cols )                const;

        ////////////////////////////////////////////////////////
        ///     \brief Search for the maximum number of nz over
        ///            all columns of the matrix in this
        ///            environment.
        ///
        ///     The maximum number of nz in a matrix column is
        ///     used to built the maximum necessary receive buffer
        ///     for the communication routines. E.g. to send
        ///     columns from one pe to another.
        ///
        ///     \param max The column of maximum length over all
        ///            pes in this environment.
        ///     \param max_nnz the maximum number of nz over all
        ///            columns in my_nbr_cols on this pe.
        ////////////////////////////////////////////////////////
        void    Dist_Max
                (   int  max,
                    int& max_nnz )                      const;

        ////////////////////////////////////////////////////////
        ///     \brief Distribute the length of all columns to
        ///            all pes in this environment.
        ///
        ///     \param len_cols Length of each column of the
        ///             work chunk on this pe
        ///     \param my_nbr_cols Nunmber of columns on this pe.
        ///     \param len_all_cols Array holding the length of
        ///             each column of the whole matrix.
        ///     \param all_nbr_cols Array which holds
        ///                         my_nbr_cols of all pes.
        ///     \param start_indices Array holding the start
        ///             index of each pe over the global matrix.
        ////////////////////////////////////////////////////////
        void    Dist_Col_Length
                (   int*        len_cols,
                    const int   my_nbr_cols,
                    int*&       len_all_cols,
                    int*        all_nbr_cols,
                    int*        start_indices )         const;

        ////////////////////////////////////////////////////////
        ///     \brief Sums up the nnz over all pes and scatters
        ///            it to all.
        ///
        ///     \param my_nnz The number of nz on the local work
        ///            chunk on this pe.
        ///     \param nnz Reveive buffer holding the sum of all
        ///            nnz on this pe in the environment.
        ////////////////////////////////////////////////////////
        void    Sum_NNZ
                (   int     my_nnz,
                    int&    nnz )                       const;

        ////////////////////////////////////////////////////////
        ///     \brief Dummy function in sequential environment.
        ///
        ///     In parallel environment only real matrices can
        ///     be solved with the HYPRE package. Check for
        ///     correct solver parameter necessary.
        ///
        ///     \param solver_param Value of requested sol
        ///                         parameter.
        ////////////////////////////////////////////////////////
        void    Check_Solver_Type
                (   const int solver_param )            const { }

        ////////////////////////////////////////////////////////
        ///     \brief Testing for correct solver parameter.
        ///
        ///     In sequential environment no other preconditioners
        ///     than FSPAI can be used in the PCG solver.
        ///     Therefore, testing for sol > 3 here, which is
        ///     not allowed in sequential environment.
        ///
        ///     \param solver_param Value of requested sol
        ///                         parameter.
        ////////////////////////////////////////////////////////
        void    Check_Solver_Param
                (   const int solver_param )            const;

	private:

	    /// Number of pes used in this environment
	    int num_procs;

	    /// Id of this pe on this environment
	    int my_id;
};

#	endif /* ENV_HANDLER_H_ */

#endif
