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

#ifndef TIMER_H_
#define TIMER_H_

//file includings
#include "ENV_Handler.h"

//C/C++ includings
#include <sys/time.h>
#include <ctime>

/// For debugging - this is the number
/// of debug timers you may want to use
const int dbg_timers = 4;

//////////////////////////////////////////////
///     \class Timer
///     \brief This class is responsible for
///            the time measurement
///
///     This class provides some useful
///     methods for time measurement. Each
///     pe will sum its time within his own
///     time arrays and after all pes
///     finished their work, the maximum
///     time for one operation will be found
///     out throgh MPI communication between
///     all pes. Furthermore this class
///     provides some debug features to
///     show the time of a specific operation
///     on each pe.
//////////////////////////////////////////////
class Timer
{
    public:

        /// Empty Constructor
        Timer();

        /// Destructor
        ~Timer() { };

        // Member variables

        /// For debugging - holding the current start
        /// time of all started time measurements
        static double   dbg_start_timers[dbg_timers];

        /// For debugging - holding the current stop
        /// time of all started time measurements
        static double   dbg_stop_timers[dbg_timers];

        /// For debugging - holding the current sum
        /// of time of all previously performed time
        /// measurements
        static double   dbg_sum_timers[dbg_timers];

        // Methods

        //////////////////////////////////////////////
        ///     \brief  Starting time for time
        ///             measurement
        ///     \param env_handler Environment interface
        ///////////////////////////////////////////////
        void    Start( ENV_Handler& env_handler );

        //////////////////////////////////////////////
        ///     \brief  Stopping time for time
        ///             measurement
        ///
        ///     Takes the start time, subtracts
        ///     stop - start and adds the result to
        ///     sum_time. Notice that with this method
        ///     no overlapping time measurement can be
        ///     done.
        ///
        ///     \param env_handler Environment interface
        ///////////////////////////////////////////////
        void    Stop( ENV_Handler& env_handler );

        //////////////////////////////////////////////
        ///     \brief  Streaming time of time
        ///             measurement to shell
        ///
        ///     To report the correct time of an
        ///     operation within MPI environment, all
        ///     pes have to exchange their sum of times.
        ///     The maxmum value is the correct time
        ///     measurement.
        ///
        ///     \param env_handler Environment interface
        ///////////////////////////////////////////////
        void    Report( ENV_Handler& env_handler );

        //////////////////////////////////////////////
        ///     \brief  For debugging - Starting time
        ///             for time measurement
        ///
        ///     To get a time measurement of a specific
        ///     operation on each pe the debug arrays
        ///     will be filled with the sum of time.
        ///     Notice that this is only for debugging
        ///     and makes the program a much more
        ///     slower.
        ///
        ///     \param id Id of operation to start the
        ///               time for
        ///     \param env_handler Environment interface
        ///////////////////////////////////////////////
        void    Dbg_Start
                (   int          id,
                    ENV_Handler& env_handler );

        //////////////////////////////////////////////
        ///     \brief  For debugging - Stopping time
        ///             for time measurement
        ///
        ///     \param id Id of operation to stop the
        ///               time for
        ///     \param env_handler Environment interface
        ///////////////////////////////////////////////
        void    Dbg_Stop
                (   int          id,
                    ENV_Handler& env_handler );

        //////////////////////////////////////////////
        ///     \brief   Streaming time of time
        ///              measurement to shell
        ///
        ///     For each pe there is the sum of time
        ///     it spent within the measured operations.
        ///
        ///     \param env_handler Environment interface
        ///////////////////////////////////////////////
        void    Dbg_Report( ENV_Handler& env_handler ) const;

    private:

        /// Current start time of one
        /// time measurement
        double  start_timer;

        /// Current stop time of one
        /// time measurement
        double  stop_timer;

        /// Current sum of time of one
        /// time measurement
        double  sum_time;

        struct  timespec tp;
};

#endif
