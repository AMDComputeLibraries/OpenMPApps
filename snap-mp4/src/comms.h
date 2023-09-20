/*******************************************************************************
Copyright (c) 2016 Advanced Micro Devices, Inc.

All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************/

#pragma once

/** \file
* \brief Communication routines (setup and sweep data transfer)
*/

#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "buffers.h"

#include "profiler.h"


/** \brief Check MPI error codes */
void check_mpi(const int err, const char *msg);

/** \brief Setup cartesian communicator and find your MPI rank */
void setup_comms(struct problem * problem, struct rankinfo * rankinfo);

/** \brief Just call MPI_Finalize */
void finish_comms(void);

void calculate_neighbours(struct problem * problem, struct rankinfo * rankinfo); 

/** \brief Receive chunk number of XY planes starting at position z_pos */
void recv_boundaries(int z_pos, const int octant, const int istep, const int jstep, const int kstep, struct problem * problem, struct rankinfo * rankinfo, struct buffers * buffers);
/** \brief Send chunk number of XY planes starting at position z_pos */
void send_boundaries(int z_pos, const int octant, const int istep, const int jstep, const int kstep, struct problem * problem, struct rankinfo * rankinfo, struct buffers * buffers);
