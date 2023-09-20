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

#include "comms.h"

void check_mpi(const int err, const char *msg)
{
}


void setup_comms(struct problem * problem, struct rankinfo * rankinfo)
{
    // Create the MPI Cartesian topology
    unsigned int dimc[] = {problem->npex, problem->npey, problem->npez};
    int *dims = (int *) dimc;
    int periods[] = {0, 0, 0};

    // Get my ranks in x, y and z

    // Note: The following assumes one tile per MPI rank
    // TODO: Change to allow for tiling

    // Calculate rankinfo sizes
    rankinfo->nx = problem->nx / problem->npex;
    rankinfo->ny = problem->ny / problem->npey;
    rankinfo->nz = problem->nz / problem->npez;

    // Calculate i,j,k lower and upper bounds in terms of problem grid
    rankinfo->ilb = rankinfo->ranks[0]     * rankinfo->nx;
    rankinfo->iub = (rankinfo->ranks[0]+1) * rankinfo->nx;
    rankinfo->jlb = rankinfo->ranks[1]     * rankinfo->ny;
    rankinfo->jub = (rankinfo->ranks[1]+1) * rankinfo->ny;
    rankinfo->klb = rankinfo->ranks[2]     * rankinfo->nz;
    rankinfo->kub = (rankinfo->ranks[2]+1) * rankinfo->nz;

    // Calculate neighbouring ranks
    calculate_neighbours(problem, rankinfo); 
}

void finish_comms(void)
{
}

void calculate_neighbours(struct problem * problem, struct rankinfo * rankinfo)
{
    int mpi_err;

    // Calculate my neighbours
    int coords[3];
    // x-dir + 1
    coords[0] = (rankinfo->ranks[0] == problem->npex - 1) ? rankinfo->ranks[0] : rankinfo->ranks[0] + 1;
    coords[1] = rankinfo->ranks[1];
    coords[2] = rankinfo->ranks[2];
    // x-dir - 1
    coords[0] = (rankinfo->ranks[0] == 0) ? rankinfo->ranks[0] : rankinfo->ranks[0] - 1;
    coords[1] = rankinfo->ranks[1];
    coords[2] = rankinfo->ranks[2];
    // y-dir + 1
    coords[0] = rankinfo->ranks[0];
    coords[1] = (rankinfo->ranks[1] == problem->npey - 1) ? rankinfo->ranks[1] : rankinfo->ranks[1] + 1;
    coords[2] = rankinfo->ranks[2];
    // y-dir - 1
    coords[0] = rankinfo->ranks[0];
    coords[1] = (rankinfo->ranks[1] == 0) ? rankinfo->ranks[1] : rankinfo->ranks[1] - 1;
    coords[2] = rankinfo->ranks[2];
    // z-dir + 1
    coords[0] = rankinfo->ranks[0];
    coords[1] = rankinfo->ranks[1];
    coords[2] = (rankinfo->ranks[2] == problem->npez - 1) ? rankinfo->ranks[2] : rankinfo->ranks[2] + 1;
    // z-dir - 1
    coords[0] = rankinfo->ranks[0];
    coords[1] = rankinfo->ranks[1];
    coords[2] = (rankinfo->ranks[2] == 0) ? rankinfo->ranks[2] : rankinfo->ranks[2] - 1;
}


void recv_boundaries(int z_pos, const int octant, const int istep, const int jstep, const int kstep,
    struct problem * problem, struct rankinfo * rankinfo,
    struct buffers * buffers)
{
    int mpi_err;

    // Check if pencil has an external boundary for this sweep direction
    // If so, set as vacuum
    size_t i_offset;
    if (kstep == -1)
    {
        // Correct XY plane position for sweep direction
        int stride = problem->nang*problem->ng*rankinfo->ny;
        i_offset = (rankinfo->nz-problem->chunk-z_pos) * stride;
    }
    else
    {
        i_offset = problem->nang*problem->ng*rankinfo->ny*z_pos;
    }

    if ( (istep == -1 && rankinfo->iub == problem->nx)
        || (istep == 1 && rankinfo->ilb == 0))
    {
        zero_buffer(buffers->flux_i, i_offset, problem->nang*problem->ng*rankinfo->ny*problem->chunk);
    }
    // Otherwise, internal boundary - get data from MPI receives
    else
    {
#if 0
        cl_err = clEnqueueWriteBuffer(context->queue, buffers->flux_i, CL_FALSE,
            sizeof(double)*i_offset,
            sizeof(double)*problem->nang*problem->ng*rankinfo->ny*problem->chunk,
            (buffers->flux_i)+i_offset, 0, NULL, &flux_i_write_event);
        check_ocl(cl_err, "Copying flux i buffer to device");
#endif
//        copy(buffers->flux_i,buffers->flux_i);
size_t nsize = problem->nang*problem->ng*rankinfo->ny*problem->chunk;
double * flux_i = buffers->flux_i;
#pragma omp target update to(flux_i[i_offset:nsize])
    }

    size_t j_offset;
    if (kstep == -1)
    {
        // Correct XY plane position for sweep direction
        int stride = problem->nang*problem->ng*rankinfo->nx;
        j_offset = (rankinfo->nz-problem->chunk-z_pos) * stride;
    }
    else
    {
        j_offset = problem->nang*problem->ng*rankinfo->nx*z_pos;
    }
    if ( (jstep == -1 && rankinfo->jub == problem->ny)
        || (jstep == 1 && rankinfo->jlb == 0))
    {
        zero_buffer(buffers->flux_j, j_offset, problem->nang*problem->ng*rankinfo->nx*problem->chunk);
    }
    else
    {
        // Copy flux_j to the device
#if 0
        cl_err = clEnqueueWriteBuffer(context->queue, buffers->flux_j, CL_FALSE,
            sizeof(double)*j_offset,
            sizeof(double)*problem->nang*problem->ng*rankinfo->nx*problem->chunk,
            (buffers->flux_j)+j_offset, 0, NULL, &flux_j_write_event);
        check_ocl(cl_err, "Copying flux j buffer to device");
#endif
//        copy(buffers->flux_j,buffers->flux_j);
size_t nsize = problem->nang*problem->ng*rankinfo->nx*problem->chunk;
double * flux_j = buffers->flux_j;
#pragma omp target update to(flux_j[j_offset:nsize])
    }
}


void send_boundaries(int z_pos, const int octant, const int istep, const int jstep, const int kstep,
    struct problem * problem, struct rankinfo * rankinfo,
    struct buffers * buffers)
{
    int mpi_err;


    // Get the edges off the device
    // I
    size_t i_offset;
    if (kstep == -1)
    {
        // Correct XY plane position for sweep direction
        int stride = problem->nang*problem->ng*rankinfo->ny;
        i_offset = (rankinfo->nz-problem->chunk-z_pos) * stride;
    }
    else
    {
        i_offset = problem->nang*problem->ng*rankinfo->ny*z_pos;
    }
#if 0
    cl_err = clEnqueueReadBuffer(context->queue, buffers->flux_i, CL_FALSE,
        sizeof(double)*i_offset,
        sizeof(double)*problem->nang*problem->ng*rankinfo->ny*problem->chunk,
        (buffers->flux_i)+i_offset, 0, NULL, &flux_i_read_event);
    check_ocl(cl_err, "Copying flux i buffer back to host");
#endif
//        copy(buffers->flux_i,buffers->flux_i);
int nsize = problem->nang*problem->ng*rankinfo->ny*problem->chunk;
double * flux_i = buffers->flux_i;
#pragma omp target update from(flux_i[i_offset:nsize])

    // J
    size_t j_offset;
    if (kstep == -1)
    {
        // Correct XY plane position for sweep direction
        int stride = problem->nang*problem->ng*rankinfo->nx;
        j_offset = (rankinfo->nz-problem->chunk-z_pos) * stride;
    }
    else
    {
        j_offset = problem->nang*problem->ng*rankinfo->nx*z_pos;
    }
#if 0
    cl_err = clEnqueueReadBuffer(context->queue, buffers->flux_j, CL_TRUE,
        sizeof(double)*j_offset,
        sizeof(double)*problem->nang*problem->ng*rankinfo->nx*problem->chunk,
        (buffers->flux_j)+j_offset, 0, NULL, &flux_j_read_event);
    check_ocl(cl_err, "Copying flux j buffer back to host");
#endif
//        copy(buffers->flux_j,buffers->flux_j);
nsize = problem->nang*problem->ng*rankinfo->nx*problem->chunk;
double * flux_j = buffers->flux_j;
#pragma omp target update from(flux_j[j_offset:nsize])

    double tick = wtime();

    // Send to neighbour with MPI_Send
    // X
    // Y
    sweep_mpi_time += wtime() - tick;
}


