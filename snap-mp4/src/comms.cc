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
MPI_Comm snap_comms;

void check_mpi(const int err, const char *msg)
{
    if (err != MPI_SUCCESS)
    {
        fprintf(stderr, "MPI Error: %d. %s\n", err, msg);
        exit(EXIT_FAILURE);
    }
}


void setup_comms(struct problem * problem, struct rankinfo * rankinfo)
{
    // Create the MPI Cartesian topology
    unsigned int dimc[] = {problem->npex, problem->npey, problem->npez};
    int *dims = (int *) dimc;
    int periods[] = {0, 0, 0};
    unsigned int mpi_err = MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &snap_comms);
    check_mpi(mpi_err, "Creating MPI Cart");

    // Get my ranks in x, y and z
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &rankinfo->rank);
    check_mpi(mpi_err, "Getting MPI rank");
    mpi_err = MPI_Cart_coords(snap_comms, rankinfo->rank, 3, rankinfo->ranks);
    check_mpi(mpi_err, "Getting Cart co-ordinates");

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
    calculate_neighbours(snap_comms, problem, rankinfo);
}

void finish_comms(void)
{
    int mpi_err = MPI_Finalize();
    check_mpi(mpi_err, "MPI_Finalize");
}

void calculate_neighbours(MPI_Comm comms,  struct problem * problem, struct rankinfo * rankinfo)
{
    int mpi_err;

    // Calculate my neighbours
    int coords[3];
    // x-dir + 1
    coords[0] = (rankinfo->ranks[0] == problem->npex - 1) ? rankinfo->ranks[0] : rankinfo->ranks[0] + 1;
    coords[1] = rankinfo->ranks[1];
    coords[2] = rankinfo->ranks[2];
    mpi_err = MPI_Cart_rank(comms, coords, &rankinfo->xup);
    check_mpi(mpi_err, "Getting x+1 rank");
    // x-dir - 1
    coords[0] = (rankinfo->ranks[0] == 0) ? rankinfo->ranks[0] : rankinfo->ranks[0] - 1;
    coords[1] = rankinfo->ranks[1];
    coords[2] = rankinfo->ranks[2];
    mpi_err = MPI_Cart_rank(comms, coords, &rankinfo->xdown);
    check_mpi(mpi_err, "Getting x-1 rank");
    // y-dir + 1
    coords[0] = rankinfo->ranks[0];
    coords[1] = (rankinfo->ranks[1] == problem->npey - 1) ? rankinfo->ranks[1] : rankinfo->ranks[1] + 1;
    coords[2] = rankinfo->ranks[2];
    mpi_err = MPI_Cart_rank(comms, coords, &rankinfo->yup);
    check_mpi(mpi_err, "Getting y+1 rank");
    // y-dir - 1
    coords[0] = rankinfo->ranks[0];
    coords[1] = (rankinfo->ranks[1] == 0) ? rankinfo->ranks[1] : rankinfo->ranks[1] - 1;
    coords[2] = rankinfo->ranks[2];
    mpi_err = MPI_Cart_rank(comms, coords, &rankinfo->ydown);
    check_mpi(mpi_err, "Getting y-1 rank");
    // z-dir + 1
    coords[0] = rankinfo->ranks[0];
    coords[1] = rankinfo->ranks[1];
    coords[2] = (rankinfo->ranks[2] == problem->npez - 1) ? rankinfo->ranks[2] : rankinfo->ranks[2] + 1;
    mpi_err = MPI_Cart_rank(comms, coords, &rankinfo->zup);
    check_mpi(mpi_err, "Getting z+1 rank");
    // z-dir - 1
    coords[0] = rankinfo->ranks[0];
    coords[1] = rankinfo->ranks[1];
    coords[2] = (rankinfo->ranks[2] == 0) ? rankinfo->ranks[2] : rankinfo->ranks[2] - 1;
    mpi_err = MPI_Cart_rank(comms, coords, &rankinfo->zdown);
    check_mpi(mpi_err, "Getting z-1 rank");
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
        if (istep == -1)
        {
            mpi_err = MPI_Recv(buffers->flux_i+i_offset, problem->nang*problem->ng*rankinfo->ny*problem->chunk, MPI_DOUBLE,
                rankinfo->xup, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            check_mpi(mpi_err, "Receiving from upward x neighbour");
        }
        else
        {
            mpi_err = MPI_Recv(buffers->flux_i+i_offset, problem->nang*problem->ng*rankinfo->ny*problem->chunk, MPI_DOUBLE,
                rankinfo->xdown, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            check_mpi(mpi_err, "Receiving from downward x neighbour");
        }
        // Copy flux_i to the device
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
        if (jstep == -1)
        {
            mpi_err = MPI_Recv(buffers->flux_j+j_offset, problem->nang*problem->ng*rankinfo->nx*problem->chunk, MPI_DOUBLE,
                rankinfo->yup, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            check_mpi(mpi_err, "Receiving from upward y neighbour");
        }
        else
        {
            mpi_err = MPI_Recv(buffers->flux_j+j_offset, problem->nang*problem->ng*rankinfo->nx*problem->chunk, MPI_DOUBLE,
                rankinfo->ydown, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            check_mpi(mpi_err, "Receiving from downward y neighbour");
        }
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
    if (istep == -1 && rankinfo->xdown != rankinfo->rank)
    {
        mpi_err = MPI_Send(buffers->flux_i+i_offset, problem->nang*problem->ng*rankinfo->ny*problem->chunk, MPI_DOUBLE,
                rankinfo->xdown, 0, MPI_COMM_WORLD);
        check_mpi(mpi_err, "Sending to downward x neighbour");
    }
    else if (istep == 1 && rankinfo->xup != rankinfo->rank)
    {
        mpi_err = MPI_Send(buffers->flux_i+i_offset, problem->nang*problem->ng*rankinfo->ny*problem->chunk, MPI_DOUBLE,
                rankinfo->xup, 0, MPI_COMM_WORLD);
        check_mpi(mpi_err, "Sending to upward x neighbour");
    }
    // Y
    if (jstep == -1 && rankinfo->ydown != rankinfo->rank)
    {
        mpi_err = MPI_Send(buffers->flux_j+j_offset, problem->nang*problem->ng*rankinfo->nx*problem->chunk, MPI_DOUBLE,
                rankinfo->ydown, 0, MPI_COMM_WORLD);
        check_mpi(mpi_err, "Sending to downward y neighbour");
    }
    else if (jstep == 1 && rankinfo->yup != rankinfo->rank)
    {
        mpi_err = MPI_Send(buffers->flux_j+j_offset, problem->nang*problem->ng*rankinfo->nx*problem->chunk, MPI_DOUBLE,
                rankinfo->yup, 0, MPI_COMM_WORLD);
        check_mpi(mpi_err, "Sending to upward y neighbour");
    }
    sweep_mpi_time += wtime() - tick;
}


