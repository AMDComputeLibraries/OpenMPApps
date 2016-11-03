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

#include "buffers.h"
#include <iostream>
#include <string>

#include <stdio.h>
#include <omp.h>

void check_device_memory_requirements(
    struct problem * problem, struct rankinfo * rankinfo)
{
#if 1
//OpenMP does not have routines to query device name and memory sizes
// just hard code this info for now
    std::wcout << "device : Quadro K4000" <<  std::endl;
    unsigned long global = 3ul*1024ul*1024ul*1024ul; // for Quadro K4000 (t0)
    std::wcout << "device memory: " << global <<  std::endl;
#else

    hc::accelerator acc;
    std::wstring device_name(acc.get_description());
    std::wcout << "device : " << device_name << std::endl;

    unsigned long global = acc.get_dedicated_memory();
//  for hsa - accelerator.get_dedicated_memory()
    printf("device memory: %ld\n",global);
    size_t tsmem = acc.get_max_tile_static_size();
    printf("tile static memory: %ld\n",tsmem);
#endif

    unsigned long total = 0;
    // Add up the memory requirements, in bytes.
    total += problem->nang*problem->ng*rankinfo->nx*rankinfo->ny*rankinfo->nz*8;
    total += problem->nang*problem->ng*rankinfo->nx*rankinfo->ny*rankinfo->nz*8;
    total += problem->nang*problem->ng*rankinfo->ny*rankinfo->nz;
    total += problem->nang*problem->ng*rankinfo->nx*rankinfo->nz;
    total += problem->nang*problem->ng*rankinfo->nx*rankinfo->ny;
    total += problem->ng*rankinfo->nx*rankinfo->ny*rankinfo->nz;
    if (problem->cmom-1 == 0)
        total += problem->ng*rankinfo->nx*rankinfo->ny*rankinfo->nz;
    else
        total += 1;
    total += problem->nang;
    total += problem->nang;
    total += problem->nang;
    total += problem->nang;
    total += problem->nang*problem->cmom*8;
    total += problem->ng;
    total += problem->ng*rankinfo->nx*rankinfo->ny*rankinfo->nz;
    total += problem->cmom*problem->ng*rankinfo->nx*rankinfo->ny*rankinfo->nz;
    total += problem->cmom*problem->ng*rankinfo->nx*rankinfo->ny*rankinfo->nz;
    total += problem->nmom*problem->ng*problem->ng;
    total += 1;
    total += problem->nang;
    total += problem->nang;
    total += problem->ng;
    total += problem->ng;
    total += problem->nang*problem->ng*rankinfo->nx*rankinfo->ny*rankinfo->nz;
    total *= sizeof(double);

    if (global < total)
    {
        fprintf(stderr,"Error: Device does not have enough global memory.\n");
        fprintf(stderr, "Required: %.1f GB\n", (double)total/(1024.0*1024.0*1024.0));
        fprintf(stderr, "Available: %.1f GB\n", (double)global/(1024.0*1024.0*1024.0));
        exit(EXIT_FAILURE);
    }
}


void zero_buffer(double *pbuffer, size_t offset, size_t size)
{
#pragma omp target data map(to: pbuffer[offset:size+offset],offset,size)
#pragma omp target 
#pragma omp parallel for
    for(int i = 0; i<size; i++)
    {
       pbuffer[i+offset] = 0.0;
    }
}

void set_buffer(double *pbuffer, double value, size_t size)
{
#pragma omp target data map(to: pbuffer[:size],value,size)
#pragma omp target teams
#pragma omp parallel for
    for(int i = 0; i<size; i++)
    {
       pbuffer[i] = value;
    }
}

void swap_angular_flux_buffers(struct buffers * buffers)
{
    double *tmp;
    for (int i=0;i<8;i++)
    {
       tmp = buffers->angular_flux_in[i];
       buffers->angular_flux_in[i] = buffers->angular_flux_out[i];
       buffers->angular_flux_out[i] = tmp;
    }
}
