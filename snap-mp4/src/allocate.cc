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

#include <stdlib.h>

#include "global.h"
#include "problem.h"
#include "buffers.h"

void allocate_buffers(struct problem * problem, struct rankinfo * rankinfo, struct buffers * buffers)
{
    // Angular flux arrays
    size_t nsize = problem->nang*problem->ng*rankinfo->nx*rankinfo->ny*rankinfo->nz;
    for (int i=0;i<8;i++)
    {
       double * angular_flux_in = (double *) malloc(sizeof(double)*nsize);
       buffers->angular_flux_in[i] = angular_flux_in;
#pragma omp target enter data map(to: angular_flux_in[:nsize])
       double * angular_flux_out = (double *) malloc(sizeof(double)*nsize);
       buffers->angular_flux_out[i] = angular_flux_out;
#pragma omp target enter data map(to: angular_flux_out[:nsize])
    }

    // Allocate edge flux arrays
    nsize = problem->nang*problem->ng*rankinfo->ny*rankinfo->nz;
    double * flux_i = (double *) malloc(sizeof(double)*nsize);
    buffers->flux_i = flux_i;
#pragma omp target enter data map(to: flux_i[:nsize])
    nsize = problem->nang*problem->ng*rankinfo->nx*rankinfo->nz;
    double * flux_j = (double *) malloc(sizeof(double)*nsize);
    buffers->flux_j = flux_j;
#pragma omp target enter data map(to: flux_j[:nsize])
    nsize = problem->nang*problem->ng*rankinfo->nx*rankinfo->ny;
    double * flux_k = (double *) malloc(sizeof(double)*nsize);
    buffers->flux_k = flux_k;
#pragma omp target enter data map(to: flux_k[:nsize])

    // Scalar flux
    // grid * ng
    nsize = rankinfo->nx*rankinfo->ny*rankinfo->nz*problem->ng;
    double * scalar_flux = (double *) malloc(sizeof(double)*nsize);
#pragma omp target enter data map(to: scalar_flux[:nsize])
    buffers->scalar_flux = scalar_flux;
    buffers->old_inner_scalar_flux = (double *) malloc(sizeof(double)*nsize);
    buffers->old_outer_scalar_flux = (double *) malloc(sizeof(double)*nsize);

    //Scalar flux moments
    if (problem->cmom-1 == 0)
        nsize = 1;
    else
        nsize = (problem->cmom-1)*problem->ng*rankinfo->nx*rankinfo->ny*rankinfo->nz;
    double * scalar_flux_moments= (double *) calloc(nsize,sizeof(double));
#pragma omp target enter data map(to: scalar_flux_moments[:nsize])
    buffers->scalar_flux_moments = scalar_flux_moments;

    // Weights and cosines
    nsize = problem->nang;
    buffers->quad_weights = (double *) malloc(sizeof(double)*nsize);
    buffers->mu = (double *) malloc(sizeof(double)*nsize);
    buffers->eta = (double *) malloc(sizeof(double)*nsize);
    buffers->xi = (double *) malloc(sizeof(double)*nsize);

    // Scattering coefficient
    nsize = problem->nang*problem->cmom*8;
    buffers->scat_coeff = (double *) malloc(sizeof(double)*nsize);

    // Material cross section
    nsize = problem->ng;
    buffers->mat_cross_section= (double *) malloc(sizeof(double)*nsize);

    // Source terms
    nsize = problem->ng*rankinfo->nx*rankinfo->ny*rankinfo->nz;
    buffers->fixed_source = (double *) malloc(sizeof(double)*nsize);
    nsize = problem->cmom*problem->ng*rankinfo->nx*rankinfo->ny*rankinfo->nz;
    double * outer_source = (double *) malloc(sizeof(double)*nsize);
    buffers->outer_source = outer_source;
#pragma omp target enter data map(to: outer_source[:nsize])
    double * inner_source = (double *) malloc(sizeof(double)*nsize);
    buffers->inner_source = inner_source;
#pragma omp target enter data map(to: inner_source[:nsize])

    // Scattering terms
    nsize = problem->nmom*problem->ng*problem->ng;
    buffers->scattering_matrix = (double *) malloc(sizeof(double)*nsize);

    // Diamond diference co-efficients
    nsize = 1;
    buffers->dd_i = (double *) malloc(sizeof(double)*nsize);
    nsize = problem->nang;
    buffers->dd_j = (double *) malloc(sizeof(double)*nsize);
    buffers->dd_k = (double *) malloc(sizeof(double)*nsize);

    // Velocities
    nsize = problem->ng;
    buffers->velocities = (double *) malloc(sizeof(double)*nsize);
    buffers->velocity_delta = (double *) malloc(sizeof(double)*nsize);

    // Denominator array
    nsize = problem->nang*problem->ng*rankinfo->nx*rankinfo->ny*rankinfo->nz;
    buffers->denominator = (double *) malloc(sizeof(double)*nsize);
}

void free_buffers(struct buffers * buffers)
{
    for (int i=0;i<8;i++)
    {
       free(buffers->angular_flux_in[i]);
       free(buffers->angular_flux_out[i]);
    }
    free(buffers->flux_i);
    free(buffers->flux_j);
    free(buffers->flux_k);
    free(buffers->scalar_flux);
    free(buffers->old_inner_scalar_flux);
    free(buffers->old_outer_scalar_flux);
    free(buffers->scalar_flux_moments);
    free(buffers->quad_weights);
    free(buffers->mu);
    free(buffers->eta);
    free(buffers->xi);
    free(buffers->scat_coeff);
    free(buffers->mat_cross_section);
    free(buffers->fixed_source);
    free(buffers->outer_source);
    free(buffers->inner_source);
    free(buffers->scattering_matrix);
    free(buffers->dd_i);
    free(buffers->dd_j);
    free(buffers->dd_k);
    free(buffers->velocities);
    free(buffers->velocity_delta);
    free(buffers->denominator);
}
