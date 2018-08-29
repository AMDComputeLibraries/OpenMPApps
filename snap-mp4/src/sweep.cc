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

#include "sweep.h"

void init_planes(struct plane** planes, unsigned int *num_planes, struct problem * problem, struct rankinfo * rankinfo)
{
    *num_planes = rankinfo->nx + rankinfo->ny + problem->chunk - 2;
    *planes = (struct plane *) malloc(sizeof(struct plane) * *num_planes);

    for (unsigned int p = 0; p < *num_planes; p++)
        (*planes)[p].num_cells = 0;

    for (unsigned int k = 0; k < problem->chunk; k++)
        for (unsigned int j = 0; j < rankinfo->ny; j++)
            for (unsigned int i = 0; i < rankinfo->nx; i++)
            {
                unsigned int p = i + j + k;
                (*planes)[p].num_cells += 1;
            }

    for (unsigned int p = 0; p < *num_planes; p++)
    {
        (*planes)[p].cell_ids = (struct cell_id *) malloc(sizeof(struct cell_id) * (*planes)[p].num_cells);
    }

    unsigned int index[*num_planes];
    for (unsigned int p = 0; p < *num_planes; p++)
        index[p] = 0;

    for (unsigned int k = 0; k < problem->chunk; k++)
        for (unsigned int j = 0; j < rankinfo->ny; j++)
            for (unsigned int i = 0; i < rankinfo->nx; i++)
            {
                unsigned int p = i + j + k;
                (*planes)[p].cell_ids[index[p]].i = i;
                (*planes)[p].cell_ids[index[p]].j = j;
                (*planes)[p].cell_ids[index[p]].k = k;
                index[p] += 1;
            }
}

#define SOURCE_INDEX(m,g,i,j,k,cmom,ng,nx,ny) ((m)+((cmom)*(g))+((cmom)*(ng)*(i))+((cmom)*(ng)*(nx)*(j))+((cmom)*(ng)*(nx)*(ny)*(k)))
#define SCAT_COEFF_INDEX(a,l,o,nang,cmom) ((a)+((nang)*(l))+((nang)*(cmom)*o))
#define FLUX_I_INDEX(a,g,j,k,nang,ng,ny) ((a)+((nang)*(g))+((nang)*(ng)*(j))+((nang)*(ng)*(ny)*(k)))
#define FLUX_J_INDEX(a,g,i,k,nang,ng,nx) ((a)+((nang)*(g))+((nang)*(ng)*(i))+((nang)*(ng)*(nx)*(k)))
#define FLUX_K_INDEX(a,g,i,j,nang,ng,nx) ((a)+((nang)*(g))+((nang)*(ng)*(i))+((nang)*(ng)*(nx)*(j)))
#define ANGULAR_FLUX_INDEX(a,g,i,j,k,nang,ng,nx,ny) ((a)+((nang)*(g))+((nang)*(ng)*(i))+((nang)*(ng)*(nx)*(j))+((nang)*(ng)*(nx)*(ny)*(k)))
#define DENOMINATOR_INDEX(a,g,i,j,k,nang,ng,nx,ny) ((a)+((nang)*(g))+((nang)*(ng)*(i))+((nang)*(ng)*(nx)*(j))+((nang)*(ng)*(nx)*(ny)*(k)))

#define source(m,g,i,j,k) source[SOURCE_INDEX((m),(g),(i),(j),(k),cmom,ng,nx,ny)]
#define scat_coeff(a,l,o) scat_coeff[SCAT_COEFF_INDEX((a),(l),(o),nang,cmom)]
#define flux_i(a,g,j,k) flux_i[FLUX_I_INDEX((a),(g),(j),(k),nang,ng,ny)]
#define flux_j(a,g,i,k) flux_j[FLUX_J_INDEX((a),(g),(i),(k),nang,ng,nx)]
#define flux_k(a,g,i,j) flux_k[FLUX_K_INDEX((a),(g),(i),(j),nang,ng,nx)]
#define angular_flux_in(a,g,i,j,k) angular_flux_in[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_out(a,g,i,j,k) angular_flux_out[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define denominator(a,g,i,j,k) denominator[DENOMINATOR_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]

void sweep_plane(
    const unsigned int z_pos,
    const int octant,
    const int istep,
    const int jstep,
    const int kstep,
    const unsigned int this_plane,
    const struct plane * planes,
    struct problem * problem,
    struct rankinfo * rankinfo,
    struct buffers * buffers
    )
{
    unsigned int nang = problem->nang;
    unsigned int nx = rankinfo->nx;
    unsigned int ny = rankinfo->ny;
    unsigned int nz = rankinfo->nz;
    unsigned int ng = problem->ng;
    unsigned int cmom = problem->cmom;

    size_t cgxyz = cmom*ng*nx*ny*nz;
    size_t agxyz = ng*nang*nx*ny*nz;

    struct cell_id * plane =  buffers->planes[this_plane];
    size_t num_cells = planes[this_plane].num_cells;
    size_t ncg = num_cells*ng;

    double * source =  buffers->inner_source;
    double * scat_coeff =  buffers->scat_coeff;
    double * dd_i =  buffers->dd_i;
    double * dd_j =  buffers->dd_j;
    double * dd_k =  buffers->dd_k;
    double * mu =  buffers->mu;

    double * velocity_delta =  buffers->velocity_delta;
    double * mat_cross_section =  buffers->mat_cross_section;
    double * denominator =  buffers->denominator;
    double * angular_flux_in =  buffers->angular_flux_in[octant];
    double * flux_i =  buffers->flux_i;
    double * flux_j =  buffers->flux_j;
    double * flux_k =  buffers->flux_k;
    double * angular_flux_out =  buffers->angular_flux_out[octant];

    // 2 dimensional kernel
    // First dimension: number of angles * number of groups
    // Second dimension: number of cells in plane

#pragma omp target data map(to: plane[:1], \
                            source[:cgxyz], \
                            scat_coeff[:nang*cmom*8], \
                            dd_i[:1], \
                            dd_j[:nang], \
                            dd_k[:nang], \
                            mu[:nang], \
                            velocity_delta[:ng], \
                            mat_cross_section[:ng], \
                            denominator[:agxyz], \
                            angular_flux_in[:agxyz], \
                            z_pos,istep,jstep,kstep, \
                            ncg,num_cells,nx,ny,nz,ng,nang,cmom,octant) \
                        map(from: angular_flux_out[:agxyz]) \
                        map(tofrom: flux_i[:ng*nang*ny*nz], \
                            flux_j[:ng*nang*nx*nz], \
                            flux_k[:ng*nang*nx*ny])

#pragma omp target teams distribute num_teams(num_cells) thread_limit(ng*nang)
    for (int cidx = 0; cidx < num_cells; cidx++)
    {
#ifndef NO_COLLAPSE
#pragma omp parallel for collapse(2)
       for(int g = 0; g<ng; g++)
       {
          for(int a = 0; a<nang; a++)
          {
#else
// The collapse clause above works fine.  This just shows how one could avoid
// the collapse clause by manually collapsing the above two loops.
#pragma omp parallel for
       for (int tid=0; tid<(ng*nang); tid++)
       {
           int g = tid / nang ;
           int a = tid % nang ;
#endif
    // Read cell index from plane buffer
    const size_t i = (istep > 0) ? plane[cidx].i         : nx - plane[cidx].i         - 1;
    const size_t j = (jstep > 0) ? plane[cidx].j         : ny - plane[cidx].j         - 1;
    const size_t k = (kstep > 0) ? plane[cidx].k + z_pos : nz - plane[cidx].k - z_pos - 1;
    //
    // Compute the angular flux (psi)
    //

    // Begin with the first scattering moment
    double source_term = source(0,g,i,j,k);

    // Add in the anisotropic scattering source moments
    for (unsigned int l = 1; l < cmom; l++)
    {
        source_term += scat_coeff(a,l,octant) * source(l,g,i,j,k);
    }

    double psi =
        source_term
        + flux_i(a,g,j,k)*mu[a]*dd_i[0]
        + flux_j(a,g,i,k)*dd_j[a]
        + flux_k(a,g,i,j)*dd_k[a];
    // Add contribution from last timestep flux if time-dependant
    if (velocity_delta[g] != 0.0)
    {
        psi += velocity_delta[g] * angular_flux_in(a,g,i,j,k);
    }

    // "Divide" by denominator
    psi *= denominator(a,g,i,j,k);

    // Compute upwind fluxes
    double tmp_flux_i = 2.0 * psi - flux_i(a,g,j,k);
    double tmp_flux_j = 2.0 * psi - flux_j(a,g,i,k);
    double tmp_flux_k = 2.0 * psi - flux_k(a,g,i,j);

    // Time difference the final flux value
    if (velocity_delta[g] != 0.0)
    {
        psi = 2.0 * psi - angular_flux_in(a,g,i,j,k);
    }

    // Fixup
    double zeros[4];
    int num_ok = 4;
    for (int fix = 0; fix < 4; fix++)
    {
        zeros[0] = (tmp_flux_i < 0.0) ? 0.0 : 1.0;
        zeros[1] = (tmp_flux_j < 0.0) ? 0.0 : 1.0;
        zeros[2] = (tmp_flux_k < 0.0) ? 0.0 : 1.0;
        zeros[3] = (psi < 0.0)        ? 0.0 : 1.0;

        if (num_ok == zeros[0] + zeros[1] + zeros[2] + zeros[3])
            continue;

        num_ok = zeros[0] + zeros[1] + zeros[2] + zeros[3];

        // Recalculate psi
        psi =
            flux_i(a,g,j,k)*mu[a]*dd_i[0]*(1.0 + zeros[0]) +
            flux_j(a,g,i,k)*dd_j[a]*(1.0 + zeros[1]) +
            flux_k(a,g,i,j)*dd_k[a]*(1.0 + zeros[2]);

        if (velocity_delta[g] != 0.0)
        {
            psi += velocity_delta[g] * angular_flux_in(a,g,i,j,k) * (1.0 + zeros[3]);
        }

        psi = 0.5 * psi + source_term;

        double new_denominator =
            mat_cross_section[g] +
            mu[a] * dd_i[0] * zeros[0] +
            dd_j[a] * zeros[1] +
            dd_k[a] * zeros[2] +
            velocity_delta[g] * zeros[3];
        if (new_denominator > 1.0E-12)
        {
            psi /= new_denominator;
        }
        else
        {
            psi = 0.0;
        }

        tmp_flux_i = 2.0 * psi - flux_i(a,g,j,k);
        tmp_flux_j = 2.0 * psi - flux_j(a,g,i,k);
        tmp_flux_k = 2.0 * psi - flux_k(a,g,i,j);

        if (velocity_delta[g] != 0.0)
        {
            psi = 2.0 * psi - angular_flux_in(a,g,i,j,k);
        }

    }

    // Write values to global memory
    flux_i(a,g,j,k) = tmp_flux_i * zeros[0];
    flux_j(a,g,i,k) = tmp_flux_j * zeros[1];
    flux_k(a,g,i,j) = tmp_flux_k * zeros[2];
    angular_flux_out(a,g,i,j,k) = psi * zeros[3];
#ifdef NO_COLLAPSE
        } // end of for each tid
#else
// end of loops for collapse clause
            } // end of for each a
        } // end of for each g
#endif
    }  // end of for each cidx

}

