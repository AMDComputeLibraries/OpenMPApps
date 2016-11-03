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

#include "scalar_flux.h"
#include <string.h>
#include <iostream>
#include <stdio.h>

#define ANGULAR_FLUX_INDEX(a,g,i,j,k,nang,ng,nx,ny) ((a)+((nang)*(g))+((nang)*(ng)*(i))+((nang)*(ng)*(nx)*(j))+((nang)*(ng)*(nx)*(ny)*(k)))
#define SCALAR_FLUX_INDEX(g,i,j,k,ng,nx,ny) ((g)+((ng)*(i))+((ng)*(nx)*(j))+((ng)*(nx)*(ny)*(k)))


#define angular_flux_in_0(a,g,i,j,k) angular_flux_in_0[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_in_1(a,g,i,j,k) angular_flux_in_1[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_in_2(a,g,i,j,k) angular_flux_in_2[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_in_3(a,g,i,j,k) angular_flux_in_3[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_in_4(a,g,i,j,k) angular_flux_in_4[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_in_5(a,g,i,j,k) angular_flux_in_5[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_in_6(a,g,i,j,k) angular_flux_in_6[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_in_7(a,g,i,j,k) angular_flux_in_7[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_out_0(a,g,i,j,k) angular_flux_out_0[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_out_1(a,g,i,j,k) angular_flux_out_1[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_out_2(a,g,i,j,k) angular_flux_out_2[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_out_3(a,g,i,j,k) angular_flux_out_3[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_out_4(a,g,i,j,k) angular_flux_out_4[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_out_5(a,g,i,j,k) angular_flux_out_5[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_out_6(a,g,i,j,k) angular_flux_out_6[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define angular_flux_out_7(a,g,i,j,k) angular_flux_out_7[ANGULAR_FLUX_INDEX((a),(g),(i),(j),(k),nang,ng,nx,ny)]
#define scalar_flux(g,i,j,k) scalar_flux[SCALAR_FLUX_INDEX((g),(i),(j),(k),ng,nx,ny)]
void compute_scalar_flux(
    struct problem * problem,
    struct rankinfo * rankinfo,
    struct buffers * buffers
    )
{

    // get closest power of 2 to nang
    size_t power = 1 << (unsigned int)ceil(log2((double)problem->nang));

    unsigned int nang = problem->nang;
    unsigned int nx = rankinfo->nx;
    unsigned int ny = rankinfo->ny;
    unsigned int nz = rankinfo->nz;
    unsigned int ng = problem->ng;

    size_t nxyz = nx*ny*nz;
    size_t agxyz = ng*nang*nx*ny*nz;
    size_t gxyz = ng*nx*ny*nz;

    double * angular_flux_out_0 =  buffers->angular_flux_out[0];
    double * angular_flux_out_1 =  buffers->angular_flux_out[1];
    double * angular_flux_out_2 =  buffers->angular_flux_out[2];
    double * angular_flux_out_3 =  buffers->angular_flux_out[3];
    double * angular_flux_out_4 =  buffers->angular_flux_out[4];
    double * angular_flux_out_5 =  buffers->angular_flux_out[5];
    double * angular_flux_out_6 =  buffers->angular_flux_out[6];
    double * angular_flux_out_7 =  buffers->angular_flux_out[7];
    double * angular_flux_in_0  =  buffers->angular_flux_in[0];
    double * angular_flux_in_1  =  buffers->angular_flux_in[1];
    double * angular_flux_in_2  =  buffers->angular_flux_in[2];
    double * angular_flux_in_3  =  buffers->angular_flux_in[3];
    double * angular_flux_in_4  =  buffers->angular_flux_in[4];
    double * angular_flux_in_5  =  buffers->angular_flux_in[5];
    double * angular_flux_in_6  =  buffers->angular_flux_in[6];
    double * angular_flux_in_7  =  buffers->angular_flux_in[7];
    double * quad_weights =  buffers->quad_weights;
    double * scalar_flux =  buffers->scalar_flux;
    double * velocity_delta =  buffers->velocity_delta;


#pragma omp target data map(to: quad_weights[:nang], \
                            velocity_delta[:ng], \
                            nxyz,nx,ny,nz,ng,nang, \
                            angular_flux_in_0[:agxyz], \
                            angular_flux_in_1[:agxyz], \
                            angular_flux_in_2[:agxyz], \
                            angular_flux_in_3[:agxyz], \
                            angular_flux_in_4[:agxyz], \
                            angular_flux_in_5[:agxyz], \
                            angular_flux_in_6[:agxyz], \
                            angular_flux_in_7[:agxyz], \
                            angular_flux_out_0[:agxyz], \
                            angular_flux_out_1[:agxyz], \
                            angular_flux_out_2[:agxyz], \
                            angular_flux_out_3[:agxyz], \
                            angular_flux_out_4[:agxyz], \
                            angular_flux_out_5[:agxyz], \
                            angular_flux_out_6[:agxyz], \
                            angular_flux_out_7[:agxyz]) \
                        map(from: scalar_flux[:gxyz])


#pragma omp target teams distribute num_teams(nxyz) thread_limit(ng)
  for (int gid = 0; gid<nxyz ; gid++)
  {
       size_t i = gid % nx;
       size_t j = (gid / nx) % ny;
       size_t k = gid / (nx*ny);
#pragma omp parallel for
   for (unsigned int g = 0; g < ng; g++)
   {

    double local_scalar = 0.0;
#pragma omp parallel for reduction(+:local_scalar)
    for (int a = 0; a<nang; a++)
    {
// We want to perform a weighted sum of angles in each cell in each energy group
// One work-group per cell per energy group, and reduce within a work-group

        const double w = quad_weights[a];
        if (velocity_delta[g] != 0.0)
        {
            local_scalar +=
                w * (0.5 * (angular_flux_out_0(a,g,i,j,k) + angular_flux_in_0(a,g,i,j,k))) +
                w * (0.5 * (angular_flux_out_1(a,g,i,j,k) + angular_flux_in_1(a,g,i,j,k))) +
                w * (0.5 * (angular_flux_out_2(a,g,i,j,k) + angular_flux_in_2(a,g,i,j,k))) +
                w * (0.5 * (angular_flux_out_3(a,g,i,j,k) + angular_flux_in_3(a,g,i,j,k))) +
                w * (0.5 * (angular_flux_out_4(a,g,i,j,k) + angular_flux_in_4(a,g,i,j,k))) +
                w * (0.5 * (angular_flux_out_5(a,g,i,j,k) + angular_flux_in_5(a,g,i,j,k))) +
                w * (0.5 * (angular_flux_out_6(a,g,i,j,k) + angular_flux_in_6(a,g,i,j,k))) +
                w * (0.5 * (angular_flux_out_7(a,g,i,j,k) + angular_flux_in_7(a,g,i,j,k)));
        }
        else
        {
            local_scalar +=
                w * angular_flux_out_0(a,g,i,j,k) +
                w * angular_flux_out_1(a,g,i,j,k) +
                w * angular_flux_out_2(a,g,i,j,k) +
                w * angular_flux_out_3(a,g,i,j,k) +
                w * angular_flux_out_4(a,g,i,j,k) +
                w * angular_flux_out_5(a,g,i,j,k) +
                w * angular_flux_out_6(a,g,i,j,k) +
                w * angular_flux_out_7(a,g,i,j,k);
        }
    }
    scalar_flux(g,i,j,k) = local_scalar;
   }
  }
}

#define SCALAR_FLUX_MOMENTS_INDEX(m,g,i,j,k,cmom,ng,nx,ny) ((m)+((cmom-1)*(g))+((cmom-1)*(ng)*(i))+((cmom-1)*(ng)*(nx)*(j))+((cmom-1)*(ng)*(nx)*(ny)*(k)))
#define SCAT_COEFF_INDEX(a,l,o,nang,cmom) ((a)+((nang)*(l))+((nang)*(cmom)*o))

#define scalar_flux_moments(l,g,i,j,k) scalar_flux_moments[SCALAR_FLUX_MOMENTS_INDEX((l),(g),(i),(j),(k),cmom,ng,nx,ny)]
#define scat_coeff(a,l,o) scat_coeff[SCAT_COEFF_INDEX((a),(l),(o),nang,cmom)]

void compute_scalar_flux_moments(
    struct problem * problem,
    struct rankinfo * rankinfo,
    struct buffers * buffers
    )
{
    // get closest power of 2 to nang
    unsigned int nang = problem->nang;
    unsigned int nx = rankinfo->nx;
    unsigned int ny = rankinfo->ny;
    unsigned int nz = rankinfo->nz;
    unsigned int ng = problem->ng;
    unsigned int cmom = problem->cmom;

    size_t nxyz = nx*ny*nz;
    size_t ngxyz = ng*nx*ny*nz;
    size_t agxyz = ng*nang*nx*ny*nz;
    size_t cgxyz = (cmom-1)*ng*nx*ny*nz;

    double * angular_flux_out_0 =  buffers->angular_flux_out[0];
    double * angular_flux_out_1 =  buffers->angular_flux_out[1];
    double * angular_flux_out_2 =  buffers->angular_flux_out[2];
    double * angular_flux_out_3 =  buffers->angular_flux_out[3];
    double * angular_flux_out_4 =  buffers->angular_flux_out[4];
    double * angular_flux_out_5 =  buffers->angular_flux_out[5];
    double * angular_flux_out_6 =  buffers->angular_flux_out[6];
    double * angular_flux_out_7 =  buffers->angular_flux_out[7];
    double * angular_flux_in_0  =  buffers->angular_flux_in[0];
    double * angular_flux_in_1  =  buffers->angular_flux_in[1];
    double * angular_flux_in_2  =  buffers->angular_flux_in[2];
    double * angular_flux_in_3  =  buffers->angular_flux_in[3];
    double * angular_flux_in_4  =  buffers->angular_flux_in[4];
    double * angular_flux_in_5  =  buffers->angular_flux_in[5];
    double * angular_flux_in_6  =  buffers->angular_flux_in[6];
    double * angular_flux_in_7  =  buffers->angular_flux_in[7];
    double * velocity_delta =  buffers->velocity_delta;
    double * quad_weights =  buffers->quad_weights;
    double * scat_coeff =  buffers->scat_coeff;
    double * scalar_flux_moments =  buffers->scalar_flux_moments;

#pragma omp target data map(to: quad_weights[:nang], \
                            scat_coeff[:nang*cmom*8], \
                            velocity_delta[:ng], \
                            nxyz,nx,ny,nz,ng,nang,cmom, \
                            angular_flux_in_0[:agxyz], \
                            angular_flux_in_1[:agxyz], \
                            angular_flux_in_2[:agxyz], \
                            angular_flux_in_3[:agxyz], \
                            angular_flux_in_4[:agxyz], \
                            angular_flux_in_5[:agxyz], \
                            angular_flux_in_6[:agxyz], \
                            angular_flux_in_7[:agxyz], \
                            angular_flux_out_0[:agxyz], \
                            angular_flux_out_1[:agxyz], \
                            angular_flux_out_2[:agxyz], \
                            angular_flux_out_3[:agxyz], \
                            angular_flux_out_4[:agxyz], \
                            angular_flux_out_5[:agxyz], \
                            angular_flux_out_6[:agxyz], \
                            angular_flux_out_7[:agxyz]) \
                        map(from: scalar_flux_moments[:cgxyz])

#pragma omp target teams distribute num_teams(nxyz) thread_limit(ng*(cmom-1))
    for (int gid = 0; gid<nxyz ; gid++)
    {
       size_t i = gid % nx;
       size_t j = (gid / nx) % ny;
       size_t k = gid / (nx*ny);
#pragma omp parallel for collapse(2)
       for (unsigned int g = 0; g < ng; g++)
       {

          for (unsigned int l = 0; l < cmom-1; l++)
          {
             double local_scalar = 0.0;
#pragma omp parallel for reduction(+:local_scalar)
             for (int a = 0; a<nang; a++)
             {
// We want to perform a weighted sum of angles in each cell in each energy group
// One work-group per cell per energy group, and reduce within a work-group

        // Load into local memory
                const double w = quad_weights[a];
                if (velocity_delta[g] != 0.0)
                {
                   local_scalar +=
                    scat_coeff(a,l+1,0) * w * (0.5 * (angular_flux_out_0(a,g,i,j,k) + angular_flux_in_0(a,g,i,j,k))) +
                    scat_coeff(a,l+1,1) * w * (0.5 * (angular_flux_out_1(a,g,i,j,k) + angular_flux_in_1(a,g,i,j,k))) +
                    scat_coeff(a,l+1,2) * w * (0.5 * (angular_flux_out_2(a,g,i,j,k) + angular_flux_in_2(a,g,i,j,k))) +
                    scat_coeff(a,l+1,3) * w * (0.5 * (angular_flux_out_3(a,g,i,j,k) + angular_flux_in_3(a,g,i,j,k))) +
                    scat_coeff(a,l+1,4) * w * (0.5 * (angular_flux_out_4(a,g,i,j,k) + angular_flux_in_4(a,g,i,j,k))) +
                    scat_coeff(a,l+1,5) * w * (0.5 * (angular_flux_out_5(a,g,i,j,k) + angular_flux_in_5(a,g,i,j,k))) +
                    scat_coeff(a,l+1,6) * w * (0.5 * (angular_flux_out_6(a,g,i,j,k) + angular_flux_in_6(a,g,i,j,k))) +
                    scat_coeff(a,l+1,7) * w * (0.5 * (angular_flux_out_7(a,g,i,j,k) + angular_flux_in_7(a,g,i,j,k)));
               }
               else
               {
                local_scalar +=
                    scat_coeff(a,l+1,0) * w * angular_flux_out_0(a,g,i,j,k) +
                    scat_coeff(a,l+1,1) * w * angular_flux_out_1(a,g,i,j,k) +
                    scat_coeff(a,l+1,2) * w * angular_flux_out_2(a,g,i,j,k) +
                    scat_coeff(a,l+1,3) * w * angular_flux_out_3(a,g,i,j,k) +
                    scat_coeff(a,l+1,4) * w * angular_flux_out_4(a,g,i,j,k) +
                    scat_coeff(a,l+1,5) * w * angular_flux_out_5(a,g,i,j,k) +
                    scat_coeff(a,l+1,6) * w * angular_flux_out_6(a,g,i,j,k) +
                    scat_coeff(a,l+1,7) * w * angular_flux_out_7(a,g,i,j,k);
               }
            }
            scalar_flux_moments(l,g,i,j,k) = local_scalar;
         }
      }
    }
}

void copy_back_scalar_flux(
    struct problem *problem,
    struct rankinfo * rankinfo,
    struct buffers * buffers,
    double * oscalar_flux,
    bool blocking
    )
{
   size_t nsize = rankinfo->nx*rankinfo->ny*rankinfo->nz*problem->ng;
   double * scalar_flux = buffers->scalar_flux;
#pragma omp target update from(scalar_flux[:nsize])
   memcpy( oscalar_flux,scalar_flux,nsize*sizeof(double));
}
