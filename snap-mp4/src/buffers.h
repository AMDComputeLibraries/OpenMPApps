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

/** \file
* \brief Manage the allocation of GPU buffers
*/

#pragma once

#include "global.h"

/** \brief Struct to contain all the GPU buffers */
struct buffers
{
    /** @{ \brief
    Angular flux - two copies for time dependence, each ocant in own buffer
    */
    double *angular_flux_in[8];
    double *angular_flux_out[8];
    /** @} */

    /** @{
    \brief Edge flux arrays */
    double *flux_i, *flux_j, *flux_k;
    /** @} */

    /** @{ \brief Scalar flux arrays */
    double *scalar_flux;
    double *scalar_flux_moments;
    double *old_inner_scalar_flux;
    double *old_outer_scalar_flux;
    /** @} */

    /** \brief Quadrature weights */
    double *quad_weights;

    /** @{ \brief Cosine coefficients */
    double * mu, * eta, * xi;
    /** @} */

    /** \brief Scattering coefficient */
    double *scat_coeff;

    /** \brief Material cross section */
    double *mat_cross_section;

    /** @{ \brief Source terms */
    double *fixed_source;
    double *outer_source;
    double *inner_source;
    /** @} */

    /** \brief Scattering terms */
    double *scattering_matrix;

    /** @{ \brief Diamond diference co-efficients */
    double *dd_i, *dd_j, *dd_k;
    /** @} */

    /** \brief Mock velocities array */
    double *velocities;

    /** \brief Time absorption coefficient */
    double *velocity_delta;

    /** \brief Transport denominator */
    double *denominator;

    /** \brief Lists of cell indicies in each plane
    Each buffer is an array of the i,j,k indicies for cells within that plane
    One buffer per plane */
    struct cell_id ** planes;
};

/** \brief Check the device has enough memory to allocate the buffers */
void check_device_memory_requirements(struct problem * problem, struct rankinfo * rankinfo);

/** \brief Launch a kernel to zero the (1D) buffer (non-blocking) */
void zero_buffer(double *buffer, size_t offset, size_t size);

/** \brief Swap the angular flux pointers around (in <-> out) */
void swap_angular_flux_buffers(struct buffers * buffers);

/** \brief Allocate the problem arrays */
void allocate_buffers(struct problem * problem, struct rankinfo * rankinfo, struct buffers * buffers);

/** \brief Free the arrays sroted in the \a mem struct */
void free_buffers(struct buffers * buffers);
