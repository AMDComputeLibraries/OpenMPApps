
#pragma once

/** \file
* \brief Source update routines
*/

#include "global.h"
#include "problem.h"

#include "buffers.h"

#include "profiler.h"

/** \ingroup MEM
* @{
* \brief Index for source arrays */
/* define macros */
#define SOURCE_INDEX(m,g,i,j,k,cmom,ng,nx,ny) ((m)+((cmom)*(g))+((cmom)*(ng)*(i))+((cmom)*(ng)*(nx)*(j))+((cmom)*(ng)*(nx)*(ny)*(k)))


#define outer_source(m,g,i,j,k) outer_source[SOURCE_INDEX((m),(g),(i),(j),(k),cmom,ng,nx,ny)]
#define inner_source(m,g,i,j,k) inner_source[SOURCE_INDEX((m),(g),(i),(j),(k),cmom,ng,nx,ny)]
#define fixed_source(g,i,j,k) fixed_source[FIXED_SOURCE_INDEX((g),(i),(j),(k),ng,nx,ny)]
#define scattering_matrix(m,g1,g2) scattering_matrix[SCATTERING_MATRIX_INDEX((m),(g1),(g2),nmom,ng)]
#define scalar_flux(g,i,j,k) scalar_flux[SCALAR_FLUX_INDEX((g),(i),(j),(k),ng,nx,ny)]
#define scalar_flux_moments(m,g,i,j,k) scalar_flux_moments[SCALAR_FLUX_MOMENTS_INDEX((m),(g),(i),(j),(k),cmom,ng,nx,ny)]

/**@}*/

/** \brief Compute the outer source on the device (non-blocking)
*
* First moment is set to fixed source. Subsequent momemnts
* use group-to-group scattering.
*/
void compute_outer_source(const struct problem * problem, const struct rankinfo * rankinfo, struct buffers * buffers);

/** \brief Compute the inner source on the device (non-blocking)
*
* Set to the outer source plus within group scattering based on scalar flux and scalar flux moments.
*/
void compute_inner_source(const struct problem * problem, const struct rankinfo * rankinfo, struct buffers * buffers);
