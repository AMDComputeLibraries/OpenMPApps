
#include "source.h"


void compute_outer_source(
    const struct problem * problem,
    const struct rankinfo * rankinfo,
    struct buffers * buffers
    )
{
    unsigned int nang = problem->nang;
    unsigned int nx = rankinfo->nx;
    unsigned int ny = rankinfo->ny;
    unsigned int nz = rankinfo->nz;
    unsigned int ng = problem->ng;
    unsigned int cmom = problem->cmom;
    unsigned int nmom = problem->nmom;

    size_t cgxyz = (cmom-1)*ng*nx*ny*nz;
    size_t gxyz = ng*nx*ny*nz;
    size_t nxyz = nx*ny*nz;

    double * outer_source =  buffers->outer_source;
    double * fixed_source =  buffers->fixed_source;
    double * scalar_flux =  buffers->scalar_flux;
    double * scattering_matrix =  buffers->scattering_matrix;
    double * scalar_flux_moments =  buffers->scalar_flux_moments;

#pragma omp target data map(to: scattering_matrix[:nmom*ng*ng], \
                            scalar_flux[:gxyz], \
                            scalar_flux_moments[:cgxyz], \
                            fixed_source[:gxyz], \
                            nxyz,nx,ny,nz,ng,nmom,cmom) \
                        map(tofrom: outer_source[:cmom*gxyz])

#if 1
#pragma omp target 
#pragma omp parallel for collapse(3)
    for (int k = 0; k<nz; k++)
    {
       for (int j = 0; j<ny; j++)
       {
          for (int i = 0; i<nx; i++)
          {

    for (unsigned int g = 0; g < ng; g++)
    {
#else
#pragma omp target teams distribute num_teams(nxyz) thread_limit(ng)
    for (int gid = 0; gid<nxyz ; gid++)
    {
       size_t i = gid % nx;
       size_t j = (gid / nx) % ny;
       size_t k = gid / nx*ny;
#pragma omp parallel for 
    for (unsigned int g = 0; g < ng; g++)
    {
#endif
        // Set first moment to the fixed source
        outer_source(0,g,i,j,k) = fixed_source(g,i,j,k);

        // Loop over groups and moments to compute out-of-group scattering
        for (unsigned int g2 = 0; g2 < ng; g2++)
        {
            if (g == g2)
                continue;
            // Compute scattering source
            outer_source(0,g,i,j,k) += scattering_matrix(0,g2,g) * scalar_flux(g2,i,j,k);
            // Other moments
            unsigned int mom = 1;
            for (unsigned int l = 1; l < nmom; l++)
            {
                for (unsigned int m = 0; m < 2*l+1; m++)
                {
                    outer_source(mom,g,i,j,k) += scattering_matrix(l,g2,g) * scalar_flux_moments(mom-1,g2,i,j,k);
                    mom += 1;
                }
            }
        }
    }
    }
    }
    }

}


void compute_inner_source(
    const struct problem * problem,
    const struct rankinfo * rankinfo,
    struct buffers * buffers
    )
{
    unsigned int nang = problem->nang;
    unsigned int nx = rankinfo->nx;
    unsigned int ny = rankinfo->ny;
    unsigned int nz = rankinfo->nz;
    unsigned int ng = problem->ng;
    unsigned int cmom = problem->cmom;
    unsigned int nmom = problem->nmom;

    size_t nxyz = nx*ny*nz;
    size_t cgxyz = (cmom-1)*ng*nx*ny*nz;
    size_t gxyz = ng*nx*ny*nz;

    double * outer_source =  buffers->outer_source;
    double * inner_source =  buffers->inner_source;
    double * scalar_flux =  buffers->scalar_flux;
    double * scattering_matrix =  buffers->scattering_matrix;
    double * scalar_flux_moments =  buffers->scalar_flux_moments;

#pragma omp target data map(to: scattering_matrix[:nmom*ng*ng], \
                            scalar_flux[:gxyz], \
                            scalar_flux_moments[:cgxyz], \
                            nxyz,nx,ny,nz,ng,nmom,cmom, \
                            outer_source[:cmom*gxyz]) \
                        map(from: inner_source[:cmom*gxyz])
#if 1
#pragma omp target 
#pragma omp parallel for collapse(3)
    for (int k = 0; k<nz; k++)
    {
       for (int j = 0; j<ny; j++)
       {
          for (int i = 0; i<nx; i++)
          {

    for (unsigned int g = 0; g < ng; g++)
    {
#else
 
#pragma omp target teams distribute num_teams(nxyz) thread_limit(ng)
    for (int gid = 0; gid<nxyz ; gid++)
    {
       size_t i = gid % nx;
       size_t j = (gid / nx) % ny;
       size_t k = gid / nx*ny;
#pragma omp parallel for 
    for (unsigned int g = 0; g < ng; g++)
    {
#endif
        // Set first moment to outer source plus scattering contribution of scalar flux
        inner_source(0,g,i,j,k) = outer_source(0,g,i,j,k) + scattering_matrix(0,g,g) * scalar_flux(g,i,j,k);

        // Set other moments similarly based on scalar flux moments
        unsigned int mom = 1;
        for (unsigned int l = 1; l < nmom; l++)
        {
            for (unsigned int m = 0; m < 2*l+1; m++)
            {
                inner_source(mom,g,i,j,k) = outer_source(mom,g,i,j,k) + scattering_matrix(l,g,g) * scalar_flux_moments(mom-1,g,i,j,k);
                mom += 1;
            }
        }
    }

          }
       }
    }
}
