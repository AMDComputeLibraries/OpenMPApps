
#include "convergence.h"

static double tolr=1.0E-12;

int inner_convergence(
    const struct problem * problem,
    const struct rankinfo * rankinfo,
    const struct buffers * buffers
    )
{
    double diffs[problem->ng];
    for (unsigned int g = 0; g < problem->ng; g++)
        diffs[g] = -DBL_MAX;

    // Calculate the maximum difference across each sub-domain for each group
    for (unsigned int k = 0; k < rankinfo->nz; k++)
        for (unsigned int j = 0; j < rankinfo->ny; j++)
            for (unsigned int i = 0; i < rankinfo->nx; i++)
                for (unsigned int g = 0; g < problem->ng; g++)
                {
                    double newsf = buffers->scalar_flux[SCALAR_FLUX_INDEX(g,i,j,k,problem->ng,rankinfo->nx,rankinfo->ny)];
                    double old = buffers->old_inner_scalar_flux[SCALAR_FLUX_INDEX(g,i,j,k,problem->ng,rankinfo->nx,rankinfo->ny)];
//printf("%d %d %d %d newsf %g old %g\n",g,i,j,k,newsf,old);
                    if (fabs(old) > tolr)
                        diffs[g] = fmax(fabs(newsf / old - 1.0), diffs[g]);
                    else
                        diffs[g] = fmax(fabs(newsf - old), diffs[g]);
                }

    // Do an AllReduce for each group to work out global maximum difference
    double recv[problem->ng];
    int result = 0;
    for (unsigned int g = 0; g < problem->ng; g++)
    {
//printf("%d diffs %g recv %g\n",g,diffs[g] , recv[g]);
        diffs[g] = 0.0;//recv[g];
        result += (diffs[g] >= problem->epsi)?1:0;
    }
    return result;
}


bool outer_convergence(
    const struct problem * problem,
    const struct rankinfo * rankinfo,
    const struct buffers * buffers,
    double * max_diff
    )
{

    *max_diff = -DBL_MAX;

    // Calculate the maximum difference across each sub-domain for each group
    for (unsigned int k = 0; k < rankinfo->nz; k++)
        for (unsigned int j = 0; j < rankinfo->ny; j++)
            for (unsigned int i = 0; i < rankinfo->nx; i++)
                for (unsigned int g = 0; g < problem->ng; g++)
                {
                    double newsf = buffers->scalar_flux[SCALAR_FLUX_INDEX(g,i,j,k,problem->ng,rankinfo->nx,rankinfo->ny)];
                    double old = buffers->old_outer_scalar_flux[SCALAR_FLUX_INDEX(g,i,j,k,problem->ng,rankinfo->nx,rankinfo->ny)];
                    if (fabs(old) > tolr)
                        *max_diff = fmax(fabs(newsf / old - 1.0), *max_diff);
                    else
                        *max_diff = fmax(fabs(newsf - old), *max_diff);
                }


    // Do an AllReduce to work out global maximum difference
    double recv =0.0;
    *max_diff = recv;
    return *max_diff <= 100.0 * problem->epsi;
}
