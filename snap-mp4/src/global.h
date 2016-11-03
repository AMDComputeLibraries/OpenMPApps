
#pragma once

/** \file
* \brief Structures containing problem sizes and MPI information
*/

#include <stdbool.h>
#include <stdlib.h>


/** \brief Problem dimensions
*
* Read in from input file or calculated from those inputs
*/
struct problem
{
    /**@{ \brief Global grid size */
    unsigned int nx, ny, nz;
    /**@}*/

    /**@{ \brief Physical grid size */
    double lx, ly, lz;
    /**@}*/

    /**@{ \brief Width of spatial cells */
    double dx, dy, dz;
    /**@}*/

    /** \brief Energy groups */
    unsigned int ng;

    /** \brief Angles per octant
        (3D assumed) */
    unsigned int nang;

    /** \brief Number of expansion moments */
    unsigned int nmom;

    /**  \brief Number of computational moments
    *
    * = nmom*nmom */
    unsigned int cmom;

    /**@{*/
    /** \brief Number of inner iterations */
    unsigned int iitm;
    /** \brief Number of outer iterations */
    unsigned int oitm;
    /**@}*/

    /**@{*/
    /** \brief Number of timesteps */
    unsigned int nsteps;
    /** \brief Total time to simulate */
    double tf;
    /** \brief Time domain stride */
    double dt;
    /**@}*/

    /** \brief Convergence criteria */
    double epsi;

    /**@{ \brief Number of MPI tasks in each direction */
    unsigned int npex, npey, npez;
    /**@}*/

    /** \brief KBA chunk size */
    unsigned int chunk;

    /** \brief Global variable to determine if there are multiple GPUs per node (if so we have to get GPUs VERY carefully) */
    bool multigpu;
};


/** \brief MPI Information
*
* Holds rankinfo information about tile size and MPI rank
*/
struct rankinfo
{
    /**  \brief My WORLD rank */
    int rank;

    /** \brief My MPI Cartesian co-ordinate ranks */
    int ranks[3];

    /**@{ \brief Local grid size */
    unsigned int nx, ny, nz;
    /**@}*/

    /**@{ \brief Global grid corners of MPI partition */
    unsigned int ilb, iub;
    unsigned int jlb, jub;
    unsigned int klb, kub;
    /**@}*/

    /**@{ \brief Global ranks of my neighbours */
    int xup, xdown;
    int yup, ydown;
    int zup, zdown;
    /**@}*/
};
