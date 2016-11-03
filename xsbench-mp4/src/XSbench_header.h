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

#ifndef __XSBENCH_HEADER_H__
#define __XSBENCH_HEADER_H__

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<strings.h>
#include<math.h>
#include<omp.h>
#include<unistd.h>
#include<sys/time.h>

// Papi Header
#ifdef PAPI
#include "papi.h"
#endif

// I/O Specifiers
#define INFO 1
#define DEBUG 1
#define SAVE 1

// Structures
typedef struct{
	double energy;
	double total_xs;
	double elastic_xs;
	double absorbtion_xs;
	double fission_xs;
	double nu_fission_xs;
} NuclideGridPoint;

typedef struct{
	double energy;
	int xs_ptrs;
} GridPoint;

typedef struct{
	int nthreads;
	long n_isotopes;
	long n_gridpoints;
	int lookups;
	char * HM;
} Inputs;

// XSutils.c function prototypes
NuclideGridPoint ** gpmatrix(size_t m, size_t n);
double *** d3tensor(size_t p, size_t q, size_t r);
void gpmatrix_free( NuclideGridPoint ** M );
int NGP_compare( const void * a, const void * b );
int binary_search( NuclideGridPoint * A, double quarry, uint n );
double rn(unsigned long * seed);
double rn_v(void);
unsigned int hash(unsigned char *str, int nbins);
size_t estimate_mem_usage(Inputs in);
void binary_dump(long n_isotopes, long n_gridpoints, NuclideGridPoint ** nuclide_grids,
     GridPoint * energy_grid, int * grid_ptrs);
void binary_read(long n_isotopes, long n_gridpoints, NuclideGridPoint ** nuclide_grids,
     GridPoint * energy_grid, int * grid_ptrs);
double timer();

// io.c funtion prototypes
void logo(int version);
void center_print(const char *s, int width);
void print_results(Inputs in, int mype, double runtime, int nprocs,
     double vval, unsigned long long vhash);
void print_inputs(Inputs in, int nprocs, int version);
void border_print(void);
void fancy_int(long a);
void print_CLI_error(void);
Inputs read_CLI( int argc, char * argv[] );

// Materials.c funtion prototypes
int * load_num_nucs(long n_isotopes);
int * load_mats_idx(int * num_nucs);
int * load_mats( int * num_nucs, int * mats_idx, int size_mats, long n_isotopes );
double * load_concs( int size_mats );
double * load_concs_v( int size_mats );
int pick_mat(double roll);

// GridInit.c funtion prototypes
void generate_grids( long n_isotopes, long n_gridpoints, 
		     NuclideGridPoint **nuclide_grids);
void generate_grids_v( long n_isotopes, long n_gridpoints, 
		       NuclideGridPoint **nuclide_grids);
void sort_nuclide_grids( long n_isotopes, long n_gridpoints, 
			 NuclideGridPoint **nuclide_grids);
int * generate_ptr_grid(int n_isotopes, int n_gridpoints);
GridPoint * generate_energy_grid( long n_isotopes, long n_gridpoints, 
     NuclideGridPoint **nuclide_grids, int * grid_ptrs);
void set_grid_ptrs( GridPoint * energy_grid, int * grid_ptrs, long n_isotopes,
		    long n_gridpoints, NuclideGridPoint **nuclide_grids);

// CalculateXS.c funtion prototypes
void calculate_micro_xs(double p_energy, int nuc, long n_isotopes, long n_gridpoints,
     GridPoint * restrict energy_grid, int * restrict grid_ptrs,
     NuclideGridPoint ** restrict nuclide_grids, int idx,
     double * restrict xs_vector );
void calculate_macro_xs( double p_energy, int mat, long n_isotopes,
    long n_gridpoints, int * restrict num_nucs, double * restrict concs,
    GridPoint * restrict energy_grid, int * restrict grid_ptrs,
    NuclideGridPoint ** restrict nuclide_grids,
    int * restrict mats, int * restrict mats_idx,
    double * restrict macro_xs_vector );
long grid_search( long n, double quarry, GridPoint * A);

// papi.c funtion prototypes
void counter_stop( int * eventset, int num_papi_events );
void counter_init( int * eventset, int * num_papi_events );

#endif
