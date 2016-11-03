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

//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
//------------------------------------------------------------------------------------------------------------------------------
void smooth(level_type * level, int x_id, int rhs_id, double a, double b){
  if(NUM_SMOOTHS&1){
    fprintf(stderr,"error - NUM_SMOOTHS must be even...\n");
    exit(0);
  }

  #ifdef USE_L1JACOBI
  double weight = 1.0;
  #else
  double weight = 2.0/3.0;
  #endif
 
  int block,s;
  for(s=0;s<NUM_SMOOTHS;s++){
    // exchange ghost zone data... Jacobi ping pongs between x_id and VECTOR_TEMP
    if((s&1)==0){exchange_boundary(level,x_id,stencil_get_shape());apply_BCs(level,x_id,stencil_get_shape());}
    else{exchange_boundary(level,VECTOR_TEMP,stencil_get_shape());apply_BCs(level,VECTOR_TEMP,stencil_get_shape());}

    // apply the smoother... Jacobi ping pongs between x_id and VECTOR_TEMP
    double _timeStart = getTime();

    blockCopy_type * my_blocks = level->my_blocks;
    box_type * my_boxes = level->my_boxes;
    double * vector_base = level->vectors[0];
    double * RedBlack_FP = level->RedBlack_FP;
    int num_my_blocks = level->num_my_blocks;
    int num_my_boxes = level->num_my_boxes;
    int num_vectors = level->numVectors;
    int box_volume = level->box_volume;
    double h = level->h;
    int box_ghosts = level->box_ghosts;

#pragma omp target							\
  map(to:my_blocks[0:num_my_blocks], my_boxes[0:num_my_boxes])		\
  map(tofrom:vector_base[0:(num_vectors*num_my_boxes*box_volume)])	\
  if(num_my_blocks >= GPU_THRESHOLD && GPU_OFFLOAD_ENABLE && GPU_ENABLE_SMOOTHER)
#pragma omp teams distribute						\
  firstprivate(num_my_blocks, num_my_boxes, num_vectors, box_volume, h,	\
	       box_ghosts, rhs_id, x_id, s)
    for(block=0;block<num_my_blocks;block++){
      const int box = my_blocks[block].read.box;
      const int ilo = my_blocks[block].read.i;
      const int jlo = my_blocks[block].read.j;
      const int klo = my_blocks[block].read.k;
      const int ihi = my_blocks[block].dim.i + ilo;
      const int jhi = my_blocks[block].dim.j + jlo;
      const int khi = my_blocks[block].dim.k + klo;

      int i,j,k;
      const double h2inv = 1.0/(h*h);
      const int ghosts =  box_ghosts;
      const int jStride = my_boxes[box].jStride;
      const int kStride = my_boxes[box].kStride;
      const int color000 = (my_boxes[box].low.i^my_boxes[box].low.j^my_boxes[box].low.k^s)&1;
      // is element 000 red or black on *THIS* sweep

      const double * __restrict__ rhs = &vector_base[rhs_id*num_my_boxes*box_volume +
						     box*box_volume +
						     ghosts*(1+jStride+kStride)];
      const double * __restrict__ alpha = &vector_base[VECTOR_ALPHA*num_my_boxes*box_volume +
						       box*box_volume +
						       ghosts*(1+jStride+kStride)];
      const double * __restrict__ beta_i = &vector_base[VECTOR_BETA_I*num_my_boxes*box_volume +
							box*box_volume +
							ghosts*(1+jStride+kStride)];
      const double * __restrict__ beta_j = &vector_base[VECTOR_BETA_J*num_my_boxes*box_volume +
							box*box_volume +
							ghosts*(1+jStride+kStride)];
      const double * __restrict__ beta_k = &vector_base[VECTOR_BETA_K*num_my_boxes*box_volume +
							box*box_volume +
							ghosts*(1+jStride+kStride)];
      #ifdef USE_L1JACOBI
      const double * __restrict__ lambda = &vector_base[VECTOR_L1INV*num_my_boxes*box_volume +
							box*box_volume +
							ghosts*(1+jStride+kStride)];
      #else
      const double * __restrict__ lambda = &vector_base[VECTOR_DINV*num_my_boxes*box_volume +
							box*box_volume +
							ghosts*(1+jStride+kStride)];
      #endif
      const double * __restrict__ x_n;
            double * __restrict__ x_np1;
      if((s&1)==0){
	x_n = &vector_base[x_id*num_my_boxes*box_volume +
			   box*box_volume +
			   ghosts*(1+jStride+kStride)];
	x_np1 = &vector_base[VECTOR_TEMP*num_my_boxes*box_volume +
			     box*box_volume +
			     ghosts*(1+jStride+kStride)];
      }
      else{
	x_n = &vector_base[VECTOR_TEMP*num_my_boxes*box_volume +
			   box*box_volume +
			   ghosts*(1+jStride+kStride)];
	x_np1 = &vector_base[x_id*num_my_boxes*box_volume +
			     box*box_volume +
			     ghosts*(1+jStride+kStride)];
      }

#pragma omp parallel for collapse(3)					\
      if(num_my_blocks >= GPU_THRESHOLD && GPU_OFFLOAD_ENABLE && GPU_ENABLE_SMOOTHER)      
      for(k=klo;k<khi;k++){
	for(j=jlo;j<jhi;j++){
	  for(i=ilo;i<ihi;i++){
	    int ijk = i + j*jStride + k*kStride;
	    double Ax_n = apply_op_ijk(x_n);
	    x_np1[ijk] = x_n[ijk] + weight*lambda[ijk]*(rhs[ijk]-Ax_n);
	  }
	}
      }

    } // box-loop
    level->timers.smooth += (double)(getTime()-_timeStart);
  } // s-loop
}

//------------------------------------------------------------------------------------------------------------------------------
