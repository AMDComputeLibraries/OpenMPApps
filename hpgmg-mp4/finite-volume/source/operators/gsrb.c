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

#if   defined(GSRB_FP)
  #warning Overriding default GSRB implementation and using pre-computed 1.0/0.0 FP array for Red-Black to facilitate vectorization...
#elif defined(GSRB_STRIDE2)
  #if defined(GSRB_OOP)
  #warning Overriding default GSRB implementation and using out-of-place and stride-2 accesses to minimize the number of flops
  #else
  #warning Overriding default GSRB implementation and using stride-2 accesses to minimize the number of flops
  #endif
#elif defined(GSRB_BRANCH)
  #if defined(GSRB_OOP)
  #warning Overriding default GSRB implementation and using out-of-place implementation with an if-then-else on loop indices...
  #else
  #warning Overriding default GSRB implementation and using if-then-else on loop indices...
  #endif
#else
#define GSRB_STRIDE2 // default implementation
#endif
//------------------------------------------------------------------------------------------------------------------------------
void smooth(level_type * level, int x_id, int rhs_id, double a, double b){
  int block,s;
  int num_my_blocks = level->num_my_blocks;
  int num_my_boxes = level->num_my_boxes;
  blockCopy_type * my_blocks = level->my_blocks;
  box_type * my_boxes = level->my_boxes;
  double * vector_base = level->vectors[0];
  
#pragma omp target data map(to:my_blocks[0:num_my_blocks], my_boxes[0:num_my_boxes])\
  if(num_my_blocks >= GPU_THRESHOLD && GPU_OFFLOAD_ENABLE && GPU_ENABLE_SMOOTHER)
  {
  for(s=0;s<2*NUM_SMOOTHS;s++){ // there are two sweeps per GSRB smooth

    // exchange the ghost zone...
    #ifdef GSRB_OOP // out-of-place GSRB ping pongs between x and VECTOR_TEMP
    if((s&1)==0){
      exchange_boundary(level,x_id,stencil_get_shape());
      apply_BCs(level,x_id,stencil_get_shape());
    }
    else{
      exchange_boundary(level,VECTOR_TEMP,stencil_get_shape());
      apply_BCs(level,VECTOR_TEMP,stencil_get_shape());
    }
    #else // in-place GSRB only operates on x
    exchange_boundary(level,x_id,stencil_get_shape());
    apply_BCs(level,x_id,stencil_get_shape());
    #endif

    // apply the smoother...
    double _timeStart = getTime();

    double * RedBlack_FP = level->RedBlack_FP;
    int num_vectors = level->numVectors;
    int box_volume = level->box_volume;
    double h = level->h;
    int box_ghosts = level->box_ghosts;

    int    xn_base_for_this_s  = ((s&1) == 0) ? x_id : VECTOR_TEMP;
    int  xnp1_base_for_this_s  = ((s&1) == 0) ? VECTOR_TEMP : x_id;
    int size = level->num_my_boxes * level->box_volume;
    int index_xn   = xn_base_for_this_s   *num_my_boxes*box_volume;
    int index_xnp1 = xnp1_base_for_this_s *num_my_boxes*box_volume;
    double *vector_xn = vector_base;
    double *vector_xnp = vector_base;

#pragma omp target\
  map(to:vector_base[0:(num_vectors*num_my_boxes*box_volume)]) \
  map(to:vector_xn[index_xn:index_xn + size]) \
  map(from:vector_xnp[index_xnp1:index_xnp1 + size]) \
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
      const double * __restrict__ Dinv = &vector_base[VECTOR_DINV*num_my_boxes*box_volume +
						      box*box_volume +
						      ghosts*(1+jStride+kStride)];
#ifdef GSRB_OOP
      const double * __restrict__ x_n;
      double * __restrict__ x_np1;
      if((s&1)==0){
	x_n = &vector_xn[x_id*num_my_boxes*box_volume +
			   box*box_volume +
			   ghosts*(1+jStride+kStride)];
	x_np1 = &vector_xnp[VECTOR_TEMP*num_my_boxes*box_volume +
			     box*box_volume +
			     ghosts*(1+jStride+kStride)];
      }
      else{
	x_n = &vector_xn[VECTOR_TEMP*num_my_boxes*box_volume +
			   box*box_volume +
			   ghosts*(1+jStride+kStride)];
	x_np1 = &vector_xnp[x_id*num_my_boxes*box_volume +
			     box*box_volume +
			     ghosts*(1+jStride+kStride)];
      }
#else
      // i.e. [0] = first non ghost zone point      
      const double * __restrict__ x_n = &vector_base[VECTOR_TEMP*num_my_boxes*box_volume +
						     box*box_volume +
						     ghosts*(1+jStride+kStride)];
      // i.e. [0] = first non ghost zone point
      double * __restrict__ x_np1 = &vector_base[x_id*num_my_boxes*box_volume +
						 box*box_volume +
						 ghosts*(1+jStride+kStride)];
#endif
	    
#if defined(GSRB_FP)
#error GSRB_FP has not been implemented 
#pragma omp parallel for						\
 if(num_my_blocks >= GPU_THRESHOLD && GPU_OFFLOAD_ENABLE && GPU_ENABLE_SMOOTHER)      
      for(k=klo;k<khi;k++){
	const double * __restrict__ RedBlack = RedBlack_FP +
	  ghosts*(1+jStride) +
	  kStride*((k^color000)&0x1);
	for(j=jlo;j<jhi;j++){
	  for(i=ilo;i<ihi;i++){
	    int ij  = i + j*jStride;
	    int ijk = i + j*jStride + k*kStride;
	    double Ax     = apply_op_ijk(x_n);
	    double lambda =     Dinv_ijk();
	    x_np1[ijk] = x_n[ijk] + RedBlack[ij]*lambda*(rhs[ijk]-Ax);
	  }
	}
      }

#elif defined(GSRB_STRIDE2)
#warning GSRB_STRIDE2 running on GPU does not give same error as CPU version, please use GSRB_BRANCH
#pragma omp parallel for					\
  if(num_my_blocks >= GPU_THRESHOLD && GPU_OFFLOAD_ENABLE && GPU_ENABLE_SMOOTHER)
	for(k=klo;k<khi;k++){
	  for(j=jlo;j<jhi;j++){
#ifdef GSRB_OOP
	    // out-of-place must copy old value...
	    for(i=ilo;i<ihi;i++){
	      int ijk = i + j*jStride + k*kStride; 
	      x_np1[ijk] = x_n[ijk];
	    }
#endif
	    for(i=ilo+((ilo^j^k^color000)&1);i<ihi;i+=2){ // stride-2 GSRB
	      int ijk = i + j*jStride + k*kStride; 
	      double Ax     = apply_op_ijk(x_n);
	      double lambda =     Dinv_ijk();
	      x_np1[ijk] = x_n[ijk] + lambda*(rhs[ijk]-Ax);
	    }
	  }
	}

#elif defined(GSRB_BRANCH)
#pragma omp parallel for						\
  collapse(3) if(num_my_blocks >= GPU_THRESHOLD && GPU_OFFLOAD_ENABLE && GPU_ENABLE_SMOOTHER)
	for(k=klo;k<khi;k++){
	  for(j=jlo;j<jhi;j++){
	    for(i=ilo;i<ihi;i++){
	      int ijk = i + j*jStride + k*kStride;
	      if((i^j^k^color000^1)&1){ // looks very clean when [0] is i,j,k=0,0,0 
		double Ax     = apply_op_ijk(x_n);
		double lambda =     Dinv_ijk();
		x_np1[ijk] = x_n[ijk] + lambda*(rhs[ijk]-Ax);
#ifdef GSRB_OOP
	      }else{
		x_np1[ijk] = x_n[ijk]; // copy old value when sweep color != cell color
#endif
	      }
	    }
	  }
	}

#else
#error no GSRB implementation was specified
#endif

    } // boxes
    level->timers.smooth += (double)(getTime()-_timeStart);
  } // s-loop
  } //target data
}


//------------------------------------------------------------------------------------------------------------------------------
