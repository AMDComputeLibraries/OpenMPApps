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
// Lu = a*alpha[]*u[] - b*divergence( beta[]*gradient(u[]) )
//------------------------------------------------------------------------------------------------------------------------------
#ifndef DEFINES_H
#define DEFINES_H
//------------------------------------------------------------------------------------------------------------------------------
#define  VECTOR_TEMP         0 // 
#define  VECTOR_E            1 // error used in residual correction FMG
#define  VECTOR_F_MINUS_AV   2 // cell centered residual (f-Av)
//------------------------------------------------------------------------------------------------------------------------------
#define  VECTOR_F            3 // original right-hand side (Au=f), cell centered
#define  VECTOR_U            4 // numerical solution
#define  VECTOR_ALPHA        5 // cell centered coefficient
#define  VECTOR_BETA_I       6 // face centered coefficient (n.b. element 0 is the left face of the ghost zone element)
#define  VECTOR_BETA_J       7 // face centered coefficient (n.b. element 0 is the back face of the ghost zone element)
#define  VECTOR_BETA_K       8 // face centered coefficient (n.b. element 0 is the bottom face of the ghost zone element)
//------------------------------------------------------------------------------------------------------------------
#define  VECTOR_DINV         9 // cell centered relaxation parameter (e.g. inverse of the diagonal)
#define  VECTOR_L1INV       10 // cell centered relaxation parameter (e.g. inverse of the L1 norm of each row)
//------------------------------------------------------------------------------------------------------------------
#define VECTORS_RESERVED    11 // total number of vectors and the starting location for any auxillary bottom solver vectors
//------------------------------------------------------------------------------------------------------------------------------
#define GPU_OFFLOAD_ENABLE 1 // General GPU offloading enabled. If disabled, it will disable smoother GPU Offloading regardless of GPU_ENABLE_SMOOTHER. 
#define GPU_ENABLE_SMOOTHER 1 // Enable GPU offloading specifically for smoother
#define GPU_THRESHOLD 16 // Number of blocks needed before using GPU offloading

#define GSRB_BRANCH // If using GSRB, the GSRB_BRANCH is the only version that gives the exact same error as CPU only version. GSRB_STRIDE2 compiles, but there error is off by about 10^-12.

#endif
