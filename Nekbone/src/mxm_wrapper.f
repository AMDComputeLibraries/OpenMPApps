      subroutine mxm(a,n1,b,n2,c,n3)

#if defined(XSMM_DISPATCH)
      USE :: LIBXSMM
#endif

#define LIBXSMM_DMM1(N, a, b, c)  LIBXSMM_DMM1_str(N, a, b, c) 
#define LIBXSMM_DMM1_str(N, a, b, c)  libxsmm_dmm_##N##x##N##_##N##_##N(a, b, c)
#define LIBXSMM_DMM2(N, a, b, c)   LIBXSMM_DMM2_str(N, a, b, c)
#define LIBXSMM_DMM2_str(N, a, b, c)  libxsmm_dmm_##N##_##N##x##N##_##N(a, b, c)
#define LIBXSMM_DMM3(N, a, b, c)  LIBXSMM_DMM3_str(N, a, b, c) 
#define LIBXSMM_DMM3_str(N, a, b, c)  libxsmm_dmm_##N##_##N##_##N(a, b, c)

c
c     Compute matrix-matrix product C = A*B
c     for contiguously packed matrices A,B, and C.
c
#if defined (MKL)
#     include  "mkl_direct_call.fi"
#endif
      real a(n1,n2),b(n2,n3),c(n1,n3)
      real alpha, beta
c
      include 'SIZE'
      include 'TOTAL'
c
      integer aligned
      integer K10_mxm
      integer init, prevn2

#if defined(XSMM_DISPATCH)
      TYPE(LIBXSMM_DMMFUNCTION) :: xmm1,xmm2,xmm3
#endif

      data init /0/, prevn2 /0/
      save init, prevn2
#if defined(XSMM_DISPATCH)
      save xsmm1, xsmm2, xsmm3
#endif

c     write(*,*) "in", init, prevn2, LOC(xsmm1), LOC(xsmm2), LOC(xsmm3)

#if defined (MKL)
      alpha = 1.0
      beta = 0.0
      call dgemm('N','N',n1,n3,n2,alpha,A,n1,B,n2,beta,C,n1)
#elif defined (BLAS_MXM)
      alpha = 1.0
      beta = 0.0
      call dgemm('N','N',n1,n3,n2,alpha,a,n1,b,n2,beta,c,n1)
#elif defined (BG)
      call bg_aligned3(a,b,c,aligned)
      if (n2.eq.2) then
         call mxm44_2(a,n1,b,n2,c,n3)
      else if ((aligned.eq.1) .and.
     $         (n1.ge.8) .and. (n2.ge.8) .and. (n3.ge.8) .and.
     $         (modulo(n1,2).eq.0) .and. (modulo(n2,2).eq.0) ) then
         if (modulo(n3,4).eq.0) then
            call bg_mxm44(a,n1,b,n2,c,n3)
         else
            call bg_mxm44_uneven(a,n1,b,n2,c,n3)
         endif
      else if((aligned.eq.1) .and.
     $        (modulo(n1,6).eq.0) .and. (modulo(n3,6).eq.0) .and.
     $        (n2.ge.4) .and. (modulo(n2,2).eq.0) ) then
         call bg_mxm3(a,n1,b,n2,c,n3)
      else
         call mxm44_0(a,n1,b,n2,c,n3)
      endif
#elif defined (K10_MXM)
      ! fow now only supported for lx1=8
      ! tuned for AMD K10
      ierr = K10_mxm(a,n1,b,n2,c,n3) 
      if (ierr.gt.0) call mxmf2(a,n1,b,n2,c,n3)
#elif defined (XSMM_DISPATCH)
      if (init == 0) then
        CALL libxsmm_init()
        init = 1
        write(*,*) "initializing libxsmm"
      end if 

      if (prevn2 /= n2) then
        prevn2 = n2

        CALL libxsmm_dispatch(xmm1, n2, n2, n2*n2, alpha=1D0, beta=0D0)
c       write(*,*) "initialized xmm1"
        IF (.NOT. libxsmm_available(xmm1)) THEN
          write(*,*) "  ** Error: unable to dispatch libxsmm call"
          STOP 
        END IF

        CALL libxsmm_dispatch(xmm2, n2, n2, n2, alpha=1D0, beta=0D0)
c       write(*,*) "initialized xmm2"
        IF (.NOT. libxsmm_available(xmm2)) THEN
          write(*,*) "  ** Error: unable to dispatch libxsmm call"
          STOP 
        END IF

        CALL libxsmm_dispatch(xmm3, n2*n2, n2, n2, alpha=1D0, beta=0D0)
c       write(*,*) "initialized xmm3"
        IF (.NOT. libxsmm_available(xmm3)) THEN
          write(*,*) "  ** Error: unable to dispatch libxsmm call"
          STOP 
        END IF
      end if

      if (n1 .eq. n2*n2) then
c       write(*,*) "call to xmm3", n1, n2, n3
        IF (.NOT. libxsmm_available(xmm3)) THEN
          write(*,*) "  ** Error: unable to dispatch libxsmm call"
          STOP 
        END IF
        call libxsmm_call(xmm3, C_LOC(a), C_LOC(b), C_LOC(c))
c       call libxsmm_dmm_256_16_16(a, b, c)
      else if (n3 .eq. n2*n2) then
c       write(*,*) "call to xmm1", n1, n2, n3
        IF (.NOT. libxsmm_available(xmm1)) THEN
          write(*,*) "  ** Error: unable to dispatch libxsmm call"
          STOP 
        END IF
        call libxsmm_call(xmm1, C_LOC(a), C_LOC(b), C_LOC(c))
c       call libxsmm_dmm_16_256_16(a, b, c)
      else
c       write(*,*) "call to xmm2", n1, n2, n3
        IF (.NOT. libxsmm_available(xmm2)) THEN
          write(*,*) "  ** Error: unable to dispatch libxsmm call"
          STOP 
        END IF
        call libxsmm_call(xmm2, C_LOC(a), C_LOC(b), C_LOC(c))
c       call libxsmm_dmm_16_16_16(a, b, c)
      end if
#elif defined (XSMM_FIXED)
      if (n2 .eq. NPOLY) then
        if (n1 .eq. n2*n2) then
          call LIBXSMM_DMM1(NPOLY, a, b, c)
        else if (n3 .eq. n2*n2) then
          call LIBXSMM_DMM2(NPOLY, a, b, c)
        else
          call LIBXSMM_DMM3(NPOLY, a, b, c)
        end if
      else
        write(*,*) "Invalid matrix size"
        stop
      end if
#elif defined (XSMM)
      alpha = 1.0
      beta = 0.0
      CALL libxsmm_dgemm('N','N',n1,n3,n2,alpha,A,n1,B,n2,beta,C,n1)
#elif defined (MXMBASIC)
      do j=1,n3
        do i=1,n1
          c(i,j) = 0.0
          do k=1,n2
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          enddo
        enddo
      enddo
#else
      call mxmf2(a,n1,b,n2,c,n3)
#endif

c     write(*,*) "out", init, prevn2, xsmm1, xsmm2, xsmm3

      return
      end
