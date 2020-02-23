#ifdef _CUDA

      attributes(global) subroutine ax_cuf2(w,u,ur,us,ut,
     &                gxyz,dxm1,dxtm1)

      include 'SIZE'

      real, intent(out) :: w(lx1,ly1,lz1,lelt)
      real u(lx1,ly1,lz1,lelt)
      real ur  (lx1,ly1,lz1,lelt)
      real us  (lx1,ly1,lz1,lelt)
      real ut  (lx1,ly1,lz1,lelt)

      real gxyz(lx1,ly1,lz1,2*ldim,lelt)

      real, intent(in) :: dxm1(lx1,lx1)
      real, intent(in) :: dxtm1(lx1,lx1)

      real rtmp,stmp,ttmp,wijke
      real, shared :: shdxm1(lx1,ly1)
      real, shared :: shdxtm1(lx1,ly1)
      integer l,e,i,j,k,kk,n,nstrides

      e = blockIdx%x
      k = threadIdx%z
      j = threadIdx%y
      i = threadIdx%x

      if (k.eq.1) then
         shdxm1(i,j) = dxm1(i,j)
         shdxtm1(i,j) = dxtm1(i,j)
      end if
      call syncthreads()

c Figure out how many strided accesses that this block needs to perform
      nstrides = lz1 / blockDim%z
      if (mod(lz1, blockDim%z) .gt. 0) then
        nstrides = nstrides + 1
      endif

c Perform the strided accesses.  Each thread in the block proceeds in
c lockstep.
      kk = k
      do n = 1, nstrides
        if (kk .le. lz1) then
          rtmp = 0.0
          stmp = 0.0
          ttmp = 0.0
          do l = 1, lx1
            rtmp = rtmp + shdxm1(i,l)  * u(l,j,kk,e)
            stmp = stmp + shdxm1(j,l)  * u(i,l,kk,e)
            ttmp = ttmp + shdxm1(kk,l) * u(i,j,l,e)
          enddo
          ur(i,j,kk,e) = gxyz(i,j,kk,1,e)*rtmp
     $                 + gxyz(i,j,kk,2,e)*stmp
     $                 + gxyz(i,j,kk,3,e)*ttmp
          us(i,j,kk,e) = gxyz(i,j,kk,2,e)*rtmp
     $                 + gxyz(i,j,kk,4,e)*stmp
     $                 + gxyz(i,j,kk,5,e)*ttmp
          ut(i,j,kk,e) = gxyz(i,j,kk,3,e)*rtmp
     $                 + gxyz(i,j,kk,5,e)*stmp
     $                 + gxyz(i,j,kk,6,e)*ttmp
        endif
        kk = kk + blockDim%z
c       rahaman 2017-03-31: The optimized kernels (e.g.nek_kernel16.cuf)
c       called synchthreads after each strided access.  I don't believe
c       this is necessary.  When I omit the thread sync, I see no
c       runtime erros and the solution matches CPU version
c       call syncthreads()
      enddo

      call syncthreads()

      kk = k
      do n = 1, nstrides
        if (kk .le. lz1) then
          wijke = 0.0
          do l = 1, lx1
            wijke = wijke + shdxtm1(i,l)  * ur(l,j,kk,e) 
     $                    + shdxtm1(j,l)  * us(i,l,kk,e)
     $                    + shdxtm1(kk,l) * ut(i,j,l,e)
          enddo
          w(i,j,kk,e) = wijke
        endif
        kk = kk + blockDim%z
c       rahaman 2017-03-31: The optimized kernels (e.g.nek_kernel16.cuf)
c       called synchthreads after each strided access.  I don't believe
c       this is necessary.  When I omit the thread sync, I see no
c       runtime erros and the solution matches CPU version
c       call syncthreads()
      enddo

      return
      end

#else

      subroutine ax_cuf2(w,u,ur,us,ut,gxyz,dxm1,dxtm1)
        call err_chk(
     $ 'ERROR: Called ax_cuf2 but did not compile with CUDA')
      return
      end

#endif
