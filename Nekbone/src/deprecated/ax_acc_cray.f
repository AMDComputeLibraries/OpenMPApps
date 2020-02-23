c-----------------------------------------------------------------------
      subroutine ax_acc(w,u,gxyz,ur,us,ut,wk,n) ! Matrix-vector product: w=A*u

      include 'SIZE'
      include 'TOTAL'

c      real w(nx1*ny1*nz1,nelt),u(nx1*ny1*nz1,nelt)
c      real gxyz(2*ldim,nx1*ny1*nz1,nelt)
c      parameter (lt=lx1*ly1*lz1*lelt)
c      real ur(lt),us(lt),ut(lt),wk(lt)
      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)

      real w(nx1,ny1,nz1,nelt)
      real u(nx1,ny1,nz1,nelt)
      real gxyz(nx1,ny1,nz1,2*ldim,lelt)

      real ur(nx1,ny1,nz1,lelt)
      real us(nx1,ny1,nz1,lelt)
      real ut(nx1,ny1,nz1,lelt)
      real wk(nx1,ny1,nz1,lelt)

      real wr,ws,wt,tmp
      integer i,j,k,l,e,n

      integer lt

      lt = nx1*ny1*nz1*nelt

!$ACC DATA PRESENT(w,u,gxyz,ur,us,ut,wk,dxm1,dxtm1)
!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR
!$acc&  private(wr,ws,wt)
!DIR NOBLOCKING
      do e = 1,nelt
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            wr = 0      ! scalar ur gets promoted to a vector register over index ‘i’
            ws = 0
            wt = 0
!$ACC LOOP SEQ
            do l=1,nx1    ! serial loop, no reduction needed
               wr = wr + dxm1(i,l)*u(l,j,k,e)
               ws = ws + dxm1(j,l)*u(i,l,k,e)
               wt = wt + dxm1(k,l)*u(i,j,l,e)
            enddo
            ur(i,j,k,e) = gxyz(i,j,k,1,e)*wr
     $                  + gxyz(i,j,k,2,e)*ws
     $                  + gxyz(i,j,k,3,e)*wt
            us(i,j,k,e) = gxyz(i,j,k,2,e)*wr
     $                  + gxyz(i,j,k,4,e)*ws
     $                  + gxyz(i,j,k,5,e)*wt
            ut(i,j,k,e) = gxyz(i,j,k,3,e)*wr
     $                  + gxyz(i,j,k,5,e)*ws
     $                  + gxyz(i,j,k,6,e)*wt
         enddo
         enddo
         enddo
      enddo
!$ACC END PARALLEL LOOP

!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR 
      do e=1,nelt
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            w(i,j,k,e) = 0.0
!$ACC LOOP SEQ
            do l=1,nx1    ! serial loop, no reduction needed
               w(i,j,k,e) = w(i,j,k,e) + dxtm1(i,l)*ur(l,j,k,e)
     $                                 + dxtm1(j,l)*us(i,l,k,e)
     $                                 + dxtm1(k,l)*ut(i,j,l,e)
            enddo
         enddo
         enddo
         enddo
      enddo
!$ACC END PARALLEL LOOP

      call dssum2_acc(w)         ! Gather-scatter operation  ! w   = QQ  w
                                                            !            L
      call add2s2_acc(w,u,.1,n)   !2n
      call maskit_acc(w,cmask,nx1,ny1,nz1)  ! Zero out Dirichlet conditions

!$ACC END DATA

      nxyz=nx1*ny1*nz1
      flop_a = flop_a + (19*nxyz+12*nx1*nxyz)*nelt

      return
      end
c-----------------------------------------------------------------------
