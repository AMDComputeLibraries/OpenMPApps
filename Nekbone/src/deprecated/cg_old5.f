c-----------------------------------------------------------------------
      subroutine cg(x,f,g,c,r,w,p,z,n,niter,flop_cg)
      include 'SIZE'

c     Solve Ax=f where A is SPD and is invoked by ax()
c
c     Output:  x - vector of length n
c
c     Input:   f - vector of length n
c     Input:   g - geometric factors for SEM operator
c     Input:   c - inverse of the counting matrix
c
c     Work arrays:   r,w,p,z  - vectors of length n
c
c     User-provided ax(w,z,n) returns  w := Az,  
c
c     User-provided solveM(z,r,n) ) returns  z := M^-1 r,  
c

      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)
      parameter (lt=lx1*ly1*lz1*lelt)
      real ur(lt),us(lt),ut(lt),wk(lt)

      real x(n),f(n),r(n),w(n),p(n),z(n),g(1),c(n)

      character*1 ans

      pap = 0.0

c     set machine tolerances
      one = 1.
      eps = 1.e-20
      if (one+eps .eq. one) eps = 1.e-14
      if (one+eps .eq. one) eps = 1.e-7

      rtz1=1.0

      call rzero(x,n)
      call copy (r,f,n)
      call maskit (r,cmask,nx1,ny1,nz1) ! Zero out Dirichlet conditions

      rnorm = sqrt(glsc3(r,c,r,n))
      iter = 0
      if (nid.eq.0)  write(6,6) iter,rnorm

      miter = niter
c     call tester(z,r,n)  
      do iter=1,miter
         call solveM(z,r,n) ! preconditioner here

         rtz2=rtz1                                                       ! OPS
         rtz1=glsc3(r,c,z,n)   ! parallel weighted inner product r^T C z ! 3n

         beta = rtz1/rtz2
         if (iter.eq.1) beta=0.0
         call add2s1(p,z,beta,n)                                         ! 2n

         call ax(w,p,g,ur,us,ut,wk,n)                                    ! flopa
         pap=glsc3(w,c,p,n)                                              ! 3n

         alpha=rtz1/pap
         alphm=-alpha
         call add2s2(x,p,alpha,n)                                        ! 2n
         call add2s2(r,w,alphm,n)                                        ! 2n
         rtr = glsc3(r,c,r,n)                                            ! 3n
         if (iter.eq.1) rlim2 = rtr*eps**2
         if (iter.eq.1) rtr0  = rtr
         rnorm = sqrt(rtr)
c        if (nid.eq.0.and.mod(iter,100).eq.0) 
c     $        write(6,6) iter,rnorm,alpha,beta,pap
        if (nid.eq.0.) 
     $        write(6,6) iter,rnorm,alpha,beta,pap

    6    format('cg:',i4,1p4e12.4)
c        if (rtr.le.rlim2) goto 1001

      enddo

 1001 continue

      if (nid.eq.0) write(6,6) iter,rnorm,alpha,beta,pap

      flop_cg = flop_cg + iter*15.*n

      return
      end
c-----------------------------------------------------------------------
      subroutine solveM(z,r,n)
      include 'INPUT'
      real z(n),r(n)

      nn = n
      call h1mg_solve(z,r,nn)

      return
      end
c-----------------------------------------------------------------------
      subroutine ax(w,u,gxyz,ur,us,ut,wk,n) ! Matrix-vector product: w=A*u

      include 'SIZE'
      include 'TOTAL'

      real w(nx1*ny1*nz1,nelt),u(nx1*ny1*nz1,nelt)
      real gxyz(nx1*ny1*nz1,2*ldim,nelt)

      parameter (lt=lx1*ly1*lz1*lelt)
      real ur(lt),us(lt),ut(lt),wk(lt)
      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)

      integer e

      do e=1,nelt                                ! ~
         call ax_e( w(1,e),u(1,e),gxyz(1,1,e)    ! w   = A  u
     $                             ,ur,us,ut,wk) !  L     L  L
      enddo                                      ! 

      call dssum(w)         ! Gather-scatter operation  ! w   = QQ  w
                                                        !            L
      call add2s2(w,u,.1,n)   !2n
      call maskit(w,cmask,nx1,ny1,nz1)  ! Zero out Dirichlet conditions

      nxyz=nx1*ny1*nz1
      flop_a = flop_a + (19*nxyz+12*nx1*nxyz)*nelt

      return
      end
c-------------------------------------------------------------------------
      subroutine ax1(w,u,n)
      include 'SIZE'
      real w(n),u(n)
      real h2i
  
      h2i = (n+1)*(n+1)  
      do i = 2,n-1
         w(i)=h2i*(2*u(i)-u(i-1)-u(i+1))
      enddo
      w(1)  = h2i*(2*u(1)-u(2  ))
      w(n)  = h2i*(2*u(n)-u(n-1))

      return
      end
c-------------------------------------------------------------------------
      subroutine ax_e(w,u,g,ur,us,ut,wk) ! Local matrix-vector product
      include 'SIZE'
      include 'TOTAL'

      parameter (lxyz=lx1*ly1*lz1)
      real ur(lxyz),us(lxyz),ut(lxyz),wk(lxyz)
      real w(nx1*ny1*nz1),u(nx1*ny1*nz1),g(nx1*ny1*nz1,2*ldim)


      nxyz = nx1*ny1*nz1
      n    = nx1-1

      call local_grad3(ur,us,ut,u,n,dxm1,dxtm1)

      do i=1,nxyz
         wr = g(i,1)*ur(i) + g(i,2)*us(i) + g(i,3)*ut(i)
         ws = g(i,2)*ur(i) + g(i,4)*us(i) + g(i,5)*ut(i)
         wt = g(i,3)*ur(i) + g(i,5)*us(i) + g(i,6)*ut(i)
         ur(i) = wr
         us(i) = ws
         ut(i) = wt
      enddo

      call local_grad3_t(w,ur,us,ut,n,dxm1,dxtm1,wk)

      return
      end
c-------------------------------------------------------------------------
      subroutine local_grad3(ur,us,ut,u,n,D,Dt)
c     Output: ur,us,ut         Input:u,n,D,Dt
      real ur(0:n,0:n,0:n),us(0:n,0:n,0:n),ut(0:n,0:n,0:n)
      real u (0:n,0:n,0:n)
      real D (0:n,0:n),Dt(0:n,0:n)
      integer e

      m1 = n+1
      m2 = m1*m1

      call mxm(D ,m1,u,m1,ur,m2)
      do k=0,n
         call mxm(u(0,0,k),m1,Dt,m1,us(0,0,k),m1)
      enddo
      call mxm(u,m2,Dt,m1,ut,m1)

      return
      end
c-----------------------------------------------------------------------
      subroutine local_grad3_t(u,ur,us,ut,N,D,Dt,w)
c     Output: ur,us,ut         Input:u,N,D,Dt
      real u (0:N,0:N,0:N)
      real ur(0:N,0:N,0:N),us(0:N,0:N,0:N),ut(0:N,0:N,0:N)
      real D (0:N,0:N),Dt(0:N,0:N)
      real w (0:N,0:N,0:N)
      integer e

      m1 = N+1
      m2 = m1*m1
      m3 = m1*m1*m1

      call mxm(Dt,m1,ur,m1,u,m2)

      do k=0,N
         call mxm(us(0,0,k),m1,D ,m1,w(0,0,k),m1)
      enddo
      call add2(u,w,m3)

      call mxm(ut,m2,D ,m1,w,m1)
      call add2(u,w,m3)

      return
      end
c-----------------------------------------------------------------------
      subroutine maskit(w,pmask,nx,ny,nz)   ! Zero out Dirichlet conditions
      include 'SIZE'
      include 'PARALLEL'

      real pmask(-1:lx1*ly1*lz1*lelt)
      real w(1)
      integer e

      nxyz = nx*ny*nz
      nxy  = nx*ny
      if(pmask(-1).lt.0) then
        j=pmask(0)
        do i = 1,j
           k = pmask(i)
           w(k)=0.0
        enddo
      else
c         Zero out Dirichlet boundaries.
c
c                      +------+     ^ Y
c                     /   3  /|     |
c               4--> /      / |     |
c                   +------+ 2 +    +----> X
c                   |   5  |  /    /
c                   |      | /    /
c                   +------+     Z   
c

        nn = 0
        do e  = 1,nelt
          call get_face(w,nx,e)
          do i = 1,nxyz
             if(w(i).eq.0) then
               nn=nn+1
               pmask(nn)=i
             endif
          enddo
        enddo     
        pmask(-1) = -1.
        pmask(0) = nn
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine masko(w)   ! Old 'mask'
      include 'SIZE'
      real w(1)

      if (nid.eq.0) w(1) = 0.  ! suitable for solvability

      return
      end
c-----------------------------------------------------------------------
      subroutine masking(w,nx,e,x0,x1,y0,y1,z0,z1)
c     Zeros out boundary
      include 'SIZE'
      integer e,x0,x1,y0,y1,z0,z1
      real w(nx,nx,nx,nelt)
      
c       write(6,*) x0,x1,y0,y1,z0,z1
      do k=z0,z1
      do j=y0,y1
      do i=x0,x1
          w(i,j,k,e)=0.0
      enddo
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine tester(z,r,n)
c     Used to test if solution to precond. is SPD
      real r(n),z(n)

      do j=1,n
         call rzero(r,n)
         r(j) = 1.0
         call solveM(z,r,n)
         do i=1,n
            write(79,*) z(i)
         enddo
      enddo
      call exitt0
      return
      end
c-----------------------------------------------------------------------
      subroutine get_face(w,nx,ie)
c     zero out all boundaries as Dirichlet
c     to change, change this routine to only zero out 
c     the nodes desired to be Dirichlet, and leave others alone.
      include 'SIZE'
      include 'PARALLEL'
      real w(1)
      integer nx,ie,nelx,nely,nelz
      integer x0,x1,y0,y1,z0,z1
      
      x0=1
      y0=1
      z0=1
      x1=nx
      y1=nx
      z1=nx
      
      nelxy=nelx*nely
      ngl = lglel(ie)        !global element number

      ir = 1+(ngl-1)/nelxy   !global z-count
      iq = mod1(ngl,nelxy)   !global y-count
      iq = 1+(iq-1)/nelx     
      ip = mod1(ngl,nelx)    !global x-count

c     write(6,1) ip,iq,ir,nelx,nely,nelz, nelt,' test it'
c  1  format(7i7,a8)

      if(mod(ip,nelx).eq.1.or.nelx.eq.1)   then  ! Face4
         x0=1
         x1=1
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(ip,nelx).eq.0)               then   ! Face2
         x0=nx
         x1=nx
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      x0=1
      x1=nx
      if(mod(iq,nely).eq.1.or.nely.eq.1) then    ! Face1
         y0=1
         y1=1
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(iq,nely).eq.0)              then    ! Face3
         y0=nx
         y1=nx
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      y0=1
      y1=nx
      if(mod(ir,nelz).eq.1.or.nelz.eq.1) then    ! Face5
         z0=1
         z1=1
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(ir,nelz).eq.0)              then    ! Face6
         z1=nx
         z0=nx
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      return
      end
c-----------------------------------------------------------------------


#ifdef _OPENACC

c-----------------------------------------------------------------------
      subroutine cg_acc(x,f,g,c,r,w,p,z,n,niter,flop_cg)
      include 'SIZE'

c     Solve Ax=f where A is SPD and is invoked by ax()
c
c     Output:  x - vector of length n
c
c     Input:   f - vector of length n
c     Input:   g - geometric factors for SEM operator
c     Input:   c - inverse of the counting matrix
c
c     Work arrays:   r,w,p,z  - vectors of length n
c
c     User-provided ax(w,z,n) returns  w := Az,  
c
c     User-provided solveM(z,r,n) ) returns  z := M^-1 r,  
c

      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)
      parameter (lt=lx1*ly1*lz1*lelt)

c      real ur(lt),us(lt),ut(lt),wk(lt)      
      common /TEMP0_ACC/ ur(lx1,lx1,lx1,lelt)
     $     ,             us(lx1,lx1,lx1,lelt)
     $     ,             ut(lx1,lx1,lx1,lelt)
     $     ,             wk(lx1,lx1,lx1,lelt)
      real ur,us,ut,wk

      real x(n),f(n),r(n),w(n),p(n),z(n),c(n)

      real g(2*ldim,lt)

      character*1 ans

      pap = 0.0

c     set machine tolerances
      one = 1.
      eps = 1.e-20
      if (one+eps .eq. one) eps = 1.e-14
      if (one+eps .eq. one) eps = 1.e-7

      rtz1=1.0

!$ACC DATA PRESENT(x,g,c,r,w,p,z,ur,us,ut,wk)
      call copy(r,f,n)
      call maskit (r,cmask,nx1,ny1,nz1) ! Zero out Dirichlet conditions

!$ACC UPDATE DEVICE(r,cmask,c,p)
      call rzero_acc(x,n)
      rnorm = sqrt(glsc3_acc(r,c,r,n))
      iter = 0
      if (nid.eq.0)  write(6,6) iter,rnorm

      miter = niter
c     call tester(z,r,n)  
      do iter=1,miter
         call solveM_acc(z,r,n)    ! preconditioner here

         rtz2=rtz1                                                       ! OPS
         rtz1=glsc3_acc(r,c,z,n)   ! parallel weighted inner product r^T C z ! 3n

         beta = rtz1/rtz2
         if (iter.eq.1) beta=0.0
         call add2s1_acc(p,z,beta,n)                                     ! 2n

         call ax_acc(w,p,g,ur,us,ut,wk,n)                                ! flopa

         pap=glsc3_acc(w,c,p,n)                                          ! 3n

         alpha=rtz1/pap
         alphm=-alpha
         call add2s2_acc(x,p,alpha,n)                                    ! 2n
         call add2s2_acc(r,w,alphm,n)                                    ! 2n
         rtr = glsc3_acc(r,c,r,n)                                        ! 3n

         if (iter.eq.1) rlim2 = rtr*eps**2
         if (iter.eq.1) rtr0  = rtr
         rnorm = sqrt(rtr)
C        if (nid.eq.0.and.mod(iter,100).eq.0) 
C    $      write(6,6) iter,rnorm,alpha,beta,pap
        if (nid.eq.0.) 
     $        write(6,6) iter,rnorm,alpha,beta,pap

    6    format('cg:',i4,1p4e12.4)
c        if (rtr.le.rlim2) goto 1001

      enddo

 1001 continue

!$ACC END DATA

      if (nid.eq.0) write(6,6) iter,rnorm,alpha,beta,pap

      flop_cg = flop_cg + iter*15.*n

      return
      end


c-----------------------------------------------------------------------
      subroutine maskit_acc(w,pmask,nx,ny,nz)   ! Zero out Dirichlet conditions
      include 'SIZE'
      include 'PARALLEL'

      real pmask(-1:lx1*ly1*lz1*lelt)
      real w(nx1*ny1*nz1*lelt)
      integer i,j,e

      nxyz = nx*ny*nz
      nxy  = nx*ny

!$ACC DATA PRESENT(w,pmask)
      if(pmask(-1).lt.0) then
        j=pmask(0)
       
!$ACC PARALLEL LOOP
        do i = 1,j
           k = pmask(i)
           w(k)=0.0
        enddo

      else
         write(*,*) "OpenACC version is not implemented yet"
         stop
c         Zero out Dirichlet boundaries.
c
c                      +------+     ^ Y
c                     /   3  /|     |
c               4--> /      / |     |
c                   +------+ 2 +    +----> X
c                   |   5  |  /    /
c                   |      | /    /
c                   +------+     Z   
c

        nn = 0
!!$ACC PARALLEL LOOP 
        do e  = 1,nelt
           call get_face(w,nx,e)
!!$ACC LOOP
          do i = 1,nxyz
             if(w(i).eq.0) then
               nn=nn+1
               pmask(nn)=i
             endif
          enddo
        enddo     
        pmask(-1) = -1.
        pmask(0) = nn
      endif
!$ACC END DATA

      return
      end

c-----------------------------------------------------------------------
      subroutine solveM_acc(z,r,n)
      include 'INPUT'
      real z(n),r(n)

      nn = n
      call h1mg_solve_acc(z,r,nn)

      return
      end


#ifdef CRAYACC
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
      real gxyz(2*ldim,nx1,ny1,nz1,lelt)

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
            ur(i,j,k,e) = gxyz(1,i,j,k,e)*wr
     $                  + gxyz(2,i,j,k,e)*ws
     $                  + gxyz(3,i,j,k,e)*wt
            us(i,j,k,e) = gxyz(2,i,j,k,e)*wr
     $                  + gxyz(4,i,j,k,e)*ws
     $                  + gxyz(5,i,j,k,e)*wt
            ut(i,j,k,e) = gxyz(3,i,j,k,e)*wr
     $                  + gxyz(5,i,j,k,e)*ws
     $                  + gxyz(6,i,j,k,e)*wt
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

      call dssum_acc(w)         ! Gather-scatter operation  ! w   = QQ  w
                                                            !            L
      call add2s2_acc(w,u,.1,n)   !2n
      call maskit_acc(w,cmask,nx1,ny1,nz1)  ! Zero out Dirichlet conditions

!$ACC END DATA

      nxyz=nx1*ny1*nz1
      flop_a = flop_a + (19*nxyz+12*nx1*nxyz)*nelt

      return
      end

#else
c-----------------------------------------------------------------------
      subroutine ax_acc(w,u,gxyz,ur,us,ut,wk,n) ! Matrix-vector product: w=A*u

#ifdef TUNED_CUF_KERNEL
      use cudafor
#endif

      include 'SIZE'
      include 'TOTAL'

#ifdef TUNED_CUF_KERNEL
      interface
      attributes(global) subroutine ax_cuf2(w,u,ur,us,ut,
     &                gxyz,dxm1,dxtm1)

      real, intent(out) :: w(nx1,ny1,nz1,nelt)
      real, intent(in)  :: u(nx1,ny1,nz1,nelt)
      real ur  (nx1,ny1,nz1,lelt)
      real us  (nx1,ny1,nz1,lelt)
      real ut  (nx1,ny1,nz1,lelt)

      real gxyz(nx1,ny1,nz1,2*ldim,lelt)

      real, intent(in) :: dxm1(nx1,nx1)
      real, intent(in) :: dxtm1(nx1,nx1)
      end subroutine
      end interface
#endif
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

#ifdef TUNED_CUF_KERNEL

!$acc host_data use_device(w,u,ur,us,ut,gxyz,dxm1,dxtm1)
       if (nx1.eq.8) then
         call ax_cuf2<<<nelt,dim3(nx1,ny1,nz1)>>>(w,u,
     $                ur,us,ut,gxyz,dxm1,dxtm1)
       else if (nx1.eq.12) then
         call ax_cuf2<<<nelt,dim3(nx1,ny1,nz1/2)>>>(w,u,
     $                ur,us,ut,gxyz,dxm1,dxtm1)
       else if (nx1.eq.16) then
         call ax_cuf2<<<nelt,dim3(nx1,ny1,nz1/4)>>>(w,u,
     $                ur,us,ut,gxyz,dxm1,dxtm1)
       endif
!$acc end host_data

#else

!$ACC KERNELS
      do e = 1,nelt
!$ACC LOOP COLLAPSE(3)
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            wr = 0      ! scalar ur gets promoted to vector register over index
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
!$ACC END KERNELS

!$ACC KERNELS
      do e=1,nelt
!$ACC LOOP COLLAPSE(3)
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
!$ACC END KERNELS

#endif
      call dssum_acc(w)         ! Gather-scatter operation  ! w   = QQ  w
                                                            !            L
      call add2s2_acc(w,u,.1,n)   !2n
      call maskit_acc(w,cmask,nx1,ny1,nz1)  ! Zero out Dirichlet conditions

!$ACC END DATA

      nxyz=nx1*ny1*nz1
      flop_a = flop_a + (19*nxyz+12*nx1*nxyz)*nelt

      return
      end

#endif

#endif 
