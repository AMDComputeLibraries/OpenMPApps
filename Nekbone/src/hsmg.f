c control boundary based on a  mask, not hardcoded.
c-----------------------------------------------------------------------
      subroutine h1mg_setup()
      include 'SIZE'
      include 'TOTAL'
      include 'HSMG'
      integer p_msk

c     call geom_reset(1)  ! Recompute g1m1 etc. with deformed only

      mg_fld = 1
      call h1mg_index_0      ! set indices to 0 

      call h1mg_setup_mg_nx
      call h1mg_setup_semhat ! SEM hat matrices for each level
      call h1mg_setup_intp   ! Interpolation operators
      call h1mg_setup_dssum  ! set direct stiffness summation handles
      call h1mg_setup_wtmask ! set restriction weight matrices 

      call h1mg_setup_fdm    ! set up fast diagonalization method

      call h1mg_setup_schwarz_wt

      l=mg_h1_lmax
      call mg_set_msk (p_msk,l)


      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_index_0 ! initialize index sets
      include 'SIZE'
      include 'HSMG'

      n = lmgn*(2)   ! max_# of MG level*2

      call izero( mg_rstr_wt_index      , n ) ! Index for each
      call izero( mg_fast_s_index       , n )
      call izero( mg_fast_d_index       , n )
      call izero( mg_schwarz_wt_index   , n )
      
      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_mg_nx()
      include 'SIZE'
      include 'HSMG'
      integer nf,nc,nr
      integer nx,ny,nz

      integer mgn2(10)
      save    mgn2
      data    mgn2 / 1, 2, 2, 2, 2, 3, 3, 5, 5, 5/
c     data    mgn2 / 1, 2, 3, 4, 5, 6, 7, 8, 9, 0

      mg_h1_lmax = 3
c     mg_h1_lmax = 4
      if (nx1.eq.4) mg_h1_lmax = 2

      mgnx1    = 1
      mg_nx(1) = mgnx1
      mg_ny(1) = mgnx1
      mg_nz(1) = mgnx1

      mgnx2 = 2*(lx2/4) + 1
      if (nx1.eq.5)  mgnx2 = 3
      if (nx1.le.10) mgnx2 = mgn2(nx1)
      if (nx1.eq.8)  mgnx2 = 3

      mgnx2 = min(3,mgnx2)  ! This choice seems best (9/24/12)

      mg_nx(2) = mgnx2
      mg_ny(2) = mgnx2
      mg_nz(2) = mgnx2

      mg_nx(3) = mgnx2+1
      mg_ny(3) = mgnx2+1
      mg_nz(3) = mgnx2+1

      mg_nx(mg_h1_lmax) = nx1-1
      mg_ny(mg_h1_lmax) = ny1-1
      mg_nz(mg_h1_lmax) = nz1-1

      if (nid.eq.0) write(*,*) 'h1_mg_nx:',(mg_nx(i),i=1,mg_h1_lmax)
c     if (nid.eq.0) write(*,*) 'h1_mg_ny:',(mg_ny(i),i=1,mg_h1_lmax)
c     if (nid.eq.0) write(*,*) 'h1_mg_nz:',(mg_nz(i),i=1,mg_h1_lmax)

      do ifld=1,ldimt1
      do l=1,mg_h1_lmax
         mg_h1_n(l,ifld)=(mg_nx(l)+1)
     $                  *(mg_ny(l)+1)
     $                  *(mg_nz(l)+1)*nelt
      enddo
      enddo
      
      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_semhat ! SEM hat matrices for each level
      include 'SIZE'
      include 'HSMG'
      include 'SEMHAT'
c     n+1 sized ah arrays for n =nx of level

      do l=1,mg_h1_lmax
         n = mg_nx(l)     ! polynomial order
         call semhat(ah,bh,ch,dh,zh,wh,n)
         call copy(mg_ah(1,l),ah,(n+1)*(n+1))
         call copy(mg_bh(1,l),bh,n+1)
         call copy(mg_dh(1,l),dh,(n+1)*(n+1))
         call transpose(mg_dht(1,l),n+1,dh,n+1)
         call copy(mg_zh(1,l),zh,n+1)

         mg_nh(l)=n+1
         mg_nhz(l)=mg_nz(l)+1

      enddo
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_intp
      include 'SIZE'
      include 'HSMG'
      integer l,nf,nc

      do l=1,mg_h1_lmax-1

         nf=mg_nh(l+1)
         nc=mg_nh(l)

!        Standard multigrid coarse-to-fine interpolation
         call h1mg_setup_intpm(
     $           mg_jh(1,l),mg_zh(1,l+1),mg_zh(1,l),nf,nc)
         call transpose(mg_jht(1,l),nc,mg_jh(1,l),nf)

      enddo
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_intpm(jh,zf,zc,nf,nc)
      include 'SIZE'
      integer nf,nc
      real jh(nf,nc),zf(1),zc(1)
      real w(2*lx1+2)

      do i=1,nf
         call fd_weights_full(zf(i),zc,nc-1,1,w)
         do j=1,nc
            jh(i,j)=w(j)
         enddo
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_dssum
      include 'SIZE'
      include 'HSMG'

      integer nx
      

      do l=1,mg_h1_lmax  ! set up direct stiffness summation for each level
         nx=mg_nh(l)
         call proxy_setupds(mg_gsh_handle(l,mg_fld),nx) !assumes nx=ny=nz
         nx=nx+2
         call proxy_setupds(mg_gsh_schwarz_handle(l,mg_fld),nx)
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_wtmask
      include 'SIZE'
      include 'HSMG'
      integer i,l

      i=0
      do l=1,mg_h1_lmax
         mg_rstr_wt_index(l,mg_fld)=i
         i=i+mg_nh(l)*mg_nhz(l)*2*ndim*nelt
         if(i .gt. lmg_rwt*2*ldim*lelt) then
            itmp = i/(2*ldim*lelt)
            write(6,*) 'parameter lmg_rwt too small',i,itmp,lmg_rwt
            call exitt0
         endif
         call h1mg_setup_rstr_wt(
     $           mg_rstr_wt(mg_rstr_wt_index(l,mg_fld))
     $          ,mg_nh(l),mg_nh(l),mg_nhz(l),l,mg_work)
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_rstr_wt(wt,nx,ny,nz,l,w)
      include 'SIZE'

c     basically just the inverse counting matrix
      integer nx,ny,nz,l
      real w(nx,ny,nz,nelt)
      real wt(nx,nz,2,ndim,nelt)
      
      integer ie
      !init border nodes to 1
      call rzero(w,nx*ny*nz*nelt)

      do ie=1,nelt
         do j=1,ny
         do i=1,nx
            w(i,j,1,ie)=1.0
            w(i,j,nz,ie)=1.0
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            w(i,1,k,ie)=1.0
            w(i,ny,k,ie)=1.0
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            w(1,j,k,ie)=1.0
            w(nx,j,k,ie)=1.0
         enddo
         enddo
      enddo

      call h1mg_dssum(w,l)

      !invert the count w to get the weight wt
      do ie=1,nelt
         do k=1,nz
         do j=1,ny
            wt(j,k,1,1,ie)=1.0/w(1,j,k,ie)
            wt(j,k,2,1,ie)=1.0/w(nx,j,k,ie)
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            wt(i,k,1,2,ie)=1.0/w(i,1,k,ie)
            wt(i,k,2,2,ie)=1.0/w(i,ny,k,ie)
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            wt(i,j,1,3,ie)=1.0/w(i,j,1,ie)
            wt(i,j,2,3,ie)=1.0/w(i,j,nz,ie)
         enddo
         enddo
      enddo
      
      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_fdm()
      include 'SIZE'
      include 'HSMG'
      
      integer l,i,j,nl

      i = mg_fast_s_index(mg_h1_lmax,mg_fld-1)
      j = mg_fast_d_index(mg_h1_lmax,mg_fld-1)
      do l=2,mg_h1_lmax
         mg_fast_s_index(l,mg_fld)=i
         nl = mg_nh(l)+2
         i=i+nl*nl*2*ndim*nelt
         if(i .gt. lmg_fasts*2*ldim*lelt) then
            itmp = i/(2*ldim*lelt)
            write(6,*) 'lmg_fasts too small',i,itmp,lmg_fasts,l
            call exitt0
         endif
         mg_fast_d_index(l,mg_fld)=j
         j=j+(nl**ndim)*nelt
         if(j .gt. lmg_fastd*lelt) then
            itmp = i/(2*ldim*lelt)
            write(6,*) 'lmg_fastd too small',i,itmp,lmg_fastd,l
            call exitt0
         endif
         call h1mg_setup_fast(
     $             mg_fast_s(mg_fast_s_index(l,mg_fld))
     $            ,mg_fast_d(mg_fast_d_index(l,mg_fld))
     $            ,mg_nh(l)+2,mg_ah(1,l),mg_bh(1,l),mg_nx(l))
      enddo
      mg_fast_s_index(l,mg_fld)=i
      mg_fast_d_index(l,mg_fld)=j
      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_fast(s,d,nl,ah,bh,n)
      include 'SIZE'
      include 'PARALLEL'
      include 'HSMG'

      real s(nl*nl,2,ndim,nelt)
      real d(nl**ndim,nelt)
      real ah(1),bh(1)

      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
      real lr ,ls ,lt
 
      real llr(lelt),lls(lelt),llt(lelt)
     $   , lmr(lelt),lms(lelt),lmt(lelt)
     $   , lrr(lelt),lrs(lelt),lrt(lelt)
      integer lbr,rbr,lbs,rbs,lbt,rbt
      
      real l(nx1,ny1,nz1,nelt)
      integer i,j,k
      integer ie,il,nr,ns,nt
      real eps,diag
      
c     zero if on edge, lmr of neigh ow, only should be done once and saved?
      do ie = 1,nelt
c         lmr(ie) = 1.0/nelx  !even element distribution.
c         lms(ie) = 1.0/nely  !each element is 1/nelx in width
c         lmt(ie) = 1.0/nelz
          lmr(ie) = 2.0  !even element distribution.
          lms(ie) = 2.0  !each element is (-1,1) in width
          lmt(ie) = 2.0
        

          l(1,  2,  2,ie)=lmr(ie)
          l(nx1,2,  2,ie)=lmr(ie)
          l(2,  1,  2,ie)=lms(ie)
          l(2,ny1,  2,ie)=lms(ie)
          l(2,  2,  1,ie)=lmt(ie)
          l(2,  2,nz1,ie)=lmt(ie)
      enddo
      call dssum(l)
      do ie = 1,nelt
         llr(ie) = l(1,  2,  2,ie)-lmr(ie)
         lrr(ie) = l(nx1,2,  2,ie)-lmr(ie)
         lls(ie) = l(2,  1,  2,ie)-lms(ie)
         lrs(ie) = l(2,ny1,  2,ie)-lms(ie)
         llt(ie) = l(2,  2,  1,ie)-lmt(ie)
         lrt(ie) = l(2,  2,nz1,ie)-lmt(ie)
c          write(6,*) ie, ' ele'
c          write(6,*) llr(ie),lmr(ie),lrr(ie)
c          write(6,*) lls(ie),lms(ie),lrs(ie)
c          write(6,*) llt(ie),lmt(ie),lrt(ie)
c          write(6,*)
      enddo
 
      do ie=1,nelt
         call get_bcs(lbr,rbr,lbs,rbs,lbt,rbt,ie)
         

         nr=nl
         ns=nl
         nt=nl

         call h1mg_setup_fast1d(s(1,1,1,ie),lr,nr,lbr,rbr
     $            ,llr(ie),lmr(ie),lrr(ie),ah,bh,n,ie)
         call h1mg_setup_fast1d(s(1,1,2,ie),ls,ns,lbs,rbs
     $            ,lls(ie),lms(ie),lrs(ie),ah,bh,n,ie)
         call h1mg_setup_fast1d(s(1,1,3,ie),lt,nt,lbt,rbt
     $                     ,llt(ie),lmt(ie),lrt(ie),ah,bh,n,ie)
         il=1
         eps = 1.e-5 * (vlmax(lr(2),nr-2)
     $               + vlmax(ls(2),ns-2) + vlmax(lt(2),nt-2))
         do k=1,nt
         do j=1,ns
         do i=1,nr
            diag = lr(i)+ls(j)+lt(k)
            if (diag.gt.eps) then
                d(il,ie) = 1.0/diag
            else
c               write(6,3) ie,'Reset Eig in h1mg setup fast:',i,j,k,l
c    $                         ,eps,diag,lr(i),ls(j),lt(k)
    3           format(i6,1x,a21,4i5,1p5e12.4)
                d(il,ie) = 0.0
            endif
            il=il+1
         enddo
         enddo
         enddo
      enddo


      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_fast1d(s,lam,nl,lbc,rbc,ll,lm,lr,ah,bh,n,ie)
      integer nl,lbc,rbc,n
      real s(nl,nl,2),lam(nl),ll,lm,lr
      real ah(0:n,0:n),bh(0:n)
      
      include 'SIZE'
      parameter(lxm=lx1+2)
      common /ctmp0/ b(2*lxm*lxm),w(2*lxm*lxm)
      
      call h1mg_setup_fast1d_a(s,lbc,rbc,ll,lm,lr,ah,n)
      call h1mg_setup_fast1d_b(b,lbc,rbc,ll,lm,lr,bh,n)
      
      call generalev(s,b,lam,nl,w)

      if(lbc.gt.0) call row_zero(s,nl,nl,1)   
      if(lbc.eq.1) call row_zero(s,nl,nl,2)   

      if(rbc.gt.0) call row_zero(s,nl,nl,nl)  
      if(rbc.eq.1) call row_zero(s,nl,nl,nl-1)
      
      call transpose(s(1,1,2),nl,s,nl)
      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_fast1d_a(a,lbc,rbc,ll,lm,lr,ah,n)
      integer lbc,rbc,n
      real a(0:n+2,0:n+2),ll,lm,lr
      real ah(0:n,0:n)
      
      real fac
      integer i,j,i0,i1

      i0=0
      if(lbc.eq.1) i0=1    
      i1=n                
      if(rbc.eq.1) i1=n-1 
      
      call rzero(a,(n+3)*(n+3))
      fac = 2.0/lm
      a(1,1)=1.0
      a(n+1,n+1)=1.0
      do j=i0,i1
         do i=i0,i1
            a(i+1,j+1)=fac*ah(i,j)
         enddo
      enddo
      if(lbc.eq.0) then
         fac = 2.0/ll
         a(0,0)=fac*ah(n-1,n-1)
         a(1,0)=fac*ah(n  ,n-1)
         a(0,1)=fac*ah(n-1,n  )
         a(1,1)=a(1,1)+fac*ah(n  ,n  )
      else
         a(0,0)=1.0
      endif
      if(rbc.eq.0) then
         fac = 2.0/lr
         a(n+1,n+1)=a(n+1,n+1)+fac*ah(0,0)
         a(n+2,n+1)=fac*ah(1,0)
         a(n+1,n+2)=fac*ah(0,1)
         a(n+2,n+2)=fac*ah(1,1)
      else
         a(n+2,n+2)=1.0
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_fast1d_b(b,lbc,rbc,ll,lm,lr,bh,n)
      integer lbc,rbc,n
      real b(0:n+2,0:n+2),ll,lm,lr
      real bh(0:n)
      
      real fac
      integer i,j,i0,i1
      i0=0
      if(lbc.eq.1) i0=1   
      i1=n                
      if(rbc.eq.1) i1=n-1 
      
      call rzero(b,(n+3)*(n+3))
      fac = 0.5*lm
      b(1,1)=1.0
      b(n+1,n+1)=1.0
      do i=i0,i1
         b(i+1,i+1)=fac*bh(i)
      enddo
      if(lbc.eq.0) then
         fac = 0.5*ll
         b(0,0)=fac*bh(n-1)
         b(1,1)=b(1,1)+fac*bh(n  )
      else
         b(0,0)=1.0
      endif
      if(rbc.eq.0) then
         fac = 0.5*lr
         b(n+1,n+1)=b(n+1,n+1)+fac*bh(0)
         b(n+2,n+2)=fac*bh(1)
      else
         b(n+2,n+2)=1.0
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine generalev(a,b,lam,n,w)
c     Solve the generalized eigenvalue problem  A x = lam B x
c
c     A -- symm.
c     B -- symm., pos. definite
c
c
      include 'SIZE'
      include 'PARALLEL'
 
      real a(n,n),b(n,n),lam(n),w(n,n)
      real aa(100),bb(100)
 
      parameter (lbw=4*lx1*ly1*lz1*lelt)
      common /bigw/ bw(lbw)
 
      lw = n*n

      call copy(aa,a,100)
      call copy(bb,b,100)
      if (ifdblas) then
         call dsygv(1,'V','U',n,a,n,b,n,lam,bw,lbw,info)
      else
         call ssygv(1,'V','U',n,a,n,b,n,lam,bw,lbw,info)
      endif
 
      if (info.ne.0) then
         if (nid.eq.0) then
            call outmat2(aa ,n,n,n,'aa  ')
            call outmat2(bb ,n,n,n,'bb  ')
            call outmat2(a  ,n,n,n,'Aeig')
            call outmat2(lam,1,n,n,'Deig')
         endif
 
         ninf = n-info
         write(6,*) 'Error in generalev, info=',info,n,ninf
         call exitt0
      endif
 
      return
      end
c-----------------------------------------------------------------------
      subroutine row_zero(a,m,n,e)
      integer m,n,e
      real a(m,n)
      do j=1,n
         a(e,j)=0.
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine get_bcs(lbr,rbr,lbs,rbs,lbt,rbt,ie)
c     based on Dirichlet boundary conditions set on all 6faces
      include 'SIZE'
      include 'PARALLEL'
      integer lbr,rbr,lbs,rbs,lbt,rbt

      lbr=0
      rbr=0
      lbs=0
      rbs=0
      lbt=0
      rbt=0

      nelxy=nelx*nely
      ngl = lglel(ie)

      ir = 1+(ngl-1)/nelxy
      iq = mod1(ngl,nelxy)
      iq = 1+(iq-1)/nelx
      ip = mod1(ngl,nelx)


      if(mod(ip,nelx).eq.1.or.nelx.eq.1) lbr=1
      if(mod(ip,nelx).eq.0)              rbr=1

      if(mod(iq,nely).eq.1.or.nely.eq.1) lbs=1
      if(mod(iq,nely).eq.0)              rbs=1

      if(mod(ir,nelz).eq.1.or.nelz.eq.1) lbt=1
      if(mod(ir,nelz).eq.0)              rbt=1 

c     write(6,1) ngl, lbr,rbr,lbs,rbs,lbt,rbt
c 1   format(i3,' element',6i3)


      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_setup_schwarz_wt
      include 'SIZE'
      include 'HSMG'
      
      integer l,i,nl,nlz

      i = mg_schwarz_wt_index(mg_h1_lmax,mg_fld-1)
      do l=2,mg_h1_lmax

         mg_schwarz_wt_index(l,mg_fld)=i
         nl  = mg_nh(l)
         nlz = mg_nhz(l)
         i   = i+nl*nlz*4*ndim*nelt

         if (i .gt. lmg_swt*4*ldim*lelt) then
            itmp = i/(4*ldim*lelt)
            write(6,*) 'lmg_swt too small',i,itmp,lmg_swt,l
            call exitt0
         endif

c        call h1mg_setup_schwarz_wt3d(
c    $          mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld))
c    $         ,nl,mg_worke)
         call h1mg_setup_schwarz_wt3d_1(
     $          mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld))
     $         ,l)

      enddo

      mg_schwarz_wt_index(l,mg_fld)=i

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_schwarz_wt3d_1(wt,l)
      include 'SIZE'
      include 'HSMG'

      real wt(1),work(1)

      integer enx,eny,enz

      zero =  0
      one  =  1
      onem = -1

      n  = mg_h1_n (l,mg_fld)

      enx=mg_nh(l)+2
      eny=mg_nh(l)+2
      enz=mg_nh(l)+2
      ns = enx*eny*enz*nelt
      i  = ns+1

      call rone(mg_work(i),ns)

c     Sum overlap region (border excluded)
      call h1mg_extrude(mg_work,0,zero,mg_work(i),0,one ,enx,eny,enz)
      call h1mg_schwarz_dssum(mg_work(i),l)
      call h1mg_extrude(mg_work(i),0,one ,mg_work,0,onem,enx,eny,enz)
      call h1mg_extrude(mg_work(i),2,one,mg_work(i),0,one,enx,eny,enz)

      call h1mg_schwarz_toreg3d(mg_work,mg_work(i),mg_nh(l))

      call h1mg_dssum(mg_work,l)                           ! sum border nodes


      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nh(l)
      nxyz = nx*ny*nz
      k    = 1
      do ie=1,nelt
c        call outmat(mg_work(k),nx,ny,'NEW WT',ie)
         call h1mg_setup_schwarz_wt3d_2(wt,ie,nx,mg_work(k))
         k = k+nxyz
      enddo


      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_schwarz_wt3d_2(wt,ie,n,work)
      include 'SIZE'
      integer n
      real wt(n,n,4,3,nelt)
      real work(n,n,n)

      integer ie,i,j,k

      
      do k=1,n
      do j=1,n
         wt(j,k,1,1,ie)=1.0/work(1,j,k)
         wt(j,k,2,1,ie)=1.0/work(2,j,k)
         wt(j,k,3,1,ie)=1.0/work(n-1,j,k)
         wt(j,k,4,1,ie)=1.0/work(n,j,k)
      enddo
      enddo
      do k=1,n
      do i=1,n
         wt(i,k,1,2,ie)=1.0/work(i,1,k)
         wt(i,k,2,2,ie)=1.0/work(i,2,k)
         wt(i,k,3,2,ie)=1.0/work(i,n-1,k)
         wt(i,k,4,2,ie)=1.0/work(i,n,k)
      enddo
      enddo
      do j=1,n
      do i=1,n
         wt(i,j,1,3,ie)=1.0/work(i,j,1)
         wt(i,j,2,3,ie)=1.0/work(i,j,2)
         wt(i,j,3,3,ie)=1.0/work(i,j,n-1)
         wt(i,j,4,3,ie)=1.0/work(i,j,n)
      enddo
      enddo
      return
      end

c----------------------------------------------------------------------
      subroutine h1mg_setup_schwarz_wt3d(wt,n,work)
      include 'SIZE'
      include 'PARALLEL'
      integer n
      real wt(n,n,4,3,nelt)
      real work(n,n,n)
      
      integer ie,i,j,k
      integer lbr,rbr,lbs,rbs,lbt,rbt

      logical ifsqrt 
 
      ifsqrt = .true.
      ifsqrt = .false.

      do ie=1,nelt
         do k=1,n
         do j=1,n
            work(1,j,k)=1.0
            work(2,j,k)=1.0
            work(n-1,j,k)=1.0
            work(n,j,k)=1.0
         enddo
         enddo
         do k=1,n
         do i=1,n
            work(i,1,k)=1.0
            work(i,2,k)=1.0
            work(i,n-1,k)=1.0
            work(i,n,k)=1.0
         enddo
         enddo
         do j=1,n
         do i=1,n
            work(i,j,1)=1.0
            work(i,j,2)=1.0
            work(i,j,n-1)=1.0
            work(i,j,n)=1.0
         enddo
         enddo

         call get_bcs(lbr,rbr,lbs,rbs,lbt,rbt,ie)
         
         if(lbr.eq.0) then
            do k=1,n
            do j=1,n
               work(1,j,k)=work(1,j,k)+1.0
               work(2,j,k)=work(2,j,k)+1.0
            enddo
            enddo
         endif
         if(rbr.eq.0) then
            do k=1,n
            do j=1,n
               work(n-1,j,k)=work(n-1,j,k)+1.0
               work(n,j,k)=work(n,j,k)+1.0
            enddo
            enddo
         endif
         if(lbs.eq.0) then
            do k=1,n
            do i=1,n
               work(i,1,k)=work(i,1,k)+1.0
               work(i,2,k)=work(i,2,k)+1.0
            enddo
            enddo
         endif
         if(rbs.eq.0) then
            do k=1,n
            do i=1,n
               work(i,n-1,k)=work(i,n-1,k)+1.0
               work(i,n,k)=work(i,n,k)+1.0
            enddo
            enddo
         endif
         if(lbt.eq.0) then
            do j=1,n
            do i=1,n
               work(i,j,1)=work(i,j,1)+1.0
               work(i,j,2)=work(i,j,2)+1.0
            enddo
            enddo
         endif
         if(rbt.eq.0) then
            do j=1,n
            do i=1,n
               work(i,j,n-1)=work(i,j,n-1)+1.0
               work(i,j,n)=work(i,j,n)+1.0
            enddo
            enddo
         endif
         do k=1,n
         do j=1,n
            wt(j,k,1,1,ie)=1.0/work(1,j,k)
            wt(j,k,2,1,ie)=1.0/work(2,j,k)
            wt(j,k,3,1,ie)=1.0/work(n-1,j,k)
            wt(j,k,4,1,ie)=1.0/work(n,j,k)
         enddo
         enddo
         do k=1,n
         do i=1,n
            wt(i,k,1,2,ie)=1.0/work(i,1,k)
            wt(i,k,2,2,ie)=1.0/work(i,2,k)
            wt(i,k,3,2,ie)=1.0/work(i,n-1,k)
            wt(i,k,4,2,ie)=1.0/work(i,n,k)
         enddo
         enddo
         do j=1,n
         do i=1,n
            wt(i,j,1,3,ie)=1.0/work(i,j,1)
            wt(i,j,2,3,ie)=1.0/work(i,j,2)
            wt(i,j,3,3,ie)=1.0/work(i,j,n-1)
            wt(i,j,4,3,ie)=1.0/work(i,j,n)
         enddo
         enddo
         if(ifsqrt) then
            do ii=1,3
            do k=1,4
            do j=1,4
            do i=1,n
               wt(i,j,k,ii,ie)=sqrt(wt(i,j,k,ii,ie))
            enddo
            enddo
            enddo
            enddo
         endif

      enddo


      return
      end
c----------------------------------------------------------------------
      subroutine mg_set_msk(p_msk ,l0)
      include 'SIZE'
      include 'HSMG'
      integer p_msk

      l                  = mg_h1_lmax
      p_mg_msk(l,mg_fld) = 0
      n                  = mg_h1_n(l,mg_fld)


      do l=mg_h1_lmax,1,-1
         nx = mg_nh  (l)
         ny = mg_nh  (l)
         nz = mg_nhz (l)

         p_msk = p_mg_msk(l,mg_fld)

         call h1mg_setup_mask
     $     (mg_imask(p_msk),nm,nx,ny,nz,nelt,l,mg_work)

         if (l.gt.1) p_mg_msk(l-1,mg_fld)=p_mg_msk(l,mg_fld)+nm

      enddo

      p_msk = p_mg_msk(l0,mg_fld)

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_mask(mask,nm,nx,ny,nz,nel,l,w)
      include 'SIZE'

      integer mask(1)        ! Pointer to Dirichlet BCs
      integer nx,ny,nz,l
      real w(nx,ny,nz,nel)
      real cmask(-1:lx1*ly1*lz1*lelt)
      
      integer e,count,ptr

      zero = 0
      nxyz = nx*ny*nz
      n    = nx*ny*nz*nel

      call rone(w,n)   ! Init everything to 1

      cmask(-1) = 0.0
      call maskit(w,cmask,nx,ny,nz)

c     call h1mg_dsprod(w,l)    ! direct stiffness multiply
c
c     Prototypical mask layout, nel=5:
c
c    e=1 ...             10
c      1  2  3  4  5 ... 10 | 11 12 13 14 | 15 | 16 |
c     11 15 16 ...          |  3 p1 p2 p3 |  0 |  0 | ...
c                              ^
c                              |
c                              |_count for e=1
c

      nm  = 1                  ! store mask
      do e=1,nel

         mask(e) = nel+nm
         count   = 0          ! # Dirchlet points on element e
         ptr     = mask(e)

         do i=1,nxyz
            if (w(i,1,1,e).eq.0) then
               nm    = nm   +1
               count = count+1
               ptr   = ptr  +1
               mask(ptr) = i + nxyz*(e-1)   ! where I mask on element e 
            endif
         enddo


         ptr       = mask(e)
         mask(ptr) = count

         nm        = nm+1     ! bump pointer to hold next count

      enddo

      nm = nel + nm-1 ! Return total number of mask pointers/counters

      return
      end
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine h1mg_dssum(u,l)
      include 'SIZE'
      include 'HSMG'

      call adelay
      call gs_op(mg_gsh_handle(l,mg_fld),u,1,1,0)

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_dsprod(u,l)
      include 'SIZE'
      include 'HSMG'
      real u(1)

      call adelay
      call gs_op(mg_gsh_handle(l,mg_fld),u,1,2,0)

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_schwarz_dssum(u,l)
      include 'SIZE'
      include 'HSMG'

      call adelay
      call gs_op(mg_gsh_schwarz_handle(l,mg_fld),u,1,1,0)
      return
      end
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine hsmg_coarse_solve(e,r)
      include 'SIZE'
      include 'HSMG'
      real e(1),r(1)

       n=mg_nh(1)
       nxyz = n*n*n*nelt
       call copy(e,r,nxyz)
       call h1mg_dssum(e,1)


      return
      end
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine h1mg_solve(z,rhs,nn)  !  Solve preconditioner: Mz=rhs
      real z(1),rhs(1)

c     Assumes that preprocessing has been completed via h1mg_setup()


      include 'SIZE'
      include 'HSMG'       ! Same array space as HSMG
      include 'TOTAL'
      
      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrmg/ e(2*lt),w(lt),r(lt)
      integer p_msk


      nel   = nelt

      op    =  1.                                     ! Coefficients for h1mg_ax
      om    = -1.
      sigma =  1.

      l     = mg_h1_lmax
      n     = mg_h1_n(l,mg_fld)
      is    = 1                                       ! solve index

      call h1mg_schwarz(z,rhs,sigma,l)                ! z := sigma W M       rhs
                                                      !               Schwarz
      call copy(r,rhs,n)                              ! r  := rhs

      do l = mg_h1_lmax-1,2,-1                        ! DOWNWARD Leg of V-cycle
         is = is + n
         n  = mg_h1_n(l,mg_fld)
                                                      !          T
         call h1mg_rstr(r,l,.true.)                   ! r   :=  J r
                                                      !  l         l+1
!        OVERLAPPING Schwarz exchange and solve:
         call h1mg_schwarz(e(is),r,sigma,l)           ! e := sigma W M       r
                                                      !  l            Schwarz l
      enddo
      
      is = is+n
                                                      !         T
      call h1mg_rstr(r,1,.false.)                     ! r  :=  J  r
                                                      !  l         l+1

      p_msk = p_mg_msk(l,mg_fld)
      call h1mg_mask(r,mg_imask(p_msk),nel)           !        -1
      call hsmg_coarse_solve ( e(is) , r )            ! e  := A   r
      call h1mg_mask(e(is),mg_imask(p_msk),nel)       !  1     1   1

c     nx = mg_nh(1)
c     call h1mg_mask(e(is),mg_imask(p_msk),nel)       !  1     1   1
c     call exitt

      do l = 2,mg_h1_lmax-1                           ! UNWIND.  No smoothing.
         im = is
         is = is - n
         n  = mg_h1_n(l,mg_fld)
         call h1mg_intp (w,e(im),l-1)                 ! w   :=  J e
         i1=is-1                                      !            l-1
         do i=1,n
            e(i1+i) = e(i1+i) + w(i)                  ! e   :=  e  + w
         enddo                                        !  l       l
      enddo

      l  = mg_h1_lmax
      n  = mg_h1_n(l,mg_fld)
      im = is  ! solve index
      call h1mg_intp(w,e(im),l-1)                     ! w   :=  J e
      do i = 1,n                                      !            l-1
         z(i) = z(i) + w(i)                           ! z := z + w
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz(e,r,sigma,l)
      include 'SIZE'
      include 'HSMG'

      real e(1),r(1)

      n = mg_h1_n(l,mg_fld)

      call h1mg_schwarz_wt    (e,l)          ! e  := W^.5* e
      call h1mg_schwarz_part1 (e,r,l)        !  l           l
      call h1mg_schwarz_wt    (e,l)          ! e  :=  e *W^.5
      call cmult              (e,sigma,n)    !  l      l

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz_part1 (e,r,l)
      include 'SIZE'
      include 'HSMG'

      real e(1),r(1)

      integer enx,eny,enz,pm

      zero =  0
      one  =  1
      onem = -1

      n  = mg_h1_n (l,mg_fld)
      pm = p_mg_msk(l,mg_fld)

      call h1mg_mask  (r,mg_imask(pm),nelt)  ! Zero Dirichlet nodes

      call h1mg_schwarz_toext3d(mg_work,r,mg_nh(l))

      enx=mg_nh(l)+2
      eny=mg_nh(l)+2
      enz=mg_nh(l)+2
      i = enx*eny*enz*nelt+1
 
c     exchange interior nodes
      call h1mg_extrude(mg_work,0,zero,mg_work,2,one,enx,eny,enz)
      call h1mg_schwarz_dssum(mg_work,l)
      call h1mg_extrude(mg_work,0,one ,mg_work,2,onem,enx,eny,enz)

      call h1mg_fdm(mg_work(i),mg_work,l) ! Do the local solves

c     Sum overlap region (border excluded)
      call h1mg_extrude(mg_work,0,zero,mg_work(i),0,one ,enx,eny,enz)
      call h1mg_schwarz_dssum(mg_work(i),l)
      call h1mg_extrude(mg_work(i),0,one ,mg_work,0,onem,enx,eny,enz)
      call h1mg_extrude(mg_work(i),2,one,mg_work(i),0,one,enx,eny,enz)

      call h1mg_schwarz_toreg3d(e,mg_work(i),mg_nh(l))

      call h1mg_dssum(e,l)                           ! sum border nodes
      call h1mg_mask (e,mg_imask(pm),nelt) ! apply mask 

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz_wt(e,l)
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      
      call h1mg_schwarz_wt3d2(
     $    e,mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),mg_nh(l))
      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz_wt3d(e,wt,n)
      include 'SIZE'
      integer n
      real e(n,n,n,nelt)
      real wt(n,n,4,3,nelt)
      
      integer ie,i,j,k

      do ie=1,nelt
         do k=1,n
         do j=1,n
            e(1  ,j,k,ie)=e(1  ,j,k,ie)*wt(j,k,1,1,ie)
            e(2  ,j,k,ie)=e(2  ,j,k,ie)*wt(j,k,2,1,ie)
            e(n-1,j,k,ie)=e(n-1,j,k,ie)*wt(j,k,3,1,ie)
            e(n  ,j,k,ie)=e(n  ,j,k,ie)*wt(j,k,4,1,ie)
         enddo
         enddo
         do k=1,n
         do i=3,n-2
            e(i,1  ,k,ie)=e(i,1  ,k,ie)*wt(i,k,1,2,ie)
            e(i,2  ,k,ie)=e(i,2  ,k,ie)*wt(i,k,2,2,ie)
            e(i,n-1,k,ie)=e(i,n-1,k,ie)*wt(i,k,3,2,ie)
            e(i,n  ,k,ie)=e(i,n  ,k,ie)*wt(i,k,4,2,ie)
         enddo
         enddo
         do j=3,n-2
         do i=3,n-2
            e(i,j,1  ,ie)=e(i,j,1  ,ie)*wt(i,j,1,3,ie)
            e(i,j,2  ,ie)=e(i,j,2  ,ie)*wt(i,j,2,3,ie)
            e(i,j,n-1,ie)=e(i,j,n-1,ie)*wt(i,j,3,3,ie)
            e(i,j,n  ,ie)=e(i,j,n  ,ie)*wt(i,j,4,3,ie)
         enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz_wt3d2(e,wt,n)
      include 'SIZE'
      integer n
      real e(n,n,n,nelt)
      real wt(n,n,4,3,nelt)
      
      integer ie,i,j,k
      do ie=1,nelt
         do k=1,n
         do j=1,n
            e(1  ,j,k,ie)=e(1  ,j,k,ie)*sqrt(wt(j,k,1,1,ie))
            e(2  ,j,k,ie)=e(2  ,j,k,ie)*sqrt(wt(j,k,2,1,ie))
            e(n-1,j,k,ie)=e(n-1,j,k,ie)*sqrt(wt(j,k,3,1,ie))
            e(n  ,j,k,ie)=e(n  ,j,k,ie)*sqrt(wt(j,k,4,1,ie))
         enddo
         enddo
         do k=1,n
         do i=3,n-2
            e(i,1  ,k,ie)=e(i,1  ,k,ie)*sqrt(wt(i,k,1,2,ie))
            e(i,2  ,k,ie)=e(i,2  ,k,ie)*sqrt(wt(i,k,2,2,ie))
            e(i,n-1,k,ie)=e(i,n-1,k,ie)*sqrt(wt(i,k,3,2,ie))
            e(i,n  ,k,ie)=e(i,n  ,k,ie)*sqrt(wt(i,k,4,2,ie))
         enddo
         enddo
         do j=3,n-2
         do i=3,n-2
            e(i,j,1  ,ie)=e(i,j,1  ,ie)*sqrt(wt(i,j,1,3,ie))
            e(i,j,2  ,ie)=e(i,j,2  ,ie)*sqrt(wt(i,j,2,3,ie))
            e(i,j,n-1,ie)=e(i,j,n-1,ie)*sqrt(wt(i,j,3,3,ie))
            e(i,j,n  ,ie)=e(i,j,n  ,ie)*sqrt(wt(i,j,4,3,ie))
         enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_mask(w,mask,nel)
      include 'SIZE'

      real    w   (1)
      integer mask(1)        ! Pointer to Dirichlet BCs
      integer e
      
      do e=1,nel
         im = mask(e)
         call mg_mask_e(w,mask(im)) ! Zero out Dirichlet conditions
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mg_mask_e(w,mask) ! Zero out Dirichlet conditions
      real w(1)
      integer mask(0:1)

      n=mask(0)
      do i=1,n
         w(mask(i)) = 0.
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz_toreg3d(b,a,n)
c     strip off ghost cell
      include 'SIZE'
      integer n
      real a(0:n+1,0:n+1,0:n+1,nelt),b(n,n,n,nelt)
      
      integer i,j,k,ie
      do ie=1,nelt
      do k=1,n
      do j=1,n
      do i=1,n
         b(i,j,k,ie)=a(i,j,k,ie)
      enddo
      enddo
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz_toext3d(a,b,n)
c     border nodes (ghost cell = 0)
      include 'SIZE'
      integer n
      real a(0:n+1,0:n+1,0:n+1,nelt),b(n,n,n,nelt)
      
      integer i,j,k,ie
      call rzero(a,(n+2)*(n+2)*(n+2)*nelt)
      do ie=1,nelt
      do k=1,n
      do j=1,n
      do i=1,n
         a(i,j,k,ie)=b(i,j,k,ie)
      enddo
      enddo
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_extrude(arr1,l1,f1,arr2,l2,f2,nx,ny,nz)
      include 'SIZE'
      integer l1,l2,nx,ny,nz
      real arr1(nx,ny,nz,nelt),arr2(nx,ny,nz,nelt)
      real f1,f2
      
      integer i,j,k,ie,i0,i1
      i0=2
      i1=nx-1
      
      do ie=1,nelt
         do k=i0,i1
         do j=i0,i1
            arr1(l1+1 ,j,k,ie) = f1*arr1(l1+1 ,j,k,ie)
     $                          +f2*arr2(l2+1 ,j,k,ie)
            arr1(nx-l1,j,k,ie) = f1*arr1(nx-l1,j,k,ie)
     $                          +f2*arr2(nx-l2,j,k,ie)
         enddo
         enddo
         do k=i0,i1
         do i=i0,i1
            arr1(i,l1+1 ,k,ie) = f1*arr1(i,l1+1 ,k,ie)
     $                          +f2*arr2(i,l2+1 ,k,ie)
            arr1(i,nx-l1,k,ie) = f1*arr1(i,nx-l1,k,ie)
     $                          +f2*arr2(i,nx-l2,k,ie)
         enddo
         enddo
         do j=i0,i1
         do i=i0,i1
            arr1(i,j,l1+1 ,ie) = f1*arr1(i,j,l1+1 ,ie)
     $                          +f2*arr2(i,j,l2+1 ,ie)
            arr1(i,j,nx-l1,ie) = f1*arr1(i,j,nx-l1,ie)
     $                          +f2*arr2(i,j,nx-l2,ie)
         enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
c     clobbers r
      subroutine h1mg_fdm(e,r,l)
      include 'SIZE'
      include 'HSMG'
      call h1mg_do_fast(e,r,
     $      mg_fast_s(mg_fast_s_index(l,mg_fld)),
     $      mg_fast_d(mg_fast_d_index(l,mg_fld)),
     $      mg_nh(l)+2)
      return
      end
c-----------------------------------------------------------------------
c     clobbers r
      subroutine h1mg_do_fast(e,r,s,d,nl)
      include 'SIZE'
      real e(nl**ndim,nelt)
      real r(nl**ndim,nelt)
      real s(nl*nl,2,ndim,nelt)
      real d(nl**ndim,nelt)
      
      integer ie,nn,i
      nn=nl**ndim

      do ie=1,nelt
         call h1mg_tnsr3d_el(e(1,ie),nl,r(1,ie),nl
     $                      ,s(1,2,1,ie),s(1,1,2,ie),s(1,1,3,ie))
         do i=1,nn
            r(i,ie)=d(i,ie)*e(i,ie)
         enddo
         call h1mg_tnsr3d_el(e(1,ie),nl,r(1,ie),nl
     $                      ,s(1,1,1,ie),s(1,2,2,ie),s(1,2,3,ie))
      enddo

      return
      end
c----------------------------------------------------------------------
c     computes
c     v = [C (x) B (x) A] u
      subroutine h1mg_tnsr3d_el(v,nv,u,nu,A,Bt,Ct)

      integer nv,nu
      real v(nv*nv*nv),u(nu*nu*nu),A(1),Bt(1),Ct(1)
      include 'SIZE'
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      integer i

      call mxm(A,nv,u,nu,work,nu*nu)
      do i=0,nu-1
         call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
      enddo
      call mxm(work2,nv*nv,Ct,nu,v,nv)

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_rstr(r,l,ifdssum) ! r =J r,   l is coarse level
      include 'SIZE'
      include 'HSMG'
      logical ifdssum

      real r(1)
      integer l

      call h1mg_do_wt(r,mg_rstr_wt(mg_rstr_wt_index(l+1,mg_fld))
     $                     ,mg_nh(l+1),mg_nh(l+1),mg_nhz(l+1))

      call h1mg_tnsr1(r,mg_nh(l),mg_nh(l+1),mg_jht(1,l),mg_jh(1,l))

      if (ifdssum) call h1mg_dssum(r,l)

      return
      end
c-----------------------------------------------------------------------
c     u = wt .* u
      subroutine h1mg_do_wt(u,wt,nx,ny,nz)
      include 'SIZE'
      integer nx,ny,nz
      real u(nx,ny,nz,nelt)
      real wt(nx,nz,2,ndim,nelt)
      
      integer e

      do ie=1,nelt
         do k=1,nz
         do j=1,ny
            u( 1,j,k,ie)=u( 1,j,k,ie)*wt(j,k,1,1,ie)
            u(nx,j,k,ie)=u(nx,j,k,ie)*wt(j,k,2,1,ie)
         enddo
         enddo
         do k=1,nz
         do i=2,nx-1
            u(i, 1,k,ie)=u(i, 1,k,ie)*wt(i,k,1,2,ie)
            u(i,ny,k,ie)=u(i,ny,k,ie)*wt(i,k,2,2,ie)
         enddo
         enddo
         do j=2,ny-1
         do i=2,nx-1
            u(i,j, 1,ie)=u(i,j, 1,ie)*wt(i,j,1,3,ie)
            u(i,j,nz,ie)=u(i,j,nz,ie)*wt(i,j,2,3,ie)
         enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_tnsr1(v,nv,nu,A,At)
c
c     v = [A (x) A (x) A] u 
c
      integer nv,nu
      real v(1),A(1),At(1)

      call h1mg_tnsr1_3d(v,nv,nu,A,At,At)

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_tnsr1_3d(v,nv,nu,A,Bt,Ct) ! v = [C (x) B (x) A] u
      integer nv,nu
      real v(1),A(1),Bt(1),Ct(1)
      include 'SIZE'
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      integer e,e0,ee,es

      e0=1
      es=1
      ee=nelt

      if (nv.gt.nu) then
         e0=nelt
         es=-1
         ee=1
      endif

      nu3 = nu**3
      nv3 = nv**3

      do e=e0,ee,es
         iu = 1 + (e-1)*nu3
         iv = 1 + (e-1)*nv3
         call mxm(A,nv,v(iu),nu,work,nu*nu)
         do i=0,nu-1
            call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
         enddo
         call mxm(work2,nv*nv,Ct,nu,v(iv),nv)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_intp(uf,uc,l) ! l is coarse level
      real uf(1),uc(1)
      integer l
      include 'SIZE'
      include 'HSMG'
      call h1mg_tnsr(uf,mg_nh(l+1),uc,mg_nh(l),mg_jh(1,l),mg_jht(1,l))
      return
      end
c-----------------------------------------------------------------------
c     computes
c     v = [A (x) A] u      or
c     v = [A (x) A (x) A] u 
      subroutine h1mg_tnsr(v,nv,u,nu,A,At)
      integer nv,nu
      real v(1),u(1),A(1),At(1)

      call h1mg_tnsr3d(v,nv,u,nu,A,At,At)

      return
      end
c-------------------------------------------------------T--------------
c     computes
c              
c     v = [C (x) B (x) A] u
      subroutine h1mg_tnsr3d(v,nv,u,nu,A,Bt,Ct)
      integer nv,nu
      real v(nv*nv*nv,nelt),u(nu*nu*nu,nelt),A(1),Bt(1),Ct(1)
      include 'SIZE'
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      integer ie, i
      do ie=1,nelt
         call mxm(A,nv,u(1,ie),nu,work,nu*nu)
         do i=0,nu-1
            call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
         enddo
         call mxm(work2,nv*nv,Ct,nu,v(1,ie),nv)
      enddo
      return
      end
c------------------------------------------   T  -----------------------
c----------------------------------------------------------------------
      subroutine outmat2(a,m,n,k,name)
      include 'SIZE'
      real a(m,n)
      character*4 name

      n2 = min(n,8)
      write(6,2) nid,name,m,n,k
      do i=1,m
         write(6,1) nid,name,(a(i,j),j=1,n2)
      enddo
c   1 format(i3,1x,a4,16f6.2)
    1 format(i3,1x,a4,1p8e14.5)
    2 format(/,'Matrix: ',i3,1x,a4,3i8)
      return
      end
c-----------------------------------------------------------------------

