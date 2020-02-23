c-----------------------------------------------------------------------
      program nekbone
      
      include 'SIZE'
      include 'TOTAL'
      include 'mpif.h'

      parameter (lxyz = lx1*ly1*lz1)
      parameter (lt=lxyz*lelt)

      real ah(lx1*lx1),bh(lx1),ch(lx1*lx1),dh(lx1*lx1)
     $    ,zpts(2*lx1),wght(2*lx1)
      
      real x(lt),f(lt),r(lt),w(lt),p(lt),z(lt),c(lt)
      real g(6,lt)
      real mfloplist(1024), avmflop
      real tstart, tstop
      integer icount

      logical ifbrick
      integer iel0,ielN,ielD   ! element range per proc.
      integer nx0,nxN,nxD      ! poly. order range
      integer npx,npy,npz      ! poly. order range
      integer mx,my,mz         ! poly. order range
      integer numthreads, omp_get_max_threads

      call iniproc(mpi_comm_world)    ! has nekmpi common block
      tstart = dnekclock()
      call read_param(ifbrick,iel0,ielN,ielD,nx0,nxN,nxD,
     +                npx,npy,npz,mx,my,mz)

      numthreads = 1
#ifdef _OPENMP
      numthreads= omp_get_max_threads()
#endif 

      if (nid.eq.0) then
        write(*,*) "Max number of threads: ", numthreads
      end if

c     GET PLATFORM CHARACTERISTICS
c     iverbose = 1
c     call platform_timer(iverbose)   ! iverbose=0 or 1

      icount = 0

#ifndef NITER 
#define NITER 100
#endif
      niter = NITER

      if (nid.eq.0) then
        write(*,*) "Number of iterations: ", niter
      end if

#ifdef LOG
#define WLOG(X) if (nid .eq. 0) write(*,*) X 
#else 
#define WLOG(X) 
#endif

c     SET UP and RUN NEKBONE
      do nx1=nx0,nxN,nxD
         WLOG("calling init_dim")
         call init_dim
         do nelt=iel0,ielN,ielD
           WLOG("calling init_mesh")
           call init_mesh(ifbrick,npx,npy,npz,mx,my,mz)
           WLOG("calling prox_setupds")
           call proxy_setupds    (gsh)     ! Has nekmpi common block
           WLOG("calling set_multiplicity")
           call set_multiplicity (c)       ! Inverse of counting matrix

           WLOG("calling proxy_setup")
           call proxy_setup(ah,bh,ch,dh,zpts,wght,g) 

           n     = nx1*ny1*nz1*nelt

           WLOG("calling set_f")
           call set_f(f,c,n)
           WLOG("calling cg")
           call cg(x,f,g,c,r,w,p,z,n,niter,flop_cg)

           WLOG("calling nekgsync")
           call nekgsync()

           WLOG("calling set_timer_flop_count")
           call set_timer_flop_cnt(0)
           WLOG("calling cg")
           call cg(x,f,g,c,r,w,p,z,n,niter,flop_cg)
           WLOG("calling set_timer_flop_count")
           call set_timer_flop_cnt(1)

           WLOG("calling gs_free")
           call gs_free(gsh)

           icount = icount + 1
           mfloplist(icount)= mflops*np
         enddo
      enddo

      avmflop = 0.0
      do i = 1, icount
        avmflop = avmflop + mfloplist(i)
      end do

      if (icount .ne. 0) then
        avmflop = avmflop/icount
      end if

      if (nid .eq. 0) then
        write(6,1) avmflop
      end if
    1 format('Av MFlops = ', 1pe12.4)

c     TEST BANDWIDTH BISECTION CAPACITY
c     call xfer(np,cr_h)

      call nekgsync()
      tstop = dnekclock()
      if (nid .eq.0) write(*,*) "Total run time = ", tstop-tstart

      call exitt0

      end
c--------------------------------------------------------------
      subroutine set_f(f,c,n)
      real f(n),c(n)
      integer i
      integer, allocatable :: seed(:)

      call RANDOM_SEED(SIZE=i)
      allocate(seed(i))
      seed = 5
      call RANDOM_SEED(PUT=seed(1:i))

      do i=1,n
         call RANDOM_NUMBER(f(i))
      enddo

      call dssum(f)
      call col2 (f,c,n)

      deallocate(seed)

      return
      end
c-----------------------------------------------------------------------
      subroutine init_dim

C     Transfer array dimensions to common

      include 'SIZE'
      include 'INPUT'
 
      ny1=nx1
      nz1=nx1
 
      ndim=ldim

      return
      end
c-----------------------------------------------------------------------
      subroutine init_mesh(ifbrick,npx,npy,npz,mx,my,mz)
      include 'SIZE'
      include 'TOTAL'
      logical ifbrick
      integer e,eg,offs
 

      if(.not.ifbrick) then   ! A 1-D array of elements of length P*lelt
  10     continue
         nelx = nelt*np
         nely = 1
         nelz = 1
   
         do e=1,nelt
            eg = e + nid*nelt
            lglel(e) = eg
         enddo
      else              ! A 3-D block of elements 
        if (npx*npy*npz .ne. np) then
          call cubic(npx,npy,npz,np)  !xyz distribution of total proc
        end if 
        if (mx*my*mz .ne. nelt) then
          call cubic(mx,my,mz,nelt)   !xyz distribution of elements per proc
        end if 
      
c       if(mx.eq.nelt) goto 10

        nelx = mx*npx
        nely = my*npy 
        nelz = mz*npz

        e = 1
        offs = (mod(nid,npx)*mx) + npx*(my*mx)*(mod(nid/npx,npy)) 
     $      + (npx*npy)*(mx*my*mz)*(nid/(npx*npy))
        do k = 0,mz-1
        do j = 0,my-1
        do i = 0,mx-1
           eg = offs+i+(j*nelx)+(k*nelx*nely)+1
           lglel(e) = eg
           e        = e+1
        enddo
        enddo
        enddo
      endif

      if (nid.eq.0) then
        write(6,*) "Processes: npx= ", npx, " npy= ", npy, " npz= ", npz
        write(6,*) "Local Elements: mx= ", mx, " my= ", my, " mz= ", mz
        write(6,*) "Elements: nelx= ", nelx, " nely= ", nely,
     &             " nelz= ", nelz
      end if

      return
      end
c-----------------------------------------------------------------------
      subroutine cubic(mx,my,mz,np)

      mx = np
      my = 1
      mz = 1
      ratio = np

      iroot3 = np**(1./3.) + 0.000001
      do i= iroot3,1,-1
        iz = i
        myx = np/iz
        nrem = np-myx*iz

        if (nrem.eq.0) then
          iroot2 = myx**(1./2.) + 0.000001
          do j=iroot2,1,-1
            iy = j
            ix = myx/iy
            nrem = myx-ix*iy
            if (nrem.eq.0) goto 20
          enddo
   20     continue

          if (ix < iy) then
            it = ix
            ix = iy
            iy = it
          end if      

          if (ix < iz) then
            it = ix
            ix = iz
            iz = it
          end if      

          if (iy < iz) then
            it = iy
            iy = iz
            iz = it
          end if      

          if ( REAL(ix)/iz < ratio) then
            ratio = REAL(ix)/iz
            mx = ix
            my = iy
            mz = iz
          end if 

        end if
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine set_multiplicity (c)       ! Inverse of counting matrix
      include 'SIZE'
      include 'TOTAL'

      real c(1)

      n = nx1*ny1*nz1*nelt

      call rone(c,n)
      call gs_op(gsh,c,1,1,0)  ! Gather-scatter operation  ! w   = QQ  w

      do i=1,n
         c(i) = 1./c(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_timer_flop_cnt(iset)
      include 'SIZE'
      include 'TOTAL'
      include 'TIMER'

      integer i, numThrd, totThd
      integer omp_get_max_threads
      real tmp1(20), tmp2(20), tmp3(20), tmp4(20)

      real time0,time1
      save time0,time1

      if (iset.eq.0) then
         flop_a  = 0
         flop_cg = 0

         do i = 1, tmax
           trzero(i) = 0
           tcopy(i) = 0
           tsolvem(i) = 0
           tglsc3a(i) = 0
           tglsc3b(i) = 0
           tglsc3c(i) = 0
           tglsc3d(i) = 0
           tadd2s1(i) = 0
           tadd2s2a(i) = 0
           tadd2s2b(i) = 0
           tadd2s2c(i) = 0
           tlocalgrad3(i) = 0
           twrwswt(i) = 0
           tlocalgrad3t(i) = 0
           tgsop(i) = 0
           tgop(1,i) = 0
           tgop(2,i) = 0
           tgop(3,i) = 0
           tgop(4,i) = 0
         end do

         time0   = dnekclock()
      else
        time1   = dnekclock()-time0
        if (time1.gt.0) mflops = (flop_a+flop_cg)/(1.e6*time1)

        if (nid.eq.0) then
          write(6,1) nelt,np,nx1, nelt*np
          write(6,2) mflops*np, mflops
          write(6,3) flop_a,flop_cg,time1
        end if

    1   format('nelt = ', i7, ', np = ', i9, ', nx1 = ', i7,
     &         ', elements =', i10 )
    2   format('Tot MFlops = ', 1pe12.5, ', MFlops = ', e12.5)
    3   format('Ax FOp = ', 1pe12.5, ', CG FOp = ', e12.5,
     &         ', Solve Time = ', e12.5)

#ifdef TIMERS
        numThrd = 1
#ifdef _OPENMP
        numThrd = omp_get_max_threads()
#endif
        totThd = numThrd*np

        do i = 1, numThrd
          tglsc3a(i) = tglsc3a(i) - tgop(1,i)
          tglsc3b(i) = tglsc3b(i) - tgop(2,i)
          tglsc3c(i) = tglsc3c(i) - tgop(3,i)
          tglsc3d(i) = tglsc3d(i) - tgop(4,i)
        end do

        do i = 1,20
          tmp1(i) = 0.0
        end do
        
        tmp1(1)= time1
        do i = 1, numThrd
          tmp1(2)= tmp1(2) + trzero(i)
          tmp1(3)= tmp1(3) + tcopy(i)
          tmp1(4)= tmp1(4) + tsolvem(i)
          tmp1(5)= tmp1(5) + tglsc3a(i)
          tmp1(6)= tmp1(6) + tglsc3b(i)
          tmp1(7)= tmp1(7) + tglsc3c(i)
          tmp1(8)= tmp1(8) + tglsc3d(i)
          tmp1(9)= tmp1(9) + tadd2s1(i)
          tmp1(10)= tmp1(10) + tadd2s2a(i)
          tmp1(11)= tmp1(11) + tadd2s2b(i)
          tmp1(12)= tmp1(12) + tadd2s2c(i)
          tmp1(13)= tmp1(13) + tlocalgrad3(i)
          tmp1(14)= tmp1(14) + twrwswt(i)
          tmp1(15)= tmp1(15) + tlocalgrad3t(i)
          tmp1(16)= tmp1(16) + tgsop(i)
          tmp1(17)= tmp1(17) + tgop(1,i)
          tmp1(18)= tmp1(18) + tgop(2,i)
          tmp1(19)= tmp1(19) + tgop(3,i)
          tmp1(20)= tmp1(20) + tgop(4,i)
        end do

        call gop(tmp1, tmp4, '+  ', 20)

        tmp2(1)= time1
        tmp2(2)= trzero(1)
        tmp2(3)= tcopy(1)
        tmp2(4)= tsolvem(1)
        tmp2(5)= tglsc3a(1)
        tmp2(6)= tglsc3b(1)
        tmp2(7)= tglsc3c(1)
        tmp2(8)= tglsc3d(1)
        tmp2(9)= tadd2s1(1)
        tmp2(10)= tadd2s2a(1)
        tmp2(11)= tadd2s2b(1)
        tmp2(12)= tadd2s2c(1)
        tmp2(13)= tlocalgrad3(1)
        tmp2(14)= twrwswt(1)
        tmp2(15)= tlocalgrad3t(1)
        tmp2(16)= tgsop(1)
        tmp2(17)= tgop(1,1)
        tmp2(18)= tgop(2,1)
        tmp2(19)= tgop(3,1)
        tmp2(20)= tgop(4,1)

        do i = 2, numThrd
          if (trzero(i) < tmp2(2)) tmp2(2)= trzero(i)
          if (tcopy(i) < tmp2(3)) tmp2(3)= tcopy(i)
          if (tsolvem(i) < tmp2(4)) tmp2(4)= tsolvem(i)
          if (tglsc3a(i) < tmp2(5)) tmp2(5)= tglsc3a(i)
          if (tglsc3b(i) < tmp2(6)) tmp2(6)= tglsc3b(i)
          if (tglsc3c(i) < tmp2(7)) tmp2(7)= tglsc3c(i)
          if (tglsc3d(i) < tmp2(8)) tmp2(8)= tglsc3d(i)
          if (tadd2s1(i) < tmp2(9)) tmp2(9)= tadd2s1(i)
          if (tadd2s2a(i) < tmp2(10)) tmp2(10)= tadd2s2a(i)
          if (tadd2s2b(i) < tmp2(11)) tmp2(11)= tadd2s2b(i)
          if (tadd2s2c(i) < tmp2(12)) tmp2(12)= tadd2s2c(i)
          if (tlocalgrad3(i) < tmp2(13)) tmp2(13)= tlocalgrad3(i)
          if (twrwswt(i) < tmp2(14)) tmp2(14)= twrwswt(i)
          if (tlocalgrad3t(i) < tmp2(15)) tmp2(15)= tlocalgrad3t(i)
          if (tgsop(i) < tmp2(16)) tmp2(16)= tgsop(i)
          if (tgop(1,i) < tmp2(17)) tmp2(17)= tgop(1,i)
          if (tgop(2,i) < tmp2(18)) tmp2(18)= tgop(2,i)
          if (tgop(3,i) < tmp2(19)) tmp2(19)= tgop(3,i)
          if (tgop(4,i) < tmp2(20)) tmp2(20)= tgop(4,i)
        end do

        call gop(tmp2, tmp4, 'm  ', 20)

        tmp3(1)= time1
        tmp3(2)= trzero(1)
        tmp3(3)= tcopy(1)
        tmp3(4)= tsolvem(1)
        tmp3(5)= tglsc3a(1)
        tmp3(6)= tglsc3b(1)
        tmp3(7)= tglsc3c(1)
        tmp3(8)= tglsc3d(1)
        tmp3(9)= tadd2s1(1)
        tmp3(10)= tadd2s2a(1)
        tmp3(11)= tadd2s2b(1)
        tmp3(12)= tadd2s2c(1)
        tmp3(13)= tlocalgrad3(1)
        tmp3(14)= twrwswt(1)
        tmp3(15)= tlocalgrad3t(1)
        tmp3(16)= tgsop(1)
        tmp3(17)= tgop(1,1)
        tmp3(18)= tgop(2,1)
        tmp3(19)= tgop(3,1)
        tmp3(20)= tgop(4,1)

        do i = 2, numThrd
          if (trzero(i) > tmp3(2)) tmp3(2)= trzero(i)
          if (tcopy(i) > tmp3(3)) tmp3(3)= tcopy(i)
          if (tsolvem(i) > tmp3(4)) tmp3(4)= tsolvem(i)
          if (tglsc3a(i) > tmp3(5)) tmp3(5)= tglsc3a(i)
          if (tglsc3b(i) > tmp3(6)) tmp3(6)= tglsc3b(i)
          if (tglsc3c(i) > tmp3(7)) tmp3(7)= tglsc3c(i)
          if (tglsc3d(i) > tmp3(8)) tmp3(8)= tglsc3d(i)
          if (tadd2s1(i) > tmp3(9)) tmp3(9)= tadd2s1(i)
          if (tadd2s2a(i) > tmp3(10)) tmp3(10)= tadd2s2a(i)
          if (tadd2s2b(i) > tmp3(11)) tmp3(11)= tadd2s2b(i)
          if (tadd2s2c(i) > tmp3(12)) tmp3(12)= tadd2s2c(i)
          if (tlocalgrad3(i) > tmp3(13)) tmp3(13)= tlocalgrad3(i)
          if (twrwswt(i) > tmp3(14)) tmp3(14)= twrwswt(i)
          if (tlocalgrad3t(i) > tmp3(15)) tmp3(15)= tlocalgrad3t(i)
          if (tgsop(i) > tmp3(16)) tmp3(16)= tgsop(i)
          if (tgop(1,i) > tmp3(17)) tmp3(17)= tgop(1,i)
          if (tgop(2,i) > tmp3(18)) tmp3(18)= tgop(2,i)
          if (tgop(3,i) > tmp3(19)) tmp3(19)= tgop(3,i)
          if (tgop(4,i) > tmp3(20)) tmp3(20)= tgop(4,i)
        end do

        call gop(tmp3, tmp4, 'M  ', 20)

        if (nid.eq.0) then
          write(6,4) "time       = ",tmp1(1)/np, tmp2(1), tmp3(1)
          write(6,4) "rzero      = ",tmp1(2)/totThd, tmp2(2), tmp3(2)
          write(6,4) "copy       = ",tmp1(3)/totThd, tmp2(3), tmp3(3)
          write(6,4) "glsc3a     = ",tmp1(5)/totThd, tmp2(5), tmp3(5)
          write(6,4) "gopa       = ",tmp1(17)/totThd, tmp2(17), tmp3(17)
          write(6,4) "solveM     = ",tmp1(4)/totThd, tmp2(4), tmp3(4)
          write(6,4) "glsc3b     = ",tmp1(6)/totThd, tmp2(6), tmp3(6)
          write(6,4) "gopb       = ",tmp1(18)/totThd, tmp2(18), tmp3(18)
          write(6,4) "add2s1     = ",tmp1(9)/totThd, tmp2(9), tmp3(9)
          write(6,4) "localgrad3 = ",tmp1(13)/totThd, tmp2(13), tmp3(13)
          write(6,4) "wrwswt     = ",tmp1(14)/totThd, tmp2(14), tmp3(14)
          write(6,4) "localgradt = ",tmp1(15)/totThd, tmp2(15), tmp3(15)
          write(6,4) "gsop       = ",tmp1(16)/totThd, tmp2(16), tmp3(16)
          write(6,4) "add2s2a    = ",tmp1(10)/totThd, tmp2(10), tmp3(10)
          write(6,4) "glsc3c     = ",tmp1(7)/totThd, tmp2(7), tmp3(7)
          write(6,4) "gopc       = ",tmp1(19)/totThd, tmp2(19), tmp3(19)
          write(6,4) "add2s2b    = ",tmp1(11)/totThd, tmp2(11), tmp3(11)
          write(6,4) "add2s2c    = ",tmp1(12)/totThd, tmp2(12), tmp3(12)
          write(6,4) "glsc3d     = ",tmp1(8)/totThd, tmp2(8), tmp3(8)
          write(6,4) "gopd       = ",tmp1(20)/totThd, tmp2(20), tmp3(20)
        endif

    4   format(A, 1pe12.4, e12.4, e12.4)

c       if (nid.eq.0) then
c         write(6,4) "av time: ", tmp2(1)/np, tmp2(2)/totThd,
c    &               tmp2(3)/totThd, tmp2(4)/totThd, tmp2(5)/totThd
c         write(6,5) "av time: ", tmp2(5)/totThd, tmp2(6)/totThd,
c    &               tmp2(7)/totThd, tmp2(8)/totThd
c       endif

c       if (nid.eq.0) then
c         write(6,4) "min time: ", tmp2(1), tmp2(2), tmp2(3),
c    &               tmp2(4), tmp2(5)
c         write(6,5) "min time: ", tmp2(5), tmp2(6), tmp2(7),
c    &               tmp2(8)
c       endif

c       if (nid.eq.0) then
c         write(6,4) "max time: ", tmp2(1), tmp2(2), tmp2(3),
c    &               tmp2(4), tmp2(5)
c         write(6,5) "max time: ", tmp2(5), tmp2(6), tmp2(7),
c    &               tmp2(8)
c       endif

c   4   format(A, ' cg= ', 1pe12.4, ', zcm= ', e12.4,
c    &         ', glsc3= ', e12.4, ', add2sx= ', e12.4,
c    &         ', ax= ', e12.4)
c   5   format(A, ' ax= ', 1pe12.4, ', add2s2= ', e12.4,
c    &         ', gsop= ', e12.4, ', axe= ', e12.4)
#endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine xfer(np,gsh)
      include 'SIZE'
      parameter(npts_max = lx1*ly1*lz1*lelt)

      real buffer(2,npts_max)
      integer ikey(npts_max)


      nbuf = 800
      npts = 1
      do itest=1,200
         npoints = npts*np

         call load_points(buffer,nppp,npoints,npts,nbuf)
         iend   = mod1(npoints,nbuf)
         istart = 1
         if(nid.ne.0)istart = iend+(nid-1)*nbuf+1
         do i = 1,nppp
            icount=istart+(i-1)
            ikey(i)=mod(icount,np)
         enddo

         call nekgsync
         time0 = dnekclock()
         do loop=1,50
            call crystal_tuple_transfer(gsh,nppp,npts_max,
     $                ikey,1,ifake,0,buffer,2,1)
         enddo
         time1 = dnekclock()
         etime = (time1-time0)/50

         if (nid.eq.0) write(6,1) np,npts,npoints,etime
   1     format(2i7,i10,1p1e12.4,' bandwidth' )
         npts = 1.02*(npts+1)
         if (npts.gt.npts_max) goto 100
      enddo
 100  continue

      return
      end
c-----------------------------------------------------------------------
      subroutine load_points(buffer,nppp,npoints,npts,nbuf)
      include 'SIZE'
      include 'PARALLEL'

      real buffer(2,nbuf)

      nppp=0
      if(nbuf.gt.npts) then
       npass = 1+npoints/nbuf

       do ipass = 1,npass
          if(nid.eq.ipass.and.ipass.ne.npass) then
            do i = 1,nbuf
             buffer(1,i)=i
             buffer(2,i)=nid
            enddo
            nppp=nbuf
          elseif (npass.eq.ipass.and.nid.eq.0) then
            mbuf=mod1(npoints,nbuf)
            do i=1,mbuf
               buffer(1,i)=i
               buffer(2,i)=nid
            enddo
            nppp=mbuf
          endif
       enddo
      else
       do i = 1,npts
          buffer(1,i)=i
          buffer(2,i)=nid
       enddo
       nppp=npts
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine read_param(ifbrick,iel0,ielN,ielD,nx0,nxN,nxD,
     +                      npx,npy,npz,mx,my,mz)
      include 'SIZE'
      logical ifbrick
      integer iel0,ielN,ielD,nx0,nxN,nxD,npx,npy,npz,mx,my,mz

      !open .rea
      if(nid.eq.0) then
         open(unit=9,file='data.rea',status='old') 
         read(9,*,err=100) ifbrick
         read(9,*,err=100) iel0,ielN,ielD
         read(9,*,err=100) nx0,nxN,nxD
         read(9,*,err=100) npx,npy,npz
         read(9,*,err=100) mx,my,mz
         close(9)
      endif
      call bcast(ifbrick,4)
      call bcast(iel0,4)
      call bcast(ielN,4)
      call bcast(ielD,4)
c     nx0=lx1
c     nxN=lx1
      call bcast(nx0,4)
      call bcast(nxN,4)
      call bcast(nxD,4)
      call bcast(npx,4)
      call bcast(npy,4)
      call bcast(npz,4)
      call bcast(mx,4)
      call bcast(my,4)
      call bcast(mz,4)
      if(iel0.gt.ielN.or.nx0.gt.nxN) goto 200

      return

  100 continue
      write(6,*) "ERROR READING data.rea....ABORT"
      call exitt0

  200 continue
      write(6,*) "ERROR data.rea :: iel0 > ielN or nx0 > nxN :: ABORT"
      call exitt0
  
      return
      end
c-----------------------------------------------------------------------
