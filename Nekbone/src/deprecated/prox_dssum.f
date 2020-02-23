c-----------------------------------------------------------------------
      subroutine dssum(f)
      include 'SIZE'
      include 'TOTAL'
      real f(1)

c     call nekgsync()
      call adelay
      call gs_op(gsh,f,1,1,0)  ! Gather-scatter operation  ! w   = QQ  w

      return
      end
c-----------------------------------------------------------------------
      subroutine proxy_setupds(gs_handle,nx)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      integer gs_handle,dof
      integer*8 glo_num(lx1*ly1*lz1*lelt)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      t0 = dnekclock()

      call set_vert_box(glo_num,nx) ! Set global-to-local map
c     call outmat_glo_num(glo_num,nx)
    

      ntot      = nx*nx*nx*nelt   ! assumes nx=ny=nz
      call gs_setup(gs_handle,glo_num,ntot,nekcomm,mp) ! Initialize gather-scatter
      dof = ntot *mp
      t1 = dnekclock() - t0
c     if (nid.eq.0) then
c        write(6,1) t1,gs_handle,nx,dof
c   1    format('   setupds time',1pe11.4,' seconds ',2i3,i12)
c     endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_vert_box(glo_num,nx)

c     Set up global numbering for elements in a box

      include 'SIZE'
      include 'PARALLEL'

      integer*8 glo_num(1),ii,kg,jg,ig ! The latter 3 for proper promotion

      integer e,ex,ey,ez,eg

      nn = nx-1  ! nn := polynomial order

      do e=1,nelt
        eg = lglel(e)                              
        call get_exyz(ex,ey,ez,eg,nelx,nely,nelz)  
        do k=0,nn
        do j=0,nn
        do i=0,nn
           kg = nn*(ez-1) + k                     
           jg = nn*(ey-1) + j                     
           ig = nn*(ex-1) + i
           ii = 1 + ig + jg*(nn*nelx+1) + kg*(nn*nelx+1)*(nn*nely+1) 
           ll = 1 + i + nx*j + nx*nx*k + nx*nx*nx*(e-1)
           glo_num(ll) = ii
        enddo
        enddo
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_exyz(ex,ey,ez,eg,nelx,nely,nelz)
      integer ex,ey,ez,eg

      nelxy = nelx*nely
 
      ez = 1 +  (eg-1)/nelxy
      ey = mod1 (eg,nelxy)
      ey = 1 +  (ey-1)/nelx
      ex = mod1 (eg,nelx)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_glo_num(glo_num,nx)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      integer*8 glo_num(lx1*ly1*lz1,lelt)

      integer e

      do e=1,nelt
         call outmat_e_i8(glo_num(1,e),e,nx)
      enddo
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_e_i8(gn,e,nx)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      integer*8 gn(lx1,ly1,lz1)

      integer e

      write(6,*)
      write(6,2) e
      write(6,*)

      do k0=3,1,-2

         k1=k0+1
            write(6,*) k0,k1
         do j=nx,1,-1
            write(6,1) ((gn(i,j,k),i=1,nx),k=k0,k1)
         enddo
         write(6,*)

      enddo
    1 format('gn:',4i5,3x,4i5)
    2 format('gn: element: ',i4)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_glo_num_general(glo_num,nx)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      integer*8 glo_num(1)

      integer e

      io = 0
      do e=1,nelt
         write(6,*)
         write(6,2) e
         write(6,*)

         do k=nx,1,-1
              write(6,*) 'k = ',k
            do j=nx,1,-1
                write(6,1) (glo_num(igo),igo=1+io,io+nx)
                io = io+nx
            enddo
         enddo

      enddo
    1 format('gn:',6i5)
    2 format('gn: element: ',i4)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_r(x,name5)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      character*5 name5

      real x(lx1*ly1*lz1,lelt)

      integer e

      do e=1,nelt
         call outmat_e_r(x(1,e),name5,e)
      enddo
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_e_r(x,name5,e)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      character*5 name5

      real x(lx1,ly1,lz1)

      integer e

      write(6,*)
      write(6,2) e,name5
      write(6,*)

      do k0=3,1,-2

         k1=k0+1
         do j=ny1,1,-1
            write(6,1) ((x(i,j,k),i=1,4),k=k0,k1)
         enddo
         write(6,*)

      enddo
    1 format('mat: ',4f8.3,3x,4f8.3)
    2 format('mat: element: ',i4,2x,a5)
 
      return
      end
c-----------------------------------------------------------------------

#ifdef _OPENACC
c-----------------------------------------------------------------------
      subroutine i8sort(a,ind,n) 
c     Sort routine for a = int*8, ind=int.
c     Uses heap sort (p 231 Num. Rec., 1st Ed.)

      integer*8 a(1),aa
      integer ind(1)

      do j=1,n
         ind(j)=j
      enddo

      if (n.le.1) return
      L=n/2+1
      ir=n
 100     continue
         if (l.gt.1) then
            l=l-1
            aa  = a  (l)
            ii  = ind(l)
         else
                 aa =   a(ir)
                 ii = ind(ir)
              a(ir) =   a( 1)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
                 a(1) = aa
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
 200              continue
         if (j.le.ir) then
            if (j.lt.ir) then
               if ( a(j).lt.a(j+1) ) j=j+1
            endif
            if (aa.lt.a(j)) then
                 a(i) = a(j)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         goto 200
         endif
           a(i) = aa
         ind(i) = ii
      goto 100
      end

c-----------------------------------------------------------------------
      subroutine ldssum(u,ug)
   
      include 'SIZE'
      include 'INPUT'
      include 'ACCNEK'

      real u  (lx1*ly1*lz1*lelt)
      real ug(lx1*ly1*lz1*lelt)

      ndssum      = ids_lgl1(0)
      nglobl      = ids_lgl2(0)

      call rzero(ug,nglobl)

      do i=1,ndssum       ! local Q^T
         il=ids_lgl1(i)
         ig=ids_lgl2(i)
         ug(ig) = ug(ig)+u(il)
      enddo

      do i=1,ndssum       ! local Q
         il=ids_lgl1(i)
         ig=ids_lgl2(i)
         u(il) = ug(ig)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine proxy_setupds_acc(gs_handle)
      include 'SIZE'
      include 'INPUT'
      include 'ACCNEK'

      integer gs_handle,dof
      integer*8 glo_num(lx1*ly1*lz1*lelt)
      integer pcount

c      integer ids_lgl(2,-1:1)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      common /nekmpi_acc/ gs_handle2
      integer gs_handle2

c     Work arrays
      integer*8 ngv,wk(lx1*ly1*lz1*lelt)
      real      u  (lx1*ly1*lz1*lelt),v(lx1*ly1*lz1*lelt)
      integer   idx(lx1*ly1*lz1*lelt),vertex(1),n_nonlocal

      t0 = dnekclock()

      call set_vert_box(glo_num,nx1) ! Set global-to-local map

      n = nx1*ny1*nz1*nelt

c      call gs_setup(gs_handle,glo_num,n,nekcomm,mp) ! Initialize gather-scatter

      do ipass=1,2
         call i8copy(wk,glo_num,n)
         call i8sort(wk,idx,n)

         ig = 0
         ip = 0
         glo_num_last = 0

         do i=1,n
            if (wk(i).ne.0) then
              ip=ip+1
              il=idx(i)
              if (ig.eq.0) then
                 ig=ig+1
                 ids_lgl1(ip) = il
                 ids_lgl2(ip) = ig
              elseif (wk(i).eq.glo_num_last) then
                 ids_lgl1(ip) = il
                 ids_lgl2(ip) = ig
              else
                 ig = ig+1
                 ids_lgl1(ip) = il
                 ids_lgl2(ip) = ig
              endif
              glo_num_last  = wk(i)
              wk(i) = abs(wk(i))
            endif
         enddo
         
         ids_lgl1(0) = ip
         ids_lgl2(0) = ig
         ndssum      = ids_lgl1(0)
         nglobl      = ids_lgl2(0)

         if (ipass.eq.1) then ! eliminate singletons
            call gs_setup(gs_handle,glo_num,n,nekcomm,mp) ! Initialize gather-scatter
            call rone    (u,n)
            call gs_op   (gs_handle,u,1,1,0)  ! 1 ==> +
c            call gs_free (gs_handle)

            call rone    (v,n)
            call ldssum  (v,wk)

            n_nonlocal = 0
            do i=1,ndssum              ! Identify singletons and nonlocals
               il=ids_lgl1(i)
               ig=ids_lgl2(i)
               if (u(il).lt.1.1) then          ! Singleton
                  glo_num(il) = 0
               elseif (u(il).ne.v(il)) then ! Nonlocal
                  n_nonlocal  = n_nonlocal+1
                  glo_num(il) = -glo_num(il)
               endif
            enddo
         else  ! ipass = 2
            igl = 0
            nnl = 0
            do i=1,n_nonlocal   ! Identify nonlocals in ug()
               ig=ids_lgl2(i)
               if (ig.ne.igl) nnl=nnl+1
               wk(nnl) = wk(i)
               igl = ig
            enddo
            n_nonlocal=nnl
           
c            if ( n_nonlocal .eq. 0 ) then
c               write(*,*) "n_nonlocal=0"
c               return
c            enddo
            call gs_setup(gsh_acc,wk,n_nonlocal,nekcomm,mp)

         endif
      enddo

      ids_lgl1(-1) = gs_handle
      ids_lgl2(-1) = n_nonlocal

C      Used in GPU calucation
      pcount = 1
      ids_ptr(1) =1
      do i=2,ndssum
         if (ids_lgl2(i) .ne. ids_lgl2(i-1)) then
            pcount = pcount+1
            ids_ptr(pcount) = i
         endif
      enddo
      ids_ptr(pcount+1) = ndssum + 1
c
c     Check  
c      write(*,*) pcount, ids_lgl1(0), ids_lgl2(0), n_nonlocal
c      if ( pcount .ne. ids_lgl2(0)) stop

      dof = n *mp
      t1 = dnekclock() - t0
c     if (nid.eq.0) then
c        write(6,1) t1,gs_handle,nx,dof
c   1    format('   setupds time',1pe11.4,' seconds ',2i3,i12)
c     endif

      return
      end

c-----------------------------------------------------------------------
      subroutine dssum_acc(u)
      include 'SIZE'
      include 'ACCNEK'

      common /nekmpi/ mp

      parameter(lt=lx1*ly1*lz1*lelt)
      common /nsmpi_acc/ ug(lt)

      real u(lt),ug(lt)

      integer ndssum, nglobl

c     call nekgsync()
      call adelay

      ndssum=ids_lgl1(0)
      nglobl=ids_lgl2(0)
      gs_hnd      = ids_lgl1(-1)
      n_nonlocal  = ids_lgl2(-1)

!$ACC DATA PRESENT(ids_ptr(1:nglobl+1),u(1:lt))
!$ACC& PRESENT(ids_lgl1)
!$ACC& PRESENT(ug(1:nglobl))
!$ACC PARALLEL LOOP COLLAPSE(1)
      do i = 1,nglobl
         ug(i) = 0.0
         ! local Q^T
!$ACC LOOP SEQ
         do j = ids_ptr(i),ids_ptr(i+1)-1
            il = ids_lgl1(j)
            ug(i) = ug(i) + u(il)
         enddo
      enddo
!$ACC END PARALLEL LOOP

      if (n_nonlocal .gt. 1) then
c         ngv = nv + n_nonlocal  ! Number that must be copied out
!$ACC UPDATE HOST(ug(1:n_nonlocal))
!!$acc wait 
         call gs_op(gsh_acc,ug(1),1,1,0) ! 1 ==> +
!$ACC UPDATE DEVICE(ug(1:n_nonlocal)) 
!!$acc wait 
      endif
	
	abc

	stop

!$ACC PARALLEL LOOP COLLAPSE(1) 
      do i = 1,nglobl
        ! local Q
!$ACC LOOP SEQ
         do j = ids_ptr(i),ids_ptr(i+1)-1
            il = ids_lgl1(j)
            u(il) = ug(i)
        enddo
      enddo
!$ACC END PARALLEL LOOP
!$ACC END DATA

      return 
      end

#endif
