#ifdef TIMERS
#define NBTIMER(a) a = dnekclock()
#define STIMER(a) a = dnekclock_sync()
#define ACCUMTIMER(b,a) b = b + (dnekclock()- a)
#else
#define NBTIMER(a)
#define STIMER(a)
#define ACCUMTIMER(a,b)
#endif


      subroutine rzeroi(a,n,start,fin)
        implicit none
  
        real a(n)
        integer n, i, start, fin

        do i = start, fin
          a(i) = 0.0
        end do 

        return
      end subroutine

c----------------------------------------------------------

      subroutine copyi(a,b,n, start, fin)
        implicit none

        real a(n),b(n)
        integer n, i, start, fin

        do i=start,fin
          a(i)=b(i)
        enddo

        return
      end subroutine

c----------------------------------------------------------

      subroutine glsc3i(val,a,b,mult,n,find,lind)
      implicit none

      include 'TIMER'

      real val,a(n),b(n),mult(n)
      real tsum,psum,work(1)
      integer n,find,lind
      integer i, tmt, thread
      integer omp_get_thread_num
CCC 
	real tsuma(100000)
        save tsuma
    
      save psum
      data psum /0.0/

      thread = 0
#ifdef _OPENMP
CCC      thread = omp_get_thread_num()
#endif
      tmt = thread + 1

      tsum = 0.0
CCC
      if (lind .gt. 100000) then
        print *, "OOPS, target region temp too small"
        stop
      end if

CCCc$OMP target map(from:tsuma)map(to:a,b,mult,find,lind,n)private(i)
CCc$OMP target map(tsuma,a,b,mult,find,lind,n)private(i)
CCc$OMP teams distribute parallel do
c$OMP target teams distribute parallel do 
c$OMP+     map(from:tsuma)map(to:a,b,mult,find,lind) private(i)
CCC   c$OMP parallel 
      do i=find, lind
CCC         tsum = tsum + a(i)*b(i)*mult(i)
         tsuma(i) = a(i)*b(i)*mult(i)
      end do
c$OMP end target teams distribute parallel do 

CCC
	do i=find, lind
	  tsum = tsum + tsuma(i)
	end do
CCC c$OMP ATOMIC update
      psum = psum + tsum
CCC c$OMP END ATOMIC

CCC c$OMP BARRIER
      NBTIMER(ttemp4)
CCC c$OMP MASTER
      call gop(psum,work,'+  ',1)
      val = psum
      psum = 0.0
CCC c$OMP END MASTER
CCC c$OMP BARRIER
      ACCUMTIMER(tgop(gopi(tmt),tmt), ttemp4)


      return
      end subroutine

c----------------------------------------------------------

      subroutine solveMi(z,r,n,start,fin)
      implicit none

      real z(n),r(n)
      integer n,start,fin

      call copyi(z,r,n,start,fin) 

      return
      end

c----------------------------------------------------------

      subroutine add2s1i(a,b,c1,n,start,fin)
      implicit none

      real a(n),b(n),c1
      integer n,start,fin
      integer i

      do i= start, fin
        a(i)=c1*a(i)+b(i)
      end do

      return
      end subroutine

c----------------------------------------------------------

      subroutine add2s2i(a,b,c1,n,start,fin)
      implicit none
 
      real a(n),b(n),c1
      integer n,start,fin
      integer i

      do i= start,fin
        a(i)=a(i)+c1*b(i)
      end do

      return
      end subroutine

