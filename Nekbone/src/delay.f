c-----------------------------------------------------------------------
      subroutine init_delay
      include 'SIZE'
      include 'mpif.h'
      include 'DELAY'
c   tmean:  20-100   0 -->
c   tavg:   10-50    0 -->  zero delay
c   drange: 40-200   0 -->
c   rms:    10-100   1 -->

      tmean   = 80
      tavg    = 40
      d_range = 80
      rms     = 30

      int_calls = 50

      if(nid.eq.0) then
         write(6,*)
         write(6,'(4e12.4)') tmean,tavg,d_range,rms
         write(6,*) " tmean, tavg, d_range, rms"
         write(6,*)
       endif



      return
      end
c-----------------------------------------------------------------------
      subroutine adelay()
      include 'SIZE'
      include 'DELAY'
      include 'TOTAL'

      parameter (imode1=0)
      parameter (imode2=0)
      parameter (imode3=0)

      save icall_sd,ncalls
      data icall_sd /0/
      data ncalls   /0/

      real*8  skew,wk,tot_time,tavg01
    
      integer*4 nms 
   
      if (ncalls.gt.int_calls) goto 102

      if (icall_sd.eq.0) then

       iseed=nid*11+1    
       call ZBQLINI(iseed)               

  102  continue
       ncalls=0  
       r01=ZBQLU01()
       t01=6.0*tavg*r01   

       f01=exp(-0.5*((t01-tmean)/rms)**2)
       fmax=1.01          
       r02=ZBQLU01()

c      Rejection method
       if ((r02*fmax).gt.f01) then
  101    continue 
         r01=ZBQLU01()
         t01=d_range*tavg*r01
         f02=exp(-0.5*((t01-tmean)/rms)**2)
         r02=ZBQLU01()
         if ((r02*fmax).gt.f02) goto 101
       endif    
       
       tot_time=t01
       t_delay =t01
       icall_sd=1
 
       call gop(tot_time,wk,'+  ',1)
       if(tot_time.eq.0.0) then
          t_delay = 0.0
       else
          tavg01 =tot_time/real(np) 
          t_delay=t01*tavg/tavg01       
       endif
c       write(6,*) "Test:",nid,t01,t_delay,tavg
      endif

      nms=int(t_delay)
c     write (6,*) "pTest ",nms,ncalls
      call sleepabit(nms)        
      ncalls=ncalls+1
c      write (6,*) "Test Done"    
      return 
      end




*******************************************************************
********    FILE: randgen.f             ***********
********    AUTHORS: Richard Chandler       ***********
********         (richard@stats.ucl.ac.uk)  ***********
********         Paul Northrop          ***********
********         (northrop@stats.ox.ac.uk)  ***********
********    LAST MODIFIED: 26/8/03          ***********
********    See file randgen.txt for details    ***********
*******************************************************************

      BLOCK DATA ZBQLBD01
*
*       Initializes seed array etc. for random number generator.
*       The values below have themselves been generated using the
*       NAG generator.
*
      COMMON /ZBQL0001/ ZBQLIX,B,C
      DOUBLE PRECISION ZBQLIX(43),B,C
      INTEGER I
      DATA (ZBQLIX(I),I=1,43) /8.001441D7,5.5321801D8,
     +1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8,
     +7.9952507D7,3.81202335D8,3.11575334D8,4.02878631D8,
     +2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8,
     +4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8,
     +2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8,
     +1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8,
     +3.20029657D8,3.26369301D8,9.441177D6,3.53244738D8,
     +2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8,
     +3.75423178D8,7.0893571D7 ,4.26059785D8,3.95854390D8,
     +2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8,
     +2.63576576D8/
      DATA B / 4.294967291D9 /
      DATA C / 0.0D0 /
      END
******************************************************************

******************************************************************
      SUBROUTINE ZBQLINI(SEED)
******************************************************************
      INTEGER LFLNO
      PARAMETER (LFLNO=80)
      INTEGER SEED,SS,MM,HH,DD,FILNO,I
      INTEGER INIT
      DOUBLE PRECISION ZBQLIX(43),B,C
      DOUBLE PRECISION TMPVAR1,DSS,DMM,DHH,DDD

      COMMON /ZBQL0001/ ZBQLIX,B,C
      SAVE INIT

*
*   Ensure we don't call this more than once in a program
*
      IF (INIT.GE.1) THEN
       IF(INIT.EQ.1) THEN
        WRITE(*,1)
        INIT = 2
       ENDIF
       RETURN
      ELSE
       INIT = 1
      ENDIF
*
*       If SEED = 0, cat the contents of the clock into a file
*       and transform to obtain ZQBLIX(1), then use a congr.
*       algorithm to set remaining elements. Otherwise take
*       specified value of SEED.
*
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>    NB FOR SYSTEMS WHICH DO NOT SUPPORT THE  >>>>>>>
*>>>>>>>    (NON-STANDARD) 'CALL SYSTEM' COMMAND,    >>>>>>>
*>>>>>>>    THIS WILL NOT WORK, AND THE FIRST CLAUSE >>>>>>>
*>>>>>>>    OF THE FOLLOWING IF BLOCK SHOULD BE  >>>>>>>
*>>>>>>>    COMMENTED OUT.               >>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (SEED.EQ.0) THEN
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>    COMMENT OUT FROM HERE IF YOU DON'T HAVE  >>>>>>>
*>>>>>>>    'CALL SYSTEM' CAPABILITY ...         >>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       temp=1.0
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<   ... TO HERE (END OF COMMENTING OUT FOR    <<<<<<<
*<<<<<<<<   USERS WITHOUT 'CALL SYSTEM' CAPABILITY    <<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ELSE
       TMPVAR1 = DMOD(DBLE(SEED),B)
      ENDIF
      ZBQLIX(1) = TMPVAR1
      DO 100 I = 2,43
       TMPVAR1 = ZBQLIX(I-1)*3.0269D4
       TMPVAR1 = DMOD(TMPVAR1,B)
       ZBQLIX(I) = TMPVAR1
 100  CONTINUE

 1    FORMAT(//5X,'****WARNING**** You have called routine ZBQLINI ',
     +'more than',/5X,'once. I''m ignoring any subsequent calls.',//)
 2    FORMAT(//5X,'**** ERROR **** In routine ZBQLINI, I couldn''t',
     +' find an',/5X,
     +'available file number. To rectify the problem, decrease the ',
     +'value of',/5X,
     +'the parameter LFLNO at the start of this routine (in file ',
     +'randgen.f)',/5X,
     +'and recompile. Any number less than 100 should work.')
      END
******************************************************************



******************************************************************
      FUNCTION ZBQLU01(DUMMY)
*
*       Returns a uniform random number between 0 & 1, using
*       a Marsaglia-Zaman type subtract-with-borrow generator.
*       Uses double precision, rather than integer, arithmetic
*       throughout because MZ's integer constants overflow
*       32-bit integer storage (which goes from -2^31 to 2^31).
*       Ideally, we would explicitly truncate all integer
*       quantities at each stage to ensure that the double
*       precision representations do not accumulate approximation
*       error; however, on some machines the use of DNINT to
*       accomplish this is *seriously* slow (run-time increased
*       by a factor of about 3). This double precision version
*       has been tested against an integer implementation that
*       uses long integers (non-standard and, again, slow) -
*       the output was identical up to the 16th decimal place
*       after 10^10 calls, so we're probably OK ...
*
      real*8 ZBQLU01,DUMMY,B,C,ZBQLIX(43),X,B2,BINV
      INTEGER CURPOS,ID22,ID43

      COMMON /ZBQL0001/ ZBQLIX,B,C
      SAVE /ZBQL0001/
      SAVE CURPOS,ID22,ID43
      DATA CURPOS,ID22,ID43 /1,22,43/

      B2 = B
      BINV = 1.0D0/B
 5    X = ZBQLIX(ID22) - ZBQLIX(ID43) - C
      IF (X.LT.0.0D0) THEN
       X = X + B
       C = 1.0D0
      ELSE
       C = 0.0D0
      ENDIF
      ZBQLIX(ID43) = X
*
*     Update array pointers. Do explicit check for bounds of each to
*     avoid expense of modular arithmetic. If one of them is 0 the others
*     won't be
*
      CURPOS = CURPOS - 1
      ID22 = ID22 - 1
      ID43 = ID43 - 1
      IF (CURPOS.EQ.0) THEN
       CURPOS=43
      ELSEIF (ID22.EQ.0) THEN
       ID22 = 43
      ELSEIF (ID43.EQ.0) THEN
       ID43 = 43
      ENDIF
*
*     The integer arithmetic there can yield X=0, which can cause
*     problems in subsequent routines (e.g. ZBQLEXP). The problem
*     is simply that X is discrete whereas U is supposed to
*     be continuous - hence if X is 0, go back and generate another
*     X and return X/B^2 (etc.), which will be uniform on (0,1/B).
*
      IF (X.LT.BINV) THEN
       B2 = B2*B
       GOTO 5
      ENDIF

      ZBQLU01 = X/B2

      END
******************************************************************

