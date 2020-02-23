#!/usr/bin/env python

size_template = \
"""C     Dimension file to be included

      parameter (ldim=3)                      ! dimension
      parameter (lx1={lx1},ly1=lx1,lz1=lx1)      ! polynomial order

      parameter (lp = {lp})                     ! max number of processors
      parameter (lelt={lelt})                    ! max number of elements, per proc

      parameter (lelg=lelt*lp)                ! max total elements in a test

      parameter (ldimt=1,ldimt1=ldimt+1)      ! used in 'include' files


      common /dimn/ nelx,nely,nelz,nelt       ! local element common block
     $            , nx1,ny1,nz1,ndim,nfield,nid
"""

rea_template = \
""".false.    = ifbrick                  ! brick or linear geometry
{nelt} {nelt} 1 = iel0,ielN(per proc),stride ! range of number of elements per proc.
{nx1} {nx1}  1 = nx0,nxN,stride           ! poly. order range for nx1
0    0  0 = npx,npy,npz              ! np distrb, if != np, nekbone handle
0    0  0 = mx,my,mz                 ! nelt distrb, if != nelt, nekbone handle
"""

def nb_make_input(lx1, lelt):

    with open('SIZE', 'w') as f:
        f.write(size_template.format(lx1=lx1, lelt=lelt))

    with open('data.rea', 'w') as f:
        f.write(rea_template.format(lx1=lx1, lelt=lelt))

if __name__ == '__main__':

    import os
    from subprocess import check_call, STDOUT
    from datetime import datetime

    # =============================================================================================
    # SETUP PHASES
    # =============================================================================================

    # Make directory for logfiles
    logdir = 'logs.{0}'.format(datetime.now().isoformat())
    if not os.path.exists(logdir):
        os.makedirs(logdir)

    # Make the SIZE file
    with open('SIZE', 'w') as f:
        f.write(size_template.format(lx1=16, lelt=2048, lp=16))

    # Print the header
    print('mode,\tnx1,\tnelt,\ttime,\tMflops')

    # =============================================================================================
    # TEST PHASES
    # =============================================================================================

    for mode in ['serial', 'mpi', 'acc', 'cuda']:

        # -----------------------------------------------------------------------------------------
        # Compilation phase
        # -----------------------------------------------------------------------------------------

        if mode == 'serial':
            makenek = './makenek.mpi'
        else:
            makenek = './makenek.{0}'.format(mode)

        makeOutFile = os.path.join(logdir, 'make.{0}.output'.format(mode))

        with open(makeOutFile, 'w') as makeOut:
            check_call([makenek, 'data', '-nocompile'], stdout=makeOut, stderr=STDOUT)
            check_call(['make', '-f', 'makefile', 'clean'], stdout=makeOut, stderr=STDOUT)
            check_call(['make', '-f', 'makefile'], stdout=makeOut, stderr=STDOUT)

        # -----------------------------------------------------------------------------------------
        # Run phases
        # -----------------------------------------------------------------------------------------

        for nx1 in [8, 16]:
            for nelt in [2**n for n in range(12)]:

                # Construct the .rea file for these nx1, nelt, np settings
                if mode == 'mpi':
                    numProcs = 16
                else:
                    numProcs = 1

                if numProcs > nelt:
                    print('{0},\t{1},\t{2},\t{3},\t{4}'.format(mode, nx1, nelt, '', ''))
                else:

                    with open('data.rea', 'w') as f:
                        f.write(rea_template.format(nx1=nx1, nelt=nelt/numProcs))

                    times = []
                    flops = []
                    numReps = 3

                    for rep in range(1, numReps+1):

                        # Name the logfiles
                        runOutFile = os.path.join(logdir, 
                                'run.{0}.nx1_{1}.nelt_{2}.rep_{3}.output'.format(mode, nx1, nelt, rep))
                        runErrFile = os.path.join(logdir,
                                'run.{0}.nx1_{1}.nelt_{2}.rep_{3}.error'.format(mode, nx1, nelt, rep))

                        # DO THE RUN!
                        with open(runOutFile, 'w') as runOut, open(runErrFile, 'w') as runErr:
                            check_call(['mpiexec', '-np', str(numProcs), '-bind-to', 'numa', './nekbone', 'data'], 
                                    stdout=runOut, stderr=runErr)

                        # Grep the output
                        with open(runOutFile, 'r') as runOut:
                            for l in runOut:
                                if 'Solve Time' in l:
                                    times.append(float(l.strip().split()[-1]))
                                elif 'Tot MFlops' in l:
                                    flops.append(float(l.strip().split()[3].strip(',')))

                    # Print results to stdout
                    print('{0},\t{1},\t{2},\t{3:0.4E},\t{4:0.4E}'.format(mode, nx1, nelt, 
                        sum(times)/numReps, sum(flops)/numReps))

