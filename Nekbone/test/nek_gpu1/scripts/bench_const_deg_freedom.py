#!/usr/bin/env python3

size_template = \
"""C     Dimension file to be included

      parameter (ldim=3)                      ! dimension
      parameter (lx1={lx1},ly1=lx1,lz1=lx1)      ! polynomial order

      parameter (lp = 10)                     ! max number of processors
      parameter (lelt={lelt})                    ! max number of elements, per proc

      parameter (lelg=lelt*lp)                ! max total elements in a test

      parameter (ldimt=1,ldimt1=ldimt+1)      ! used in 'include' files


      common /dimn/ nelx,nely,nelz,nelt       ! local element common block
     $            , nx1,ny1,nz1,ndim,nfield,nid
"""

rea_template = \
""".true.    = ifbrick                  ! brick or linear geometry
{lelt} {lelt} 1 = iel0,ielN(per proc),stride ! range of number of elements per proc.
{lx1} {lx1}  1 = nx0,nxN,stride           ! poly. order range for nx1
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

    logdir = 'logs'
    os.makedirs(logdir, exist_ok=True)

    print('mode,\tnx1,\tnelt,\ttime,')

    for mode in ['mpi', 'acc', 'cuda']:

        for lx1 in range(2, 17, 2):

            lelt = 256 * 16**3 // lx1**3

            nb_make_input(lx1=lx1, lelt=lelt)

            makenek = './makenek.{0}'.format(mode)
            makeOutFile = os.path.join(logdir,
                'make.{0}.nx1_{1}.nelt_{2}.output'.format(mode, lx1, lelt))

            with open(makeOutFile, 'w') as makeOut:
                check_call([makenek, 'data', '-nocompile'], stdout=makeOut, stderr=STDOUT)
                check_call(['make', '-f', 'makefile', 'clean'], stdout=makeOut, stderr=STDOUT)
                check_call(['make', '-f', 'makefile'], stdout=makeOut, stderr=STDOUT)

            times = []
            numReps = 3

            for rep in range(1, numReps+1):

                runOutFile = os.path.join(logdir, 
                    'run.{0}.nx1_{1}.nelt_{2}.rep_{3}.output'.format(mode, lx1, lelt, rep))
                runErrFile = os.path.join(logdir,
                    'run.{0}.nx1_{1}.nelt_{2}.rep_{3}.error'.format(mode, lx1, lelt, rep))

                with open(runOutFile, 'w') as runOut, open(runErrFile, 'w') as runErr:
                    check_call(['mpiexec', '-n', '1', './nekbone', 'data'], stdout=runOut, stderr=runErr)

                with open(runOutFile, 'r') as runOut:
                    for l in runOut:
                        if 'Solve Time' in l:
                            times.append(float(l.strip().split()[-1]))

            print('{0},\t{1},\t{2},\t{3},'.format(mode, lx1, lelt, sum(times)/numReps))

