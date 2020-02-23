import os, stat, sys, unittest
from subprocess import check_call, STDOUT

testRoot = os.path.dirname(os.path.abspath(__file__))

makenekTemplate = \
"""#!/bin/bash
# Nek5000 build config file
# (c) 2008,2009,2010 UCHICAGO ARGONNE, LLC

# source path
SOURCE_ROOT="../../src"

# Fortran compiler
F77="{f77}"

# C compiler
CC="{cc}"

# pre-processor symbol list
# (set PPLIST=? to get a list of available symbols)
# NEKCOMM, NEKDLAY, BG, MGRID
#PPLIST="?"


# OPTIONAL SETTINGS
# -----------------

# enable MPI (default true)
IFMPI="{ifmpi}"

# auxilliary files to compile
# NOTE: source files have to located in the same directory as makenek
#       a makefile_usr.inc has to be provided containing the build rules
#USR="foo.o"

# linking flags
USR_LFLAGS="{usr_lflags}"

# generic compiler flags
G="{g}"

# optimization flags
#OPT_FLAGS_STD=""
#OPT_FLAGS_MAG=""

###############################################################################
# DONT'T TOUCH WHAT FOLLOWS !!!
###############################################################################
# assign version tag
mver=1
# overwrite source path with optional 2nd argument
if [ -d $2 ] && [ $# -eq 2 ]; then
  SOURCE_ROOT="$2"
  echo "change source code directory to: ", $SOURCE_ROOT
fi
# do some checks and create makefile
source $SOURCE_ROOT/makenek.inc
# compile
make -f makefile 2>&1 | tee compiler.out
exit 0
"""


class NekboneTestCase():
    testDir = None

    @staticmethod
    def makeNekbone(mode, cwd=None, f77='pgfortran', cc='pgcc', ifmpi='false', usr_lflags='', g=''):

        makenek = 'makenek.{0}'.format(mode)
        makeOutFile = 'make.{0}.output'.format(mode)

        if cwd:
            makenek = os.path.join(cwd, makenek)
            makeOutFile = os.path.join(cwd, makeOutFile)

        with open(makenek, 'w') as f:
            f.write(makenekTemplate.format(f77=f77, cc=cc, ifmpi=ifmpi, usr_lflags=usr_lflags, g=g))
        os.chmod(makenek,
                 stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH |
                 stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH |
                 stat.S_IWUSR)

        with open(makeOutFile, 'w') as makeOut:
            check_call([makenek, 'data', '-nocompile'], cwd=cwd, stdout=makeOut)
            check_call(['make', '-f', 'makefile', 'clean'], cwd=cwd, stdout=makeOut)
            check_call(['make', '-f', 'makefile'], cwd=cwd, stdout=makeOut)

    @staticmethod
    def runNekboneSerial(mode, cwd=None, reaFile='data'):

        runOutFile = 'run.{0}.output'.format(mode)
        if cwd:
            runOutFile = os.path.join(cwd, runOutFile)

        with open(runOutFile, 'w') as runOut:
            check_call(['./nekbone', reaFile], cwd=cwd, stdout=runOut)

    @staticmethod
    def runNekboneParallel(mode, cwd=None, nProcs=2, reaFile='data'):

        runOutFile = 'run.{0}.output'.format(mode)
        if cwd:
            runOutFile = os.path.join(cwd, runOutFile)

        with open(runOutFile, 'w') as runOut:
            check_call(['mpirun', '-n', nProcs, './nekbone', reaFile], cwd=cwd, stdout=runOut)

    @classmethod
    def setUpClass(cls):

        cls.makeNekbone(mode='serial', cwd=cls.testDir, f77='pgfortran', cc='pgcc', ifmpi='false', usr_lflags='', g='')
        cls.runNekboneSerial(mode='serial', cwd=cls.testDir, reaFile='data')

    def assertOutputEqual(self, testMode, refMode='serial'):
        cwd = self.__class__.testDir

        testOutFile = os.path.join(cwd, 'run.{0}.output'.format(testMode))
        with open(testOutFile, 'r') as f:
            testOutList = [l.strip() for l in f if l.startswith('cg:')]

        refOutFile = os.path.join(cwd, 'run.{0}.output'.format(refMode))
        with open(refOutFile, 'r') as f:
            refOutList = [l.strip() for l in f if l.startswith('cg:')]

        for testLine, refLine in zip(testOutList, refOutList):
            errMsg = '\n'.join([
                'Results from {refMode} and {testMode} did not match',
                '{refMode}:', refLine,
                '{testMode}:', testLine,
            ])
            self.assertEqual(testLine, refLine, msg=errMsg)

    def test_AccDevice(self):
        testMode = 'acc.device'

        self.makeNekbone(mode=testMode, cwd=self.__class__.testDir,
                         f77='pgfortran', cc='pgcc', ifmpi='false', usr_lflags='-ta=nvidia:cc50',
                         g='-acc -Minfo=accel -ta=nvidia:cc50')

        self.runNekboneSerial(mode=testMode, cwd=self.__class__.testDir, reaFile='data')

        self.assertOutputEqual(testMode=testMode, refMode='serial')

    def test_AccHost(self):
        testMode = 'acc.host'

        self.makeNekbone(mode=testMode, cwd=self.__class__.testDir,
                         f77='pgfortran', cc='pgcc', ifmpi='false', usr_lflags='-ta=host',
                         g='-acc -Minfo=accel -ta=host')

        self.runNekboneSerial(mode=testMode, cwd=self.__class__.testDir, reaFile='data')

        self.assertOutputEqual(testMode=testMode, refMode='serial')

    def test_CudaDevice(self):
        testMode = 'cuda.device'

        self.makeNekbone(mode=testMode, cwd=self.__class__.testDir,
                         f77='pgfortran', cc='pgcc', ifmpi='false', usr_lflags='-Mcuda=cc50 -ta=nvidia:cc50',
                         g='-acc -Minfo=accel -Mcuda=cc50 -ta=nvidia:cc50')

        self.runNekboneSerial(mode=testMode, cwd=self.__class__.testDir, reaFile='data')

        self.assertOutputEqual(testMode=testMode, refMode='serial')


class NekboneGpu1(NekboneTestCase, unittest.TestCase):
    testDir = os.path.join(testRoot, 'nek_gpu1')


class NekboneExample1(NekboneTestCase, unittest.TestCase):
    testDir = os.path.join(testRoot, 'example1')


class NekboneExample2(NekboneTestCase, unittest.TestCase):
    testDir = os.path.join(testRoot, 'example2')


class NekboneExample3(NekboneTestCase, unittest.TestCase):
    testDir = os.path.join(testRoot, 'example3')


#class NekboneMgrid(NekboneTestCase):
#    testDir = os.path.join(testRoot, 'nek_mgrid')


if __name__ == '__main__':
    unittest.main()
