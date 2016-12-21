# Default build suggestion of MPI + OPENMP with Clang on IBM (Power 8) + NVIDIA GPU machines.
# You might have to change the compiler name and flags.

SHELL = /bin/sh
.SUFFIXES: .cc .o

LULESH_EXEC = lulesh2.0

MPI_INC = /opt/local/include/openmpi
MPI_LIB = /opt/local/lib

# Point your mpicc to Clang
CXX = $(CLANG)

SOURCES2.0 = \
	lulesh.cc \
	lulesh-comm.cc \
	lulesh-viz.cc \
	lulesh-util.cc \
	lulesh-init.cc
OBJECTS2.0 = $(SOURCES2.0:.cc=.o)

teams = 
ifdef TEAMS
teams = -DTEAMS=$(TEAMS)
endif

threads = 
ifdef THREADS
threads = -DTHREADS=$(THREADS)
endif

gpu = 
ifdef USE_GPU
gpu = -DUSE_GPU=$(USE_GPU)
endif

mpi = -DUSE_MPI=0
ifdef USE_MPI
mpi = -DUSE_MPI=$(USE_MPI)
endif

# Tuning flags for Power 8
CXXFLAGS = -fopenmp -O3 -fopenmp-targets=nvptx64-nvidia-cuda $(shared) $(mpi) $(teams) $(threads) $(gpu) -I$(LIBOMP_LIB) -target x86_64-pc-linux-gnu

LDFLAGS = -L$(LIBOMP_LIB) -L$(CUDA)/lib64  -lcuda -lcudart -lelf -lffi -target x86_64-pc-linux-gnu

.cc.o: lulesh.h
	@echo "Building $<"
	$(CXX) -c $(CXXFLAGS) -o $@  $<

all: $(LULESH_EXEC)

lulesh2.0: $(OBJECTS2.0)
	@echo "Linking"
	$(CXX) -fopenmp=libomp $(CXXFLAGS) $(OBJECTS2.0) $(LDFLAGS) -lomp -lomptarget -lstdc++ -lm -o $@

clean:
	/bin/rm -f *.o *~ *.tgt* $(OBJECTS) $(LULESH_EXEC)
	/bin/rm -rf *.dSYM

tar: clean
	cd .. ; tar cvf lulesh-2.0.tar LULESH-2.0 ; mv lulesh-2.0.tar LULESH-2.0

