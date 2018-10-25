# Default build suggestion of MPI + OPENMP with Clang on IBM (Power 8) + NVIDIA GPU machines.
# You might have to change the compiler name and flags.

SHELL = /bin/sh
.SUFFIXES: .cc .o

HCC2 ?= /opt/rocm/hcc2
HCC2_GPU ?= gfx803
HCC2_GPUTARGET = amdgcn-amd-amdhsa
# Default is Radeon fiji/gfx803, for nvidia "export HCC2_GPU=sm_35"
ifeq (sm_,$(findstring sm_,$(HCC2_GPU)))
  HCC2_GPUTARGET = nvptx64-nvidia-cuda
endif
HCC2_CPUTARGET ?= x86_64-pc-linux-gnu
DEBUG_LEVEL ?= 0

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
CXXFLAGS = -O3 $(shared) $(mpi) $(teams) $(threads) $(gpu)
CXXFLAGS += -target $(HCC2_CPUTARGET) -fopenmp -fopenmp-targets=$(HCC2_GPUTARGET) -Xopenmp-target=$(HCC2_GPUTARGET) -march=$(HCC2_GPU)
CXX = $(HCC2)/bin/clang++
LDFLAGS = -fopenmp -fopenmp-targets=$(HCC2_GPUTARGET) -target $(HCC2_CPUTARGET) -lm
ifeq (nvptx,$(findstring nvptx,$(HCC2_GPUTARGET)))
   CUDA ?= /usr/local/cuda
   LDFLAGS += -L$(CUDA)/lib64 -lcuda -lcudart -lelf -lffi
endif

ifeq ($(DEBUG_LEVEL),0)
  SETDEBUGLIB =
  SETDEBUGRUNENV =
else
# if DEBUG_LEVEL != 0, build with the debug libraries which are slow.
  SETDEBUGLIB = LIBRARY_PATH=$(HCC2)/lib-debug
# Once you build with debug libs, you have different levels of debug
# by setting these environment variables.
# LIBOMPTARGET_DEBUG: Host device runtime debug messages
  SETDEBUGRUNENV = LIBOMPTARGET_DEBUG=1
# LIBOMPTARGET_DEVICE_RTL_DEBUG: Device runtime debug messages.
#SETDEBUGRUNENV = LIBOMPTARGET_DEVICE_RTL_DEBUG=-1
# ATMI_DEBUG: Print ATMI debug messages
#SETDEBUGRUNENV = ATMI_DEBUG=1
endif

# NOTE: You do NOT need debug libs to turn kernel tracing on.
# Set LIBOMPTARGET_KERNEL_DEBUG anytime to get kernel launch trace.
# Uncomment next line to activate kernel launch trace in this makefile.
#SETDEBUGRUNENV += LIBOMPTARGET_KERNEL_DEBUG=1

.cc.o: lulesh.h
	@echo "Building $<"
	$(CXX) -c $(CXXFLAGS) -o $@  $<

all: $(LULESH_EXEC)

lulesh2.0: $(OBJECTS2.0)
	@echo "Linking"
	$(SETDEBUGLIB) $(CXX) $(CXXFLAGS) $(OBJECTS2.0) -o $@ $(LDFLAGS)

run:	lulesh2.0
	@echo "Running lulesh2.0"
	$(SETDEBUGRUNENV) ./lulesh2.0

clean:
	/bin/rm -f *.o *~ *.tgt* $(OBJECTS) $(LULESH_EXEC)
	/bin/rm -rf *.dSYM

tar: clean
	cd .. ; tar cvf lulesh-2.0.tar LULESH-2.0 ; mv lulesh-2.0.tar LULESH-2.0

