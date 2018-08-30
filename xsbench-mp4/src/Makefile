#===============================================================================
# User Options
#===============================================================================

COMPILER    ?= hcc2
OPTIMIZE    = no
DEBUG_LEVEL ?= 0
PROFILE     = no
MPI         = no
VEC_INFO    = no
VERIFY      = no
BINARY_DUMP = no
BINARY_READ = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = XSBench

source = \
Main.c \
io.c \
CalculateXS.c \
GridInit.c \
XSutils.c \
Materials.c

obj = $(source:.c=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Linker Flags
LDFLAGS = -lm

ifeq ($(COMPILER),mp4)
  OMPTARGET=nvptx64-nvidia-cuda
  CPUTARGET=x86_64-pc-linux-gnu
  CFLAGS = -O3 -fopenmp -fopenmp-targets=$(OMPTARGET) -target $(CPUTARGET)
  CC = $(CLANG)
  CFLAGS += -D__GPU__ -I$(LIBOMP_LIB)
  LDFLAGS = -fopenmp -fopenmp-targets=$(OMPTARGET) -target $(CPUTARGET) -L$(LIBOMP_LIB) -L$(CUDA)/lib64 -lcuda -lcudart -lelf -lffi
endif

# Default is Radeon fiji/gfx803, for nvidia "export HCC2_GPU=sm_35"
ifeq ($(COMPILER),hcc2)
  HCC2 ?= /opt/rocm/hcc2
  HCC2_GPU ?= gfx803
  HCC2_GPUTARGET = amdgcn-amd-amdhsa
  ifeq (sm_,$(findstring sm_,$(HCC2_GPU)))
    HCC2_GPUTARGET = nvptx64-nvidia-cuda
  endif
  HCC2_CPUTARGET ?= x86_64-pc-linux-gnu
  CFLAGS = -O3 -target $(HCC2_CPUTARGET) -fopenmp -fopenmp-targets=$(HCC2_GPUTARGET) -Xopenmp-target=$(HCC2_GPUTARGET) -march=$(HCC2_GPU) -D__GPU__
  CC = $(HCC2)/bin/clang
  LDFLAGS = -fopenmp -fopenmp-targets=$(HCC2_GPUTARGET) -target $(HCC2_CPUTARGET) -lm
  ifeq (nvptx,$(findstring nvptx,$(HCC2_GPUTARGET)))
    CUDA ?= /usr/local/cuda
    LDFLAGS += -L$(CUDA)/lib64 -lcuda -lcudart -lelf -lffi
  endif
endif

# GNU Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
  CFLAGS += -std=gnu99 -Wall -fopenmp
  # Optimization
  ifeq ($(OPTIMIZE),yes)
    CFLAGS += -O3
    # Compiler Vectorization (needs -O3 flag) information
    ifeq ($(VEC_INFO),yes)
      CFLAGS += -ftree-vectorizer-verbose=6
    endif
  endif
endif

# Intel Compiler
ifeq ($(COMPILER),intel)
  CC = icc
  CFLAGS += -std=gnu99 -w1 -openmp
  # Optimization
  ifeq ($(OPTIMIZE),yes)
    CFLAGS += -O3
    # Compiler Vectorization (needs -O3 flag) information
    ifeq ($(VEC_INFO),yes)
      CFLAGS += -ftree-vectorizer-verbose=6
    endif
  endif
endif


# BG/Q gcc Cross-Compiler
ifeq ($(MACHINE),bluegene)
  CC = mpicc
endif

# Debug Flags
ifneq ($(word 1,$(DEBUG_LEVEL)),0)

ifeq ($(COMPILER),hcc2)
  SETDEBUGLIB = LIBRARY_PATH=$(HCC2)/lib-debug
  SETDEBUGRUNENV = LIBOMPTARGET_DEBUG=1
#  once yiou build with the debug libs you habve lots of switches you can turn on  
#  with these environment variablse
# LIBOMPTARGET_DEBUG=1
# LIBOMPTARGET_DEVICE_RTL_DEBUG=-1
# ATMI_DEBUG=1
else
  CFLAGS += -g
  LDFLAGS  += -g
endif

else
  SETDEBUGLIB = 
  SETDEBUGRUNENV =
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  CFLAGS += -pg
  LDFLAGS  += -pg
endif

# PAPI source (you may need to provide -I and -L pointing
# to PAPI depending on your installation
ifeq ($(PAPI),yes)
  source += papi.c
  CFLAGS += -DPAPI
  #CFLAGS += -I/soft/apps/packages/papi/papi-5.1.1/include
  #LDFLAGS += -L/soft/apps/packages/papi/papi-5.1.1/lib -lpapi
  LDFLAGS += -lpapi
endif

# MPI
ifeq ($(MPI),yes)
  CC = mpicc
  CFLAGS += -DMPI
endif

# Verification of results mode
ifeq ($(VERIFY),yes)
  CFLAGS += -DVERIFICATION
endif

# Adds outer 'benchmarking' loop to do multiple trials for
# 1 < threads <= max_threads
ifeq ($(BENCHMARK),yes)
  CFLAGS += -DBENCHMARK
endif

# Binary dump for file I/O based initialization
ifeq ($(BINARY_DUMP),yes)
  CFLAGS += -DBINARY_DUMP
endif

# Binary read for file I/O based initialization
ifeq ($(BINARY_READ),yes)
  CFLAGS += -DBINARY_READ
endif


#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) XSbench_header.h
	$(SETDEBUGLIB) $(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj) *.o.*

edit:
	vim -p $(source) papi.c XSbench_header.h

run: $(program)
	$(SETDEBUGRUNENV) ./$(program) -s small

bgqrun: XSBench
	qsub -t 10 -n 1 -O test XSBench
