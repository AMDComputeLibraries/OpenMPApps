#===============================================================================
# User Options
#===============================================================================

COMPILER    ?= aomp
OPTIMIZE    = no
PROFILE     = no
MPI         ?= no
VEC_INFO    = no
VERIFY      = yes
BINARY_DUMP = no
BINARY_READ = no

ifeq ($(DEBUG),no)
  DEBUG_LEVEL = 0
else
ifeq ($(DEBUG),yes)
  DEBUG_LEVEL ?= 1
endif
endif
DEBUG_LEVEL ?= 0

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

# Default is Radeon vega/gfx900, for nvidia "export AOMP_GPU=sm_35"
ifeq ($(COMPILER),aomp)
# --- Standard Makefile check for AOMP installation ---
ifeq ("$(wildcard $(AOMP))","")
  ifneq ($(AOMP),)
    $(warning AOMP not found at $(AOMP))
  endif
  AOMP = $(HOME)/rocm/aomp
  ifeq ("$(wildcard $(AOMP))","")
    $(warning AOMP not found at $(AOMP))
    AOMP = /opt/rocm/aomp
    ifeq ("$(wildcard $(AOMP))","")
      $(warning AOMP not found at $(AOMP))
      $(error Please install AOMP or correctly set env-var AOMP)
    endif
  endif
endif
# --- End Standard Makefile check for AOMP installation ---
  INSTALLED_GPU  = $(shell $(AOMP)/bin/mygpu -d gfx900)
  AOMP_GPU ?= $(INSTALLED_GPU)
  AOMP_GPUTARGET = amdgcn-amd-amdhsa
  ifeq (sm_,$(findstring sm_,$(AOMP_GPU)))
    AOMP_GPUTARGET = nvptx64-nvidia-cuda
  endif
  AOMP_CPUTARGET ?= x86_64-pc-linux-gnu
  CFLAGS = -O3 -target $(AOMP_CPUTARGET) -fopenmp -fopenmp-targets=$(AOMP_GPUTARGET) -Xopenmp-target=$(AOMP_GPUTARGET) -march=$(AOMP_GPU) -D__GPU__
  CC = $(AOMP)/bin/clang
  LDFLAGS = -fopenmp -fopenmp-targets=$(AOMP_GPUTARGET) -target $(AOMP_CPUTARGET) -lm
  ifeq (nvptx,$(findstring nvptx,$(AOMP_GPUTARGET)))
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
ifneq ($(DEBUG_LEVEL),0)

ifeq ($(COMPILER),aomp)
  SETDEBUGLIB = LIBRARY_PATH=$(AOMP)/lib-debug
# Once you build with debug libs, you have different levels of debug
# by setting these environment variables.
# LIBOMPTARGET_DEBUG: Host device runtime debug messages
  SETDEBUGRUNENV = LIBOMPTARGET_DEBUG=1
# LIBOMPTARGET_DEVICE_RTL_DEBUG: Device runtime debug messages.
#SETDEBUGRUNENV = LIBOMPTARGET_DEVICE_RTL_DEBUG=-1
# ATMI_DEBUG: Print ATMI debug messages
#SETDEBUGRUNENV = ATMI_DEBUG=1
else
  CFLAGS += -g
  LDFLAGS  += -g
endif

else
  SETDEBUGLIB =
  SETDEBUGRUNENV =
endif
# NOTE: You do NOT need debug libs to turn kernel tracing on.
# Set LIBOMPTARGET_KERNEL_TRACE anytime to get kernel launch trace.
# Uncomment next line to activate kernel launch trace in this makefile.
#SETDEBUGRUNENV += LIBOMPTARGET_KERNEL_TRACE=1

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
  MPI_INC = $(HOME)/openmpi-4.0.3-install/include
  MPI_LIB = $(HOME)/openmpi-4.0.3-install/lib
  MPI_RUN = $(HOME)/openmpi-4.0.3-install/bin/mpirun
  CC      = $(AOMP)/bin/clang
  CFLAGS  += -DMPI -I$(MPI_INC)
  LDFLAGS += -L$(MPI_LIB) -lmpi
  SHIM = $(AOMP)/bin/openmpi_set_cu_mask
  RANKS ?= 1
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

CFLAGS += $(host_debug_flag)

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

mpirun: $(program)
	@echo "Running mpirun -np $(RANKS)  $(program)"
	LD_LIBRARY_PATH=$(MPI_LIB) $(SETDEBUGRUNENV) $(MPI_RUN) -np $(RANKS) $(SHIM) $(PWD)/$(program)

bgqrun: XSBench
	qsub -t 10 -n 1 -O test XSBench
