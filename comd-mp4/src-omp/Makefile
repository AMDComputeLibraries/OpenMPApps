# Makefile for src-openmp

SHELL = /bin/bash
host_debug_flag ?=

### CoMD can be built in either double or single precision and with or
### without MPI.  Select desired precision and MPI here.

# double precision (ON/OFF)
DOUBLE_PRECISION = ON
# MPI for parallel (ON/OFF)
DO_MPI = OFF

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
# Default is Radeon vega/gfx900, for nvidia "export AOMP_GPU=sm_35"
ifeq (sm_,$(findstring sm_,$(AOMP_GPU)))
  AOMP_GPUTARGET = nvptx64-nvidia-cuda
endif
AOMP_CPUTARGET ?= x86_64-pc-linux-gnu
DEBUG_LEVEL ?= 0

### Set your desired C compiler and any necessary flags.  Note that CoMD
### uses some c99 features.  You can also set flags for optimization and
### specify paths to include files that the compiler can't find on its
### own.  If you need any -L or -l switches to get C standard libraries
### (such as -lm for the math library) put them in C_LIB.
#CC = mpicc
CC = $(AOMP)/bin/clang
FLAGS = -O3 -target $(AOMP_CPUTARGET) -fopenmp -fopenmp-targets=$(AOMP_GPUTARGET) -Xopenmp-target=$(AOMP_GPUTARGET) -march=$(AOMP_GPU) $(host_debug_flag)

CFLAGS = -O3 -std=c99
FLAGS += -D__GPU__
OPTFLAGS = 
INCLUDES = 
C_LIB = -lm
LFLAGS = -fopenmp -fopenmp-targets=$(AOMP_GPUTARGET) -target $(AOMP_CPUTARGET)
ifeq (nvptx,$(findstring nvptx,$(AOMP_GPUTARGET)))
  CUDA ?= /usr/local/cuda
  LFLAGS += -L$(CUDA)/lib64 -lcuda -lcudart -lelf -lffi
endif

### If you need to specify include paths, library paths, or link flags
### for MPI, put them here.  Put both -L and -l switches into MPI_LIB.
MPI_LIB =
MPI_INCLUDE =

### A place to specify any other include or library switches your
### platform requires.
OTHER_LIB =
OTHER_INCLUDE =


ifeq ($(DEBUG_LEVEL),0)
  SETDEBUGLIB =
  SETDEBUGRUNENV =
else
# if DEBUG_LEVEL != 0, build with the debug libraries which are slow.
  SETDEBUGLIB = LIBRARY_PATH=$(AOMP)/lib-debug
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
# Set LIBOMPTARGET_KERNEL_TRACE anytime to get kernel launch trace.
# Uncomment next line to activate kernel launch trace in this makefile.
#SETDEBUGRUNENV += LIBOMPTARGET_KERNEL_TRACE=1


#########################################
### Below here, it is pitch black.  
### You are likely to be eaten by a grue.
##########################################

# clear all suffixes
.SUFFIXES:
# list only those that we use 
.SUFFIXES: .c .o

.PHONY: DEFAULT clean distclean depend

BIN_DIR=../bin

# Check for double precision
ifeq ($(DOUBLE_PRECISION), ON)
CFLAGS += -DDOUBLE
else
CFLAGS += -DSINGLE
endif

# Set executable name and add includes & libraries for MPI if needed.
CoMD_VARIANT = CoMD-openmp
ifeq ($(DO_MPI), ON)
CoMD_VARIANT = CoMD-openmp-mpi
INCLUDES += ${MPI_INCLUDE}
CFLAGS += -DDO_MPI
LDFLAGS += ${MPI_LIB}
endif
CoMD_EXE = ${BIN_DIR}/${CoMD_VARIANT}

LDFLAGS += ${C_LIB} ${OTHER_LIB}
CFLAGS  += ${OPTFLAGS} ${INCLUDES} ${OTHER_INCLUDE}  


SOURCES=$(wildcard *.c)
OBJECTS=$(SOURCES:.c=.o)


DEFAULT: ${CoMD_EXE}

%.o: %.c
	${CC} ${CFLAGS} ${FLAGS} -c $< -o $@

ljForce.o: ljForce.c
	${CC} ${CFLAGS} $(FLAGS) -c $< -o $@

${CoMD_EXE}: ${BIN_DIR} CoMD_info.h ${OBJECTS} 
	$(SETDEBUGLIB) ${CC} ${CFLAGS} $(FLAGS) -o ${CoMD_EXE} ${OBJECTS} $(LFLAGS) ${LDFLAGS}

run: ${CoMD_EXE}
	$(SETDEBUGRUNENV) ${CoMD_EXE}

CoMD_info.h: Makefile
	./generate_info_header ${CoMD_VARIANT} "$(CC)" "$(CFLAGS)" "$(LDFLAGS)"


${BIN_DIR}:
	@if [ ! -d ${BIN_DIR} ]; then mkdir -p ${BIN_DIR} ; fi

clean:
	rm -f ${CoMD_EXE} .depend.bak
	rm -f *.o CoMD_info.h .depend

distclean: clean
	rm -f ${CoMD_EXE} .depend.bak
	rm -rf html latex

.depend: $(SOURCES)
	@touch .depend
	@$(MAKE) --no-print-directory depend

depend:
	@echo "Rebuilding dependencies..."
	@$(MAKE) CoMD_info.h
	@makedepend -f .depend -Y. --$(CFLAGS)-- $(SOURCES) 2> /dev/null ; true


-include .depend
