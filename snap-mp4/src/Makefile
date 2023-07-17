MP4BUILD = 1
#host_debug_flag = -g
host_debug_flag ?=

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

ifdef MP4BUILD 
CLANG=$(AOMP)/bin/clang++

OMPIDIR ?= $(HOME)/local/openmpi
ifeq "$(wildcard $(OMPIDIR))" ""
  OMPIDIR = /opt/openmpi-4.1.5
  ifeq "$(wildcard $(OMPIDIR))" ""
    OMPIDIR = /opt/openmpi-4.1.4
    ifeq "$(wildcard $(OMPIDIR))" ""
      OMPIDIR = /usr/local/openmpi
      ifeq "$(wildcard $(OMPIDIR))" ""
        OMPIDIR = /usr/lib/x86_64-linux-gnu/openmpi
        ifeq "$(wildcard $(OMPIDIR)/lib/libmpi*)" ""
          OMPIDIR = /usr/lib/openmpi
          ifeq "$(wildcard $(OMPIDIR)/lib/libmpi*)" ""
            OMPIDIR = /usr/lib64/openmpi
            ifeq "$(wildcard $(OMPIDIR)/lib/libmpi*)" ""
              OMPIDIR = /usr/local
              ifeq "$(wildcard $(OMPIDIR)/lib/libmpi*)" ""
                OMPIDIR = /usr/local
                MPI64 = 64
                ifeq "$(wildcard $(OMPIDIR)/lib$(MPI64)/libmpi*)" ""
                  $(error Openmpi not found on the system and is required. Set OMPIDIR to openmpi installation)
                endif
              endif
            endif
          endif
        endif
      endif
    endif
  endif
endif

$(info -- openmpi found at $(OMPIDIR))
MPILIBS ?= -L${OMPIDIR}/lib$(MPI64) -lmpi

CXXFLAGS = $(host_debug_flag) -O3 -I${OMPIDIR} -I${OMPIDIR}/include -target $(AOMP_CPUTARGET) -fopenmp -fopenmp-targets=$(AOMP_GPUTARGET) -Xopenmp-target=$(AOMP_GPUTARGET) -march=$(AOMP_GPU)
CXX = $(AOMP)/bin/clang++
LFLAGS = -fopenmp -fopenmp-targets=$(AOMP_GPUTARGET) -target $(AOMP_CPUTARGET)
ifeq (nvptx,$(findstring nvptx,$(AOMP_GPUTARGET)))
   CUDA ?= /usr/local/cuda
   LFLAGS += -L$(CUDA)/lib64 -lcuda -lcudart -lelf -lffi
endif
CXXLNKFLAGS = $(host_debug_flag)

else
CXXFLAGS = $(host_debug_flag) -O3 -I${OMPIDIR/include
CXXLNKFLAGS = $(host_debug_flag)
LFLAGS = 
LIBS = -lm 
endif

MPILIBS ?= -L${OMPIDIR}/lib -lmpi

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

SRCS = 	snap_main.cc \
	input.cc \
	allocate.cc \
	comms.cc \
	problem.cc \
	ssc.cc \
	convergence.cc \
	population.cc \
	profiler.cc


OBJS = $(SRCS:%.cc=objs/%.o)
all: objs snap

ssc.cc: buffers.cc sweep.cc source.cc scalar_flux.cc 
	cat buffers.cc source.cc scalar_flux.cc sweep.cc > ssc.cc
objs:
	mkdir objs

snap: $(OBJS)
	@echo
	@echo linking snap
	@echo
	$(SETDEBUGLIB) $(CLANG) $(CXXFLAGS) $(OBJS) $(CXXLNKFLAGS) $(MPILIBS) -o $@ $(LIBS) $(LFLAGS)

run: snap
	#$(SETDEBUGRUNENV) OMP_MAX_TEAMS_DEFAULT=128 ./snap snap_input
	LD_LIBRARY_PATH=$(OMPIDIR)/lib:$(LD_LIBRARY_PATH) $(SETDEBUGRUNENV) ./snap snap_input
#$(CLANG) $(CXXLNKFLAGS) -o $@ $(MPILIBS) $(LIBS) $(OBJS) $(TFLAGS)

objs/snap_main.o: snap_main.cc mp4_planes.h
	$(CLANG) $(CXXFLAGS) -c snap_main.cc -o $@

objs/%.o: %.cc
	$(CLANG) $(CXXFLAGS) -c $< -o $@


.PHONY: clean

clean:
	rm -f snap objs/*o ssc.cc
