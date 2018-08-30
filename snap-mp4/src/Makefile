MP4BUILD = 1
#host_debug_flag = -g
host_debug_flag = 

HCC2 ?= /opt/rocm/hcc2
HCC2_GPU ?= gfx803
HCC2_GPUTARGET = amdgcn-amd-amdhsa
# Default is Radeon fiji/gfx803, for nvidia "export HCC2_GPU=sm_35"
ifeq (sm_,$(findstring sm_,$(HCC2_GPU)))
  HCC2_GPUTARGET = nvptx64-nvidia-cuda
endif
HCC2_CPUTARGET ?= x86_64-pc-linux-gnu
DEBUG_LEVEL ?= 0

ifdef MP4BUILD 
CLANG=$(HCC2)/bin/clang++

CXXFLAGS = $(host_debug_flag) -O3 -I/usr/lib/openmpi/include -target $(HCC2_CPUTARGET) -fopenmp -fopenmp-targets=$(HCC2_GPUTARGET) -Xopenmp-target=$(HCC2_GPUTARGET) -march=$(HCC2_GPU)
CXX = $(HCC2)/bin/clang++
LFLAGS = -fopenmp -fopenmp-targets=$(HCC2_GPUTARGET) -target $(HCC2_CPUTARGET)
ifeq (nvptx,$(findstring nvptx,$(HCC2_GPUTARGET)))
   CUDA ?= /usr/local/cuda
   LFLAGS += -L$(CUDA)/lib64 -lcuda -lcudart -lelf -lffi
endif
CXXLNKFLAGS = $(host_debug_flag)

else
CXXFLAGS = $(host_debug_flag) -O3 -I/usr/lib/openmpi/include 
CXXLNKFLAGS = $(host_debug_flag)
LFLAGS = 
LIBS = -lm 
endif

MPILIBS = -L/usr/lib/openmpi/lib -lmpi

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
	$(SETDEBUGRUNENV) ./snap snap_input
#$(CLANG) $(CXXLNKFLAGS) -o $@ $(MPILIBS) $(LIBS) $(OBJS) $(TFLAGS)

objs/snap_main.o: snap_main.cc mp4_planes.h
	$(CLANG) $(CXXFLAGS) -c snap_main.cc -o $@

objs/%.o: %.cc
	$(CLANG) $(CXXFLAGS) -c $< -o $@


.PHONY: clean

clean:
	rm -f snap objs/*o ssc.cc
