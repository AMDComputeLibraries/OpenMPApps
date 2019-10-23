HPGMG_ARCH := $(if $(PETSC_ARCH),$(PETSC_ARCH),build)

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
ifeq (sm_,$(findstring sm_,$(AOMP_GPU)))
   AOMP_GPUTARGET = nvptx64-nvidia-cuda
endif
AOMP_CPUTARGET ?= x86_64-pc-linux-gnu
ifeq (nvptx,$(findstring nvptx,$(AOMP_GPUTARGET)))
  CUDA ?= /usr/local/cuda
  CUDA_LIBS=-L$(CUDA)/lib64 -lcuda -lcudart -lelf -lffi
endif

DEBUG_LEVEL ?= 0
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

all: build


hpgmgconf.pyc:
	@echo "Running config for directory $(HPGMG_ARCH).  "
	./configure --arch=$(HPGMG_ARCH) --no-fe --no-fv-mpi --CC=$(AOMP)/bin/clang --CFLAGS="-O3 $(host_debug_flags) -target $(AOMP_CPUTARGET) -fopenmp -fopenmp-targets=$(AOMP_GPUTARGET) -Xopenmp-target=$(AOMP_GPUTARGET) -march=$(AOMP_GPU)" --LDFLAGS="-lm -fopenmp -fopenmp-targets=$(AOMP_GPUTARGET) -target $(AOMP_CPUTARGET)" --LDLIBS="$(CUDA_LIBS)"

build:	hpgmgconf.pyc
	$(SETDEBUGLIB) $(MAKE) -C $(HPGMG_ARCH)
	@echo "Build complete in $(HPGMG_ARCH).  Use make run to test"

run:
	$(SETDEBUGRUNENV) ./$(HPGMG_ARCH)/bin/hpgmg-fv 4 128

#./$(HPGMG_ARCH)/bin/hpgmg-fv sample -local 50,10000 -repeat 3 -op_type poisson2

test : all
	$(SETDEBUGLIB) $(MAKE) -C $(HPGMG_ARCH) test

clean :
	$(MAKE) -C $(HPGMG_ARCH) clean ; rm -rf build ; rm hpgmgconf.pyc

.PHONY: all test clean
