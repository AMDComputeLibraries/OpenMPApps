HPGMG_ARCH := $(if $(PETSC_ARCH),$(PETSC_ARCH),build)

HCC2 ?= /opt/rocm/hcc2
HCC2_GPU ?= gfx803
HCC2_GPUTARGET = amdgcn-amd-amdhsa
ifeq (sm_,$(findstring sm_,$(HCC2_GPU)))
   HCC2_GPUTARGET = nvptx64-nvidia-cuda
endif
HCC2_CPUTARGET ?= x86_64-pc-linux-gnu
ifeq (nvptx,$(findstring nvptx,$(HCC2_GPUTARGET)))
  CUDA ?= /usr/local/cuda
  CUDA_LIBS=-L$(CUDA)/lib64 -lcuda -lcudart -lelf -lffi
endif

DEBUG_LEVEL ?= 0
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

all: build


hpgmgconf.pyc:
	@echo "Running config for directory $(HPGMG_ARCH).  "
	./configure --arch=$(HPGMG_ARCH) --no-fe --no-fv-mpi --CC=$(HCC2)/bin/clang --CFLAGS="-O3 -target $(HCC2_CPUTARGET) -fopenmp -fopenmp-targets=$(HCC2_GPUTARGET) -Xopenmp-target=$(HCC2_GPUTARGET) -march=$(HCC2_GPU)" --LDFLAGS="-lm -fopenmp -fopenmp-targets=$(HCC2_GPUTARGET) -target $(HCC2_CPUTARGET)" --LDLIBS="$(CUDA_LIBS)"

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
