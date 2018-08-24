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

all :
	./configure --arch=$(HPGMG_ARCH) --no-fe --no-fv-mpi --CC=$(HCC2)/bin/clang --CFLAGS="-O3 -target $(HCC2_CPUTARGET) -fopenmp -fopenmp-targets=$(HCC2_GPUTARGET) -Xopenmp-target=$(HCC2_GPUTARGET) -march=$(HCC2_GPU)" --LDFLAGS="-fopenmp -fopenmp-targets=$(HCC2_GPUTARGET) -target $(HCC2_CPUTARGET)" --LDLIBS="$(CUDA_LIBS)"
	$(MAKE) -C $(HPGMG_ARCH)
	@echo "Build complete in $(HPGMG_ARCH).  Use make -C $(HPGMG_ARCH) test to test."

run:
	./$(HPGMG_ARCH)/bin/hpgmg-fv 4 128

#./$(HPGMG_ARCH)/bin/hpgmg-fv sample -local 50,10000 -repeat 3 -op_type poisson2

test : all
	$(MAKE) -C $(HPGMG_ARCH) test

clean :
	$(MAKE) -C $(HPGMG_ARCH) clean

.PHONY: all test clean
