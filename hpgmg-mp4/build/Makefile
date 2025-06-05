HPGMG_ARCH = build
HPGMG_CC = /COD/LATEST/aomp/llvm/bin/clang
HPGMG_CFLAGS = -O3  -target x86_64-unknown-linux-gnu -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a
HPGMG_CPPFLAGS = 
HPGMG_LDFLAGS = -lm -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -target x86_64-unknown-linux-gnu
HPGMG_LDLIBS = 
PETSC_DIR = 
PETSC_ARCH = 
PYTHON = /usr/bin/python
SRCDIR = /work/estewart/git/aomp-test/openmpapps/hpgmg-mp4
CONFIG_FV = y
CONFIG_FV_CPPFLAGS = -DUSE_BICGSTAB=1 -DUSE_SUBCOMM=1 -DUSE_FCYCLES=1 -DUSE_GSRB=1
include $(SRCDIR)/base.mk
