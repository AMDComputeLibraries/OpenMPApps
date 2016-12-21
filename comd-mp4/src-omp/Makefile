# Makefile for src-openmp

SHELL = /bin/bash

### CoMD can be built in either double or single precision and with or
### without MPI.  Select desired precision and MPI here.

# double precision (ON/OFF)
DOUBLE_PRECISION = ON
# MPI for parallel (ON/OFF)
DO_MPI = OFF

OMPTARGET=nvptx64-nvidia-cuda

### Set your desired C compiler and any necessary flags.  Note that CoMD
### uses some c99 features.  You can also set flags for optimization and
### specify paths to include files that the compiler can't find on its
### own.  If you need any -L or -l switches to get C standard libraries
### (such as -lm for the math library) put them in C_LIB.
#CC = mpicc
CC = $(CLANG)
FLAGS = -O3
FLAGS = -O3 -fopenmp -fopenmp-targets=$(OMPTARGET) -target x86_64-pc-linux-gnu

CFLAGS = -O3 -std=c99
FLAGS += -D__GPU__ -I$(LIBOMP_LIB)
OPTFLAGS = 
INCLUDES = 
C_LIB = -lm
LFLAGS = -L$(LIBOMP_LIB) -L$(CUDA)/lib64  -lcuda -lcudart -lelf -lffi -target x86_64-pc-linux-gnu


### If you need to specify include paths, library paths, or link flags
### for MPI, put them here.  Put both -L and -l switches into MPI_LIB.
MPI_LIB =
MPI_INCLUDE =

### A place to specify any other include or library switches your
### platform requires.
OTHER_LIB =
OTHER_INCLUDE =




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
	${CC} ${CFLAGS} $(FLAGS) -o ${CoMD_EXE} ${OBJECTS} $(LFLAGS) ${LDFLAGS}

CoMD_info.h: Makefile
	./generate_info_header ${CoMD_VARIANT} "$(CC)" "$(CFLAGS)" "$(LDFLAGS)"


${BIN_DIR}:
	@if [ ! -d ${BIN_DIR} ]; then mkdir -p ${BIN_DIR} ; fi

clean:
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
	@makedepend -f .depend -Y. --$(CFLAGS)-- $(SOURCES) 2> /dev/null


-include .depend
