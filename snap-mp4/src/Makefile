MP4BUILD = 1

CLANG=/usr/local/llvm_openmp/bin/clang++

OMPTARGET=nvptx64sm_50-linux-gnu
CPUTARGET=x86_64-pc-linux-gnu

ifdef MP4BUILD 
CXXFLAGS = -g -O3 -fopenmp=libomp -omptargets=$(OMPTARGET) -I/usr/lib/openmpi/include -I$(LIBOMP_LIB)
CXXLNKFLAGS = -g -fopenmp=libomp -omptargets=$(OMPTARGET)
LFLAGS = -target nvptx64sm_50-linux-gnu
LIBS = -L$(LIBOMP_LIB) -L$(CUDA)/lib64 -lm -lcuda -lcudart -lelf -lffi 
else
CXXFLAGS = -g -O3 -I/usr/lib/openmpi/include 
CXXLNKFLAGS = -g 
LFLAGS = 
LIBS = -lm 
endif

MPILIBS = -L/usr/lib/openmpi/lib -lmpi


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
	$(CLANG) $(CXXLNKFLAGS) $^ -o $@ $(LIBS) $(MPILIBS)

objs/snap_main.o: snap_main.cc mp4_planes.h
	$(CLANG) $(CXXFLAGS) -c snap_main.cc -o $@

objs/%.o: %.cc
	$(CLANG) $(CXXFLAGS) -c $< -o $@


.PHONY: clean

clean:
	rm -f snap objs/*o ssc.cc
