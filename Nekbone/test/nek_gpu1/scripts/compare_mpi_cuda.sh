#!/usr/bin/env bash

./makenek.mpi  clean
rm -f compile.mpi.output compile.mpi.error run.mpi.output run.mpi.error
./makenek.mpi  > compile.mpi.output  2> compile.mpi.error  || cat compile.mpi.error
mpiexec -n 1 ./nekbone data  > run.mpi.output  2> run.mpi.error  || cat run.mpi.error

./makenek.cuda  clean
rm -f compile.cuda.output compile.cuda.error run.cuda.output run.cuda.error
./makenek.cuda  > compile.cuda.output  2> compile.cuda.error  || cat compile.cuda.error
mpiexec -n 1 ./nekbone data  > run.cuda.output  2> run.cuda.error  || cat run.cuda.error

vimdiff run.mpi.output run.cuda.output

