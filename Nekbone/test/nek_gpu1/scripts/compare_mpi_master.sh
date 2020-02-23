#!/usr/bin/env bash

./makenek.mpi  clean
rm -f compile.mpi.output compile.mpi.error run.mpi.output run.mpi.error
./makenek.mpi  > compile.mpi.output  2> compile.mpi.error  || cat compile.mpi.error
mpiexec -n 1 ./nekbone data  > run.mpi.output  2> run.mpi.error  || cat run.mpi.error

 vimdiff scripts/master_branch.output run.mpi.output 

