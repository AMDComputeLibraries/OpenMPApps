#!/usr/bin/env bash

./makenek.mpi  clean
rm -f compile.mpi.output compile.mpi.error run.mpi.output run.mpi.error
./makenek.mpi  > compile.mpi.output  2> compile.mpi.error  || cat compile.mpi.error
mpiexec -n 1 ./nekbone data  > run.mpi.output  2> run.mpi.error  || cat run.mpi.error

./makenek.acc  clean
rm -f compile.acc.output compile.acc.error run.acc.output run.acc.error
./makenek.acc  > compile.acc.output  2> compile.acc.error  || cat compile.acc.error
mpiexec -n 1 ./nekbone data  > run.acc.output  2> run.acc.error  || cat run.acc.error

vimdiff run.mpi.output run.acc.output
