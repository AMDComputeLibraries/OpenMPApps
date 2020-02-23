#!/bin/bash
#BSUB -P csc235
#BSUB -J nekbone.mpi.pgi
#BSUB -o nekbone.mpi.pgi.%J.output
#BSUB -e nekbone.mpi.pgi.%J.error
#BSUB -W 10
#BSUB -n 10

module load pgi

cd "$(pwd)"
mpirun -n 10 --map-by node ./nekbone data
