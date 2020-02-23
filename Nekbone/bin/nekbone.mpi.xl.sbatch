#!/bin/bash
#BSUB -P csc235
#BSUB -J nekbone.mpi.xl
#BSUB -o nekbone.mpi.xl.%J.output
#BSUB -e nekbone.mpi.xl.%J.error
#BSUB -W 10
#BSUB -n 10

module load xl

cd "$(pwd)"
mpirun -n 10 --map-by node ./nekbone data
