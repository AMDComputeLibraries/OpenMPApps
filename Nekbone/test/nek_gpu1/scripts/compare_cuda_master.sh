#!/usr/bin/env bash

./makenek.cuda  clean
rm -f compile.cuda.output compile.cuda.error run.cuda.output run.cuda.error
./makenek.cuda  > compile.cuda.output  2> compile.cuda.error  || cat compile.cuda.error
mpiexec -n 1 ./nekbone data  > run.cuda.output  2> run.cuda.error  || cat run.cuda.error

vimdiff scripts/master_branch.output run.cuda.output

