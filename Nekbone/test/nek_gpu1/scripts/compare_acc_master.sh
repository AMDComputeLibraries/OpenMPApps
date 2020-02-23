#!/usr/bin/env bash

./makenek.acc  clean
rm -f compile.acc.output compile.acc.error run.acc.output run.acc.error
./makenek.acc  > compile.acc.output  2> compile.acc.error  || cat compile.acc.error
mpiexec -n 1 ./nekbone data  > run.acc.output  2> run.acc.error  || cat run.acc.error

 vimdiff scripts/master_branch.output run.acc.output

