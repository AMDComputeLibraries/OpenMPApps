#!/bin/bash -x
#COBALT -t 10
#COBALT -q default
#COBALT -A Comp_Perf_Workshop
#COBALT -n 2


echo "COBALT_JOBID" $COBALT_JOBID                                                                   
echo "COBALT_PARTNAME" $COBALT_PARTNAME
echo "COBALT_JOBSIZE" $COBALT_JOBSIZE

# user defines below as well as "-n"
rpn=2  
thr=1
case=data        

# -n  number of nodes  
# -N  --pes-per-node  (ranks per node)
# -d  --cpus-per-pe   depth   cpus (hyperthreads) per rank: nek w/o OpnMP thoreads (d=1)
# -cc --cpu-binding depth
# -j  cpus (hyperthreads) per compute unit (core)   i.e. 1,2,3,4  (j MPI processes per core)


echo $case     >  SESSION.NAME
echo `pwd`'/' >>  SESSION.NAME
touch $case.rea
rm -f ioinfo
mv -f $case.his $case.his1
mv -f $case.sch $case.sch1

aprun -n $((COBALT_JOBSIZE*rpn)) -N $rpn -d $thr -cc depth -j 1 ./nekbone      

exit $?

