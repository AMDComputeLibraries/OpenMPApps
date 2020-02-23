#include <unistd.h>

#ifndef FNAME_H
#define FNAME_H

/*
   FORTRAN naming convention
     default      cpgs_setup, etc.
     -DUPCASE     CPGS_SETUP, etc.
     -DUNDERSCORE cpgs_setup_, etc.
*/

#ifdef UPCASE
#  define FORTRAN_NAME(low,up) up
#else
#ifdef UNDERSCORE
#  define FORTRAN_NAME(low,up) low##_
#else
#  define FORTRAN_NAME(low,up) low
#endif
#endif

#endif

#define sleepabit FORTRAN_NAME(sleepabit, SLEEPABIT)

void sleepabit(long int *nn)
{

unsigned int usecs;

usecs=*nn;
//printf("test %d %u", *nn, usecs);
usleep(usecs);  
}
