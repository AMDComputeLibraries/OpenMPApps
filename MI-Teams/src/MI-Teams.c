#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>

const int numthreads = 256;
const int numteams = 64;
#ifndef NUMTHREADS
#define NUMTHREADS numthreads
#endif
#ifndef NUMTEAMS
#define NUMTEAMS numteams
#endif

    typedef float float16;

    const uint64_t MB =64*1024;   // global batch size = 64*1024
    const uint64_t E = 128;       // embedding dimension = 128
    const uint64_t M = 20000000;  // table size = 20 million
    const uint64_t P = 256;       // average segment (bag) len
    const uint64_t PT = 1024;// MB ;//* P;  // sum of all segments, = MB * ~P

// tensors:
    float16 g[MB*E];   // gradients, size = MB * E, f16
    float   h[M];      // gradient history, size = M, f32
    int lengths[MB];
    int indices[PT];

    float lr = 0.12345;
    float epsilon = 0.0000033333;
   
int main(int argc, char **argv) {
  struct timespec t0,t1,t2,t3;
  clock_gettime(CLOCK_REALTIME, &t0);
  float *myW = malloc(E*M*4);
// initialize on host
  fprintf(stderr, "Starting MI-Teams\n");
  // time the data copy separate from gflops
  clock_gettime(CLOCK_REALTIME, &t1);
  #pragma omp target data map(tofrom: g[:MB*E], h[:M], lengths[:MB], indices[:PT]) map(from:myW[:M*E])
  {
    // leaving num_threads and teams off, to allow easier env-var control
    #pragma omp target teams num_teams(NUMTEAMS) thread_limit(NUMTHREADS)
  {
    #pragma omp distribute parallel for 
    for (uint64_t i=0; i < MB*E; i++) {
      g[i] = (float) i;
    }
    //#pragma omp target teams num_teams(NUMTEAMS) thread_limit(NUMTHREADS)
    #pragma omp distribute parallel for 
    for (uint64_t i=0; i < MB; i++) {
      lengths[i] = 100 + i % E * 10;
    } 
    //#pragma omp target teams num_teams(NUMTEAMS) thread_limit(NUMTHREADS)
    #pragma omp distribute parallel for 
    for (uint64_t i=0; i < PT; i++) {
      indices[i] = i % M;
    }
   // #pragma omp target teams num_teams(NUMTEAMS) thread_limit(NUMTHREADS)
    #pragma omp distribute parallel for 
    for (uint64_t i=0; i < M; i++) {
      h[i] = i % E + 10;
    }
  }
  clock_gettime(CLOCK_REALTIME, &t2);

    int current = 0;
    float final_sum;
    #pragma omp target teams num_teams(NUMTEAMS) thread_limit(NUMTHREADS)
     //  distribute needs to distribute a loop
      #pragma omp distribute parallel for 
      for (uint64_t sample=0; sample <MB; sample++) {
        int partial_sum = 0;
    //  #pragma omp parallel for // reduction (+:partial_sum)
        for (int e=0; e<E; e++) {
	  float embedding = g[sample*e];
          partial_sum += (embedding * embedding);
	}
        final_sum += partial_sum / E;
        int len = lengths[sample];
        for (uint64_t i=0; i < len; i++) {
	  current++;
	  if (current>=PT) current = 0;
          int idx = indices[current] % M;
          float hi = h[idx] + final_sum;
	  if (idx >= M) printf("Died\n");
          h[idx] = hi ; // store to embedding	     
          float float_step = lr / hi + epsilon;
          //float float_step = lr / sqrt(hi) + epsilon;
          for (uint64_t e=0; e<E; e++)
            myW[idx*E+e] += g[sample*E+e] * float_step;//   # update weights
        }
      }
    clock_gettime(CLOCK_REALTIME, &t3);
  }
  double m = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)/1e9;
  fprintf(stderr, "Time %f for startup\n", m);
  double t = (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec)/1e9;
  fprintf(stderr, "Time %f gpu init\n", t);
  double c = (t2.tv_sec - t0.tv_sec) + (t2.tv_nsec - t0.tv_nsec)/1e9;
  fprintf(stderr, "Time %f for compute\n", c);
  double d = (t2.tv_sec - t0.tv_sec) + (t2.tv_nsec - t0.tv_nsec)/1e9;
  fprintf(stderr, "Time %f for startup + init + compute\n", d);

  fprintf(stderr, "Passed\n");
  return 0;
}
